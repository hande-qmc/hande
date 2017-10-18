module blocking

use qmc_data, only: dt_numerator, dt_denominator, dt_shift, dt_proj_energy
use const, only: p, depsilon

! Module for performing reblocking on the fly.

! References
! ----------
!   [1] "Monte Carlo errors with less errors", U. Wolff, Comput. Phys.
!   Commun. 156, 143 (2004) and arXiv:hep-lat/0306017.
!   [2] "Strategies for improving the efficiency of quantum Monte Carlo
!   calculations", R. M. Lee, G. J. Conduit, N. Nemec, P. Lopez Rios, and N.
!   D.  Drummond, Phys. Rev. E. 83, 066706 (2011).
!   [3] "Error estimates on averages of correlated data", H. Flyvbjerg and
!   H.G. Petersen, J. Chem. Phys. 91, 461 (1989).

implicit none

contains

    subroutine allocate_blocking(qmc_in, blocking_in, bl)

        ! Allocate different arrays in blocking_t data according to the input
        ! options.

        ! In:
        !   qmc_in: Input options relating to QMC methods.
        !   blocking_in: input options for blocking on the fly.
        ! In/Out:
        !   bl: Information needed to peform blocking on the fly. Maximum block
        !       size, frequency at which the data for changing start position
        !       and number of possible start positions are determined

        use qmc_data, only: qmc_in_t, blocking_t, blocking_in_t
        use checking, only: check_allocate

        type(qmc_in_t), intent(in) :: qmc_in
        type(blocking_t), intent(inout) :: bl
        type(blocking_in_t), intent(in) :: blocking_in
        integer :: ierr

        ! Set 2^(lg_max) to a value that is approximately 4 times larger than nreport.
        bl%lg_max = nint(log(real(qmc_in%nreport))/log(2.0)) + 2

        ! Set save_fq to 2^(lg_max - 10) when start_save_frequency in blocking_in
        ! is negative (Default). If otherwise, set save_fq to 2^(start_save_frequency).
        if (blocking_in%start_save_frequency < 0) then
            bl%save_fq = 2 ** (bl%lg_max - 10)
        else
            bl%save_fq = 2 ** (blocking_in%start_save_frequency)
        end if

        ! number of saved startpoints is enough to cover the nreport in qmc_in
        ! (default) unless specified in blocking_in.

        if (blocking_in%start_point_number < 0) then
            bl%n_saved_startpoints = int(qmc_in%nreport/(bl%save_fq))
        else
            bl%n_saved_startpoints = blocking_in%start_point_number
        end if

        ! Each dimension and what they mean are listed in qmc_data.f90
        ! blocking_t and reblock_data_t.

        allocate(bl%reblock_data(3, 0:bl%lg_max), stat=ierr)
        call check_allocate('bl%reblock_data',(bl%lg_max+1)*3,ierr)
        allocate(bl%data_product(0:bl%lg_max), stat=ierr)
        call check_allocate('bl%data_product',(bl%lg_max+1),ierr)
        allocate(bl%reblock_data_2(3, 0:bl%lg_max), stat=ierr)
        call check_allocate('bl%reblock_data_2', (bl%lg_max+1)*3,ierr)
        allocate(bl%data_product_2(0:bl%lg_max), stat=ierr)
        call check_allocate('bl%data_product_2',(bl%lg_max+1),ierr)

        ! Mean, standard deviation and covariance of Proj. energy and reference
        ! population is calculated for each block size and added to block_mean,
        ! block_std and block_cov. Only mean and standard deviation is
        ! calculated for the shift.

        allocate(bl%block_mean(3, 0:bl%lg_max), stat=ierr)
        call check_allocate('bl%block_mean', 3*(bl%lg_max+1), ierr)
        allocate(bl%block_std(3, 0:bl%lg_max), stat=ierr)
        call check_allocate('bl%block_std', 3*(bl%lg_max+1), ierr)
        allocate(bl%block_cov(0:bl%lg_max), stat=ierr)
        call check_allocate('bl%block_cov', bl%lg_max+1, ierr)

        ! Array in which the reblock_data and data_product are saved every
        ! start_fq.

        allocate(bl%reblock_save(3, 0:bl%lg_max, 0:bl%n_saved_startpoints) , stat=ierr)
        call check_allocate('bl%reblock_save',(bl%n_saved_startpoints+1)*3*(bl%lg_max+1), ierr)
        allocate(bl%product_save(0:bl%lg_max, 0:bl%n_saved_startpoints),stat=ierr)
        call check_allocate('bl%product_save',(bl%n_saved_startpoints+1)*(bl%lg_max+1),ierr)

        ! 1/(sqrt(number of data)) * fractional error for both \sum H_0j N_jand
        ! reference population is saved for comparison.

        allocate(bl%err_compare(3, 0:bl%n_saved_startpoints))
        call check_allocate('bl%err_compare', (bl%n_saved_startpoints+1)*3, ierr)

        bl%data_product = 0
        bl%data_product_2 = 0
        bl%block_mean = 0
        bl%block_std = 0
        bl%block_cov = 0
        bl%product_save = 0
        bl%err_compare = 0

    end subroutine allocate_blocking

    subroutine deallocate_blocking(bl)

        ! Deallocate different arrays in blocking_t data.

        ! In/Out:
        !   bl: Information needed to peform blocking on the fly. It is
        !   unchanged in the output.

        use qmc_data, only: blocking_t
        use checking, only: check_deallocate

        integer :: ierr
        type(blocking_t), intent(inout) :: bl

        deallocate(bl%reblock_data, stat=ierr)
        call check_deallocate('bl%reblock_data', ierr)
        deallocate(bl%reblock_data_2, stat=ierr)
        call check_deallocate('bl%reblock_data_2', ierr)
        deallocate(bl%data_product, stat=ierr)
        call check_deallocate('bl%data_product', ierr)
        deallocate(bl%data_product_2, stat=ierr)
        call check_deallocate('bl%data_product_2', ierr)
        deallocate(bl%block_mean, stat=ierr)
        call check_deallocate('bl%block_mean', ierr)
        deallocate(bl%block_std, stat=ierr)
        call check_deallocate('bl%block_std', ierr)
        deallocate(bl%block_cov, stat=ierr)
        call check_deallocate('bl%block_cov', ierr)
        deallocate(bl%reblock_save, stat=ierr)
        call check_deallocate('bl%reblock_save', ierr)
        deallocate(bl%product_save, stat=ierr)
        call check_deallocate('bl%product_save', ierr)
        deallocate(bl%err_compare, stat=ierr)
        call check_deallocate('bl%err_compare', ierr)

    end subroutine deallocate_blocking

    subroutine write_blocking_report_header(iunit)

        ! Write the block analysis report to a file identified by iunit.

        ! In:
        !   iunit: identifies the file to which the header is written out to.

        integer, intent(in) :: iunit

        write(iunit, '(1X, "#iterations")', advance = 'no')
        write(iunit, '(1X, "Start")', advance = 'no')
        write(iunit, '(2X, "Mean \sum H_0j N_j")', advance = 'no')
        write(iunit, '(1X, "Std \sum H_0 N_j")', advance = 'no')
        write(iunit, '(4X, "Error in error")', advance = 'no')
        write(iunit, '(3X, "Mean N_0")', advance = 'no')
        write(iunit, '(9X, "Std N_0")', advance = 'no')
        write(iunit, '(10 X, "Error in error")', advance = 'no')
        write(iunit, '(3X, "Mean Shift")', advance = 'no')
        write(iunit, '(7X, "Std Shift")', advance = 'no')
        write(iunit, '(8X, "Error in error")', advance = 'no')
        write(iunit, '(2X, "Mean Proj. energy")', advance = 'no')
        write(iunit, '(1X, "Std Proj. energy")', advance = 'no')
        write(iunit, '(1X, "Error in error")')
        flush(iunit)

    end subroutine write_blocking_report_header

    subroutine collect_data(qmc_in, qs, bl, ireport)

        ! Collecting data for the reblock analysis

        ! In:
        !   qmc_in: input options relating to QMC methods.
        !   qs: qmc_state where the data for current iteration is taken
        !   ireport: number of reports
        ! In/Out:
        !   bl: Information needed to peform blocking on the fly. The data is
        !       collected every report and reblock_data and data_product is updated.

        use qmc_data, only: qmc_in_t, qmc_state_t, blocking_t
        use const, only: int_64

        type(qmc_in_t), intent(in) :: qmc_in
        type(qmc_state_t), intent(in) :: qs
        integer, intent(in) :: ireport
        type(blocking_t), intent(inout) :: bl
        integer :: i
        integer(int_64) :: reblock_size

        ! \sum H_0j N_j, reference population and shift are added to
        ! data_accumulator of every block size.

        bl%reblock_data(dt_numerator,:)%data_accumulator = bl%reblock_data(dt_numerator,:)%data_accumulator + &
                                                                                        qs%estimators(1)%proj_energy
        bl%reblock_data(dt_denominator,:)%data_accumulator = bl%reblock_data(dt_denominator,:)%data_accumulator + &
                                                                                        qs%estimators(1)%D0_population
        bl%reblock_data(dt_shift,:)%data_accumulator = bl%reblock_data(dt_shift,:)%data_accumulator + qs%shift(1)

        ! Everytime enough data is collected for each block size, the
        ! data_accumulator is divied by the block size and added to
        ! sum_of_blocks and squared and added to sum_of_block_squares.

        do i = 0, bl%lg_max
            if (mod(bl%n_reports_blocked, 2_int_64 ** i) == 0) then
                reblock_size = 2_int_64 ** i

                bl%reblock_data(:,i)%n_blocks = int(bl%n_reports_blocked/reblock_size)

                bl%reblock_data(:,i)%sum_of_blocks = bl%reblock_data(:,i)%sum_of_blocks + bl%reblock_data(:,i)%data_accumulator/&
                                                                                                                     reblock_size

                bl%reblock_data(:,i)%sum_of_block_squares = bl%reblock_data(:,i)%sum_of_block_squares + &
                                                                (bl%reblock_data(:,i)%data_accumulator/reblock_size) ** 2

                bl%data_product(i) = bl%data_product(i) + (bl%reblock_data(dt_numerator,i)%data_accumulator/reblock_size) * &
                                                          (bl%reblock_data(dt_denominator,i)%data_accumulator/reblock_size)
                bl%reblock_data(:,i)%data_accumulator = 0
            end if
        end do
    end subroutine collect_data

    subroutine mean_std_cov(bl, energy_estimate_dist)

        ! Mean, standard deviation and covariance for each block size is
        ! calculated from the collected data.

        ! In/Out:
        !   bl: Information needed to peform blocking on the fly. block_mean,
        !       block_std and block_cov for each block size is calculated. 0 is
        !       returned if there aren't sufficient blocks to calculate them.
        ! Out (optional):
        !   sd_shift_dist: standard distribution of the shift distribution.
        !       Obtained as sqrt(<S^2>-<S>^2) for instantaneous shift values.
        !   sd_proj_numer_dist: standard distribution of the projected energy
        !       numerator distribution.
        !   sd_proj_denom_dist: stardard deviation of the instantaneous reference
        !       population distribution.

        use qmc_data, only: blocking_t
        use const, only: p

        type(blocking_t), intent(inout) :: bl
        integer :: i
        real(p), intent(out), optional :: energy_estimate_dist(dt_numerator:dt_shift,2)

        do i = 0, bl%lg_max

            if (bl%reblock_data_2(dt_numerator,i)%n_blocks > 0) then
                bl%block_mean(:,i) = bl%reblock_data_2(:,i)%sum_of_blocks/bl%reblock_data_2(:,i)%n_blocks

            else
                bl%block_mean(:,i) = 0.0_p
            end if
            if (bl%reblock_data_2(dt_numerator,i)%n_blocks > 1) then
                bl%block_std(:,i) = (bl%reblock_data_2(:,i)%sum_of_block_squares/bl%reblock_data_2(:,i)%n_blocks - &
                                                                                                bl%block_mean(:,i) ** 2)
                if (i==0 .and. present(energy_estimate_dist)) then
                    energy_estimate_dist(dt_numerator,1) = bl%block_mean(dt_numerator,i)
                    energy_estimate_dist(dt_numerator,2) = sqrt(bl%block_std(dt_numerator,i))

                    energy_estimate_dist(dt_denominator,1) = bl%block_mean(dt_denominator,i)
                    energy_estimate_dist(dt_denominator,2) = sqrt(bl%block_std(dt_denominator,i))

                    energy_estimate_dist(dt_shift,1) = bl%block_mean(dt_shift,i)
                    energy_estimate_dist(dt_shift,2) = sqrt(bl%block_std(dt_shift,i))
                end if
                bl%block_std(:,i) = sqrt(bl%block_std(:,i)/(bl%reblock_data_2(:,i)%n_blocks - 1))

                bl%block_cov(i) = (bl%data_product_2(i) - (bl%block_mean(1,i)*bl%block_mean(2,i)  * &
                                                                    bl%reblock_data_2(dt_numerator,i)%n_blocks))/&
                                                                    (bl%reblock_data_2(dt_numerator,i)%n_blocks - 1)
            else
                bl%block_std(:,i) = 0.0_p
                bl%block_cov(i) = 0.0_p
                if (i==0 .and. present(energy_estimate_dist)) then
                    energy_estimate_dist = 0.0_p
                end if
            end if
        end do

    end subroutine mean_std_cov

    pure function fraction_error(mean_1, mean_2, data_number, std_1, std_2, cov_in) result(error_est)

        ! Function to calculate the error of a fraction when the error and mean
        ! of the denominator and numerator is known.

        ! In:
        !   mean_1: Mean of the numerator.
        !   mean_2: Mean of the denominator.
        !   data_number: Number of data for which the mean and standard
        !       deviation is calculated from.
        !   std_1: Standard deviation of the numerator.
        !   std_2: Standard deviation of the denominator.
        !   cov_in: Covariance between the numerator and the denominator.
        ! Returns:
        !   error_est: The estimated standard deviation of the fraction.
        use const

        real(p) :: error_est
        real(p), intent(in) :: mean_1
        real(p), intent(in) :: mean_2
        integer, intent(in) :: data_number
        real(p), intent(in) :: std_1
        real(p), intent(in) :: std_2
        real(p), intent(in) :: cov_in
        real(p) :: mean_cur

        mean_cur = mean_1/mean_2
        error_est = abs(mean_cur*sqrt((std_1/mean_1)**2 + (std_2/mean_2)**2 - 2*cov_in/(data_number*mean_1*mean_2)))

    end function fraction_error

    subroutine find_optimal_block(bl)

        ! Finds the optimal mean and the standard deviation satisfying the
        ! condition B^3 > 2 * (B * (number of blocks)) * (std(B)/std(0))^4
        ! as suggested by Wolff [1] and Lee et al [2].
        ! Where B is the block size and std(B) is the standard deviation
        ! calculated at block size B.
        ! Returns 0 if none of the block sizes satisfy the condition.
        ! Optimal block size for the fraction between \sum H_0j N_j and reference
        ! population is the larger optimal block size between the optimal block
        ! size of \sum H_0j N_j and reference population.

        ! References
        ! ----------
        !   [1] "Monte Carlo errors with less errors", U. Wolff, Comput. Phys.
        !   Commun. 156, 143 (2004) and arXiv:hep-lat/0306017.
        !   [2] "Strategies for improving the efficiency of quantum Monte Carlo
        !   calculations", R. M. Lee, G. J. Conduit, N. Nemec, P. Lopez Rios, and N.
        !   D.  Drummond, Phys. Rev. E. 83, 066706 (2011).

        ! In/Out:
        !   bl: Information needed to peform blocking on the fly. Suggested
        !       optimal mean and standard deviation is returned.

        use qmc_data, only: blocking_t
        use const, only: p

        integer :: i, j, B
        type(blocking_t), intent(inout) :: bl

        ! Smallest block size satisfying the condition is found.
        do i = dt_numerator, dt_shift
            do j = 0, (bl%lg_max)
                B = 2**(j)
                if (B > (2*bl%reblock_data_2(i,j)%n_blocks*real(B)*((bl%block_std(i,j) / bl%block_std(i,0))**4))**(1.0/3.0)) then
                    bl%optimal_size(i) = j
                    exit
                end if
            end do
            if (bl%optimal_size(i) == 1) then
                bl%optimal_mean(i) = 0.0_p
                bl%optimal_std(i) = 0.0_p
            else
                bl%optimal_std(i) = bl%block_std(i, bl%optimal_size(i))
                if (abs(bl%optimal_std(i)) < depsilon) then
                    bl%optimal_mean(i) = 0.0_p
                else
                    bl%optimal_mean(i) = bl%block_mean(i, bl%optimal_size(i))
                end if

            end if

            if (abs(bl%optimal_std(i)) < depsilon) then
                bl%optimal_err(i) = 0.0_p
            ! calculated assuming the normal distribution following central
            ! limit theorem.
            else
                bl%optimal_err(i) = bl%optimal_std(i)/ sqrt(real(2*(bl%reblock_data_2(i, bl%optimal_size(i))%n_blocks - 1)))
            end if
        end do
        ! Larger optimal block size between the two is used.
        if (bl%optimal_size(dt_numerator) > bl%optimal_size(dt_denominator)) then
            bl%opt_bl_size = bl%optimal_size(dt_numerator)

        else
            bl%opt_bl_size = bl%optimal_size(dt_denominator)
        end if

        if (abs(bl%optimal_mean(dt_numerator)) < depsilon .or. abs(bl%optimal_mean(dt_denominator)) < depsilon) then
            bl%optimal_mean(dt_proj_energy) = 0.0_p
            bl%optimal_std(dt_proj_energy) = 0.0_p
            bl%optimal_err(dt_proj_energy) = 0.0_p
        else
            bl%optimal_mean(dt_proj_energy) = bl%block_mean(dt_numerator, bl%opt_bl_size) / &
                                                                    bl%block_mean(dt_denominator,bl%opt_bl_size)

            bl%optimal_std(dt_proj_energy) = fraction_error(bl%block_mean(dt_numerator,bl%opt_bl_size), &
                                                        bl%block_mean(dt_denominator,bl%opt_bl_size), &
                                                        bl%reblock_data_2(dt_numerator,bl%opt_bl_size)%n_blocks,&
                                                        bl%block_std(dt_numerator,bl%opt_bl_size), &
                                                        bl%block_std(dt_denominator, bl%opt_bl_size), bl%block_cov(bl%opt_bl_size))

            bl%optimal_err(dt_proj_energy) = sqrt((bl%optimal_err(dt_numerator)/bl%optimal_std(dt_numerator)) ** 2 + &
                                                 (bl%optimal_err(dt_denominator)/bl%optimal_std(dt_denominator)) ** 2) * &
                                                                                                 bl%optimal_std(dt_proj_energy)

        end if

    end subroutine find_optimal_block

    subroutine copy_block(bl)

        ! Copying the reblock_data and data_product at each potential start
        ! points to be able to change the point at which the reblocking analysis
        ! starts from can be changed.

        ! In/out:
        !   bl: Information needed to peform blocking on the fly. reblock_save
        !       product_save is returned with the reblock_data and product_data
        !       copied.

        use qmc_data, only: blocking_t

        type(blocking_t), intent(inout) :: bl
        if ((mod(bl%n_reports_blocked,(bl%save_fq * bl%n_saved)) == 0) .and. (bl%n_saved_startpoints >= bl%n_saved)) then
            bl%reblock_save(:,:,bl%n_saved) = bl%reblock_data(:,:)
            bl%product_save(:,bl%n_saved) = bl%data_product(:)
            bl%n_saved = bl%n_saved + 1
        end if

    end subroutine copy_block

    subroutine change_start(bl, ireport, restart)

        ! Return the reblock_data_2 and data_product_2 from different start
        ! position. If the wanted start point falls within the block, the entire
        ! block is removed. Not all block sizes will have the same start
        ! position when it is modified.

        ! In:
        !   ireport: Number of reports.
        !   restart: restart * save_fq is the modified point from which the
        !       reblock analysis is carried out from
        ! In/Out:
        !   bl: Information needed to peform blocking on the fly. reblock_save_2
        !       and data_product_2 are updated to contain information for
        !       reblock analysis starting from restart point.

        use qmc_data, only: blocking_t
        use const, only: p

        type(blocking_t), intent(inout) :: bl
        integer :: i, k
        integer, intent(in) :: ireport, restart
        logical :: switch

        bl%reblock_data_2 = bl%reblock_data
        bl%data_product_2 = bl%data_product

        if (bl%n_reports_blocked > bl%save_fq * restart) then

            do i = 0, bl%lg_max
                if (bl%reblock_data_2(dt_numerator, i)%n_blocks == 1) then
                    bl%reblock_data_2(:,i)%n_blocks = 0
                    bl%reblock_data_2(:,i)%data_accumulator = 0.0_p
                    bl%reblock_data_2(:,i)%sum_of_blocks = 0.0_p
                    bl%reblock_data_2(:,i)%sum_of_block_squares = 0.0_p
                    bl%data_product_2(i) = 0.0_p
                end if
            end do

            do i = 0, bl%lg_max
                switch = .true.
                do k = 0, bl%n_saved_startpoints
                    if (abs(bl%reblock_save(dt_numerator,i,k)%data_accumulator) < depsilon .and. &
                            bl%reblock_save(dt_numerator,0,k)%n_blocks >= restart*bl%save_fq) then
                        bl%reblock_data_2(:,i)%n_blocks = bl%reblock_data_2(:,i)%n_blocks - bl%reblock_save(:,i,k)%n_blocks
                        bl%reblock_data_2(:,i)%data_accumulator = bl%reblock_data_2(:,i)%data_accumulator - &
                                                                                     bl%reblock_save(:,i,k)%data_accumulator
                        bl%reblock_data_2(:,i)%sum_of_blocks = bl%reblock_data_2(:,i)%sum_of_blocks - &
                                                                                     bl%reblock_save(:,i,k)%sum_of_blocks
                        bl%reblock_data_2(:,i)%sum_of_block_squares = bl%reblock_data_2(:,i)%sum_of_block_squares - &
                                                                                    bl%reblock_save(:,i,k)%sum_of_block_squares

                        bl%data_product_2(i) = bl%data_product_2(i) - bl%product_save(i,k)

                        switch = .false.

                        exit
                    end if
                end do

                if (switch .eqv. .true.) then
                    bl%reblock_data_2(:,i)%n_blocks = 0
                    bl%reblock_data_2(:,i)%data_accumulator = 0
                    bl%reblock_data_2(:,i)%sum_of_blocks = 0
                    bl%reblock_data_2(:,i)%sum_of_block_squares = 0
                    bl%data_product_2(i) = 0.0_p
                end if
            end do
        end if

    end subroutine change_start

    subroutine err_comparison(bl, ireport)

        ! Compares the fractional error weighted by 1/sqrt(number of data) at
        ! each of the different possible start points and returns the optimal
        ! start position.

        ! In:
        !   ireport: Number of reports.
        ! In/Out:
        !   bl: Information needed to peform blocking on the fly. start_point
        !       that is the most optimal is returned.

        use qmc_data, only: blocking_t
        use const, only: p

        type(blocking_t), intent(inout) :: bl
        integer, intent(in) :: ireport
        integer :: i, j
        integer ::  minimum(3)=0

        do i = 0, bl%n_saved_startpoints

            if (bl%n_reports_blocked > bl%save_fq * i) then
                call change_start(bl, ireport, i)
                call mean_std_cov(bl)
                call find_optimal_block(bl)

                do j = 1, 3
                    if (abs(bl%optimal_std(j)) < depsilon) then
                        bl%err_compare(j,i) = 0.0_p
                    else
                        if (bl%optimal_size(j) - 1 < log(real(bl%save_fq))/ log(2.0)) then
                            bl%err_compare(j,i) = bl%optimal_err(j)/ (bl%optimal_std(j) * &
                                                                    sqrt(real((int(real(bl%n_reports_blocked &
                                                                    - i * bl%save_fq)/bl%optimal_size(j))) * bl%optimal_size(j))))
                        else
                            bl%err_compare(j,i) = bl%optimal_err(j)/(bl%optimal_std(j) * &
                                                                    sqrt(real((int(real(bl%n_reports_blocked - &
                                                                    i * bl%save_fq)/bl%optimal_size(j))-1) * bl%optimal_size(j))))
                        end if
                    end if
                end do
            else
                bl%err_compare(:,i) = 0.0_p
            end if

        end do

        do i = 1, 3
            minimum(i) = minloc(bl%err_compare(i,:), dim = 1, mask = (bl%err_compare(i,:)>0.0_p)) - 1
        end do

        bl%start_point = maxval(minimum)

        if (bl%start_point == -1) bl%start_point = 0

    end subroutine err_comparison

    subroutine do_blocking(bl, qs, qmc_in, ireport, iter, iunit, blocking_in)

        ! Carries out blocking on the fly and changes the start point to
        ! minimise the fractional error in standard deviation weighted by the
        ! number of data points.

        ! In:
        !   qmc_in: input options relating to QMC methods.
        !   ireport: Number of reports.
        !   iter: Number of iterations.
        !   iunit: Number identifying the file to which the blocking report is
        !   written out to.
        !   blocking_in: input options for blocking on the fly.
        ! In/Out:
        !   bl: Information needed to peform blocking on the fly.
        !   qs: qmc_state where the data for current iteration is taken.

        use qmc_data, only: blocking_t, qmc_state_t, qmc_in_t, blocking_in_t

        type(blocking_t), intent(inout) :: bl
        type(qmc_state_t), intent(inout) :: qs
        type(qmc_in_t), intent(in) :: qmc_in
        integer, intent(in) :: ireport, iter
        integer, intent(in) :: iunit
        type(blocking_in_t), intent(in) :: blocking_in
        ! Array containing info on distributions of energy estimators (proj
        ! energy components and shift).
        ! energy_estimate_dist(:,1) contains means, energy_estimate_dist(:,2)
        !   contains standard deviations.
        ! energy_estimate_dist(i,:) gives information on parameter
        ! corresponding to enum within qmc_data.
        real(p) :: energy_estimate_dist(dt_numerator:dt_shift,2)

        if (bl%start_ireport == -1 .and. blocking_in%start_point<0 .and. qs%vary_shift(1)) then
            bl%start_ireport = ireport
        else if (blocking_in%start_point>=0) then
            bl%start_ireport = nint(real(blocking_in%start_point)/qmc_in%ncycles)
        end if

        ! Once the shift is varied the data needed for reblocking is
        ! collected.

        if (ireport >= bl%start_ireport .and. bl%start_ireport>=0) then
            bl%n_reports_blocked = ireport - bl%start_ireport + 1
            call collect_data(qmc_in, qs, bl, ireport)
            call copy_block(bl)
        end if

        ! For every 50 reports, the optimal mean and standard deviation
        ! and the optimal error in error is calculated and printed.
        if (mod(ireport,50) ==0 .and. ireport >= bl%start_ireport .and. &
                    bl%start_ireport>=0) then
            call change_start(bl, ireport, bl%start_point)
            call mean_std_cov(bl, energy_estimate_dist)
            call find_optimal_block(bl)
            call write_blocking(bl, qmc_in, ireport, iter, iunit)
            call check_error(bl, qs, blocking_in, ireport)
            if (blocking_in%auto_shift_damping) then
                call update_shift_damping(qs, bl, energy_estimate_dist)
            end if
        end if

        if (mod(bl%n_reports_blocked,bl%save_fq) == 0 .and. bl%n_reports_blocked > 0) then
            call err_comparison(bl, ireport)
        end if

    end subroutine do_blocking

    subroutine write_blocking(bl, qmc_in, ireport, iter, iunit)

        ! Writes the blocking report the a separate output file.

        ! In:
        !   bl: Information needed to peform blocking on the fly.
        !   qmc_in: input options relating to QMC methods.
        !   ireport: Number of reports.
        !   iter: Number of iterations.
        !   iunit: Number identifying the file to which the blocking report is
        !       written out to.

        use qmc_data, only: blocking_t, qmc_in_t

        integer, intent(in) :: iunit
        type(blocking_t), intent(in) :: bl
        type(qmc_in_t), intent(in) :: qmc_in
        integer, intent(in) :: ireport, iter
        integer :: i

        write(iunit, '(1X, I8)', advance = 'no') iter
        ! Prints the point from which reblock analysis is being
        ! carried out in terms of iterations.
        write(iunit, '(1X, I8)', advance = 'no') &
             ((bl%start_point*bl%save_fq+bl%start_ireport)*qmc_in%ncycles)
        ! Prints the mean, standard deviation and the error in error for \sum
        ! H_0j N_j and N_0 and mean and standard deviation of projected
        ! energy. Returns 0 if there are insufficient data.
        write(iunit, '(4X, ES16.7)', advance = 'no') (bl%optimal_mean(dt_numerator))
        write(iunit, '(1X, ES16.7)', advance = 'no') (bl%optimal_std(dt_numerator))
        write(iunit, '(1X, ES16.7)', advance = 'no') (bl%optimal_err(dt_numerator))

        do i = 2,3
            write(iunit, '(1X, ES16.7)', advance = 'no') (bl%optimal_mean(i))
            write(iunit, '(1X, ES16.7)', advance = 'no') (bl%optimal_std(i))
            write(iunit, '(1X, ES16.7)', advance = 'no') (bl%optimal_err(i))
        end do

        write(iunit, '(1X, ES16.7)', advance = 'no') (bl%optimal_mean(dt_proj_energy))
        write(iunit, '(1X, ES16.7)', advance = 'no') (bl%optimal_std(dt_proj_energy))
        write(iunit, '(1X, ES16.7)') (bl%optimal_err(dt_proj_energy))

        flush(iunit)

    end subroutine write_blocking

    subroutine check_error(bl, qs, blocking_in, ireport)

        ! The maximum error estimate and the inverse of estimated fractional
        ! error in the projected energy is compared to the limit specified by the
        ! user. If the condition is satisfied, reblock_done = true is returned
        ! which modifies soft_exit to also be true.

        ! In:
        !   bl: Information needed to peform blocking on the fly.
        !   blocking_in: input options for blocking on the fly.
        !   ireport: Number of reports.
        ! In/Out:
        !   qs: qmc_state where the data for current iteration is taken.

        use qmc_data, only: blocking_t, blocking_in_t, qmc_state_t
        use const

        type(qmc_state_t), intent(inout) :: qs
        type(blocking_t), intent(in) :: bl
        type(blocking_in_t), intent(in) :: blocking_in
        integer, intent(in) :: ireport
        real(p) :: error_sum = 50, number_of_blocks = 0

        if (bl%optimal_std(dt_proj_energy) > 0.0_p) then

            error_sum = bl%optimal_err(dt_proj_energy) + bl%optimal_std(dt_proj_energy)

            number_of_blocks =int((ireport-(bl%start_point*bl%save_fq+bl%start_ireport))/(2.0**bl%opt_bl_size))

        end if

        if ((error_sum < blocking_in%error_limit .and. number_of_blocks > blocking_in%min_blocks_used) &
                                    .or. number_of_blocks > blocking_in%blocks_used) then
            qs%reblock_done = .true.
        end if

    end subroutine check_error

    subroutine write_blocking_report(bl, qs)

        ! Once the calculation is finished after satisfying the user-specified
        ! target error and fractional error in projected energy, the total
        ! energy, correlation energy, errer in correlation energy and reference
        ! energy is printed out at the bottom of the HANDE output file.

        ! In:
        !   bl: Information needed to peform blocking on the fly.
        !   qs: qmc_state where the data for current iteration is taken.

        use qmc_data, only: blocking_t, qmc_state_t
        use const

        type(blocking_t), intent(in) :: bl
        type(qmc_state_t), intent(in) :: qs

        write (6, '(1X,a12,/,1X,12("-"),/)') 'Total Energy'
        write (6, '(1X, "Correlation energy:",15X,es13.6)') bl%optimal_mean(dt_proj_energy)
        write (6, '(1X, "Reference energy:",17X,es13.6)') qs%ref%H00
        write (6, '(1X, "Total energy:",21X,es13.6)') bl%optimal_mean(dt_proj_energy) + qs%ref%H00
        write (6, '(1X, "Error in correlation energy:",6X,es13.6)') bl%optimal_std(dt_proj_energy)
        write (6, '()')

    end subroutine write_blocking_report

    subroutine reset_blocking_info(bl)

        ! Subroutine to clear all stored blocking information to it's original state.
        ! This is of use when using an initial short period of calculation to optimise
        ! various parameters, before using optimal parameters for the rest of the
        ! calculation.

        ! In/Out:
        !   bl: Information needed to peform blocking on the fly.

        use qmc_data, only: blocking_t

        type(blocking_t), intent(inout) :: bl

        bl%n_reports_blocked = 0
        bl%start_ireport = -1
        bl%start_point = 0
        bl%n_saved = 1
        bl%data_product = 0
        bl%data_product_2 = 0
        bl%block_mean = 0
        bl%block_std = 0
        bl%block_cov = 0
        bl%product_save = 0
        bl%err_compare = 0

    end subroutine reset_blocking_info

    subroutine update_shift_damping(qs, bl, energy_estimate_dist)

        use qmc_data, only: qmc_state_t, blocking_t

        ! Array containing info on distributions of energy estimators (proj
        ! energy components and shift).
        ! energy_estimate_dist(:,1) contains means, energy_estimate_dist(:,2)
        !   contains standard deviations.
        ! energy_estimate_dist(i,:) gives information on parameter
        ! corresponding to enum within qmc_data.

        use parallel

        real(p), intent(in) :: energy_estimate_dist(dt_numerator:dt_shift,2)
        type(qmc_state_t), intent(inout) :: qs
        type(blocking_t), intent(inout) :: bl
        real(p) :: sd_proj_energy_dist
#ifdef PARALLEL
        integer :: ierr
#endif
        select case (bl%shift_damping_status)
        case(0)
            ! Need to wait one iteration to ensure enough information to properly calculate update.
            bl%shift_damping_status = 2
#ifdef PARALLEL
            call mpi_bcast(bl%shift_damping_status, 1, mpi_integer, root, MPI_COMM_WORLD, ierr)
#endif
        case(1)
            ! Need to use currently collected information to update shift damping.
            sd_proj_energy_dist = fraction_error(energy_estimate_dist(dt_numerator,1), &
                    energy_estimate_dist(dt_denominator,1), 1, energy_estimate_dist(dt_numerator,2), &
                    energy_estimate_dist(dt_denominator,2), bl%block_cov(0))

            ! Use linear relationship between sqrt(<S^2>-<S>^2) and shift damping
            ! parameter to set shift damping such that the standard deviation of
            ! the shift is approximately 1.5 times that of the projected energy.
            qs%shift_damping = qs%shift_damping * 1.5_p * sd_proj_energy_dist / energy_estimate_dist(dt_shift,2)
            write (6, '("# Shift damping changed to",1X,es17.10)') qs%shift_damping
            qs%shift= qs%estimators(1)%proj_energy/qs%estimators(1)%D0_population
            bl%shift_damping_status = 2
            call reset_blocking_info(bl)
#ifdef PARALLEL
            call mpi_bcast(bl%shift_damping_status, 1, mpi_integer, root, MPI_COMM_WORLD, ierr)
            call mpi_bcast(qs%shift_damping, 1, mpi_preal, root, MPI_COMM_WORLD, ierr)
            call mpi_bcast(qs%shift, size(qs%shift), mpi_preal, root, MPI_COMM_WORLD, ierr)
#endif
        case(2)
            ! Need to check that standard distribution after change is within acceptable range.
            sd_proj_energy_dist = fraction_error(energy_estimate_dist(dt_numerator,1), &
                    energy_estimate_dist(dt_denominator,1), 1, energy_estimate_dist(dt_numerator,2), &
                    energy_estimate_dist(dt_denominator,2), bl%block_cov(0))

            if (sd_proj_energy_dist > 2 * energy_estimate_dist(dt_shift,2) .or. &
                        sd_proj_energy_dist < energy_estimate_dist(dt_shift,2)) then
                ! Outside acceptable range- retry optimisation.
                bl%shift_damping_status = 1
                write (6, '("# Reattemping shift damping optimisation.")')
            else
                bl%shift_damping_status = 3
                write (6, '("# Shift damping optimised; variance within acceptable range.")')
            end if
#ifdef PARALLEL
            call mpi_bcast(bl%shift_damping_status, 1, mpi_integer, root, MPI_COMM_WORLD, ierr)
#endif
        case default
            ! On a settled value of shift damping; don't need to do anything.
            continue
        end select

    end subroutine update_shift_damping

    subroutine receive_shift_updates(ireport, start_ireport, shift_damping_status, shift_damping, shift)

        use parallel

        integer, intent(in) :: ireport, start_ireport
        integer, intent(inout) :: shift_damping_status
        real(p), intent(inout) :: shift_damping
        real(p), intent(inout), allocatable :: shift(:)
#ifdef PARALLEL
        integer :: ierr


        if (mod(ireport,50) ==0 .and. ireport >= start_ireport .and. &
                    start_ireport>=0) then
            select case (shift_damping_status)
            case(0)
                call mpi_bcast(shift_damping_status, 1, mpi_integer, root, MPI_COMM_WORLD, ierr)
            case(1)
                call mpi_bcast(shift_damping_status, 1, mpi_integer, root, MPI_COMM_WORLD, ierr)
                call mpi_bcast(shift_damping, 1, mpi_preal, root, MPI_COMM_WORLD, ierr)
                call mpi_bcast(shift, size(shift), mpi_preal, root, MPI_COMM_WORLD, ierr)
            case(2)
                call mpi_bcast(shift_damping_status, 1, mpi_integer, root, MPI_COMM_WORLD, ierr)
            case default
                continue
            end select
        end if
#endif
    end subroutine receive_shift_updates

end module blocking

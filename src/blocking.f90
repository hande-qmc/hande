module blocking

! [review] - JSS: cite blocking (and DDDA?) paper(s).
! Module for performing reblocking on the fly.

! References
! ----------
!   [1] "Monte Carlo errors with less errors", U. Wolff, Comput. Phys.
!   Commun. 156, 143 (2004) and arXiv:hep-lat/0306017.
!   [2] "Strategies for improving the efficiency of quantum Monte Carlo
!   calculations", R. M. Lee, G. J. Conduit, N. Nemec, P. Lopez Rios, and N.
!   D.  Drummond, Phys. Rev. E. 83, 066706 (2011).

implicit none

contains

    subroutine allocate_blocking(qmc_in, blocking_in, bl)

        ! Allocate different arrays in blocking_t data according to the input
        ! options.

        ! In:
        !   qmc_in: input options relating to QMC methods.
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

        ! 2^(lg_max) is approximately 4 times larger than the nreport.
        bl%lg_max = nint(log(real(qmc_in%nreport))/log(2.0)) + 2

        ! Save frequency is approximately nreport/256.
        ! Calculated by 2^(lg_max - 10) when start_save_frequency in blocking_in
        ! is negative (Default). If it is specified, 2^(start_save_frequency) is
        ! used.
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

        allocate(bl%reblock_data(3, 0:bl%lg_max, 4), stat=ierr)
        call check_allocate('bl%reblock_data',(bl%lg_max+1)*4*3,ierr)
        allocate(bl%data_product(0:bl%lg_max), stat=ierr)
        call check_allocate('bl%data_product',(bl%lg_max+1),ierr)
        allocate(bl%reblock_data_2(3, 0:bl%lg_max, 4), stat=ierr)
        call check_allocate('bl%reblock_data_2', (bl%lg_max+1)*4*3,ierr)
        allocate(bl%data_product_2(0:bl%lg_max), stat=ierr)
        call check_allocate('bl%data_product_2',(bl%lg_max+1),ierr)
        ! Mean, standard deviation and covariance of Proj. energy and reference
        ! population is calculated for each block size and added to ilock_mean,
        ! block_std and block_cov
        allocate(bl%block_mean(3, 0:bl%lg_max), stat=ierr)
        call check_allocate('bl%block_mean', 3*(bl%lg_max+1), ierr)
        allocate(bl%block_std(3, 0:bl%lg_max), stat=ierr)
        call check_allocate('bl%block_std', 3*(bl%lg_max+1), ierr)
        allocate(bl%block_cov(0:bl%lg_max), stat=ierr)
        call check_allocate('bl%block_cov', bl%lg_max+1, ierr)
        ! Array in which the reblock_data and data_product are saved every
        ! start_fq.
        allocate(bl%reblock_save(3, 0:bl%lg_max, 4, 0:bl%n_saved_startpoints) , stat=ierr)
        call check_allocate('bl%reblock_save',(bl%n_saved_startpoints+1)*4*3*(bl%lg_max+1), ierr)
        allocate(bl%product_save(0:bl%lg_max, 0:bl%n_saved_startpoints),stat=ierr)
        call check_allocate('bl%product_save',(bl%n_saved_startpoints+1)*(bl%lg_max+1),ierr)
        ! 1/(sqrt(number of data)) * fractional error for both \sum H_0j N_jand
        ! reference population is saved for comparison.
        allocate(bl%err_comp(3, 0:bl%n_saved_startpoints))
        call check_allocate('bl%err_comp', (bl%n_saved_startpoints+1)*3, ierr)

        bl%reblock_data = 0
        bl%data_product = 0
        bl%reblock_data_2 = 0
        bl%data_product_2 = 0
        bl%block_mean = 0
        bl%block_std = 0
        bl%block_cov = 0
        bl%reblock_save = 0
        bl%product_save = 0
        bl%err_comp = 0

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
        deallocate(bl%err_comp, stat=ierr)
        call check_deallocate('bl%err_comp', ierr)

    end subroutine deallocate_blocking

    subroutine write_blocking_report_header(iunit)

        ! Write the block analysis report to a file identified by iunit.

        ! In:
        !   iunit: identifies the file to which the header is written out to.

        integer, intent(in) :: iunit

        ! [review] - JSS: simpler is: write (iunit, '(1X,"#iterations")') etc.
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

        type(qmc_in_t), intent(in) :: qmc_in
        type(qmc_state_t), intent(in) :: qs
        integer, intent(in) :: ireport
        type(blocking_t), intent(inout) :: bl
        integer :: i, reblock_size

        ! \sum H_0j N_j, reference population and shift are added to column 2 of every
        ! block size.

        bl%reblock_data(1,:,2) = bl%reblock_data(1,:,2) + qs%estimators%proj_energy
        bl%reblock_data(2,:,2) = bl%reblock_data(2,:,2) + qs%estimators%D0_population
        bl%reblock_data(3,:,2) = bl%reblock_data(3,:,2) + qs%shift(1)

        ! Everytime enough data is collected for each block size, the sum in
        ! column 2 is divied by the block size and added to column 3 and squared
        ! and added to column 4.
        ! Column 1 is the number of blocks that is added in column 3 and 4.
        ! data_product contains the product of \sum H_0j N_j and reference
        ! population.

        do i = 0, (bl%lg_max)
            if (mod(bl%n_reports_blocked, 2 ** i) == 0) then
                reblock_size = 2 ** i

                bl%reblock_data(:,i,1) = (bl%n_reports_blocked)/reblock_size

                bl%reblock_data(:,i,3) = bl%reblock_data(:,i,3) + bl%reblock_data(:,i,2)/reblock_size

                bl%reblock_data(:,i,4) = bl%reblock_data(:,i,4) + (bl%reblock_data(:,i,2)/reblock_size) ** 2

                bl%data_product(i) = bl%data_product(i) + (bl%reblock_data(1,i,2)/reblock_size) * &
                                                          (bl%reblock_data(2,i,2)/reblock_size)
                bl%reblock_data(:,i,2) = 0
            end if
        end do
    end subroutine collect_data

    subroutine mean_std_cov(bl)

        ! Mean, standard deviation and covariance for each block size is
        ! calculated from the collected data.

        ! In/Out:
        !   bl: Information needed to peform blocking on the fly. block_mean,
        !       block_std and block_cov for each block size is calculated. 0 is
        !       returned if there aren't sufficient blocks to calculate them.

        use qmc_data, only: blocking_t

        type(blocking_t), intent(inout) :: bl
        integer :: i

        do i = 0, bl%lg_max

            if (bl%reblock_data_2(1,i,1) > 0.0) then
                bl%block_mean(:,i) = bl%reblock_data_2(:,i,3)/bl%reblock_data_2(:,i,1)

            else
                bl%block_mean(:,i) = 0
            end if

            if (bl%reblock_data_2(1,i,1) > 1.0) then
                bl%block_std(:,i) = (bl%reblock_data_2(:,i,4)/bl%reblock_data_2(:,i,1) - bl%block_mean(:,i) ** 2)

                bl%block_std(:,i) = sqrt(bl%block_std(:,i)/(bl%reblock_data_2(:,i,1) - 1))

                bl%block_cov(i) = (bl%data_product_2(i) - (bl%block_mean(1,i)*bl%block_mean(2,i) &
                                                        *  bl%reblock_data_2(1,i,1)))/(bl%reblock_data_2(1,i,1) - 1)
            else
                bl%block_std(:,i) = 0
                bl%block_cov(i) = 0
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
        real(p), intent(in) :: data_number
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

        integer :: i, j, B, size_e
        type(blocking_t), intent(inout) :: bl

        ! Smallest block size satisfying the condition is found.
        do i = 1, 3
            do j = 0, (bl%lg_max)
                B = 2**(j)
! [review] - CJCS: From what you said above shouldn't this be B**3 >?
                if (B > (2*bl%reblock_data_2(i,j,1)*real(B)*((bl%block_std(i,j) / bl%block_std(i,0))**4))**(1.0/3.0)) then
                    bl%optimal_size(i) = j
                    exit
                end if
            end do
            if (bl%optimal_size(i) == 1) then
                bl%optimal_mean(i) = 0
                bl%optimal_std(i) = 0
            else
                bl%optimal_std(i) = bl%block_std(i, bl%optimal_size(i))
                if (bl%optimal_std(i) == 0) then
                    bl%optimal_mean(i) = 0
                else
                    bl%optimal_mean(i) = bl%block_mean(i, bl%optimal_size(i))
                end if

            end if

            if (bl%optimal_std(i) == 0) then
                bl%optimal_err(i) = 0
            ! calculated assuming the normal distribution following central
            ! limit theorem.
            else
                bl%optimal_err(i) = bl%optimal_std(i)/ sqrt(2*(bl%reblock_data_2(i, bl%optimal_size(i), 1) - 1))
            end if
        end do
        ! Larger optimal block size between the two is used.
        if (bl%optimal_size(1) > bl%optimal_size(2)) then
            size_e = bl%optimal_size(1)

        else
            size_e = bl%optimal_size(2)
        end if

        if (bl%optimal_mean(1) == 0 .or. bl%optimal_mean(2) == 0) then
            bl%optimal_mean(4) = 0
            bl%optimal_std(4) = 0
            bl%optimal_err(4) = 0
        else
            bl%optimal_mean(4) = bl%block_mean(1, size_e) / bl%block_mean(2,size_e)

            bl%optimal_std(4) = fraction_error(bl%block_mean(1,size_e), bl%block_mean(2,size_e), bl%reblock_data_2(1,size_e,1), &
                                               bl%block_std(1,size_e), bl%block_std(2, size_e), bl%block_cov(size_e))
            bl%optimal_err(4) = sqrt((bl%optimal_err(1)/bl%optimal_std(1)) ** 2 + (bl%optimal_err(2)/bl%optimal_std(2)) ** 2) * &
                                                                                                            bl%optimal_std(4)

        end if

    end subroutine find_optimal_block

    subroutine copy_block(bl, ireport)

        ! Copying the reblock_data and data_product at each potential start
        ! points to be able to change the point at which the reblocking analysis
        ! starts from can be changed.

        ! In:
        !   ireport: Number of reports.
        ! In/out:
        !   bl: Information needed to peform blocking on the fly. reblock_save
        !       product_save is returned with the reblock_data and product_data
        !       copied.

        use qmc_data, only: blocking_t

        type(blocking_t), intent(inout) :: bl
        integer :: i
        integer, intent(in) :: ireport

        do i = 1, bl%n_saved_startpoints
            if (mod(bl%n_reports_blocked,(bl%save_fq * i)) == 0 .and. &
                        all(bl%reblock_save(:,:,:,i) == 0)) then
                bl%reblock_save(:,:,:,i) = bl%reblock_data(:,:,:)
                bl%product_save(:,i) = bl%data_product(:)
            end if
        end do
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

        type(blocking_t), intent(inout) :: bl
        integer :: i, k
        integer, intent(in) :: ireport, restart
        logical :: switch

        bl%reblock_data_2 = bl%reblock_data
        bl%data_product_2 = bl%data_product

        if (bl%n_reports_blocked > bl%save_fq * restart) then

            do i = 0, bl%lg_max
                if (bl%reblock_data_2(1, i, 1) == 1) then
                    bl%reblock_data_2(:,i,:) = 0
                    bl%data_product_2(i) = 0
                end if
            end do

            do i = 0, bl%lg_max
                switch = .true.
                do k = 0, bl%n_saved_startpoints
                    if (bl%reblock_save(1,i,2,k) == 0 .and. bl%reblock_save(1,0,1,k) >= restart*bl%save_fq) then
                        bl%reblock_data_2(:,i,:) = bl%reblock_data_2(:,i,:) - bl%reblock_save(:,i,:,k)

                        bl%data_product_2(i) = bl%data_product_2(i) - bl%product_save(i,k)

                        switch = .false.

                        exit
                    end if
                end do

                if (switch .eqv. .true.) then
                    bl%reblock_data_2(:,i,:) = 0
                    bl%data_product_2(i) = 0
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
                    if (bl%optimal_std(j) == 0.0) then
                        bl%err_comp(j,i) = 0.0
                    else
                        if (bl%optimal_size(j) - 1 < log(real(bl%save_fq))/ log(2.0)) then
                            bl%err_comp(j,i) = bl%optimal_err(j)/ (bl%optimal_std(j) * sqrt(real((int(real(bl%n_reports_blocked &
                                                                    - i * bl%save_fq)/bl%optimal_size(j))) * bl%optimal_size(j))))
                        else
                            ! [review] - JSS: such a long line is rather hard to parse.
                            bl%err_comp(j,i) = bl%optimal_err(j)/&
                                (bl%optimal_std(j) * sqrt(real((int(real(bl%n_reports_blocked &
                                - i * bl%save_fq)/bl%optimal_size(j))-1) &
                                * bl%optimal_size(j))))
                        end if
                    end if
                end do
            else
                bl%err_comp(:,i) = 0
            end if

        end do

        do i = 1, 3
            minimum(i) = minloc(bl%err_comp(i,:), dim = 1, mask = (bl%err_comp(i,:)>0)) - 1
        end do

        bl%start_point = maxval(minimum)

        if (bl%start_point == -1) bl%start_point = 0

    end subroutine err_comparison

    subroutine do_blocking(bl, qs, qmc_in, ireport, iter, iunit, blocking_in, soft_exit)

        ! Carries out blocking on the fly and changes the start point to
        ! minimise the fractional error in standard deviation weighted by the
        ! number of data points.

        ! In:
        !   qs: qmc_state where the data for current iteration is taken.
        !   qmc_in: input options relating to QMC methods.
        !   ireport: Number of reports.
        !   iter: Number of iterations.
        !   iunit: Number identifying the file to which the blocking report is
        !   written out to.
        !   blocking_in: input options for blocking on the fly.
        ! In/Out:
        !   bl: Information needed to peform blocking on the fly.
        !   soft_exit: if true, the calculation is terminated.

        use qmc_data, only: blocking_t, qmc_state_t, qmc_in_t, blocking_in_t

        type(blocking_t), intent(inout) :: bl
        type(qmc_state_t), intent(inout) :: qs
        type(qmc_in_t), intent(in) :: qmc_in
        integer, intent(in) :: ireport, iter
        integer, intent(in) :: iunit
        type(blocking_in_t), intent(in) :: blocking_in
        logical, intent(inout) :: soft_exit

        if (bl%start_ireport == -1 .and. blocking_in%start_point<0 .and. qs%vary_shift(1)) then
            bl%start_ireport = ireport
        else if (blocking_in%start_point>0) then
            bl%start_ireport = blocking_in%start_point/qmc_in%ncycles
        end if

        ! Once the shift is varied the data needed for reblocking is
        ! collected.

! [review] - CJCS: Can just write 'if (qs%vary_shift(1)) then'
        if (ireport >= bl%start_ireport .and. bl%start_ireport>0) then
            bl%n_reports_blocked = ireport - bl%start_ireport + 1
            call collect_data(qmc_in, qs, bl, ireport)
            call copy_block(bl, ireport)
        end if

        ! For every 50 reports, the optimal mean and standard deviation
        ! and the optimal error in error is calculated and printed.
        if (mod(ireport,50) ==0 .and. ireport >= bl%start_ireport .and. &
                    bl%start_ireport>0) then
            call change_start(bl, ireport, bl%start_point)
            call mean_std_cov(bl)
            call find_optimal_block(bl)
            call write_blocking(bl, qmc_in, ireport, iter, iunit)
            call check_error(bl, qs, blocking_in, soft_exit)
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
        use utils, only: get_free_unit

        integer, intent(in) :: iunit
        type(blocking_t), intent(in) :: bl
        type(qmc_in_t), intent(in) :: qmc_in
        integer, intent(in) :: ireport, iter
        integer :: i

        ! [review] - JSS: inconsistent indentation.
        write(iunit, '(1X, I8)', advance = 'no') iter
        ! Prints the point from which reblock analysis is being
        ! carried out in terms of iterations.
        write(iunit, '(1X, I8)', advance = 'no') &
             ((bl%start_point*bl%save_fq+bl%start_ireport)*qmc_in%ncycles)
        ! Prints the mean, standard deviation and the error in error for \sum
        ! H_0j N_j and N_0 and mean and standard deviation of projected
        ! energy. Returns 0 if there are insufficient data.
        ! [review] - JSS: do this in a loop over i=1,4.
        write(iunit, '(4X, ES16.7)', advance = 'no') (bl%optimal_mean(1))
        write(iunit, '(1X, ES16.7)', advance = 'no') (bl%optimal_std(1))
        write(iunit, '(1X, ES16.7)', advance = 'no') (bl%optimal_err(1))

        do i = 2,3
            write(iunit, '(1X, ES16.7)', advance = 'no') (bl%optimal_mean(i))
            write(iunit, '(1X, ES16.7)', advance = 'no') (bl%optimal_std(i))
            write(iunit, '(1X, ES16.7)', advance = 'no') (bl%optimal_err(i))
        end do

        write(iunit, '(1X, ES16.7)', advance = 'no') (bl%optimal_mean(4))
        write(iunit, '(1X, ES16.7)', advance = 'no') (bl%optimal_std(4))
        write(iunit, '(1X, ES16.7)') (bl%optimal_err(4))
        ! [review] - JSS: optimal_err(4)?

        flush(iunit)

        end subroutine write_blocking

        subroutine check_error(bl, qs, blocking_in, soft_exit)

            ! checks the sum ofstandard deviation and error in error of
            ! projected energy
            ! and if the value is smaller than the user specified value and if
            ! enough blocks
            ! are used to estimate the error in error, calculation is
            ! terminated.

            ! In:
            !   bl: Information needed to peform blocking on the fly.
            !   blocking_in: projected energy. (\sum H_0j Nj/N_0).
            ! In/Out:
            !   soft_exit: if true, the calculation is terminated.

            use qmc_data, only: blocking_t, blocking_in_t, qmc_state_t
            use const

            logical, intent(inout) :: soft_exit
            type(qmc_state_t), intent(inout) :: qs
            type(blocking_t), intent(in) :: bl
            type(blocking_in_t), intent(in) :: blocking_in
            real(p) :: error_sum = 50, error_ratio = 0

            if (bl%optimal_std(4) > 0) then

                error_sum = bl%optimal_err(4) + bl%optimal_std(4)

                error_ratio = bl%optimal_std(4)/bl%optimal_err(4)

            end if

            if (error_sum < blocking_in%error_limit .and. error_ratio > &
                                        blocking_in%min_ratio) then
                qs%reblock_done = .true.
            end if

        end subroutine check_error

end module blocking

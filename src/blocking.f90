module blocking

! Module for performing reblocking on the fly. 

implicit none

contains
    
    subroutine allocate_arrays(qmc_in, bl)
        
        ! Allocate different arrays in blocking_t data according to the input
        ! options.

        ! In:
        !   qmc_in: input options relating to QMC methods.
        ! In/Out:
        !   bl: Information needed to peform blocking on the fly. Maximum block
        !       size, frequency at which the data for changing start position
        !       and number of possible start positions are determined

        use qmc_data, only: qmc_in_t, blocking_t

        type(qmc_in_t), intent(in) :: qmc_in
        type(blocking_t), intent(inout) :: bl

        ! 2^(max_2n) is approximately 4 times larger than the nreport. 
        bl%max_2n = nint(log(real(qmc_in%nreport))/log(2.0)) + 2 

        ! Save frequency is approximately nreport/128. 
        ! Calculated by 2^(max_2n - 9)  
        bl%save_fq = 2 ** (bl%max_2n - 10)

        bl%save_number = int(qmc_in%nreport/(bl%save_fq))
        
        ! 1st column of reblock_data and reblock_data_2 is the number of blocks
        ! for different block sizes where the rows represent the block size from
        ! 0 to 2^(max_2n). 2nd column is the sum of numbers, once the block size
        ! is reached the sum is divided by the block size and copied to column 3
        ! and squared and copied to column 4. The product of proj. energy and
        ! referecne population is copied to data_product.
        allocate(bl%reblock_data(4, bl%max_2n, 2), bl%data_product(bl%max_2n))
        allocate(bl%reblock_data_2(4, bl%max_2n, 2), &
        bl%data_product_2(bl%max_2n))
        ! Mean, standard deviation and covariance of Proj. energy and reference
        ! population is calculated for each block size and added to block_mean,
        ! block_std and block_cov
        allocate(bl%block_mean(bl%max_2n, 2), bl%block_std(bl%max_2n, 2), &
        bl%block_cov(bl%max_2n))
        ! Array in which the reblock_data and data_product are saved every
        ! start_fq.
        allocate(bl%reblock_save(0:bl%save_number, 4, bl%max_2n, 2), &
        bl%product_save(0:bl%save_number, bl%max_2n))
        ! 1/(sqrt(number of data)) * fractional error for both proj. energy and
        ! reference population is saved for comparison.  
        allocate(bl%err_comp(0:bl%save_number, 2))

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

    end subroutine allocate_arrays


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
        integer :: i, j, reblock_size
        
        ! proj. energy and reference population are added to column 2 of every
        ! block size.

        bl%reblock_data(2,:,1) = bl%reblock_data(2,:,1) + qs%estimators%proj_energy
        bl%reblock_data(2,:,2) = bl%reblock_data(2,:,2) + qs%estimators%D0_population
        
        bl%report = ireport - bl%start_ireport + 1     

        ! Everytime enough data is collected for each block size, the sum in
        ! column 2 is divied by the block size and added to column 3 and squared
        ! and added to column 4.
        ! Column 1 is the number of blocks that is added in column 3 and 4.
        ! data_product contains the product of proj. energy and reference
        ! population.

        do i = 0, (bl%max_2n-1)
            if (mod(bl%report, 2 ** i) == 0) then
                reblock_size = 2 ** i
        
                bl%reblock_data(1,i+1,:) = (bl%report)/reblock_size
        
                bl%reblock_data(3,i+1,:) = bl%reblock_data(3,i+1,:) &
                + bl%reblock_data(2,i+1,:)/reblock_size 
        
                bl%reblock_data(4,i+1,:) = bl%reblock_data(4,i+1,:) &
                + (bl%reblock_data(2,i+1,:)/reblock_size) ** 2
        
                bl%data_product(i+1) = bl%data_product(i+1) + (bl%reblock_data(2,i+1,1)/reblock_size) &
                * (bl%reblock_data(2,i+1,2)/reblock_size)
        
                bl%reblock_data(2,i+1,:) = 0
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
        
        do i = 1, bl%max_2n

            if (bl%reblock_data_2(1,i,1) > 0.0) then 
                bl%block_mean(i,:) = bl%reblock_data_2(3,i,:)/bl%reblock_data_2(1,i,:)

            else
                bl%block_mean(i,:) = 0
            end if

            if (bl%reblock_data_2(1,i,1) > 1.0) then
                bl%block_std(i,:) = (bl%reblock_data_2(4,i,:)/bl%reblock_data_2(1,i,:) &
                - bl%block_mean(i,:) ** 2)
                
                bl%block_std(i,:) = sqrt(bl%block_std(i,:)/(bl%reblock_data_2(1,i,:) &
                - 1))

        
                 bl%block_cov(i) = (bl%data_product_2(i) &
                - bl%block_mean(i,1)*bl%block_mean(i,2) &
                *bl%reblock_data_2(1,i,1))/(bl%reblock_data_2(1,i,1) - 1)

            else 
                bl%block_std(i,:) = 0
                bl%block_cov(i) = 0
            end if
        end do

    end subroutine mean_std_cov

    function fraction_error(mean_1, mean_2, data_number, std_1, std_2, cov_in) &
    result(error_est)
       
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
        error_est = abs(mean_cur*sqrt((std_1/mean_1)**2 + (std_2/mean_2)**2 &
        - 2*cov_in/(data_number*mean_1*mean_2))) 
        
    end function fraction_error

    subroutine find_optimal_block(bl) 
        
        ! Finds the optimal mean and the standard deviation satisfying the
        ! condition B^3 > 2 * (B * (number of blocks)) * (std(B)/std(0))^4.
        ! Where B is the block size and std(B) is the standard deviation
        ! calculated at block size B.
        ! Returns 0 if none of the block sizes satisfy the condition.
        ! Optimal block size for the fraction between proj. energy and reference
        ! population is the larger optimal block size between the optimal block
        ! size of proj. energy and reference population.

        ! In/Out:
        !   bl: Information needed to peform blocking on the fly. Suggested
        !       optimal mean and standard deviation is returned.

        use qmc_data, only: blocking_t

        integer :: i, j, B, size_e
        type(blocking_t), intent(inout) :: bl
        
        ! Smallest block size satisfying the condition is found.
        do i = 1, 2
            do j = 1, (bl%max_2n)
                B = 2**(j-1)
                if (B > (2*bl%reblock_data_2(1,j,i)&
                    *real(B)*((bl%block_std(j,i)/&
                     bl%block_std(1,i))**4))**(1.0/3.0)) then
                    bl%optimal_size(i) = j
                    exit

                end if
            end do
            if (bl%optimal_size(i) == 1) then
                bl%optimal_mean(i) = 0
                bl%optimal_std(i) = 0
            else
                bl%optimal_std(i) = bl%block_std(bl%optimal_size(i),i)
                if (bl%optimal_std(i) == 0) then 
                    bl%optimal_mean(i) = 0
                else
                    bl%optimal_mean(i) = bl%block_mean(bl%optimal_size(i),i)
                end if
            
            end if
            
            if (bl%optimal_std(i) == 0) then
                bl%optimal_err(i) = 0
            ! calculated assuming the normal distribution following central
            ! limit theorem.
            else 
                bl%optimal_err(i) = bl%optimal_std(i)/ &
                sqrt(2*(bl%reblock_data_2(1, bl%optimal_size(i), i) - 1))
            end if
        end do
        ! Larger optimal block size between the two is used.
        if (bl%optimal_size(1) > bl%optimal_size(2)) then
            size_e = bl%optimal_size(1)

        else 
            size_e = bl%optimal_size(2)
        end if
        
        if (bl%optimal_mean(1) == 0 .or. bl%optimal_mean(2) == 0) then
            bl%optimal_mean(3) = 0
            bl%optimal_std(3) = 0
        else    
            bl%optimal_mean(3) = bl%block_mean(size_e, 1)/bl%block_mean(size_e,2)

            bl%optimal_std(3) = fraction_error(bl%block_mean(size_e,1), bl%block_mean(size_e, 2), & 
            bl%reblock_data_2(1, size_e,1), bl%block_std(size_e,1), bl%block_std(size_e, 2), &
            bl%block_cov(size_e))
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

        do i = 1, bl%save_number
            if (mod(bl%report,(bl%save_fq * i)) == 0 .and. &
            all(bl%reblock_save(i,:,:,:) == 0)) then
                bl%reblock_save(i,:,:,:) = bl%reblock_data(:,:,:)
                bl%product_save(i,:) = bl%data_product(:)
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

        if (bl%report > bl%save_fq * restart) then

            do i = 1, bl%max_2n
                if (bl%reblock_data_2(1, i, 1) == 1) then
                    bl%reblock_data_2(:,i,:) = 0
                    bl%data_product_2(i) = 0
                end if
            end do
        
            do i = 1, bl%max_2n
                switch = .true.
                do k = 0, bl%save_number
                    if (bl%reblock_save(k,2,i,1) == 0 .and. &
                    bl%reblock_save(k,1,1,1) >= restart*bl%save_fq) then
                        bl%reblock_data_2(:,i,:) = bl%reblock_data_2(:,i,:) - &
                        bl%reblock_save(k,:,i,:)
                    
                        bl%data_product_2(i) = bl%data_product_2(i) - bl%product_save(k,i)

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
        integer ::  minimum(2)=0 

        do i = 0, bl%save_number

            if (bl%report > bl%save_fq * i) then
                call change_start(bl, ireport, i)
                call mean_std_cov(bl)
                call find_optimal_block(bl)
            
                do j = 1, 2
                    if (bl%optimal_std(j) == 0.0) then
                        bl%err_comp(i,j) = 0.0
                    else
                        if (bl%optimal_size(j) - 1 < log(real(bl%save_fq))/ &
                        log(2.0)) then
                            bl%err_comp(i,j) = bl%optimal_err(j)/&
                            (bl%optimal_std(j) * sqrt(real((int(real(bl%report & 
                            - i * bl%save_fq)/bl%optimal_size(j))) &
                            * bl%optimal_size(j))))
                        else
                            bl%err_comp(i,j) = bl%optimal_err(j)/&
                            (bl%optimal_std(j) * sqrt(real((int(real(bl%report &
                            - i * bl%save_fq)/bl%optimal_size(j))-1) &
                            * bl%optimal_size(j))))
                        end if
                       ! bl%err_comp(i,j) = bl%optimal_err(j)
                    end if
                end do
            else 
                bl%err_comp(i,:) = 0
            end if

        end do


        do i = 1, 2
            minimum(i) = minloc(bl%err_comp(:,i), dim = 1, mask = &
            (bl%err_comp(:,i)>0)) - 1
        end do

        if (minimum(1) > minimum(2)) then
            bl%start_point = minimum(1)
        else 
            bl%start_point = minimum(2)
        end if
    
   end subroutine err_comparison 

end module blocking        
              

         


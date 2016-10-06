module blocking

implicit none

contains
    
    subroutine allocate_arrays(qmc_in, bl)

        use qmc_data, only: qmc_in_t, blocking_t

        type(qmc_in_t), intent(in) :: qmc_in
        type(blocking_t), intent(inout) :: bl

        bl%max_2n = nint(log(real(qmc_in%nreport))/log(2.0)) + 2 
        bl%save_fq = 2 ** (bl%max_2n - 9)
        bl%save_number = int(qmc_in%nreport/(bl%save_fq))

        allocate(bl%reblock_data(4, bl%max_2n, 2), bl%data_product(bl%max_2n))
        allocate(bl%reblock_data_2(4, bl%max_2n, 2), &
        bl%data_product_2(bl%max_2n))
        allocate(bl%block_mean(bl%max_2n, 2), bl%block_std(bl%max_2n, 2), &
        bl%block_cov(bl%max_2n))
        allocate(bl%reblock_save(0:bl%save_number, 4, bl%max_2n, 2), &
        bl%product_save(0:bl%save_number, bl%max_2n))
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
                            
   
        use qmc_data, only: qmc_in_t, qmc_state_t, blocking_t     
  
        type(qmc_in_t), intent(in) :: qmc_in
        type(qmc_state_t), intent(in) :: qs 
        integer, parameter ::  data_number = 2, statistic = 4
        integer, intent(in) :: ireport
        type(blocking_t), intent(inout) :: bl
        integer :: i, j, reblock_size
        

        bl%reblock_data(2,:,1) = bl%reblock_data(2,:,1) + qs%estimators%proj_energy
        bl%reblock_data(2,:,2) = bl%reblock_data(2,:,2) + qs%estimators%D0_population
        
        bl%report = ireport - bl%start_ireport + 1     

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
                            

        use qmc_data, only: blocking_t

        type(blocking_t), intent(inout) :: bl
        integer, parameter :: data_number = 2
        integer :: i
        
        do i = 1, bl%max_2n

            if (bl%reblock_data_2(1,i,1) > 0) then 
                bl%block_mean(i,:) = bl%reblock_data_2(3,i,:)/bl%reblock_data_2(1,i,:)
            end if

            if (bl%reblock_data_2(1,i,1) > 1) then
                bl%block_std(i,:) = (bl%reblock_data_2(4,i,:)/bl%reblock_data_2(1,i,:) &
                - bl%block_mean(i,:) ** 2)
                
                bl%block_std(i,:) = sqrt(bl%block_std(i,:)/(bl%reblock_data_2(1,i,:) &
                - 1))

        
                 bl%block_cov(i) = (bl%data_product_2(i) &
                - bl%block_mean(i,1)*bl%block_mean(i,2) &
                *bl%reblock_data_2(1,i,1))/(bl%reblock_data_2(1,i,1) - 1)
            end if
        end do

    end subroutine mean_std_cov

    function fraction_error(mean_1, mean_2, data_number, std_1, std_2, cov_in) &
    result(error_est)
        
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

        use qmc_data, only: blocking_t

        integer :: optimal_size(2)
        integer :: i, j, B, size_e
        type(blocking_t), intent(inout) :: bl

        do i = 1, 2
            do j = (bl%max_2n), 1, -1
                B = 2**(j-1)
                if (B**3 > &
                    2*bl%reblock_data_2(1,j,i)*B*(bl%block_std(j,i)/bl%block_std(1,i))**4) then
                    optimal_size(i) = j                
                end if
            end do
            if (optimal_size(i) == 1) then
                bl%optimal_mean(i) = 0
                bl%optimal_std(i) = 0
            else
                bl%optimal_std(i) = bl%block_std(optimal_size(i),i)
                if (bl%optimal_std(i) == 0) then 
                    bl%optimal_mean(i) = 0
                else
                    bl%optimal_mean(i) = bl%block_mean(optimal_size(i),i)
                end if
            
            end if
            
            if (bl%optimal_std(i) == 0) then
                bl%optimal_err(i) = 0
            else 
                bl%optimal_err(i) = bl%optimal_std(i)/ &
                sqrt(2*(bl%reblock_data_2(1, optimal_size(i), i) - 1))
            end if
        end do
       
        if (optimal_size(1) > optimal_size(2)) then
            size_e = optimal_size(1)

        else 
            size_e = optimal_size(2)
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

        use qmc_data, only: blocking_t

        type(blocking_t), intent(inout) :: bl
        integer, intent(in) :: ireport
        integer :: i, j
        integer ::  minimum(2) = 0

        do i = 0, bl%save_number

            if (bl%report > bl%save_fq * i) then
                call change_start(bl, ireport, i)
                call mean_std_cov(bl)
                call find_optimal_block(bl)
            
                do j = 1, 2
                    if (bl%optimal_std(j) == 0) then
                        bl%err_comp(i,j) = 0
                    else
                        bl%err_comp(i,j) = bl%optimal_err(j)/(bl%optimal_std(j) * sqrt(real(bl%report - i * bl%save_fq))) 
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
              

         


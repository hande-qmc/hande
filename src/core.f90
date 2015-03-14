program hubbard_fciqmc

    use hande_top_level
    use system, only: sys_t
    use qmc_data, only: reference_t

    implicit none

    real :: start_cpu_time
    integer :: start_wall_time

    type(sys_t) :: sys
    type(reference_t) :: reference

    call init_hande(start_cpu_time, start_wall_time)
    call init_calc(sys, reference)

    call run_calc(sys, reference)

    call end_calc(sys, reference)
    call end_hande(start_cpu_time, start_wall_time)

end program hubbard_fciqmc

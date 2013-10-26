program hubbard_fciqmc

    use hande_top_level
    use system, only: sys_global

    implicit none

    real :: start_cpu_time
    integer :: start_wall_time

    call init_calc(sys_global, start_cpu_time, start_wall_time)

    call run_calc(sys_global)

    call end_calc(sys_global, start_cpu_time, start_wall_time)

end program hubbard_fciqmc

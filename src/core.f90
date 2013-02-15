program hubbard_fciqmc

    use hande_top_level

    implicit none

    real :: start_cpu_time
    integer :: start_wall_time

    call init_calc(start_cpu_time, start_wall_time)

    call run_calc()

    call end_calc(start_cpu_time, start_wall_time)

end program hubbard_fciqmc

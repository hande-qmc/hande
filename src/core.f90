program hubbard_fciqmc

    use hande_top_level
    use system, only: sys_t

    implicit none

    real :: start_cpu_time
    integer :: start_wall_time

    type(sys_t) :: sys

    call init_calc(sys, start_cpu_time, start_wall_time)

    call run_calc(sys)

    call end_calc(sys, start_cpu_time, start_wall_time)

end program hubbard_fciqmc

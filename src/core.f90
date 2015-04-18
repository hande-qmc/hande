program hande

    use hande_top_level, only: init_hande, end_hande
    use lua_hande, only: run_lua_hande
    use system, only: sys_t
    use qmc_data, only: reference_t

    implicit none

    real :: start_cpu_time
    integer :: start_wall_time

    type(sys_t) :: sys
    type(reference_t) :: reference

    call init_hande(start_cpu_time, start_wall_time)

    call run_lua_hande()

    call end_hande(start_cpu_time, start_wall_time)

end program hande

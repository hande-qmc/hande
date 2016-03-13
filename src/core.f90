program hande

    use hande_top_level, only: init_hande, end_hande
    use lua_hande, only: run_lua_hande

    implicit none

    real :: start_cpu_time
    integer :: start_wall_time

    call init_hande(start_cpu_time, start_wall_time)

    call run_lua_hande()

    call end_hande(start_cpu_time, start_wall_time)

end program hande

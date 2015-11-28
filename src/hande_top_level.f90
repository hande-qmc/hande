module hande_top_level

! A very coarse interface to HANDE.

implicit none

contains

    subroutine init_hande(start_cpu_time, start_wall_time)

        ! Initialise HANDE (minimum global state).
        ! Print out information about the compiled executable.

        ! Out:
        !     start_cpu_time: cpu_time at the start of the calculation.
        !     start_wall_time: system_clock at the start of the calculation.

        use report, only: environment_report, comm_global_uuid, VCS_VERSION, GLOBAL_UUID
        use calc, only: init_calc_defaults
        use parallel, only: init_parallel, parallel_report, nprocs, nthreads, parent

        real, intent(out) :: start_cpu_time
        integer, intent(out) :: start_wall_time

        call init_parallel()

        call cpu_time(start_cpu_time)
        call system_clock(start_wall_time)

        if (parent) then
            write (6,'(/,a8,/)') 'HANDE'
            call environment_report()
        end if
        call comm_global_uuid()

        call init_calc_defaults(VCS_VERSION, GLOBAL_UUID)

        if ((nprocs > 1 .or. nthreads > 1) .and. parent) call parallel_report()

    end subroutine init_hande

    subroutine end_hande(start_cpu_time, start_wall_time)

        ! Terminate HANDE: print timing report and terminate MPI stack...

        ! In:
        !     start_cpu_time: cpu_time at the start of the calculation.
        !     start_wall_time: system_clock at the start of the calculation.

        use parallel, only: parent, end_parallel
        use report, only: end_report

        real, intent(in) :: start_cpu_time
        integer, intent(in) :: start_wall_time
        real :: end_cpu_time, wall_time
        integer :: end_wall_time, count_rate, count_max

        ! Calculation time.
        call cpu_time(end_cpu_time)
        call system_clock(end_wall_time, count_rate, count_max)
        if (end_wall_time < start_wall_time) then
            ! system_clock returns the time modulo count_max.
            ! Have ticked over to the next "block" (assume only one as this
            ! happens roughly once every 1 2/3 years with gfortran!)
            end_wall_time = end_wall_time + count_max
        end if
        wall_time = real(end_wall_time-start_wall_time)/count_rate
        if (parent) call end_report(wall_time, end_cpu_time-start_cpu_time)

        call end_parallel()

    end subroutine end_hande

end module hande_top_level

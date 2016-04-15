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

        use report, only: environment_report, comm_global_uuid, HANDE_VCS_VERSION, GLOBAL_UUID
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

        call init_calc_defaults(HANDE_VCS_VERSION, GLOBAL_UUID)

        if ((nprocs > 1 .or. nthreads > 1) .and. parent) call parallel_report()

    end subroutine init_hande

    subroutine end_hande(start_cpu_time, start_wall_time)

        ! Terminate HANDE: print timing report and terminate MPI stack...

        ! In:
        !     start_cpu_time: cpu_time at the start of the calculation.
        !     start_wall_time: system_clock at the start of the calculation.

        use parallel, only: parent, end_parallel
        use report, only: wrapper_end_report

        real, intent(in) :: start_cpu_time
        integer, intent(in) :: start_wall_time

        if (parent) call wrapper_end_report("Finished running on", start_wall_time, start_cpu_time, .true.)

        call end_parallel()

    end subroutine end_hande

end module hande_top_level

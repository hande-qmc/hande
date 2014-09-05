module non_blocking_comm_m

    ! Core routines for initialisation and termination for non-blocking
    ! communication QMC algorithm.

    implicit none

contains

    subroutine init_non_blocking_comm(spawn, request, send_counts, restart_list, restart)

        ! Deal with initial send of data when using non-blocking communications.

        ! In/Out:
        !    spawn: spawned array we we be sending.
        !    request: array of requests for non-blocking communications
        !    send_counts: number of elements to send from each processor
        !    restart_list: walkers spawned from final iteration of a restart
        !                  calculation. These need to be evolved and merged.
        ! In:
        !    restart: doing a restart calculation or not.

        use spawn_data, only: spawn_t, non_blocking_send
        use parallel, only: iproc, nthreads

        type(spawn_t), intent(inout) :: spawn
        integer, intent(inout) :: request(0:), send_counts(0:)
        type(spawn_t), intent(inout) :: restart_list
        logical, intent(in) :: restart

        integer :: start, nspawn

        send_counts = 0

        if (restart) then
            start = spawn%head_start(0,iproc) + nthreads
            nspawn = restart_list%head(0,0)
            send_counts(iproc) = restart_list%head(0,0)
            spawn%sdata(:,start:start+nspawn) = restart_list%sdata(:,:nspawn)
        else
            send_counts = 0
        end if

        call non_blocking_send(spawn, send_counts, request)

    end subroutine init_non_blocking_comm

    subroutine end_non_blocking_comm(sys, rng, tinitiator, ireport, spawn, request_s, request_rep, report_time, &
                                     ntot_particles, shift)

        ! Subroutine dealing with the last iteration when using non-blocking communications.

        ! Need to deal with:
        ! 1. Final send of walkers from time step T_final.
        ! 2. Final send of report loop quantites.

        ! In:
        !    sys: system being studied.
        !    tinitiator: true if the initiator approximation is being used.
        !    ireport: index of current report loop.
        ! In/Out:
        !    rng: random number generator.
        !    spawn: spawn_t object containing spawned walkers from final
        !        iteration.
        !    request_s: array of requests for completing this final send of
        !        walkers.
        !    request_rep: array of requests for completing communication of
        !        report loop quantitites.
        !    report_time: time at the start of the current report loop.  Returns
        !        the current time (ie the time for the start of the next report
        !        loop.
        !    ntot_particles: total number of particles in main walker list.
        !    shift: current shift value.


        use annihilation, only: annihilate_main_list_wrapper
        use spawn_data, only: spawn_t, non_blocking_send, receive_spawned_walkers, &
                              annihilate_wrapper_non_blocking_spawn
        use energy_evaluation, only: update_energy_estimators_recv
        use system, only: sys_t
        use qmc_common, only: write_fciqmc_report
        use parallel, only: parent
        use fciqmc_data, only: dump_restart_file, sampling_size

        use const, only: p, dp
        use dSFMT_interface, only: dSFMT_t

        type(sys_t), intent(in) :: sys
        type(dSFMT_t), intent(inout) :: rng
        logical, intent(in) :: tinitiator
        integer, intent(in) :: ireport
        type(spawn_t), intent(inout) :: spawn
        integer, intent(inout) :: request_s(:), request_rep(:)
        real, intent(inout) :: report_time
        real(dp), intent(inout) :: ntot_particles(sampling_size)
        real(p), intent(inout) :: shift

        real :: curr_time
        real(p) :: shift_save
        real(dp) :: ntot_particles_save(sampling_size)

        call cpu_time(curr_time)

        ! Need to receive walkers sent from final iteration and merge into main list.
        call receive_spawned_walkers(spawn, request_s)
        if (.not. dump_restart_file) then
            call annihilate_wrapper_non_blocking_spawn(spawn, tinitiator)
            call annihilate_main_list_wrapper(sys, rng, tinitiator, spawn)
            ! Receive final send of report loop quantities.
        end if
        ntot_particles_save = ntot_particles
        shift_save = shift
        call update_energy_estimators_recv(request_rep, ntot_particles)
        if (parent) call write_fciqmc_report(ireport, ntot_particles, curr_time-report_time, .false.)
        ! The call to update_energy_estimators updates the shift and ntot_particles.
        ! When restarting a calculation we actually need the old (before the call)
        ! values of these quantites to be written to the restart file, so reset
        ! them in this case.
        if (dump_restart_file) then
            ntot_particles = ntot_particles_save
            shift = shift_save
        end if

    end subroutine end_non_blocking_comm

end module non_blocking_comm_m

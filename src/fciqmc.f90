module fciqmc

! Module for performing optimised (hopefully!) full configuration interaction
! quantum monte carlo (FCIQMC) calculations.

use qmc_io
use proc_pointers
implicit none

contains

    subroutine do_fciqmc(sys, qmc_in, fciqmc_in, semi_stoch_in, restart_in, load_bal_in, io_unit, &
                         reference_in, logging_in, blocking_in, qs, qmc_state_restart)

        ! Run the FCIQMC or initiator-FCIQMC algorithm starting from the initial walker
        ! distribution using the timestep algorithm.

        ! See notes about the implementation of this using function pointers
        ! in fciqmc_main.

        ! In:
        !    sys: system being studied.
        !    semi_stoch_in: Input options for the semi-stochastic adaptation.
        !    fciqmc_in: input options relating to FCIQMC.
        !    restart_in: input options for HDF5 restart files.
        !    reference_in: current reference determinant.  If not set (ie
        !       components allocated) then a best guess is made based upon the
        !       desired spin/symmetry.
        !    qmc_in: input options relating to QMC methods.
        !    load_bal_in: input options for load balancing.
        !    io_unit: io unit to write all calculation output to.
        ! In/Out:
        !    qmc_state_restart (optional): if present, restart from a previous fciqmc calculation.
        !       Deallocated on exit.
        ! Out:
        !    qs: qmc_state for use if restarting the calculation

        use parallel
        use checking, only: check_allocate, check_deallocate
        use json_out

        use const, only: debug
        use errors, only: stop_all

        use bloom_handler, only: init_bloom_stats_t, bloom_mode_fixedn, bloom_stats_warning, &
                                 bloom_stats_t, accumulate_bloom_stats, write_bloom_report
        use determinants, only: alloc_det_info_t, dealloc_det_info_t, sum_fock_values_occ_list
        use determinant_data, only: det_info_t
        use excitations, only: excit_t, create_excited_det, get_excitation
        use annihilation, only: direct_annihilation, direct_annihilation_received_list, &
                                direct_annihilation_spawned_list, deterministic_annihilation
        use calc, only: doing_calc
        use death, only: stochastic_death
        use importance_sampling, only: importance_sampling_weight
        use ifciqmc
        use non_blocking_comm_m, only: init_non_blocking_comm, end_non_blocking_comm
        use spawning, only: create_spawned_particle_initiator
        use qmc, only: init_qmc
        use qmc_common
        use dSFMT_interface, only: dSFMT_t, dSFMT_init, dSFMT_end, dSFMT_state_t_to_dSFMT_t, dSFMT_t_to_dSFMT_state_t, &
                                   free_dSFMT_state_t
        use semi_stoch, only: semi_stoch_t, check_if_determ, determ_projection
        use semi_stoch, only: dealloc_semi_stoch_t, init_semi_stoch_t, init_semi_stoch_t_flags, set_determ_info
        use system, only: sys_t, sys_t_json, read_in
        use restart_hdf5, only: init_restart_info_t, restart_info_t, dump_restart_hdf5, dump_restart_file_wrapper
        use spawn_data, only: receive_spawned_walkers, annihilate_wrapper_non_blocking_spawn, write_memcheck_report
        use qmc_data, only: qmc_in_t, fciqmc_in_t, semi_stoch_in_t, restart_in_t, load_bal_in_t, empty_determ_space, &
                            qmc_state_t, annihilation_flags_t, semi_stoch_separate_annihilation, qmc_in_t_json,      &
                            fciqmc_in_t_json, semi_stoch_in_t_json, restart_in_t_json, load_bal_in_t_json, &
                            blocking_t, blocking_in_t, blocking_in_t_json
        use reference_determinant, only: reference_t, reference_t_json
        use check_input, only: check_qmc_opts, check_fciqmc_opts, check_load_bal_opts, check_blocking_opts
        use hamiltonian_data
        use energy_evaluation, only: get_sanitized_projected_energy

        use logging, only: init_logging, end_logging, prep_logging_mc_cycle, write_logging_calc_fciqmc, &
                            logging_in_t, logging_t, logging_in_t_json, logging_t_json
        use blocking, only: write_blocking_report_header, init_blocking, do_blocking, deallocate_blocking, &
                            write_blocking_report, update_shift_damping
        use report, only: write_date_time_close

        type(sys_t), intent(in) :: sys
        type(qmc_in_t), intent(in) :: qmc_in
        type(fciqmc_in_t), intent(in) :: fciqmc_in
        type(semi_stoch_in_t), intent(in) :: semi_stoch_in
        type(restart_in_t), intent(in) :: restart_in
        type(load_bal_in_t), intent(in) :: load_bal_in
        type(reference_t), intent(in) :: reference_in
        type(qmc_state_t), intent(inout), optional :: qmc_state_restart
        type(qmc_state_t), intent(out), target :: qs
        type(blocking_in_t), intent(in) :: blocking_in

        type(logging_in_t), intent(in) :: logging_in
        integer, intent(in) :: io_unit
        type(logging_t) :: logging_info

        type(det_info_t) :: cdet
        type(dSFMT_t) :: rng
        type(bloom_stats_t) :: bloom_stats
        type(semi_stoch_t) :: determ
        type(json_out_t) :: js
        type(qmc_in_t) :: qmc_in_loc

        integer :: idet, ireport, icycle, ideterm, ierr, ispace
        integer :: iter, semi_stoch_iter
        integer(int_64) :: nattempts
        real(dp), allocatable :: nparticles_old(:)

        integer(int_p) :: ndeath
        ! [todo] - Should some of these be int_p?
        integer :: nspawn_events, nattempts_current_det_ispace
        type(excit_t) :: connection
        type(hmatel_t) :: hmatel
        real(p), allocatable :: real_population(:), weighted_population(:)
        integer :: send_counts(0:nprocs-1), req_data_s(0:nprocs-1)
        type(annihilation_flags_t) :: annihilation_flags
        type(restart_info_t) :: ri, ri_shift
        character(36) :: uuid_restart

        logical :: soft_exit, write_restart_shift, error
        logical :: determ_parent, restart_proj_est

        real :: t1, t2
        logical :: update_tau, restarting, imag

        type(blocking_t) :: bl
        integer :: iunit, restart_version_restart
        integer :: date_values(8)
        character(:), allocatable :: err_msg

        if (parent) then
            write (io_unit,'(1X,"FCIQMC")')
            write (io_unit,'(1X,"------",/)')
        end if

        restarting = present(qmc_state_restart) .or. restart_in%read_restart
        if (parent) then
            ! Check input options.
            if (present(qmc_state_restart)) then
                call check_qmc_opts(qmc_in, sys, .not.present(qmc_state_restart), restarting, &
                    qmc_state_restart=qmc_state_restart, fciqmc_in=fciqmc_in)
            else
                call check_qmc_opts(qmc_in, sys, .not.present(qmc_state_restart), restarting, fciqmc_in=fciqmc_in)
            end if
            call check_fciqmc_opts(sys, fciqmc_in, blocking_in)
            call check_load_bal_opts(load_bal_in)
            call check_blocking_opts(sys, blocking_in, restart_in)
        end if

        ! Initialise data.
        call init_qmc(sys, qmc_in, restart_in, load_bal_in, reference_in, io_unit, annihilation_flags, qs, uuid_restart, &
                      restart_version_restart, fciqmc_in=fciqmc_in, qmc_state_restart=qmc_state_restart)

        if (debug) call init_logging(logging_in, logging_info, 0)

        if (parent) then
            call json_object_init(js, tag=.true., io=io_unit)
            call sys_t_json(js, sys)
            ! The default values of pattempt_* are not in qmc_in
            qmc_in_loc = qmc_in
            qmc_in_loc%pattempt_single = qs%excit_gen_data%pattempt_single
            qmc_in_loc%pattempt_double = qs%excit_gen_data%pattempt_double
            qmc_in_loc%shift_damping = qs%shift_damping
            qmc_in_loc%pattempt_parallel = qs%excit_gen_data%pattempt_parallel
            qmc_in_loc%quasi_newton_threshold = qs%propagator%quasi_newton_threshold
            qmc_in_loc%quasi_newton_value = qs%propagator%quasi_newton_value
            qmc_in_loc%quasi_newton_pop_control = qs%propagator%quasi_newton_pop_control
            call qmc_in_t_json(js, qmc_in_loc)
            call fciqmc_in_t_json(js, fciqmc_in)
            call semi_stoch_in_t_json(js, semi_stoch_in)
            call restart_in_t_json(js, restart_in, uuid_restart)
            call blocking_in_t_json(js, blocking_in)
            call load_bal_in_t_json(js, load_bal_in)
            call reference_t_json(js, qs%ref, sys)
            call logging_in_t_json(js, logging_in)
            call logging_t_json(js, logging_info, terminal=.true.)
            call json_object_end(js, terminal=.true., tag=.true.)
            write (js%io, '()')
        end if

        allocate(nparticles_old(qs%psip_list%nspaces), stat=ierr)
        call check_allocate('nparticles_old', qs%psip_list%nspaces, ierr)
        allocate(real_population(qs%psip_list%nspaces), stat=ierr)
        call check_allocate('real_population', qs%psip_list%nspaces, ierr)
        allocate(weighted_population(qs%psip_list%nspaces), stat=ierr)
        call check_allocate('weighted_population', qs%psip_list%nspaces, ierr)

        call dSFMT_init(qmc_in%seed+iproc, 50000, rng)
        if (restart_in%restart_rng .and. allocated(qs%rng_state%dsfmt_state)) then
            call dSFMT_state_t_to_dSFMT_t(rng, qs%rng_state, err_msg=err_msg)
            if (allocated(err_msg)) call stop_all('do_fciqmc', 'Failed to reset RNG state: '//err_msg)
            call free_dSFMT_state_t(qs%rng_state)
        end if

        ! Initialise bloom_stats components to the following parameters.
        call init_bloom_stats_t(bloom_stats, mode=bloom_mode_fixedn, encoding_factor=qs%psip_list%pop_real_factor)

        ! Allocate det_info_t components.
        call alloc_det_info_t(sys, cdet, .false.)

        ! Should we dump a restart file just before the shift is turned on?
        write_restart_shift = restart_in%write_restart_shift
        call init_restart_info_t(ri, write_id=restart_in%write_id)
        call init_restart_info_t(ri_shift, write_id=restart_in%write_shift_id)

        ! In case this is not set.
        nspawn_events = 0

        ! Some initial semi-stochastic parameters.
        ! Turn semi-stochastic on immediately unless asked otherwise.
        semi_stoch_iter = max(semi_stoch_in%start_iter, qs%mc_cycles_done+1)
        if (all(qs%vary_shift) .and. semi_stoch_in%shift_iter /= -1) then
            ! User wanted the shift to start shift_iter iterations after the
            ! shift was enabled. We don't know when this happened in the
            ! previous calculation so just start semi-stochastic now.
            semi_stoch_iter = qs%mc_cycles_done+1
        end if

        call init_semi_stoch_t_flags(determ, size(qs%psip_list%states, dim=2))

        ! from restart
        nparticles_old = qs%psip_list%tot_nparticles

        ! Main fciqmc loop.
        if (parent) call write_qmc_report_header(qs%psip_list%nspaces, cmplx_est=sys%read_in%comp, io_unit=io_unit)
        restart_proj_est = present(qmc_state_restart) .or. (restart_in%read_restart .and. restart_version_restart >= 2)
        if (.not.restart_proj_est) call initial_ci_projected_energy(sys, qs, fciqmc_in%non_blocking_comm, nparticles_old)
        if (fciqmc_in%non_blocking_comm) then
            call init_non_blocking_comm(qs, req_data_s, send_counts, qmc_in%ncycles, restart_in%read_restart, restart_proj_est)
        else
            call initial_qmc_status(sys, qmc_in, qs, nparticles_old, .false., io_unit)
        end if

        ! Initialise timer.
        call cpu_time(t1)
        ! Allocate arrays needed for reblock analysis
        if (blocking_in%blocking_on_the_fly) call init_blocking(qmc_in, blocking_in, bl, qs%shift_damping_status)

        if (parent .and. blocking_in%blocking_on_the_fly) then
            open(newunit=iunit, file=blocking_in%filename, status='unknown')
            call write_blocking_report_header(iunit, sys%read_in%comp)
        end if

        do ireport = 1, qmc_in%nreport

            qs%estimators%proj_energy_old = get_sanitized_projected_energy(qs)

            ! Zero report cycle quantities.
            call init_report_loop(qs, bloom_stats)

            do icycle = 1, qmc_in%ncycles

                iter = qs%mc_cycles_done + (ireport-1)*qmc_in%ncycles + icycle

                if (debug) call prep_logging_mc_cycle(iter, logging_in, logging_info, sys%read_in%comp)

                ! Should we turn semi-stochastic on now?
                if (iter == semi_stoch_iter .and. semi_stoch_in%space_type /= empty_determ_space) then
                    determ%doing_semi_stoch = .true.
                    call init_semi_stoch_t(determ, semi_stoch_in, sys, qs%propagator, qs%psip_list, qs%ref, annihilation_flags, &
                                           qs%spawn_store%spawn, qmc_in%use_mpi_barriers, io_unit)
                end if

                call init_mc_cycle(qs%psip_list, qs%spawn_store%spawn, nattempts, ndeath, &
                                            complx = sys%read_in%comp)
                call load_balancing_wrapper(sys, qs%propagator, qs%ref, load_bal_in, annihilation_flags, &
                                            fciqmc_in%non_blocking_comm, io_unit, rng, qs%psip_list, qs%spawn_store%spawn, &
                                            qs%par_info, determ)
                if (fciqmc_in%non_blocking_comm) qs%spawn_store%spawn_recv%proc_map = qs%par_info%load%proc_map
                ideterm = 0

                do idet = 1, qs%psip_list%nstates ! loop over walkers/dets

                    cdet%f => qs%psip_list%states(:,idet)
                    cdet%data => qs%psip_list%dat(:,idet)

                    call decoder_ptr(sys, cdet%f, cdet, qs%excit_gen_data)
                    
                    if (qs%propagator%quasi_newton) &
                        cdet%fock_sum = sum_fock_values_occ_list(sys, qs%propagator%sp_fock, cdet%occ_list) - qs%ref%fock_sum

                    do ispace = 1, qs%psip_list%nspaces
                        ! Extract the real sign from the encoded sign.
                        real_population(ispace) = real(qs%psip_list%pops(ispace,idet),p)/qs%psip_list%pop_real_factor
                        weighted_population(ispace) = importance_sampling_weight(qs%trial, cdet, real_population(ispace))
                    end do

                    ! If this is a deterministic state then copy its population
                    ! across to the determ%vector array. (Both replicas use the
                    ! same deterministic space.)
                    call set_determ_info(idet, real_population, ideterm, determ, determ_parent)

                    ! Is this determinant an initiator?
                    call set_parent_flag(real_population, qmc_in%initiator_pop, determ%flags(idet), &
                                         fciqmc_in%quadrature_initiator, cdet%initiator_flag)

                    do ispace = 1, qs%psip_list%nspaces

                        imag = sys%read_in%comp .and. mod(ispace,2) == 0

                        ! It is much easier to evaluate the projected energy at the
                        ! start of the i-FCIQMC cycle than at the end, as we're
                        ! already looping over the determinants.
                        connection = get_excitation(sys%nel, sys%basis, cdet%f, qs%ref%f0)
                        if (.not. sys%read_in%comp) then
                            call update_proj_energy_ptr(sys, qs%ref%f0, qs%trial%wfn_dat, cdet, [weighted_population(ispace)], &
                                                        qs%estimators(ispace), connection, hmatel)
                        else if (.not. imag) then
                            call update_proj_energy_ptr(sys, qs%ref%f0, qs%trial%wfn_dat, cdet, &
                                                        weighted_population(ispace:ispace+1), &
                                                        qs%estimators(ispace), connection, hmatel)
                        end if

                        nattempts_current_det_ispace = decide_nattempts(rng, real_population(ispace))

                        call do_fciqmc_spawning_attempt(rng, qs%spawn_store%spawn, bloom_stats, sys, qs, &
                                                        nattempts_current_det_ispace, &
                                                        cdet, determ, determ_parent, qs%psip_list%pops(ispace, idet), &
                                                        sys%read_in%comp .and. modulo(ispace,2)==0, &
                                                        ispace, logging_info)

                        ! Clone or die.
                        if (.not. determ_parent) then
                            call stochastic_death(rng, qs, cdet%fock_sum, qs%psip_list%dat(1,idet), &
                                            qs%shift(ispace), qs%estimators(ispace)%proj_energy_old, logging_info, &
                                            qs%psip_list%pops(ispace,idet), qs%psip_list%nparticles(ispace), ndeath)
                        end if
                    end do
                end do

                associate(pl=>qs%psip_list, spawn=>qs%spawn_store%spawn, spawn_recv=>qs%spawn_store%spawn_recv)
                    if (fciqmc_in%non_blocking_comm) then
                        call receive_spawned_walkers(spawn_recv, req_data_s)
                        call evolve_spawned_walkers(sys, qmc_in, qs, spawn_recv, spawn, cdet, rng, ndeath, &
                                                    fciqmc_in%quadrature_initiator, logging_info, bloom_stats)
                        call direct_annihilation_received_list(sys, rng, qs%ref, annihilation_flags, pl, spawn_recv)
                        ! Need to add walkers which have potentially moved processor to the spawned walker list.
                        if (qs%par_info%load%needed) then
                            call redistribute_particles(pl%states, pl%pop_real_factor,  pl%pops, pl%nstates,  pl%nparticles, spawn)
                            qs%par_info%load%needed = .false.
                        end if
                        call direct_annihilation_spawned_list(sys, rng, qs%ref, annihilation_flags, pl, spawn, send_counts, &
                                                              req_data_s, qs%par_info%report_comm%nb_spawn, nspawn_events)
                        call end_mc_cycle(qs%par_info%report_comm%nb_spawn(1), ndeath, pl%pop_real_factor, nattempts, &
                                          qs%spawn_store%rspawn)
                    else
                        ! If using semi-stochastic then perform the deterministic projection step.
                        if (determ%doing_semi_stoch) call determ_projection(rng, qmc_in, qs, spawn, determ)

                        call direct_annihilation(sys, rng, qs%ref, annihilation_flags, pl, spawn, nspawn_events, determ)
                        if (debug) call write_logging_calc_fciqmc(logging_info, iter, nspawn_events, ndeath, nattempts)
                        call end_mc_cycle(nspawn_events, ndeath, pl%pop_real_factor, nattempts, qs%spawn_store%rspawn)
                    end if
                end associate

            end do

            update_tau = bloom_stats%nblooms_curr > 0

            error = qs%spawn_store%spawn%error .or. qs%psip_list%error

            call end_report_loop(io_unit, qmc_in, iter, update_tau, qs, nparticles_old, &
                                 nspawn_events, semi_stoch_in%shift_iter, semi_stoch_iter, soft_exit, &
                                 load_bal_in, bloom_stats=bloom_stats, doing_lb=fciqmc_in%doing_load_balancing, &
                                 nb_comm=fciqmc_in%non_blocking_comm, comp=sys%read_in%comp, &
                                 error=error)
            if (error) exit

            if (update_tau) call rescale_tau(qs%tau)

            call cpu_time(t2)
            if (parent) then
                if (bloom_stats%nblooms_curr > 0) call bloom_stats_warning(bloom_stats, io_unit=io_unit)
                call write_qmc_report(qmc_in, qs, ireport, nparticles_old, t2-t1, .false., &
                                        fciqmc_in%non_blocking_comm, io_unit=io_unit, cmplx_est=sys%read_in%comp)
                if (blocking_in%blocking_on_the_fly) then
                    call do_blocking(bl, qs, qmc_in, ireport, iter, iunit, blocking_in, sys%read_in%comp)
                end if

            end if
            if (blocking_in%auto_shift_damping) call update_shift_damping(qs, bl, ireport, sys%read_in%comp)

            ! Update the time for the start of the next iteration.
            t1 = t2

            call dump_restart_file_wrapper(qs, write_restart_shift, restart_in%write_freq, nparticles_old, ireport, &
                                           qmc_in%ncycles, sys%basis%nbasis, ri, ri_shift, fciqmc_in%non_blocking_comm, &
                                           sys%basis%info_string_len, rng)

            qs%psip_list%tot_nparticles = nparticles_old

            if (soft_exit) exit

            ! Should we try and update the reference determinant now?
            if (mod(ireport, fciqmc_in%select_ref_det_every_nreports) == 0) &
                    call select_ref_det(sys, fciqmc_in%ref_det_factor, qs)

        end do

        call dSFMT_t_to_dSFMT_state_t(rng, qs%rng_state)

        if (blocking_in%blocking_on_the_fly) call deallocate_blocking(bl)
        if (blocking_in%blocking_on_the_fly .and. parent) call date_and_time(VALUES=date_values)
        if (blocking_in%blocking_on_the_fly .and. parent) call write_date_time_close(iunit, date_values)

        if (fciqmc_in%non_blocking_comm) call end_non_blocking_comm(sys, rng, qmc_in, annihilation_flags, ireport, &
                                                                    qs, qs%spawn_store%spawn_recv,  req_data_s,  &
                                                                    qs%par_info%report_comm%request, t1, nparticles_old, &
                                                                    qs%shift(1), restart_in%write_restart, load_bal_in)

        if (parent) write (io_unit,'()')
        call write_bloom_report(bloom_stats, io_unit=io_unit)
        associate(pl=>qs%psip_list, spawn=>qs%spawn_store%spawn)
            if (determ%doing_semi_stoch .and. determ%projection_mode == semi_stoch_separate_annihilation) then
                call load_balancing_report(pl%nparticles, pl%nstates, qmc_in%use_mpi_barriers, spawn%mpi_time, &
                                           determ%mpi_time, io_unit=io_unit)
            else
                call load_balancing_report(pl%nparticles, pl%nstates, qmc_in%use_mpi_barriers, spawn%mpi_time, io_unit=io_unit)
            end if
        if (parent .and. blocking_in%blocking_on_the_fly .and. soft_exit) call write_blocking_report(bl, qs)
        end associate
        call write_memcheck_report(qs%spawn_store%spawn, io_unit=io_unit)

        if (soft_exit .or. error) then
            qs%mc_cycles_done = qs%mc_cycles_done + qmc_in%ncycles*ireport
        else
            qs%mc_cycles_done = qs%mc_cycles_done + qmc_in%ncycles*qmc_in%nreport
        end if

        if (qs%restart_in%write_restart) then
            call dump_restart_hdf5(qs, qs%mc_cycles_done, nparticles_old, sys%basis%nbasis, fciqmc_in%non_blocking_comm, &
                                    sys%basis%info_string_len)
            if (parent) write (io_unit,'()')
        end if

        if (determ%doing_semi_stoch) call dealloc_semi_stoch_t(determ, .false.)
        if (debug) call end_logging(logging_info)

        call dealloc_det_info_t(cdet, .false.)
        
        call dSFMT_end(rng)

    end subroutine do_fciqmc

    subroutine evolve_spawned_walkers(sys, qmc_in, qs, spawn_recv, spawn_to_send, cdet, rng, ndeath, &
                                    quadrature_initiator, logging_info, bloom_stats)

        ! Evolve spawned list of walkers one time step.
        ! Used for non-blocking communications.

        ! In:
        !   sys: system being studied.
        !   qmc_in: input options relating to QMC methods.
        !   quadrature_initiator: how to apply initiator approximation in complex systems.
        !   logging_info: information on level of logging to use within calculation.
        ! In/Out:
        !   qs: qmc_state_t containing information about the reference det and estimators.
        !   spawn: spawn_t object containing walkers spawned onto this processor during previous time step.
        !   cdet: type containing information about determinant. (easier to take this in as it is allocated
        !        / deallocated in do_fciqmc).
        !   rng: random number generator.
        !   ndeath: running total of number of particles which have died or been cloned.
        !   bloom_stats: information on blooming within calculation.

        use proc_pointers, only: sc0_ptr
        use death, only: stochastic_death
        use determinants, only: sum_fock_values_occ_list
        use determinant_data, only: det_info_t
        use dSFMT_interface, only: dSFMT_t
        use excitations, only: excit_t, get_excitation
        use ifciqmc
        use qmc_data, only: qmc_in_t, qmc_state_t
        use logging, only: logging_t
        use spawn_data, only: spawn_t
        use system, only: sys_t
        use qmc_common, only: decide_nattempts
        use hamiltonian_data
        use bloom_handler, only: bloom_stats_t
        use semi_stoch, only: semi_stoch_t

        type(sys_t), intent(in) :: sys
        type(qmc_in_t), intent(in) :: qmc_in
        type(qmc_state_t), intent(inout) :: qs
        type(spawn_t), intent(inout) :: spawn_recv, spawn_to_send
        type(dSFMT_t), intent(inout) :: rng
        type(det_info_t), intent(inout) :: cdet
        integer(int_p), intent(inout) :: ndeath
        logical, intent(in) :: quadrature_initiator

        type(excit_t) :: connection
        type(hmatel_t) :: hmatel
        integer :: idet, ispace
        integer :: nattempts_current_det
        integer(int_p) :: int_pop(spawn_recv%ntypes)
        real(p) :: real_pop(spawn_recv%ntypes)
        real(dp) :: list_pop
        logical :: imag

        type(logging_t), intent(in) :: logging_info
        type(bloom_stats_t), intent(inout) :: bloom_stats
        type(semi_stoch_t) :: determ

        allocate(cdet%f(sys%basis%tensor_label_len))
        allocate(cdet%data(1))

        do idet = 1, spawn_recv%head(0,0) ! loop over walkers/dets

            int_pop = int(spawn_recv%sdata(spawn_recv%bit_str_len+1:spawn_recv%bit_str_len+spawn_recv%ntypes, idet), int_p)

            real_pop = real(int_pop,p) / qs%psip_list%pop_real_factor
            cdet%f = int(spawn_recv%sdata(:sys%basis%tensor_label_len,idet),i0)
            ! Need to generate spawned walker data to perform evolution.
            cdet%data(1) = sc0_ptr(sys, cdet%f) - qs%ref%H00

            call decoder_ptr(sys, cdet%f, cdet, qs%excit_gen_data)
            if (qs%propagator%quasi_newton) cdet%fock_sum = sum_fock_values_occ_list(sys, qs%propagator%sp_fock, cdet%occ_list) &
                - qs%ref%fock_sum

            ! Is this determinant an initiator?
            ! [todo] - pass determ_flag rather than 1.
            call set_parent_flag(real_pop, qmc_in%initiator_pop, 1, quadrature_initiator, cdet%initiator_flag)


            do ispace = 1, qs%psip_list%nspaces

                nattempts_current_det = decide_nattempts(rng, real_pop(ispace))
                imag = sys%read_in%comp .and. mod(ispace,2) == 0

                ! It is much easier to evaluate the projected energy at the
                ! start of the i-FCIQMC cycle than at the end, as we're
                ! already looping over the determinants.
                connection = get_excitation(sys%nel, sys%basis, cdet%f, qs%ref%f0)
                call update_proj_energy_ptr(sys, qs%ref%f0, qs%trial%wfn_dat, cdet, real_pop, qs%estimators(ispace), &
                                            connection, hmatel)

                call do_fciqmc_spawning_attempt(rng, spawn_to_send, bloom_stats, sys, qs, nattempts_current_det, &
                                            cdet, determ, .false., int_pop(ispace), &
                                            sys%read_in%comp .and. modulo(ispace,2) == 0, &
                                            ispace, logging_info)

                ! Clone or die.
                ! list_pop is meaningless as particle_t%nparticles is updated upon annihilation.
                call stochastic_death(rng, qs, cdet%fock_sum, cdet%data(1), qs%shift(ispace), &
                                      qs%estimators(ispace)%proj_energy, logging_info, int_pop(ispace), list_pop, ndeath)

                ! Update population of walkers on current determinant.
                spawn_recv%sdata(spawn_recv%bit_str_len+1:spawn_recv%bit_str_len+spawn_recv%ntypes, idet) = int_pop
            end do

        end do

        deallocate(cdet%f, cdet%data)

    end subroutine evolve_spawned_walkers

    subroutine do_fciqmc_spawning_attempt(rng, spawn, bloom_stats, sys, qs, nattempts_current_det, &
                                          cdet, determ, determ_parent, pop, imag_parent, ispace, &
                                          logging_info)

        ! Perform spawning from a given determinant in a given space.

        ! In:
        !   sys: information on system under consideration.
        !   logging_info: information on current logging
        !       settings.
        !   nattempts_current_det: total number of spawning attempts
        !       to make on this determinant.
        !   ispace: space currently under consideration.
        !   determ_parent: true if parent determinant is within the
        !       semistochastic space, otherwise false.
        !   imag_parent: true if spawning from psips within an imaginary
        !       space.
        !   pop: population of given determinant in given space.
        !   determ: derived type containing information on semistochastic
        !       space within propogation.
        ! In/Out:
        !   rng: random number generator.
        !   bloom_stats: information on blooms during calculation.
        !   qs: qmc_state_t derived type with information on
        !       current calculation.
        !   spawn: stored information on spawning.
        !   cdet: determinant spawning is originating from.

        use dSFMT_interface, only: dSFMT_t
        use system, only: sys_t
        use qmc_data, only: qmc_state_t
        use logging, only: logging_t
        use determinant_data, only: det_info_t
        use bloom_handler, only: bloom_stats_t, accumulate_bloom_stats
        use semi_stoch, only: semi_stoch_t, check_if_determ

        use excitations, only: excit_t, create_excited_det
        use death, only: stochastic_death
        use spawn_data, only: spawn_t

        type(sys_t), intent(in) :: sys
        type(qmc_state_t), intent(inout) :: qs
        type(logging_t), intent(in) :: logging_info
        integer, intent(in) :: nattempts_current_det, ispace
        type(det_info_t), intent(inout) :: cdet
        logical, intent(in) :: determ_parent, imag_parent
        integer(int_p), intent(in) :: pop

        type(dSFMT_t), intent(inout) :: rng
        type(bloom_stats_t), intent(inout) :: bloom_stats
        type(spawn_t), intent(inout) :: spawn
        integer(int_p) :: nspawned, nspawned_im
        type(semi_stoch_t), intent(in) :: determ

        type(excit_t) :: connection
        integer :: iparticle, space_real, space_imag
        integer(i0) :: f_child(sys%basis%tot_string_len)
        logical :: determ_child
        integer(int_p) :: scratch

        ! First, determine the particle types possibly created by spawning.
        ! If we have a more sophisticated approach to multiple spaces this will
        ! need to be changed.
        ! NB this implicitly assumes the spaces are ordered (real, imaginary)

        if (imag_parent) then
            space_imag = ispace
            space_real = ispace - 1
        else
            space_real = ispace
            space_imag = ispace + 1
        end if

        ! Attempt spawning and death for a single determinant in a single
        ! space.
        do iparticle = 1, nattempts_current_det

            ! Attempt to spawn.
            call spawner_ptr(rng, sys, qs, qs%spawn_store%spawn%cutoff, qs%psip_list%pop_real_factor, &
                            cdet, pop, gen_excit_ptr, qs%trial%wfn_dat, &
                            logging_info, nspawned, nspawned_im, connection)
            if (imag_parent) then
                ! If imaginary parent have to factor into resulting signs/reality.
                scratch = nspawned_im
                nspawned_im = nspawned
                nspawned = -scratch
            end if

            ! Spawn if attempt was successful.
            if (nspawned /= 0_int_p) then
                if (determ_parent) then
                    call create_excited_det(sys%basis, cdet%f, connection, f_child)
                    determ_child = check_if_determ(determ%hash_table, determ%dets, f_child)
                    ! If the spawning is both from and to the deterministic space, cancel it.
                    if (.not. determ_child) then
                        call create_spawned_particle_ptr(sys%basis, qs%ref, cdet, connection, nspawned, &
                                                         space_real, spawn, f_child)
                    else
                        nspawned = 0_int_p
                    end if
                else
                    call create_spawned_particle_ptr(sys%basis, qs%ref, cdet, connection, nspawned, space_real, &
                                                     spawn)
                end if
                call accumulate_bloom_stats(bloom_stats, nspawned)
            end if
            if (nspawned_im /= 0_int_p) then
                call create_spawned_particle_ptr(sys%basis, qs%ref, cdet, connection, nspawned_im, space_imag, &
                                                     spawn)
                call accumulate_bloom_stats(bloom_stats, nspawned_im)
            end if
        end do

    end subroutine do_fciqmc_spawning_attempt

end module fciqmc

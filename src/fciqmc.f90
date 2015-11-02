module fciqmc

! Module for performing optimised (hopefully!) full configuration interaction
! quantum monte carlo (FCIQMC) calculations.

use fciqmc_data
use proc_pointers
implicit none

contains

    subroutine do_fciqmc(sys, qmc_in, fciqmc_in, semi_stoch_in, restart_in, load_bal_in, reference_in)

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
        ! In/Out:
        !    qmc_in: input options relating to QMC methods.
        !    load_bal_in: input options for load balancing.

        use parallel
        use checking, only: check_allocate
        use json_out

        use bloom_handler, only: init_bloom_stats_t, bloom_mode_fixedn, bloom_stats_warning, &
                                 bloom_stats_t, accumulate_bloom_stats, write_bloom_report
        use determinants, only: det_info_t, alloc_det_info_t, dealloc_det_info_t
        use excitations, only: excit_t, create_excited_det, get_excitation
        use annihilation, only: direct_annihilation, direct_annihilation_received_list, &
                                direct_annihilation_spawned_list, deterministic_annihilation
        use calc, only: doing_calc
        use death, only: stochastic_death
        use ifciqmc, only: set_parent_flag
        use importance_sampling, only: importance_sampling_weight
        use non_blocking_comm_m, only: init_non_blocking_comm, end_non_blocking_comm
        use spawning, only: create_spawned_particle_initiator
        use qmc, only: init_qmc
        use qmc_common
        use dSFMT_interface, only: dSFMT_t, dSFMT_init
        use semi_stoch, only: semi_stoch_t, check_if_determ, determ_projection
        use semi_stoch, only: dealloc_semi_stoch_t, init_semi_stoch_t, init_semi_stoch_t_flags, set_determ_info
        use system, only: sys_t, sys_t_json
        use restart_hdf5, only: init_restart_info_t, restart_info_t, dump_restart_hdf5
        use spawn_data, only: receive_spawned_walkers, non_blocking_send, annihilate_wrapper_non_blocking_spawn, &
                              write_memcheck_report

        use qmc_data, only: qmc_in_t, fciqmc_in_t, semi_stoch_in_t, restart_in_t, load_bal_in_t, empty_determ_space, &
                            qmc_state_t, annihilation_flags_t, reference_t, semi_stoch_separate_annihilation,        &
                            qmc_in_t_json, fciqmc_in_t_json, semi_stoch_in_t_json, restart_in_t_json, load_bal_in_t_json, &
                            reference_t_json
        use check_input, only: check_qmc_opts, check_fciqmc_opts, check_load_bal_opts

        type(sys_t), intent(in) :: sys
        type(qmc_in_t), intent(in) :: qmc_in
        type(fciqmc_in_t), intent(in) :: fciqmc_in
        type(semi_stoch_in_t), intent(in) :: semi_stoch_in
        type(restart_in_t), intent(in) :: restart_in
        type(load_bal_in_t), intent(in) :: load_bal_in
        type(reference_t), intent(in) :: reference_in

        type(det_info_t) :: cdet
        type(dSFMT_t) :: rng
        type(bloom_stats_t) :: bloom_stats
        type(semi_stoch_t) :: determ
        type(json_out_t) :: js
        type(qmc_in_t) :: qmc_in_loc

        integer :: idet, ireport, icycle, iparticle, ideterm, ierr
        integer :: iter, semi_stoch_iter
        integer(int_64) :: nattempts
        real(p), allocatable :: nparticles_old(:)

        integer(i0) :: f_child(sys%basis%string_len)
        integer(int_p) :: nspawned, ndeath
        integer :: nattempts_current_det, nspawn_events
        type(excit_t) :: connection
        real(p) :: hmatel
        real(p) :: real_population, weighted_population
        integer :: send_counts(0:nprocs-1), req_data_s(0:nprocs-1)
        type(qmc_state_t), target :: qs
        type(annihilation_flags_t) :: annihilation_flags
        type(restart_info_t) :: ri, ri_shift

        logical :: soft_exit, write_restart_shift, error
        logical :: determ_parent, determ_child

        real :: t1, t2

        logical :: update_tau

        if (parent) then
            write (6,'(1X,"FCIQMC")')
            write (6,'(1X,"------",/)')
        end if

        if (parent) then
            ! Check input options.
            call check_qmc_opts(qmc_in, .false.)
            call check_fciqmc_opts(sys, fciqmc_in)
            call check_load_bal_opts(load_bal_in)
        end if

        ! Initialise data.
        call init_qmc(sys, qmc_in, restart_in, load_bal_in, reference_in, annihilation_flags, qs, fciqmc_in=fciqmc_in)

        if (parent) then
            call json_object_init(js, tag=.true.)
            call sys_t_json(js, sys)
            ! The default values of pattempt_* are not in qmc_in
            qmc_in_loc = qmc_in
            qmc_in_loc%pattempt_single = qs%pattempt_single
            qmc_in_loc%pattempt_double = qs%pattempt_double
            call qmc_in_t_json(js, qmc_in_loc)
            call fciqmc_in_t_json(js, fciqmc_in)
            call semi_stoch_in_t_json(js, semi_stoch_in)
            call restart_in_t_json(js, restart_in)
            call load_bal_in_t_json(js, load_bal_in)
            call reference_t_json(js, qs%ref, .true.)
            call json_object_end(js, terminal=.true., tag=.true.)
            write (js%io, '()')
        end if

        allocate(nparticles_old(qs%psip_list%nspaces), stat=ierr)
        call check_allocate('nparticles_old', size(nparticles_old), ierr)

        call dSFMT_init(qmc_in%seed+iproc, 50000, rng)

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

        call init_semi_stoch_t_flags(determ, size(qs%psip_list%states, dim=2))

        ! from restart
        nparticles_old = qs%psip_list%tot_nparticles

        ! Main fciqmc loop.
        if (parent) call write_fciqmc_report_header(qs%psip_list%nspaces)

        if (fciqmc_in%non_blocking_comm) then
            call init_non_blocking_comm(qs%spawn_store%spawn, req_data_s, send_counts, qs%spawn_store%spawn_recv, &
                                        restart_in%read_restart)
            call initial_fciqmc_status(sys, qmc_in, qs, .true., send_counts(iproc)/qs%spawn_store%spawn_recv%element_len)
        else
            call initial_fciqmc_status(sys, qmc_in, qs)
        end if
        ! Initialise timer.
        call cpu_time(t1)

        do ireport = 1, qmc_in%nreport

            ! Zero report cycle quantities.
            call init_report_loop(qs, bloom_stats)

            do icycle = 1, qmc_in%ncycles

                iter = qs%mc_cycles_done + (ireport-1)*qmc_in%ncycles + icycle

                ! Should we turn semi-stochastic on now?
                if (iter == semi_stoch_iter .and. semi_stoch_in%space_type /= empty_determ_space) then
                    determ%doing_semi_stoch = .true.
                    call init_semi_stoch_t(determ, semi_stoch_in, sys, qs%psip_list, qs%ref, annihilation_flags, &
                                           qs%spawn_store%spawn, qmc_in%use_mpi_barriers)
                end if

                call init_mc_cycle(qs%psip_list, qs%spawn_store%spawn, nattempts, ndeath)
                call load_balancing_wrapper(sys, qmc_in, qs%ref, load_bal_in, annihilation_flags, fciqmc_in%non_blocking_comm, &
                                            rng, qs%psip_list, qs%spawn_store%spawn, qs%par_info, determ)
                if (fciqmc_in%non_blocking_comm) qs%spawn_store%spawn_recv%proc_map = qs%par_info%load%proc_map
                ideterm = 0

                do idet = 1, qs%psip_list%nstates ! loop over walkers/dets

                    cdet%f => qs%psip_list%states(:,idet)
                    cdet%data => qs%psip_list%dat(:,idet)

                    call decoder_ptr(sys, cdet%f, cdet)

                    ! Extract the real sign from the encoded sign.
                    real_population = real(qs%psip_list%pops(1,idet),p)/qs%psip_list%pop_real_factor
                    weighted_population = importance_sampling_weight(qs%trial, cdet, real_population)

                    ! If this is a deterministic state then copy its population
                    ! across to the determ%vector array.
                    call set_determ_info(idet, real_population, ideterm, determ, determ_parent)

                    ! It is much easier to evaluate the projected energy at the
                    ! start of the i-FCIQMC cycle than at the end, as we're
                    ! already looping over the determinants.
                    connection = get_excitation(sys%nel, sys%basis, cdet%f, qs%ref%f0)
                    call update_proj_energy_ptr(sys, qs%ref%f0, qs%trial%wfn_dat, cdet, weighted_population, &
                                                qs%estimators%D0_population, qs%estimators%proj_energy, connection, hmatel)

                    ! Is this determinant an initiator?
                    call set_parent_flag(real_population, qmc_in%initiator_pop, cdet%f, determ%flags(idet), cdet%initiator_flag)

                    nattempts_current_det = decide_nattempts(rng, real_population)

                    do iparticle = 1, nattempts_current_det

                        ! Attempt to spawn.
                        call spawner_ptr(rng, sys, qs, qs%spawn_store%spawn%cutoff, qs%psip_list%pop_real_factor, &
                                        cdet, qs%psip_list%pops(1,idet), gen_excit_ptr, qs%trial%wfn_dat, nspawned, connection)

                        ! Spawn if attempt was successful.
                        if (nspawned /= 0_int_p) then
                            if (determ_parent) then
                                call create_excited_det(sys%basis, cdet%f, connection, f_child)
                                determ_child = check_if_determ(determ%hash_table, determ%dets, f_child)
                                ! If the spawning is both from and to the
                                ! deterministic space, cancel it.
                                if (.not. determ_child) then
                                    call create_spawned_particle_ptr(sys%basis, qs%ref, cdet, connection, nspawned, &
                                                                     1, qs%spawn_store%spawn, f_child)
                                else
                                    nspawned = 0_int_p
                                end if
                            else
                                call create_spawned_particle_ptr(sys%basis, qs%ref, cdet, connection, nspawned, 1, &
                                                                 qs%spawn_store%spawn)
                            end if
                            if (abs(nspawned) >= bloom_stats%nparticles_encoded) &
                                call accumulate_bloom_stats(bloom_stats, nspawned)
                        end if

                    end do

                    ! Clone or die.
                    if (.not. determ_parent) call stochastic_death(rng, qs, qs%psip_list%dat(1,idet), qs%shift(1), &
                                                        qs%psip_list%pops(1,idet), qs%psip_list%nparticles(1), ndeath)

                end do

                associate(pl=>qs%psip_list, spawn=>qs%spawn_store%spawn, spawn_recv=>qs%spawn_store%spawn_recv)
                    if (fciqmc_in%non_blocking_comm) then
                        call receive_spawned_walkers(spawn_recv, req_data_s)
                        call evolve_spawned_walkers(sys, qmc_in, qs, spawn_recv, spawn, cdet, rng, ndeath, &
                                                    load_bal_in%nslots)
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

                        call end_mc_cycle(nspawn_events, ndeath, pl%pop_real_factor, nattempts, qs%spawn_store%rspawn)
                    end if
                end associate

            end do

            update_tau = bloom_stats%nblooms_curr > 0

            error = qs%spawn_store%spawn%error .or. qs%psip_list%error

            call end_report_loop(sys, qmc_in, iter, update_tau, qs, nparticles_old, &
                                 nspawn_events, semi_stoch_in%shift_iter, semi_stoch_iter, soft_exit, &
                                 load_bal_in, bloom_stats=bloom_stats, doing_lb=fciqmc_in%doing_load_balancing, &
                                 nb_comm=fciqmc_in%non_blocking_comm, error=error)
            if (error) exit

            if (update_tau) call rescale_tau(qs%tau)

            call cpu_time(t2)
            if (parent) then
                if (bloom_stats%nblooms_curr > 0) call bloom_stats_warning(bloom_stats)
                call write_fciqmc_report(qmc_in, qs, ireport, nparticles_old, t2-t1, .false., &
                                         fciqmc_in%non_blocking_comm)
            end if

            ! Update the time for the start of the next iteration.
            t1 = t2

            call dump_restart_file_wrapper(qs, write_restart_shift, restart_in%write_freq, nparticles_old, ireport, &
                                           qmc_in%ncycles, sys%basis%nbasis, ri, ri_shift, fciqmc_in%non_blocking_comm)

            if (soft_exit) exit

            ! Should we try and update the reference determinant now?
            if (mod(ireport, fciqmc_in%select_ref_det_every_nreports) == 0) &
                    call select_ref_det(sys, fciqmc_in%ref_det_factor, qs)

        end do

        if (fciqmc_in%non_blocking_comm) call end_non_blocking_comm(sys, rng, qmc_in, annihilation_flags, ireport, &
                                                                    qs, qs%spawn_store%spawn_recv,  req_data_s,  &
                                                                    qs%par_info%report_comm%request, t1, nparticles_old, &
                                                                    qs%shift(1), restart_in%write_restart, load_bal_in)

        if (parent) write (6,'()')
        call write_bloom_report(bloom_stats)
        associate(pl=>qs%psip_list, spawn=>qs%spawn_store%spawn)
            if (determ%doing_semi_stoch .and. determ%projection_mode == semi_stoch_separate_annihilation) then
                call load_balancing_report(pl%nparticles, pl%nstates, qmc_in%use_mpi_barriers, spawn%mpi_time, determ%mpi_time)
            else
                call load_balancing_report(pl%nparticles, pl%nstates, qmc_in%use_mpi_barriers, spawn%mpi_time)
            end if
        end associate
        call write_memcheck_report(qs%spawn_store%spawn)

        if (soft_exit .or. error) then
            qs%mc_cycles_done = qs%mc_cycles_done + qmc_in%ncycles*ireport
        else
            qs%mc_cycles_done = qs%mc_cycles_done + qmc_in%ncycles*qmc_in%nreport
        end if

        if (restart_in%write_restart) then
            call dump_restart_hdf5(ri, qs, qs%mc_cycles_done, nparticles_old, sys%basis%nbasis, fciqmc_in%non_blocking_comm)
            if (parent) write (6,'()')
        end if

        if (determ%doing_semi_stoch) call dealloc_semi_stoch_t(determ, .false.)

        call dealloc_det_info_t(cdet, .false.)

    end subroutine do_fciqmc

    subroutine evolve_spawned_walkers(sys, qmc_in, qs, spawn_recv, spawn_to_send, cdet, rng, ndeath, nload_slots)

        ! Evolve spawned list of walkers one time step.
        ! Used for non-blocking communications.

        ! In:
        !   sys: system being studied.
        !   qmc_in: input options relating to QMC methods.
        !   nload_slots: number of load balancing slots (per processor).
        ! In/Out:
        !   qs: qmc_state_t containing information about the reference det and estimators.
        !   spawn: spawn_t object containing walkers spawned onto this processor during previous time step.
        !   cdet: type containing information about determinant. (easier to take this in as it is allocated
        !        / deallocated in do_fciqmc).
        !   rng: random number generator.
        !   ndeath: running total of number of particles which have died or been cloned.

        use proc_pointers, only: sc0_ptr
        use death, only: stochastic_death
        use determinants, only: det_info_t
        use dSFMT_interface, only: dSFMT_t
        use excitations, only: excit_t, get_excitation
        use ifciqmc, only: set_parent_flag
        use qmc_data, only: qmc_in_t, qmc_state_t
        use system, only: sys_t
        use qmc_common, only: decide_nattempts

        type(sys_t), intent(in) :: sys
        type(qmc_in_t), intent(in) :: qmc_in
        type(qmc_state_t), intent(inout) :: qs
        type(spawn_t), intent(inout) :: spawn_recv, spawn_to_send
        type(dSFMT_t), intent(inout) :: rng
        type(det_info_t), intent(inout) :: cdet
        integer(int_p), intent(inout) :: ndeath
        integer, intent(in) :: nload_slots

        real(p), target :: tmp_data(1)
        type(excit_t) :: connection
        real(p) :: hmatel
        integer :: idet, iparticle, nattempts_current_det
        integer(int_p) :: nspawned
        integer(int_p) :: int_pop(spawn_recv%ntypes)
        real(p) :: real_pop, weighted_pop
        real(p) :: list_pop
        integer(i0), target :: ftmp(sys%basis%tensor_label_len)

        cdet%f => ftmp

        do idet = 1, spawn_recv%head(0,0) ! loop over walkers/dets

            int_pop = int(spawn_recv%sdata(spawn_recv%bit_str_len+1:spawn_recv%bit_str_len+spawn_recv%ntypes, idet), int_p)
            real_pop = real(int_pop(1),p) / qs%psip_list%pop_real_factor
            ftmp = int(spawn_recv%sdata(:sys%basis%tensor_label_len,idet),i0)
            ! Need to generate spawned walker data to perform evolution.
            tmp_data(1) = sc0_ptr(sys, cdet%f) - qs%ref%H00
            cdet%data => tmp_data

            call decoder_ptr(sys, cdet%f, cdet)

            ! It is much easier to evaluate the projected energy at the
            ! start of the i-FCIQMC cycle than at the end, as we're
            ! already looping over the determinants.
            connection = get_excitation(sys%nel, sys%basis, cdet%f, qs%ref%f0)
            call update_proj_energy_ptr(sys, qs%ref%f0, qs%trial%wfn_dat, cdet, real_pop, qs%estimators%D0_population, &
                                        qs%estimators%proj_energy, connection, hmatel)

            ! Is this determinant an initiator?
            ! [todo] - pass determ_flag rather than 1.
            call set_parent_flag(real_pop, qmc_in%initiator_pop, cdet%f, 1, cdet%initiator_flag)

            nattempts_current_det = decide_nattempts(rng, real_pop)

            ! Possibly redundant if only one walker spawned at each spawning event.
            do iparticle = 1, nattempts_current_det

                ! Attempt to spawn.
                call spawner_ptr(rng, sys, qs, spawn_to_send%cutoff, qs%psip_list%pop_real_factor, cdet, int_pop(1), &
                                 gen_excit_ptr, qs%trial%wfn_dat, nspawned, connection)

                ! Spawn if attempt was successful.
                if (nspawned /= 0) then
                    call create_spawned_particle_ptr(sys%basis, qs%ref, cdet, connection, nspawned, 1, spawn_to_send)
                end if

            end do

            ! Clone or die.
            ! list_pop is meaningless as particle_t%nparticles is updated upon annihilation.
            call stochastic_death(rng, qs, tmp_data(1), qs%shift(1), int_pop(1), list_pop, ndeath)
            ! Update population of walkers on current determinant.
            spawn_recv%sdata(spawn_recv%bit_str_len+1:spawn_recv%bit_str_len+spawn_recv%ntypes, idet) = int_pop

        end do

    end subroutine evolve_spawned_walkers

end module fciqmc

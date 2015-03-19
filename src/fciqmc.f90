module fciqmc

! Module for performing optimised (hopefully!) full configuration interaction
! quantum monte carlo (FCIQMC) calculations.

use fciqmc_data
use proc_pointers
implicit none

contains

    subroutine do_fciqmc(sys, qmc_in, fciqmc_in, semi_stoch_in, restart_in, &
                         load_bal_in, reference, annihilation_flags)

        ! Run the FCIQMC or initiator-FCIQMC algorithm starting from the initial walker
        ! distribution using the timestep algorithm.

        ! See notes about the implementation of this using function pointers
        ! in fciqmc_main.

        ! In:
        !    sys: system being studied.
        !    semi_stoch_in: Input options for the semi-stochastic adaptation.
        !    fciqmc_in: input options relating to FCIQMC.
        !    restart_in: input options for HDF5 restart files.
        !    load_bal_in: input options for load balancing.
        !    annihilation_flags: calculation specific annihilation flags.
        ! In/Out:
        !    qmc_in: input options relating to QMC methods.
        !    reference: reference determinant. May change during the calculation.

        use parallel

        use bloom_handler, only: init_bloom_stats_t, bloom_mode_fixedn, &
                                 bloom_stats_t, accumulate_bloom_stats, write_bloom_report
        use determinants, only: det_info_t, alloc_det_info_t, dealloc_det_info_t
        use excitations, only: excit_t, create_excited_det, get_excitation
        use annihilation, only: direct_annihilation, direct_annihilation_received_list, &
                                direct_annihilation_spawned_list, deterministic_annihilation
        use calc, only: doing_calc, use_mpi_barriers
        use death, only: stochastic_death
        use non_blocking_comm_m, only: init_non_blocking_comm, end_non_blocking_comm
        use spawning, only: create_spawned_particle_initiator
        use qmc_common
        use dSFMT_interface, only: dSFMT_t, dSFMT_init
        use utils, only: rng_init_info
        use semi_stoch, only: semi_stoch_t, check_if_determ, determ_projection, determ_projection_separate_annihil
        use semi_stoch, only: dealloc_semi_stoch_t, init_semi_stoch_t
        use system, only: sys_t
        use restart_hdf5, only: restart_info_global, dump_restart_hdf5
        use spawn_data, only: receive_spawned_walkers, non_blocking_send, annihilate_wrapper_non_blocking_spawn

        use qmc_data, only: qmc_in_t, fciqmc_in_t, semi_stoch_in_t, restart_in_t, load_bal_in_t
        use qmc_data, only: reference_t, empty_determ_space, walker_global, annihilation_flags_t

        type(sys_t), intent(in) :: sys
        type(qmc_in_t), intent(inout) :: qmc_in
        type(fciqmc_in_t), intent(inout) :: fciqmc_in
        type(semi_stoch_in_t), intent(in) :: semi_stoch_in
        type(restart_in_t), intent(in) :: restart_in
        type(load_bal_in_t), intent(in) :: load_bal_in
        type(reference_t), intent(inout) :: reference
        type(annihilation_flags_t), intent(in) :: annihilation_flags

        type(det_info_t) :: cdet
        type(dSFMT_t) :: rng
        type(bloom_stats_t) :: bloom_stats
        type(semi_stoch_t) :: determ

        integer :: idet, ireport, icycle, iparticle, ideterm
        integer :: iter, semi_stoch_iter
        integer(int_64) :: nattempts
        real(p) :: nparticles_old(walker_global%sampling_size)

        integer(i0) :: f_child(sys%basis%string_len)
        integer(int_p) :: nspawned, ndeath
        integer :: nattempts_current_det, nspawn_events
        type(excit_t) :: connection
        real(p) :: hmatel
        real(p) :: real_population
        integer :: send_counts(0:nprocs-1), req_data_s(0:nprocs-1)

        logical :: soft_exit, dump_restart_file_shift
        logical :: semi_stochastic, determ_parent, determ_child

        real :: t1

        logical :: update_tau

        if (parent) call rng_init_info(qmc_in%seed+iproc)
        call dSFMT_init(qmc_in%seed+iproc, 50000, rng)

        ! Initialise bloom_stats components to the following parameters.
        call init_bloom_stats_t(bloom_stats, mode=bloom_mode_fixedn, encoding_factor=real_factor)

        ! Allocate det_info_t components.
        call alloc_det_info_t(sys, cdet, .false.)

        ! The iteration on which to start performing semi-stochastic.
        semi_stoch_iter = semi_stoch_in%start_iter

        ! Should we dump a restart file just before the shift is turned on?
        dump_restart_file_shift = restart_in%dump_restart_file_shift

        ! Create the semi_stoch_t object, determ.
        ! If the user has asked to use semi-stochastic from the first iteration
        ! then turn it on now. Otherwise, use an empty deterministic space.
        if (semi_stoch_iter == 0) then
            call init_semi_stoch_t(determ, sys, reference, annihilation_flags, qmc_spawn, &
                                   semi_stoch_in%determ_space_type, semi_stoch_in%target_size, &
                                   semi_stoch_in%separate_annihil, use_mpi_barriers, semi_stoch_in%write_determ_space, &
                                   load_bal_in%nslots)
        else
            call init_semi_stoch_t(determ, sys, reference, annihilation_flags, qmc_spawn, empty_determ_space, 0, &
                                   .false., .false., .false., load_bal_in%nslots)
        end if

        ! In case this is not set.
        nspawn_events = 0

        ! Are we using a non-empty semi-stochastic space?
        semi_stochastic = .not. (determ%space_type == empty_determ_space)

        ! from restart
        nparticles_old = walker_global%tot_nparticles

        ! Main fciqmc loop.
        if (parent) call write_fciqmc_report_header()

        if (fciqmc_in%non_blocking_comm) then
            call init_non_blocking_comm(qmc_spawn, req_data_s, send_counts, received_list, restart_in%read_restart)
            call initial_fciqmc_status(sys, qmc_in, reference, .true., par_info%report_comm, &
                                       send_counts(iproc)/received_list%element_len)
        else
            call initial_fciqmc_status(sys, qmc_in, reference)
        end if
        ! Initialise timer.
        call cpu_time(t1)

        do ireport = 1, qmc_in%nreport

            ! Zero report cycle quantities.
            call init_report_loop(bloom_stats)

            do icycle = 1, qmc_in%ncycles

                iter = mc_cycles_done + (ireport-1)*qmc_in%ncycles + icycle

                ! Should we turn semi-stochastic on now?
                if (iter == semi_stoch_iter) then
                    call dealloc_semi_stoch_t(determ, .false.)
                    call init_semi_stoch_t(determ, sys, reference, annihilation_flags, qmc_spawn, semi_stoch_in%determ_space_type, &
                                           semi_stoch_in%target_size, semi_stoch_in%separate_annihil, &
                                           use_mpi_barriers, semi_stoch_in%write_determ_space, &
                                           load_bal_in%nslots)
                    semi_stochastic = .true.
                end if

                call init_mc_cycle(rng, sys, qmc_in, reference, load_bal_in, annihilation_flags, real_factor, &
                                   nattempts, ndeath, doing_lb=fciqmc_in%doing_load_balancing, &
                                   nb_comm=fciqmc_in%non_blocking_comm, determ=determ)
                ideterm = 0

                do idet = 1, walker_global%tot_walkers ! loop over walkers/dets

                    cdet%f => walker_global%walker_dets(:,idet)
                    cdet%data => walker_global%walker_data(:,idet)

                    call decoder_ptr(sys, cdet%f, cdet)

                    ! Extract the real sign from the encoded sign.
                    real_population = real(walker_global%walker_population(1,idet),p)/real_factor

                    ! If this is a deterministic state then copy its population
                    ! across to the determ%vector array.
                    if (determ%flags(idet) == 0) then
                        ideterm = ideterm + 1
                        determ%vector(ideterm) = real_population
                        if (determ%separate_annihilation) determ%indices(ideterm) = idet
                        determ_parent = .true.
                    else
                        determ_parent = .false.
                    end if

                    ! It is much easier to evaluate the projected energy at the
                    ! start of the i-FCIQMC cycle than at the end, as we're
                    ! already looping over the determinants.
                    connection = get_excitation(sys%nel, sys%basis, cdet%f, reference%f0)
                    call update_proj_energy_ptr(sys, reference%f0, cdet, real_population, D0_population, &
                                                proj_energy, connection, hmatel)

                    ! Is this determinant an initiator?
                    call set_parent_flag_ptr(real_population, qmc_in%initiator_pop, cdet%f, determ%flags(idet), cdet%initiator_flag)

                    nattempts_current_det = decide_nattempts(rng, real_population)

                    do iparticle = 1, nattempts_current_det

                        ! Attempt to spawn.
                        call spawner_ptr(rng, sys, qmc_in, qmc_spawn%cutoff, real_factor, cdet, walker_global%walker_population(1,idet), &
                                         gen_excit_ptr, nspawned, connection)

                        ! Spawn if attempt was successful.
                        if (nspawned /= 0_int_p) then
                            if (determ_parent) then
                                call create_excited_det(sys%basis, cdet%f, connection, f_child)
                                determ_child = check_if_determ(determ%hash_table, determ%dets, f_child)
                                ! If the spawning is both from and to the
                                ! deterministic space, cancel it.
                                if (.not. determ_child) then
                                    call create_spawned_particle_ptr(sys%basis, reference, cdet, connection, nspawned, &
                                                                     1, qmc_spawn, load_bal_in%nslots, f_child)
                                else
                                    nspawned = 0_int_p
                                end if
                            else
                                call create_spawned_particle_ptr(sys%basis, reference, cdet, connection, nspawned, 1, qmc_spawn, &
                                                                 load_bal_in%nslots)
                            end if
                            if (abs(nspawned) >= bloom_stats%nparticles_encoded) &
                                call accumulate_bloom_stats(bloom_stats, nspawned)
                        end if

                    end do

                    ! Clone or die.
                    if (.not. determ_parent) call stochastic_death(rng, qmc_in%tau, walker_global%walker_data(1,idet), shift(1), &
                                                            walker_global%walker_population(1,idet), walker_global%nparticles(1), ndeath)

                end do

                if (fciqmc_in%non_blocking_comm) then
                    call receive_spawned_walkers(received_list, req_data_s)
                    call evolve_spawned_walkers(sys, qmc_in, reference, received_list, cdet, rng, ndeath, load_bal_in%nslots)
                    call direct_annihilation_received_list(sys, rng, qmc_in, reference, annihilation_flags, walker_global)
                    ! Need to add walkers which have potentially moved processor to the spawned walker list.
                    if (par_info%load%needed) then
                        call redistribute_particles(walker_global%walker_dets, real_factor,  walker_global%walker_population, walker_global%tot_walkers, &
                                                    walker_global%nparticles, qmc_spawn, load_bal_in%nslots)
                        par_info%load%needed = .false.
                    end if
                    call direct_annihilation_spawned_list(sys, rng, qmc_in, reference, annihilation_flags, walker_global, send_counts, &
                                                          req_data_s, par_info%report_comm%nb_spawn, nspawn_events)
                    call end_mc_cycle(par_info%report_comm%nb_spawn(1), ndeath, nattempts)
                else
                    ! If using semi-stochastic then perform the deterministic
                    ! projection step.
                    if (semi_stochastic) then
                        if (determ%separate_annihilation) then
                            call determ_projection_separate_annihil(determ, qmc_in%tau)
                            call deterministic_annihilation(sys, rng, walker_global, determ)
                        else
                            call determ_projection(rng, qmc_in, qmc_spawn, determ)
                        end if
                    end if

                    if (semi_stochastic) then
                        call direct_annihilation(sys, rng, qmc_in, reference, annihilation_flags, walker_global, nspawn_events, determ)
                    else
                        call direct_annihilation(sys, rng, qmc_in, reference, annihilation_flags, walker_global, nspawn_events)
                    end if
                    call end_mc_cycle(nspawn_events, ndeath, nattempts)
                end if

            end do

            update_tau = bloom_stats%nblooms_curr > 0

            call end_report_loop(sys, qmc_in, reference, ireport, iter, update_tau, nparticles_old, nspawn_events, t1, &
                                 semi_stoch_in%shift_iter, semi_stoch_iter, soft_exit, dump_restart_file_shift, &
                                 load_bal_in, bloom_stats=bloom_stats, doing_lb=fciqmc_in%doing_load_balancing, &
                                 nb_comm=fciqmc_in%non_blocking_comm, rep_comm=par_info%report_comm)

            if (soft_exit) exit

            ! Should we try and update the reference determinant now?
            if (mod(ireport, fciqmc_in%select_ref_det_every_nreports) == 0) &
                    call select_ref_det(sys, fciqmc_in%ref_det_factor, reference)

        end do

        if (fciqmc_in%non_blocking_comm) call end_non_blocking_comm(sys, rng, qmc_in, reference, annihilation_flags, ireport, &
                                                                    received_list,  req_data_s, par_info%report_comm%request, &
                                                                    t1, nparticles_old, shift(1), restart_in%dump_restart, &
                                                                    load_bal_in)

        if (parent) write (6,'()')
        call write_bloom_report(bloom_stats)
        if (semi_stochastic .and. determ%separate_annihilation) then
            call load_balancing_report(qmc_spawn%mpi_time, determ%mpi_time)
        else
            call load_balancing_report(qmc_spawn%mpi_time)
        end if

        if (soft_exit) then
            mc_cycles_done = mc_cycles_done + qmc_in%ncycles*ireport
        else
            mc_cycles_done = mc_cycles_done + qmc_in%ncycles*qmc_in%nreport
        end if

        if (restart_in%dump_restart) then
            call dump_restart_hdf5(restart_info_global, reference, mc_cycles_done, nparticles_old, fciqmc_in%non_blocking_comm)
            if (parent) write (6,'()')
        end if

        call dealloc_semi_stoch_t(determ, .false.)

        call dealloc_det_info_t(cdet, .false.)

    end subroutine do_fciqmc

    subroutine evolve_spawned_walkers(sys, qmc_in, reference, spawn, cdet, rng, ndeath, nload_slots)

        ! Evolve spawned list of walkers one time step.
        ! Used for non-blocking communications.

        ! In:
        !   sys: system being studied.
        !   qmc_in: input options relating to QMC methods.
        !   reference: current reference determinant.
        !   nload_slots: number of load balancing slots (per processor).
        ! In/Out:
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
        use qmc_data, only: qmc_in_t, reference_t, walker_global
        use system, only: sys_t
        use qmc_common, only: decide_nattempts

        type(sys_t), intent(in) :: sys
        type(qmc_in_t), intent(in) :: qmc_in
        type(reference_t), intent(in) :: reference
        type(spawn_t), intent(inout) :: spawn
        type(dSFMT_t), intent(inout) :: rng
        type(det_info_t), intent(inout) :: cdet
        integer(int_p), intent(inout) :: ndeath
        integer, intent(in) :: nload_slots

        real(p), target :: tmp_data(walker_global%sampling_size)
        type(excit_t) :: connection
        real(p) :: hmatel
        integer :: idet, iparticle, nattempts_current_det
        integer(int_p) :: nspawned
        integer(int_p) :: int_pop(spawn%ntypes)
        real(p) :: real_pop
        real(p) :: list_pop
        integer(i0), target :: ftmp(sys%basis%tensor_label_len)

        cdet%f => ftmp

        do idet = 1, spawn%head(0,0) ! loop over walkers/dets

            int_pop = int(spawn%sdata(spawn%bit_str_len+1:spawn%bit_str_len+spawn%ntypes, idet), int_p)
            real_pop = real(int_pop(1),p) / real_factor
            ftmp = int(spawn%sdata(:sys%basis%tensor_label_len,idet),i0)
            ! Need to generate spawned walker data to perform evolution.
            tmp_data(1) = sc0_ptr(sys, cdet%f) - reference%H00
            cdet%data => tmp_data

            call decoder_ptr(sys, cdet%f, cdet)

            ! It is much easier to evaluate the projected energy at the
            ! start of the i-FCIQMC cycle than at the end, as we're
            ! already looping over the determinants.
            connection = get_excitation(sys%nel, sys%basis, cdet%f, reference%f0)
            call update_proj_energy_ptr(sys, reference%f0, cdet, real_pop, D0_population, &
                                        proj_energy, connection, hmatel)

            ! Is this determinant an initiator?
            ! [todo] - pass determ_flag rather than 1.
            call set_parent_flag_ptr(real_pop, qmc_in%initiator_pop, cdet%f, 1, cdet%initiator_flag)

            nattempts_current_det = decide_nattempts(rng, real_pop)

            ! Possibly redundant if only one walker spawned at each spawning event.
            do iparticle = 1, nattempts_current_det

                ! Attempt to spawn.
                call spawner_ptr(rng, sys, qmc_in, qmc_spawn%cutoff, real_factor, cdet, int_pop(1), gen_excit_ptr, &
                                 nspawned, connection)

                ! Spawn if attempt was successful.
                if (nspawned /= 0) then
                    call create_spawned_particle_ptr(sys%basis, reference, cdet, connection, nspawned, 1, qmc_spawn, &
                                                     nload_slots)
                end if

            end do

            ! Clone or die.
            ! list_pop is meaningless as walker_global%nparticles is updated upon annihilation.
            call stochastic_death(rng, qmc_in%tau, tmp_data(1), shift(1), int_pop(1), list_pop, ndeath)
            ! Update population of walkers on current determinant.
            spawn%sdata(spawn%bit_str_len+1:spawn%bit_str_len+spawn%ntypes, idet) = int_pop

        end do

    end subroutine evolve_spawned_walkers

end module fciqmc

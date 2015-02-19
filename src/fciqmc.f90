module fciqmc

! Module for performing optimised (hopefully!) full configuration interaction
! quantum monte carlo (FCIQMC) calculations.

use fciqmc_data
use proc_pointers
implicit none

contains

    subroutine do_fciqmc(sys)

        ! Run the FCIQMC or initiator-FCIQMC algorithm starting from the initial walker
        ! distribution using the timestep algorithm.

        ! See notes about the implementation of this using function pointers
        ! in fciqmc_main.

        ! In:
        !    sys: system being studied.

        use parallel

        use bloom_handler, only: init_bloom_stats_t, bloom_mode_fixedn, &
                                 bloom_stats_t, accumulate_bloom_stats, write_bloom_report
        use determinants, only: det_info_t, alloc_det_info_t, dealloc_det_info_t
        use excitations, only: excit_t, create_excited_det
        use annihilation, only: direct_annihilation, direct_annihilation_received_list, &
                                direct_annihilation_spawned_list
        use calc, only: folded_spectrum, doing_calc, seed, initiator_approximation, non_blocking_comm, &
                        doing_load_balancing
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
        use energy_evaluation, only: update_energy_estimators_recv

        type(sys_t), intent(in) :: sys

        type(det_info_t) :: cdet
        type(dSFMT_t) :: rng
        type(bloom_stats_t) :: bloom_stats
        type(semi_stoch_t) :: determ

        integer :: idet, ireport, icycle, iparticle, ideterm
        integer :: iter
        integer(int_64) :: nattempts
        real(dp) :: nparticles_old(sampling_size)

        integer(i0) :: f_child(sys%basis%string_len)
        integer(int_p) :: nspawned, ndeath
        integer :: nattempts_current_det, nspawn_events
        type(excit_t) :: connection
        real(p) :: hmatel
        real(dp) :: real_population
        integer :: send_counts(0:nprocs-1), req_data_s(0:nprocs-1)

        logical :: soft_exit
        logical :: semi_stochastic, determ_parent, determ_child

        real :: t1

        logical :: update_tau

        if (parent) call rng_init_info(seed+iproc)
        call dSFMT_init(seed+iproc, 50000, rng)

        ! Initialise bloom_stats components to the following parameters.
        call init_bloom_stats_t(bloom_stats, mode=bloom_mode_fixedn, encoding_factor=real_factor)

        ! Allocate det_info_t components.
        call alloc_det_info_t(sys, cdet, .false.)

        ! Create the semi_stoch_t object, determ.
        ! If the user has asked to use semi-stochastic from the first iteration
        ! then turn it on now. Otherwise, use an empty deterministic space.
        if (semi_stoch_start_iter == 0) then
            call init_semi_stoch_t(determ, sys, qmc_spawn, determ_space_type, determ_target_size, &
                                    separate_determ_annihil, use_mpi_barriers, write_determ_space)
        else
            call init_semi_stoch_t(determ, sys, qmc_spawn, empty_determ_space, 0, .false., .false., .false.)
        end if

        ! Are we using a non-empty semi-stochastic space?
        semi_stochastic = .not. (determ%space_type == empty_determ_space)

        ! from restart
        nparticles_old = tot_nparticles

        ! Main fciqmc loop.
        if (parent) call write_fciqmc_report_header()

        if (non_blocking_comm) then
            call init_non_blocking_comm(qmc_spawn, req_data_s, send_counts, received_list, restart)
            call initial_fciqmc_status(sys, par_info%report_comm, send_counts(iproc)/received_list%element_len)
        else
            call initial_fciqmc_status(sys)
        end if
        ! Initialise timer.
        call cpu_time(t1)

        do ireport = 1, nreport

            ! Zero report cycle quantities.
            call init_report_loop(bloom_stats)

            do icycle = 1, ncycles

                iter = mc_cycles_done + (ireport-1)*ncycles + icycle

                ! Should we turn semi-stochastic on now?
                if (iter == semi_stoch_start_iter) then
                    call dealloc_semi_stoch_t(determ)
                    call init_semi_stoch_t(determ, sys, qmc_spawn, determ_space_type, determ_target_size, &
                                            separate_determ_annihil, use_mpi_barriers, write_determ_space)
                    semi_stochastic = .true.
                end if

                call init_mc_cycle(real_factor, nattempts, ndeath)
                ideterm = 0

                do idet = 1, tot_walkers ! loop over walkers/dets

                    cdet%f => walker_dets(:,idet)
                    cdet%data => walker_data(:,idet)

                    call decoder_ptr(sys, cdet%f, cdet)

                    ! Extract the real sign from the encoded sign.
                    real_population = real(walker_population(1,idet),dp)/real_factor

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
                    call update_proj_energy_ptr(sys, f0, cdet, real_population, D0_population, &
                                                proj_energy, connection, hmatel)

                    ! Is this determinant an initiator?
                    call set_parent_flag_ptr(real_population, cdet%f, determ%flags(idet), cdet%initiator_flag)

                    nattempts_current_det = decide_nattempts(rng, real_population)

                    do iparticle = 1, nattempts_current_det

                        ! Attempt to spawn.
                        call spawner_ptr(rng, sys, qmc_spawn%cutoff, real_factor, cdet, walker_population(1,idet), &
                                         gen_excit_ptr, nspawned, connection)

                        ! Spawn if attempt was successful.
                        if (nspawned /= 0_int_p) then
                            if (determ_parent) then
                                call create_excited_det(sys%basis, cdet%f, connection, f_child)
                                determ_child = check_if_determ(determ%hash_table, determ%dets, f_child)
                                ! If the spawning is both from and to the
                                ! deterministic space, cancel it.
                                if (.not. determ_child) then
                                    call create_spawned_particle_ptr(sys%basis, cdet, connection, nspawned, 1, qmc_spawn, f_child)
                                else
                                    nspawned = 0_int_p
                                end if
                            else
                                call create_spawned_particle_ptr(sys%basis, cdet, connection, nspawned, 1, qmc_spawn)
                            end if
                            if (abs(nspawned) >= bloom_stats%nparticles_encoded) &
                                call accumulate_bloom_stats(bloom_stats, nspawned)
                        end if

                    end do

                    ! Clone or die.
                    if (.not. determ_parent) call death_ptr(rng, walker_data(1,idet), shift(1), &
                                                            walker_population(1,idet), nparticles(1), ndeath)

                end do

                if (non_blocking_comm) then
                    call receive_spawned_walkers(received_list, req_data_s)
                    call evolve_spawned_walkers(sys, received_list, cdet, rng, ndeath)
                    call direct_annihilation_received_list(sys, rng, initiator_approximation)
                    ! Need to add walkers which have potentially moved processor to the spawned walker list.
                    if (doing_load_balancing) call redistribute_load_balancing_dets(walker_dets, real_factor, walker_population, &
                                                                        tot_walkers, nparticles, qmc_spawn, par_info%load%needed)
                    call direct_annihilation_spawned_list(sys, rng, initiator_approximation, send_counts, req_data_s, &
                                                          par_info%report_comm%nb_spawn)
                    call end_mc_cycle(par_info%report_comm%nb_spawn(1), ndeath, nattempts)
                else
                    if (doing_load_balancing) call redistribute_load_balancing_dets(walker_dets, real_factor, walker_population, &
                                                                        tot_walkers, nparticles, qmc_spawn, par_info%load%needed)
                    if (semi_stochastic) then
                        if (determ%separate_annihilation) then
                            call determ_projection_separate_annihil(determ)
                        else
                            call determ_projection(rng, qmc_spawn, determ)
                        end if
                        call direct_annihilation(sys, rng, initiator_approximation, nspawn_events, determ)
                    else
                        call direct_annihilation(sys, rng, initiator_approximation, nspawn_events)
                    end if
                    call end_mc_cycle(nspawn_events, ndeath, nattempts)
                end if

            end do

            update_tau = bloom_stats%nblooms_curr > 0

            call end_report_loop(sys, ireport, update_tau, nparticles_old, t1, soft_exit, bloom_stats=bloom_stats, &
                                 rep_comm=par_info%report_comm)

            if (soft_exit) exit

        end do

        if (non_blocking_comm) call end_non_blocking_comm(sys, rng, initiator_approximation, ireport, received_list, req_data_s,  &
                                                          par_info%report_comm%request, t1, nparticles_old, shift(1))

        if (parent) write (6,'()')
        call write_bloom_report(bloom_stats)
        if (semi_stochastic .and. determ%separate_annihilation) then
            call load_balancing_report(qmc_spawn%mpi_time, determ%mpi_time)
        else
            call load_balancing_report(qmc_spawn%mpi_time)
        end if

        if (soft_exit) then
            mc_cycles_done = mc_cycles_done + ncycles*ireport
        else
            mc_cycles_done = mc_cycles_done + ncycles*nreport
        end if

        if (dump_restart_file) then
            call dump_restart_hdf5(restart_info_global, mc_cycles_done, nparticles_old)
            if (parent) write (6,'()')
        end if

        call dealloc_semi_stoch_t(determ)

        call dealloc_det_info_t(cdet, .false.)

    end subroutine do_fciqmc

    subroutine evolve_spawned_walkers(sys, spawn, cdet, rng, ndeath)

        ! Evolve spawned list of walkers one time step.
        ! Used for non-blocking communications.

        ! In:
        !   sys: system being studied.
        ! In/Out:
        !   spawn: spawn_t object containing walkers spawned onto this processor during previous time step.
        !   cdet: type containing information about determinant. (easier to take this in as it is allocated
        !        / deallocated in do_fciqmc).
        !   rng: random number generator.
        !   ndeath: running total of number of particles which have died or been cloned.

        use proc_pointers, only: sc0_ptr
        use determinants, only: det_info_t
        use dSFMT_interface, only: dSFMT_t
        use excitations, only: excit_t
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        type(spawn_t), intent(inout) :: spawn
        type(dSFMT_t), intent(inout) :: rng
        type(det_info_t), intent(inout) :: cdet
        integer(int_p), intent(inout) :: ndeath

        ! [todo] - Check types with Nick's real coefficient work (which will probably be
        ! [todo] - merged before this).
        real(p), target :: tmp_data(sampling_size)
        type(excit_t) :: connection
        real(p) :: hmatel
        integer :: idet, iparticle
        integer(int_p) :: nspawned
        integer(int_p) :: pop(spawn%ntypes)
        real(dp) :: list_pop
        integer(i0), target :: ftmp(sys%basis%tensor_label_len)

        cdet%f => ftmp

        do idet = 1, spawn%head(0,0) ! loop over walkers/dets

            pop = int(spawn%sdata(spawn%bit_str_len+1:spawn%bit_str_len+spawn%ntypes, idet), int_p)
            ftmp = int(spawn%sdata(:sys%basis%tensor_label_len,idet),i0)
            ! Need to generate spawned walker data to perform evolution.
            tmp_data(1) = sc0_ptr(sys, cdet%f) - H00
            cdet%data => tmp_data

            ! [todo] - Population encoding and decoding when merged with the real coefficient work.

            call decoder_ptr(sys, cdet%f, cdet)

            ! It is much easier to evaluate the projected energy at the
            ! start of the i-FCIQMC cycle than at the end, as we're
            ! already looping over the determinants.
            call update_proj_energy_ptr(sys, f0, cdet, real(pop(1),p), D0_population, &
                                        proj_energy, connection, hmatel)

            ! Is this determinant an initiator?
            ! [todo] - pass determ_flag rather than 1.
            call set_parent_flag_ptr(real(pop(1),dp), cdet%f, 1, cdet%initiator_flag)

            ! Possibly redundant if only one walker spawned at each spawning event.
            do iparticle = 1, abs(pop(1))

                ! Attempt to spawn.
                call spawner_ptr(rng, sys, qmc_spawn%cutoff, real_factor, cdet, pop(1), gen_excit_ptr, nspawned, connection)

                ! Spawn if attempt was successful.
                if (nspawned /= 0) then
                    call create_spawned_particle_ptr(sys%basis, cdet, connection, nspawned, 1, qmc_spawn)
                end if

            end do

            ! Clone or die.
            ! list_pop is meaningless as nparticles is updated upon annihilation.
            call death_ptr(rng, tmp_data(1), shift(1), pop(1), list_pop, ndeath)
            ! Update population of walkers on current determinant.
            spawn%sdata(spawn%bit_str_len+1:spawn%bit_str_len+spawn%ntypes, idet) = pop

        end do

    end subroutine evolve_spawned_walkers

end module fciqmc

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

        use annihilation, only: direct_annihilation, direct_annihilation_non_blocking
        use basis, only: basis_length, nbasis
        use bloom_handler, only: bloom_stats_t, accumulate_bloom_stats
        use calc, only: folded_spectrum, doing_calc, seed, initiator_approximation, non_blocking_comm, &
                        nb_rep_t
        use determinants, only: det_info, alloc_det_info, dealloc_det_info
        use excitations, only: excit
        use spawning, only: create_spawned_particle_initiator
        use qmc_common
        use ifciqmc, only: set_parent_flag
        use folded_spectrum_utils, only: cdet_excit
        use dSFMT_interface, only: dSFMT_t, dSFMT_init
        use utils, only: rng_init_info
        use system, only: sys_t
        use restart_hdf5, only: restart_info_global, dump_restart_hdf5
        use spawn_data, only: receive_spawned_walkers, non_blocking_send, annihilate_wrapper_non_blocking_spawn
        use energy_evaluation, only: update_energy_estimators_recv

        type(sys_t), intent(in) :: sys

        integer :: idet, ireport, icycle, iparticle
        integer(lint) :: nattempts, nparticles_old(sampling_size)
        type(det_info) :: cdet
        type(dSFMT_t) :: rng
        type(bloom_stats_t) :: bloom_stats

        integer :: nspawned, ndeath
        type(excit) :: connection
        real(p) :: hmatel
        ! [review] - JSS: ir (and similarly named quantities) are really badly named.
        ! [review] - JSS: This is my fault but ir really isn't meaningful.
        ! [review] - JSS: Could it (and similarly named quantities) be given a better name?
        ! [reply] - FM: report_quant, report_vals, report_info or report_loop_info?
        integer :: send_counts(0:nprocs-1)
        ! [review] - JSS: derived type for the non-blocking energy evaluation variables?
        ! [reply] - FM: I thought about this, hopefully there aren't any weird
        ! [reply] - FM : mpi exceptions. Will investigate.
        ! [reply] - JSS: The simplest thing is to do MPI calls on the components of the derived type, to avoid pain with MPI interfaces.
        integer :: req_data_s(0:nprocs-1)

        logical :: soft_exit

        real :: t1

        logical :: update_tau

        if (parent) call rng_init_info(seed+iproc)
        call dSFMT_init(seed+iproc, 50000, rng)

        ! Allocate det_info components.
        call alloc_det_info(sys, cdet, .false.)
        ! Folded spectrum *needs* the bit strings to be allocated as it needs
        ! be able to manipulate the bit string to create excited states.
        if (doing_calc(folded_spectrum)) call alloc_det_info(sys, cdet_excit)

        ! from restart
        nparticles_old = tot_nparticles

        ! Main fciqmc loop.
        if (parent) call write_fciqmc_report_header()
        if (non_blocking_comm) then
            ! For non-blocking communications we need to initially send zero walkers
            ! to all processors this is so we don't have to deal with a special case
            ! in the receive step later on.
            send_counts = 0
            call non_blocking_send(qmc_spawn, send_counts, req_data_s)
            call initial_fciqmc_status(sys, report_comm)
        else
            call initial_fciqmc_status(sys)
        end if
        ! Initialise timer.
        call cpu_time(t1)

        do ireport = 1, nreport

            ! Zero report cycle quantities.
            call init_report_loop(bloom_stats)

            do icycle = 1, ncycles

                call init_mc_cycle(nattempts, ndeath)

                do idet = 1, tot_walkers ! loop over walkers/dets

                    cdet%f => walker_dets(:,idet)
                    cdet%data => walker_data(:,idet)

                    call decoder_ptr(sys, cdet%f, cdet)

                    ! It is much easier to evaluate the projected energy at the
                    ! start of the i-FCIQMC cycle than at the end, as we're
                    ! already looping over the determinants.
                    call update_proj_energy_ptr(sys, f0, cdet, real(walker_population(1,idet),p), D0_population_cycle, &
                                                proj_energy, connection, hmatel)

                    ! Is this determinant an initiator?
                    call set_parent_flag_ptr(walker_population(1,idet), cdet%f, cdet%initiator_flag)

                    do iparticle = 1, abs(walker_population(1,idet))

                        ! Attempt to spawn.
                        call spawner_ptr(rng, sys, cdet, walker_population(1,idet), gen_excit_ptr, nspawned, connection)

                        ! Spawn if attempt was successful.
                        if (nspawned /= 0) then
                            call create_spawned_particle_ptr(cdet, connection, nspawned, 1, qmc_spawn)
                            if (abs(nspawned) >= bloom_stats%n_bloom) &
                                call accumulate_bloom_stats(bloom_stats, nspawned)
                        end if

                    end do

                    ! Clone or die.
                    call death_ptr(rng, walker_data(1,idet), shift(1), walker_population(1,idet), nparticles(1), ndeath)

                end do

                if (non_blocking_comm) then
                    call receive_spawned_walkers(received_list, req_data_s)
                    call evolve_spawned_walkers(sys, received_list, cdet, rng, ndeath)
                    call direct_annihilation_non_blocking(sys, initiator_approximation, send_counts, req_data_s, report_comm%nb_spawn)
                    call end_mc_cycle(ndeath, nattempts, report_comm%nb_spawn)
                else
                    call direct_annihilation(sys, initiator_approximation)
                    call end_mc_cycle(ndeath, nattempts)
                end if

            end do

            update_tau = bloom_stats%nwarnings_curr > 0

            call end_report_loop(ireport, update_tau, nparticles_old, t1, soft_exit, report_comm)

            if (soft_exit) exit

        end do

        if (non_blocking_comm) call end_non_blocking_comm(sys, initiator_approximation, ireport, received_list, &
                                                          req_data_s, report_comm%request, t1, nparticles_old)

        if (parent) then
            call write_fciqmc_final(ireport)
            write (6,'()')
        end if

        call load_balancing_report()

        if (soft_exit) then
            mc_cycles_done = mc_cycles_done + ncycles*ireport
        else
            mc_cycles_done = mc_cycles_done + ncycles*nreport
        end if

        if (dump_restart_file) then
            call dump_restart_hdf5(restart_info_global, mc_cycles_done, nparticles_old)
            if (parent) write (6,'()')
        end if

        call dealloc_det_info(cdet, .false.)
        if (doing_calc(folded_spectrum)) call dealloc_det_info(cdet_excit)

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
        use determinants, only: det_info
        use dSFMT_interface, only: dSFMT_t
        use excitations, only: excit
        use system, only: sys_t
        use basis, only: total_basis_length

        type(sys_t), intent(in) :: sys
        type(spawn_t), intent(inout) :: spawn
        type(dSFMT_t), intent(inout) :: rng
        type(det_info), intent(inout) :: cdet
        integer, intent(inout) :: ndeath

        ! [todo] - Check types with Nick's real coefficient work (which will probably be
        ! [todo] - merged before this).
        real(p), target :: tmp_data(sampling_size)
        type(excit) :: connection
        real(p) :: hmatel
        integer :: idet, iparticle, nspawned
        integer :: pop(spawn%ntypes)
        integer(lint) :: list_pop

        do idet = 1, spawn%head(0,0) ! loop over walkers/dets

            pop = spawn%sdata(spawn%bit_str_len+1:spawn%bit_str_len+spawn%ntypes, idet)
            cdet%f => spawn%sdata(:total_basis_length,idet)
            ! Need to generate spawned walker data to perform evolution.
            tmp_data(1) = sc0_ptr(sys, cdet%f) - H00
            cdet%data => tmp_data

            ! [todo] - Population encoding and decoding when merged with the real coefficient work.

            call decoder_ptr(sys, cdet%f, cdet)

            ! It is much easier to evaluate the projected energy at the
            ! start of the i-FCIQMC cycle than at the end, as we're
            ! already looping over the determinants.
            call update_proj_energy_ptr(sys, f0, cdet, real(pop(1),p), D0_population_cycle, &
                                        proj_energy, connection, hmatel)

            ! Is this determinant an initiator?
            call set_parent_flag_ptr(pop(1), cdet%f, cdet%initiator_flag)

            ! Possibly redundant if only one walker spawned at each spawning event.
            do iparticle = 1, abs(pop(1))

                ! Attempt to spawn.
                call spawner_ptr(rng, sys, cdet, pop(1), gen_excit_ptr, nspawned, connection)

                ! Spawn if attempt was successful.
                if (nspawned /= 0) then
                    call create_spawned_particle_ptr(cdet, connection, nspawned, 1, qmc_spawn)
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

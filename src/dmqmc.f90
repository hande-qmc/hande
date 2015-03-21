module dmqmc

! Main loop for performing DMQMC calculations.

use fciqmc_data
use proc_pointers
implicit none

contains

    subroutine do_dmqmc(sys, qmc_in, dmqmc_in, restart_in, reference, load_bal_in, annihilation_flags)

        ! Run DMQMC calculation. We run from a beta=0 to a value of beta
        ! specified by the user and then repeat this main loop beta_loops
        ! times, to accumulate statistics for each value for beta.

        ! In:
        !    restart_in: input options for HDF5 restart files.
        !    reference: reference determinant.
        !    load_bal_in: input options for load balancing.
        !    annihilation_flags: calculation specific annihilation flags.
        ! In/Out:
        !    sys: system being studied.  NOTE: if modified inside a procedure,
        !         it should be returned in its original (ie unmodified state)
        !         at the end of the procedure.
        !    qmc_in: input options relating to QMC methods.
        !    dmqmc_in: input options relating to DMQMC.

        use parallel
        use annihilation, only: direct_annihilation
        use bit_utils, only: count_set_bits
        use bloom_handler, only: init_bloom_stats_t, bloom_mode_fixedn, &
                                 bloom_stats_t, accumulate_bloom_stats, write_bloom_report
        use death, only: stochastic_death
        use determinants, only: det_info_t, alloc_det_info_t, dealloc_det_info_t
        use dmqmc_estimators
        use dmqmc_procedures
        use excitations, only: excit_t
        use qmc_common
        use restart_hdf5, only: restart_info_global, dump_restart_hdf5
        use system
        use dSFMT_interface, only: dSFMT_t
        use utils, only: rng_init_info
        use qmc_data, only: qmc_in_t, restart_in_t, reference_t, load_bal_in_t, annihilation_flags_t, walker_global
        use dmqmc_data, only: dmqmc_in_t

        type(sys_t), intent(inout) :: sys
        type(qmc_in_t), intent(inout) :: qmc_in
        type(restart_in_t), intent(in) :: restart_in
        type(reference_t), intent(in) :: reference
        type(load_bal_in_t), intent(in) :: load_bal_in
        type(annihilation_flags_t), intent(in) :: annihilation_flags
        type(dmqmc_in_t), intent(inout) :: dmqmc_in

        integer :: idet, ireport, icycle, iparticle, iteration, ireplica
        integer :: beta_cycle
        integer :: unused_int_1 = -1, unused_int_2 = 0
        integer(int_64) :: init_tot_nparticles
        real(p) :: tot_nparticles_old(walker_global%sampling_size)
        real(p) :: real_population(walker_global%sampling_size)
        integer(int_64) :: nattempts
        integer :: nel_temp, nattempts_current_det
        type(det_info_t) :: cdet1, cdet2
        integer(int_p) :: nspawned, ndeath
        type(excit_t) :: connection
        integer :: spawning_end, nspawn_events
        logical :: soft_exit, dump_restart_file_shift
        real :: t1, t2
        type(dSFMT_t) :: rng
        type(bloom_stats_t) :: bloom_stats

        ! Allocate det_info_t components. We need two cdet objects for each 'end'
        ! which may be spawned from in the DMQMC algorithm.
        call alloc_det_info_t(sys, cdet1, .false.)
        call alloc_det_info_t(sys, cdet2, .false.)

        ! Initialise bloom_stats components to the following parameters.
        call init_bloom_stats_t(bloom_stats, mode=bloom_mode_fixedn, encoding_factor=real_factor)

        ! Main DMQMC loop.
        if (parent) then
            call rng_init_info(qmc_in%seed+iproc)
            call write_fciqmc_report_header(dmqmc_in)
        end if
        ! Initialise timer.
        call cpu_time(t1)

        ! When we accumulate data throughout a run, we are actually accumulating
        ! results from the psips distribution from the previous iteration.
        ! For example, in the first iteration, the trace calculated will be that
        ! of the initial distribution, which corresponds to beta=0. Hence, in the
        ! output we subtract one from the iteration number, and run for one more
        ! report loop, asimplemented in the line of code below.
        qmc_in%nreport = qmc_in%nreport+1

        if (dmqmc_in%all_spin_sectors) nel_temp = sys%nel
        init_tot_nparticles = nint(qmc_in%D0_population, int_64)

        ! Should we dump a restart file just before the shift is turned on?
        dump_restart_file_shift = restart_in%dump_restart_file_shift

        do beta_cycle = 1, dmqmc_in%beta_loops

            call init_dmqmc_beta_loop(rng, qmc_in, dmqmc_in, beta_cycle)

            ! Distribute psips uniformly along the diagonal of the density
            ! matrix.
            call create_initial_density_matrix(rng, sys, qmc_in, dmqmc_in, reference, annihilation_flags, &
                                               init_tot_nparticles, walker_global%tot_nparticles, load_bal_in%nslots)

            ! Allow the shift to vary from the very start of the beta loop, if
            ! this condition is met.
            vary_shift = walker_global%tot_nparticles >= qmc_in%target_particles

            do ireport = 1, qmc_in%nreport

                call init_report_loop(bloom_stats)
                tot_nparticles_old = walker_global%tot_nparticles

                do icycle = 1, qmc_in%ncycles

                    call init_mc_cycle(rng, sys, qmc_in, reference, load_bal_in, annihilation_flags, real_factor, &
                                       qmc_spawn, nattempts, ndeath)

                    iteration = (ireport-1)*qmc_in%ncycles + icycle

                    do idet = 1, walker_global%tot_walkers ! loop over walkers/dets

                        ! f points to the bitstring that is spawning, f2 to the
                        ! other bit string.
                        cdet1%f => walker_global%walker_dets(:sys%basis%string_len,idet)
                        cdet1%f2 => walker_global%walker_dets((sys%basis%string_len+1):(2*sys%basis%string_len),idet)
                        cdet2%f => walker_global%walker_dets((sys%basis%string_len+1):(2*sys%basis%string_len),idet)
                        cdet2%f2 => walker_global%walker_dets(:sys%basis%string_len,idet)

                        ! If using multiple symmetry sectors then find the
                        ! symmetry labels of this particular det.
                        if (dmqmc_in%all_spin_sectors) then
                            sys%nel = sum(count_set_bits(cdet1%f))
                            sys%nvirt = sys%lattice%nsites - sys%nel
                        end if

                        ! Decode and store the the relevant information for
                        ! both bitstrings. Both of these bitstrings are required
                        ! to refer to the correct element in the density matrix.
                        call decoder_ptr(sys, cdet1%f, cdet1)
                        call decoder_ptr(sys, cdet2%f, cdet2)

                        ! Extract the real signs from the encoded signs.
                        real_population = real(walker_global%walker_population(:,idet),p)/real_factor

                        ! Call wrapper function which calls routines to update
                        ! all estimators being calculated, and also always
                        ! updates the trace separately.
                        ! Note DMQMC averages over multiple loops over
                        ! temperature/imaginary time so only get data from one
                        ! temperature value per ncycles.
                        if (icycle == 1) then
                            call update_dmqmc_estimators(sys, dmqmc_in, idet, iteration, cdet1, &
                                                         reference%H00, load_bal_in%nslots)
                        end if

                        do ireplica = 1, walker_global%sampling_size

                            ! If this condition is met then there will only be
                            ! one det in this symmetry sector, so don't attempt
                            ! to spawn.
                            if (.not. ((sys%nel == 0) .or. (sys%nel == sys%basis%nbasis))) then
                                nattempts_current_det = decide_nattempts(rng, real_population(ireplica))
                                do iparticle = 1, nattempts_current_det
                                    ! When using importance sampling in DMQMC we
                                    ! only spawn from one end of a density
                                    ! matrix element.
                                    ! Spawn from the first end.
                                    spawning_end = 1
                                    ! Attempt to spawn.
                                    call spawner_ptr(rng, sys, qmc_in, qmc_spawn%cutoff, real_factor, cdet1, &
                                                     walker_global%walker_population(ireplica,idet), gen_excit_ptr, nspawned, connection)
                                    ! Spawn if attempt was successful.
                                    if (nspawned /= 0_int_p) then
                                        call create_spawned_particle_dm_ptr(sys%basis, cdet1%f, cdet2%f, connection, nspawned, &
                                                                            spawning_end, ireplica, qmc_spawn, load_bal_in%nslots)

                                        if (abs(nspawned) >= bloom_stats%nparticles_encoded) &
                                            call accumulate_bloom_stats(bloom_stats, nspawned)
                                    end if

                                    ! Now attempt to spawn from the second end.
                                    if (.not. dmqmc_in%propagate_to_beta) then
                                        spawning_end = 2
                                        call spawner_ptr(rng, sys, qmc_in, qmc_spawn%cutoff, real_factor, cdet2, &
                                                         walker_global%walker_population(ireplica,idet), gen_excit_ptr, nspawned, connection)
                                        if (nspawned /= 0_int_p) then
                                            call create_spawned_particle_dm_ptr(sys%basis, cdet2%f, cdet1%f, connection, nspawned, &
                                                                                spawning_end, ireplica, qmc_spawn, &
                                                                                load_bal_in%nslots)

                                            if (abs(nspawned) >= bloom_stats%nparticles_encoded) &
                                                call accumulate_bloom_stats(bloom_stats, nspawned)
                                        end if
                                    end if
                                end do
                            end if

                            ! Clone or die.
                            ! We have contributions to the clone/death step from
                            ! both ends of the current walker. We do both of
                            ! these at once by using walker_global%walker_data(:,idet) which,
                            ! when running a DMQMC algorithm, stores the average
                            ! of the two diagonal elements corresponding to the
                            ! two indicies of the density matrix.
                            call stochastic_death(rng, qmc_in%tau, walker_global%walker_data(ireplica,idet), shift(ireplica), &
                                           walker_global%walker_population(ireplica,idet), walker_global%nparticles(ireplica), ndeath)
                        end do
                    end do

                    ! Now we have finished looping over all determinants, set
                    ! the symmetry labels back to their default value, if
                    ! necessary.
                    if (dmqmc_in%all_spin_sectors) then
                        sys%nel = nel_temp
                        sys%nvirt = sys%lattice%nsites - sys%nel
                    end if

                    ! Perform the annihilation step where the spawned walker
                    ! list is merged with the main walker list, and walkers of
                    ! opposite sign on the same sites are annihilated.
                    call direct_annihilation(sys, rng, qmc_in, reference, annihilation_flags, walker_global, qmc_spawn, nspawn_events)

                    call end_mc_cycle(nspawn_events, ndeath, nattempts)

                    ! If doing importance sampling *and* varying the weights of
                    ! the trial function, call a routine to update these weights
                    ! and alter the number of psips on each excitation level
                    ! accordingly.
                    if (dmqmc_in%vary_weights .and. iteration <= dmqmc_in%finish_varying_weights) &
                        call update_sampling_weights(rng, sys%basis, qmc_in, walker_global)

                end do

                ! Sum all quantities being considered across all MPI processes.
                call dmqmc_estimate_comms(dmqmc_in, nspawn_events, sys%max_number_excitations, qmc_in%ncycles)

                call update_shift_dmqmc(qmc_in, walker_global%tot_nparticles, tot_nparticles_old, ireport)

                ! Forcibly disable update_tau as need to average over multiple loops over beta
                ! and hence want to use the same timestep throughout.
                call end_report_loop(sys, qmc_in, reference, ireport, iteration, .false., tot_nparticles_old, nspawn_events, t1, &
                                     unused_int_1, unused_int_2, soft_exit, dump_restart_file_shift, load_bal_in, &
                                     .false., bloom_stats=bloom_stats, dmqmc_in=dmqmc_in)

                if (soft_exit) exit

            end do

            if (soft_exit) exit

            ! Calculate and output all requested estimators based on the reduced
            ! density matrix. This is for ground-state RDMs only.
            if (calc_ground_rdm) call call_ground_rdm_procedures(beta_cycle)
            ! Calculate and output new weights based on the psip distirubtion in
            ! the previous loop.
            if (dmqmc_in%find_weights) call output_and_alter_weights(dmqmc_in, sys%max_number_excitations)

        end do

        if (parent) write (6,'()')
        call write_bloom_report(bloom_stats)
        call load_balancing_report(qmc_spawn%mpi_time)

        if (soft_exit) then
            mc_cycles_done = mc_cycles_done + qmc_in%ncycles*ireport
        else
            mc_cycles_done = mc_cycles_done + qmc_in%ncycles*qmc_in%nreport
        end if

        if (restart_in%dump_restart) then
            call dump_restart_hdf5(restart_info_global, walker_global, reference, mc_cycles_done, walker_global%tot_nparticles, .false.)
            if (parent) write (6,'()')
        end if

        call dealloc_det_info_t(cdet1, .false.)
        call dealloc_det_info_t(cdet2, .false.)

    end subroutine do_dmqmc

    subroutine init_dmqmc_beta_loop(rng, qmc_in, dmqmc_in, beta_cycle)

        ! Initialise/reset DMQMC data for a new run over the temperature range.

        ! In/Out:
        !    rng: random number generator.
        ! In:
        !    initial_shift: the initial shift used for population control.
        !    beta_cycle: The index of the beta loop about to be started.
        !    dmqmc_in: input options for DMQMC.

        use dSFMT_interface, only: dSFMT_t, dSFMT_init
        use parallel
        use qmc_data, only: qmc_in_t, walker_global
        use dmqmc_data, only: dmqmc_in_t
        use utils, only: int_fmt

        type(dSFMT_t) :: rng
        type(qmc_in_t), intent(in) :: qmc_in
        type(dmqmc_in_t), intent(in) :: dmqmc_in
        integer, intent(in) :: beta_cycle
        integer :: new_seed

        ! Reset the current position in the spawning array to be the slot
        ! preceding the first slot.
        qmc_spawn%head = qmc_spawn%head_start

        ! Set all quantities back to their starting values.
        walker_global%tot_walkers = 0
        shift = qmc_in%initial_shift
        walker_global%nparticles = 0.0_dp
        if (allocated(reduced_density_matrix)) reduced_density_matrix = 0.0_p
        if (dmqmc_in%vary_weights) accumulated_probs = 1.0_p
        if (dmqmc_in%find_weights) excit_dist = 0.0_p

        new_seed = qmc_in%seed+iproc+(beta_cycle-1)*nprocs

        if (beta_cycle /= 1 .and. parent) then
            write (6,'(a32,'//int_fmt(beta_cycle,1)//')') " # Resetting beta... Beta loop =", beta_cycle
            write (6,'(a52,'//int_fmt(new_seed,1)//',a1)') " # Resetting random number generator with a seed of:", new_seed, "."
        end if

        ! Reset the random number generator with new_seed = old_seed +
        ! nprocs (each beta loop)
        call dSFMT_init(new_seed, 50000, rng)

    end subroutine init_dmqmc_beta_loop

end module dmqmc

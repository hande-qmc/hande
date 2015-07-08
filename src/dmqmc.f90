module dmqmc

! Main loop for performing DMQMC calculations.

use fciqmc_data
use proc_pointers
implicit none

contains

    subroutine do_dmqmc(sys, qmc_in, dmqmc_in, dmqmc_estimates, restart_in, load_bal_in, reference_in)

        ! Run DMQMC calculation. We run from a beta=0 to a value of beta
        ! specified by the user and then repeat this main loop beta_loops
        ! times, to accumulate statistics for each value for beta.

        ! In:
        !    sys: system being studied.
        !    restart_in: input options for HDF5 restart files.
        !    load_bal_in: input options for load balancing.
        !    reference_in: current reference determinant.  If not set (ie
        !       components allocated) then a best guess is made based upon the
        !       desired spin/symmetry.
        ! In/Out:
        !    qmc_in: input options relating to QMC methods.
        !    dmqmc_in: input options relating to DMQMC.
        !    dmqmc_estimates: type containing all DMQMC estimates.

        use parallel
        use checking, only: check_allocate
        use annihilation, only: direct_annihilation
        use bit_utils, only: count_set_bits
        use bloom_handler, only: init_bloom_stats_t, bloom_mode_fixedn, bloom_stats_warning, &
                                 bloom_stats_t, accumulate_bloom_stats, write_bloom_report
        use death, only: stochastic_death
        use determinants, only: det_info_t, alloc_det_info_t, dealloc_det_info_t, decode_det, update_sys_spin_info
        use dmqmc_estimators
        use dmqmc_procedures
        use dmqmc_initialise_dm, only: create_initial_density_matrix
        use excitations, only: excit_t, connection_exists
        use qmc, only: init_qmc
        use qmc_common
        use restart_hdf5, only: restart_info_t, dump_restart_hdf5, init_restart_info_t
        use system
        use dSFMT_interface, only: dSFMT_t
        use utils, only: rng_init_info
        use qmc_data, only: qmc_in_t, restart_in_t, reference_t, load_bal_in_t, annihilation_flags_t, qmc_state_t
        use dmqmc_data, only: dmqmc_in_t, dmqmc_estimates_t, dmqmc_weighted_sampling_t
        use check_input, only: check_qmc_opts, check_dmqmc_opts
        use spawn_data, only: write_memcheck_report

        type(sys_t), intent(inout) :: sys
        type(qmc_in_t), intent(inout) :: qmc_in
        type(dmqmc_in_t), intent(inout) :: dmqmc_in
        type(dmqmc_estimates_t), intent(inout) :: dmqmc_estimates
        type(restart_in_t), intent(in) :: restart_in
        type(load_bal_in_t), intent(inout) :: load_bal_in
        type(reference_t), intent(in) :: reference_in

        integer :: idet, ireport, icycle, iparticle, iteration, ireplica, ierr
        integer :: beta_cycle
        integer :: unused_int_1 = -1, unused_int_2 = 0
        integer(int_64) :: init_tot_nparticles
        real(p), allocatable :: tot_nparticles_old(:), real_population(:)
        integer(int_64) :: nattempts
        integer :: nel_temp, nattempts_current_det
        type(det_info_t) :: cdet1, cdet2
        integer(int_p) :: nspawned, ndeath
        type(excit_t) :: connection
        integer :: spawning_end, nspawn_events
        logical :: soft_exit, write_restart_shift, update_tau
        logical :: error, rdm_error, attempt_spawning
        real :: t1, t2
        type(dSFMT_t) :: rng
        type(bloom_stats_t) :: bloom_stats
        type(qmc_state_t), target :: qs
        type(annihilation_flags_t) :: annihilation_flags
        type(dmqmc_weighted_sampling_t) :: weighted_sampling
        type(restart_info_t) :: ri
        integer :: ms, occ_list(sys%nel)
        type(sys_t) :: sys_copy

        if (parent) then
            write (6,'(1X,"DMQMC")')
            write (6,'(1X,"-----",/)')
        end if

        if (parent) then
            call check_qmc_opts(qmc_in, .false.)
            call check_dmqmc_opts(sys, dmqmc_in)
        end if

        ! Initialise data.
        call init_qmc(sys, qmc_in, restart_in, load_bal_in, reference_in, annihilation_flags, qs, dmqmc_in=dmqmc_in)

        allocate(tot_nparticles_old(qs%psip_list%nspaces), stat=ierr)
        call check_allocate('tot_nparticles_old', size(tot_nparticles_old), ierr)
        allocate(real_population(qs%psip_list%nspaces), stat=ierr)
        call check_allocate('real_population', size(real_population), ierr)

        ! Initialise all the required arrays, ie to store thermal quantities,
        ! and to initalise reduced density matrix quantities if necessary.
        call init_dmqmc(sys, qmc_in, dmqmc_in, qs%psip_list%nspaces, qs, dmqmc_estimates, weighted_sampling)

        ! Allocate det_info_t components. We need two cdet objects for each 'end'
        ! which may be spawned from in the DMQMC algorithm.
        call alloc_det_info_t(sys, cdet1, .false.)
        call alloc_det_info_t(sys, cdet2, .false.)

        ! Initialise bloom_stats components to the following parameters.
        call init_bloom_stats_t(bloom_stats, mode=bloom_mode_fixedn, encoding_factor=real_factor)

        ! Main DMQMC loop.
        if (parent) then
            call rng_init_info(qmc_in%seed+iproc)
            call write_fciqmc_report_header(qs%psip_list%nspaces, dmqmc_in, sys%max_number_excitations)
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
        write_restart_shift = restart_in%write_restart_shift
        call init_restart_info_t(ri, write_id=restart_in%write_id)

        ! If averaging over spin then spin polarisation will change in general.
        ! Copy sys to sys_copy so that sys can be set back to its input value
        ! at the end of this routine.
        call copy_sys_spin_info(sys, sys_copy)

        rdm_error = .false.

        outer_loop: do beta_cycle = 1, dmqmc_in%beta_loops

            call init_dmqmc_beta_loop(rng, qmc_in, dmqmc_in, dmqmc_estimates, qs, beta_cycle, qs%psip_list%nstates, &
                                      qs%psip_list%nparticles, qs%spawn_store%spawn, weighted_sampling%probs)

            ! Distribute psips uniformly along the diagonal of the density
            ! matrix.
            call create_initial_density_matrix(rng, sys, qmc_in, dmqmc_in, qs%ref, annihilation_flags, &
                                               init_tot_nparticles, qs%psip_list, qs%spawn_store%spawn)

            ! Allow the shift to vary from the very start of the beta loop, if
            ! this condition is met.
            qs%vary_shift = qs%psip_list%tot_nparticles >= qmc_in%target_particles

            do ireport = 1, qmc_in%nreport

                call init_dmqmc_report_loop(bloom_stats, dmqmc_estimates, qs%spawn_store%rspawn)
                tot_nparticles_old = qs%psip_list%tot_nparticles

                do icycle = 1, qmc_in%ncycles

                    call init_mc_cycle(qs%psip_list, qs%spawn_store%spawn, nattempts, ndeath)

                    iteration = (ireport-1)*qmc_in%ncycles + icycle

                    do idet = 1, qs%psip_list%nstates ! loop over walkers/dets

                        ! f points to the bitstring that is spawning, f2 to the
                        ! other bit string.
                        cdet1%f => qs%psip_list%states(:sys%basis%string_len,idet)
                        cdet1%f2 => qs%psip_list%states((sys%basis%string_len+1):(2*sys%basis%string_len),idet)
                        cdet1%data => qs%psip_list%dat(:,idet)
                        cdet2%f => qs%psip_list%states((sys%basis%string_len+1):(2*sys%basis%string_len),idet)
                        cdet2%f2 => qs%psip_list%states(:sys%basis%string_len,idet)

                        ! Decode and store the the relevant information for
                        ! both bitstrings. Both of these bitstrings are required
                        ! to refer to the correct element in the density matrix.
                        call decoder_ptr(sys, cdet1%f, cdet1)
                        call decoder_ptr(sys, cdet2%f, cdet2)

                        ! If using multiple symmetry sectors then find the
                        ! symmetry labels of this particular det.
                        if (dmqmc_in%all_spin_sectors) call update_sys_spin_info(cdet1, sys)

                        ! Extract the real signs from the encoded signs.
                        real_population = real(qs%psip_list%pops(:,idet),p)/real_factor

                        ! Call wrapper function which calls routines to update
                        ! all estimators being calculated, and also always
                        ! updates the trace separately.
                        ! Note DMQMC averages over multiple loops over
                        ! temperature/imaginary time so only get data from one
                        ! temperature value per ncycles.
                        if (icycle == 1) then
                            call update_dmqmc_estimators(sys, dmqmc_in, idet, iteration, cdet1, qs%ref%H00, &
                                                         qs%psip_list, dmqmc_estimates, weighted_sampling, rdm_error)
                        end if

                        ! Only attempt spawning if a valid connection exists.
                        attempt_spawning = connection_exists(sys)

                        do ireplica = 1, qs%psip_list%nspaces

                            ! Only attempt spawning if a valid excitation exists.
                            if (attempt_spawning) then
                                nattempts_current_det = decide_nattempts(rng, real_population(ireplica))
                                do iparticle = 1, nattempts_current_det
                                    ! When using importance sampling in DMQMC we
                                    ! only spawn from one end of a density
                                    ! matrix element.
                                    ! Spawn from the first end.
                                    spawning_end = 1
                                    ! Attempt to spawn.
                                    call spawner_ptr(rng, sys, qmc_in, qs%tau, qs%spawn_store%spawn%cutoff, real_factor, &
                                                     cdet1, qs%psip_list%pops(ireplica,idet), gen_excit_ptr, &
                                                     weighted_sampling%probs, nspawned, connection)
                                    ! Spawn if attempt was successful.
                                    if (nspawned /= 0_int_p) then
                                        call create_spawned_particle_dm_ptr(sys%basis, qs%ref, cdet1%f, cdet2%f, connection, &
                                                                            nspawned, spawning_end, ireplica, &
                                                                            qs%spawn_store%spawn)

                                        if (abs(nspawned) >= bloom_stats%nparticles_encoded) &
                                            call accumulate_bloom_stats(bloom_stats, nspawned)
                                    end if

                                    ! Now attempt to spawn from the second end.
                                    if (.not. dmqmc_in%propagate_to_beta) then
                                        spawning_end = 2
                                        call spawner_ptr(rng, sys, qmc_in, qs%tau, qs%spawn_store%spawn%cutoff, real_factor, &
                                                         cdet2, qs%psip_list%pops(ireplica,idet), gen_excit_ptr, &
                                                         weighted_sampling%probs, nspawned, connection)
                                        if (nspawned /= 0_int_p) then
                                            call create_spawned_particle_dm_ptr(sys%basis, qs%ref, cdet2%f, cdet1%f, &
                                                                                connection, nspawned, spawning_end, ireplica, &
                                                                                qs%spawn_store%spawn)

                                            if (abs(nspawned) >= bloom_stats%nparticles_encoded) &
                                                call accumulate_bloom_stats(bloom_stats, nspawned)
                                        end if
                                    end if
                                end do
                            end if

                            ! Clone or die.
                            ! We have contributions to the clone/death step from
                            ! both ends of the current walker. We do both of
                            ! these at once by using qs%psip_list%dat(:,idet) which,
                            ! when running a DMQMC algorithm, stores the average
                            ! of the two diagonal elements corresponding to the
                            ! two indicies of the density matrix.
                            call stochastic_death(rng, qs, qs%psip_list%dat(ireplica,idet), qs%shift(ireplica), &
                                           qs%psip_list%pops(ireplica,idet), qs%psip_list%nparticles(ireplica), ndeath)
                        end do
                    end do

                    ! Perform the annihilation step where the spawned walker
                    ! list is merged with the main walker list, and walkers of
                    ! opposite sign on the same sites are annihilated.
                    call direct_annihilation(sys, rng, qmc_in, qs%ref, annihilation_flags, qs%psip_list, &
                                             qs%spawn_store%spawn, nspawn_events)

                    call end_mc_cycle(nspawn_events, ndeath, nattempts, qs%spawn_store%rspawn)

                    ! If doing importance sampling *and* varying the weights of
                    ! the trial function, call a routine to update these weights
                    ! and alter the number of psips on each excitation level
                    ! accordingly.
                    if (dmqmc_in%vary_weights .and. iteration <= dmqmc_in%finish_varying_weights) &
                        call update_sampling_weights(rng, sys%basis, qmc_in, qs%psip_list, weighted_sampling)

                end do

                ! Sum all quantities being considered across all MPI processes.
                error = qs%spawn_store%spawn%error .or. qs%psip_list%error .or. rdm_error
                call dmqmc_estimate_comms(dmqmc_in, error, nspawn_events, sys%max_number_excitations, qmc_in%ncycles, &
                                          qs%psip_list, qs, weighted_sampling%probs_old, dmqmc_estimates)
                if (error) exit outer_loop

                call update_shift_dmqmc(qmc_in, qs, qs%psip_list%tot_nparticles, tot_nparticles_old)

                ! Forcibly disable update_tau as need to average over multiple loops over beta
                ! and hence want to use the same timestep throughout.
                update_tau = .false.
                call end_report_loop(sys, qmc_in, iteration, update_tau, qs, tot_nparticles_old, &
                                     nspawn_events, unused_int_1, unused_int_2, soft_exit, &
                                     load_bal_in, .false., bloom_stats=bloom_stats)

                call cpu_time(t2)
                if (parent) then
                    if (bloom_stats%nblooms_curr > 0) call bloom_stats_warning(bloom_stats)
                    call write_fciqmc_report(qmc_in, qs, ireport, tot_nparticles_old, t2-t1, .false., &
                                             .false., dmqmc_in, dmqmc_estimates)
                end if

                ! Update the time for the start of the next iteration.
                t1 = t2

                if (soft_exit) exit outer_loop

            end do

            ! Calculate and output all requested estimators based on the reduced
            ! density matrix. This is for ground-state RDMs only.
            if (dmqmc_in%rdm%calc_ground_rdm) call call_ground_rdm_procedures(dmqmc_estimates, beta_cycle, dmqmc_in%rdm)
            ! Calculate and output new weights based on the psip distirubtion in
            ! the previous loop.
            if (dmqmc_in%find_weights) call output_and_alter_weights(dmqmc_in, sys%max_number_excitations, &
                                                                     dmqmc_estimates%excit_dist, weighted_sampling)

        end do outer_loop

        if (parent) write (6,'()')
        call write_bloom_report(bloom_stats)
        call load_balancing_report(qs%psip_list%nparticles, qs%psip_list%nstates, qmc_in%use_mpi_barriers, &
                                   qs%spawn_store%spawn%mpi_time)
        call write_memcheck_report(qs%spawn_store%spawn)

        if (soft_exit .or. error) then
            qs%mc_cycles_done = qs%mc_cycles_done + qmc_in%ncycles*ireport
        else
            qs%mc_cycles_done = qs%mc_cycles_done + qmc_in%ncycles*qmc_in%nreport
        end if

        if (restart_in%write_restart) then
            call dump_restart_hdf5(ri, qs, qs%mc_cycles_done, qs%psip_list%tot_nparticles, .false.)
            if (parent) write (6,'()')
        end if

        call copy_sys_spin_info(sys_copy, sys)
        call dealloc_det_info_t(cdet1, .false.)
        call dealloc_det_info_t(cdet2, .false.)

    end subroutine do_dmqmc

    subroutine init_dmqmc_beta_loop(rng, qmc_in, dmqmc_in, dmqmc_estimates, qs, beta_cycle, nstates_active, &
                                    nparticles, spawn, accumulated_probs)

        ! Initialise/reset DMQMC data for a new run over the temperature range.

        ! In/Out:
        !    rng: random number generator.
        !    spawn: spawn_t object.  Reset on exit.
        !    dmqmc_estimates: type containing dmqmc estimates.
        !    qs: state of QMC calculation. Shift is reset on exit.
        ! In:
        !    qmc_in: input options relating to QMC calculations.
        !    beta_cycle: The index of the beta loop about to be started.
        !    dmqmc_in: input options for DMQMC.
        ! Out:
        !    nparticles: number of particles in each space/of each type on
        !       processor.  Set to 0.
        !    nstates_active: number of occupied density matrix elements on
        !       processor.  Set to 0.
        !    accumulated_probs:  This holds the factors by which the populations
        !        on each excitation level (from 0 to max_number_excitations) are
        !        reduced, relative to DMQMC without any importance sampling.

        use dSFMT_interface, only: dSFMT_t, dSFMT_init
        use parallel
        use qmc_data, only: qmc_in_t, qmc_state_t
        use dmqmc_data, only: dmqmc_in_t, dmqmc_estimates_t
        use spawn_data, only: spawn_t
        use utils, only: int_fmt

        type(dSFMT_t), intent(inout) :: rng
        type(spawn_t), intent(inout) :: spawn
        type(qmc_in_t), intent(in) :: qmc_in
        type(dmqmc_in_t), intent(in) :: dmqmc_in
        type(dmqmc_estimates_t), intent(inout) :: dmqmc_estimates
        integer, intent(in) :: beta_cycle
        type(qmc_state_t), intent(inout) :: qs
        integer, intent(out) :: nstates_active
        real(p), intent(out) :: nparticles(:)
        real(p), intent(out) :: accumulated_probs(:)
        integer :: new_seed

        ! Reset the current position in the spawning array to be the slot
        ! preceding the first slot.
        spawn%head = spawn%head_start

        ! Set all quantities back to their starting values.
        nstates_active = 0
        qs%shift = qmc_in%initial_shift
        nparticles = 0.0_dp
        if (allocated(dmqmc_estimates%ground_rdm%rdm)) dmqmc_estimates%ground_rdm%rdm = 0.0_p
        if (dmqmc_in%vary_weights) accumulated_probs = 1.0_p
        if (dmqmc_in%find_weights) dmqmc_estimates%excit_dist = 0.0_p

        new_seed = qmc_in%seed+iproc+(beta_cycle-1)*nprocs

        if (beta_cycle /= 1 .and. parent) then
            write (6,'(a32,'//int_fmt(beta_cycle,1)//')') " # Resetting beta... Beta loop =", beta_cycle
            write (6,'(a52,'//int_fmt(new_seed,1)//',a1)') " # Resetting random number generator with a seed of:", new_seed, "."
        end if

        ! Reset the random number generator with new_seed = old_seed +
        ! nprocs (each beta loop)
        call dSFMT_init(new_seed, 50000, rng)

    end subroutine init_dmqmc_beta_loop

    subroutine init_dmqmc_report_loop(bloom_stats, dmqmc_estimates, rspawn)

        ! Initialise a report loop (basically zero quantities accumulated over
        ! a report loop).

        ! In/Out:
        !    bloom_stats: type containing information regarding blooming events.
        !    dmqmc_estimates: type containing estimates of observables.
        ! Out:
        !    rspawn: spawning rate.

        use bloom_handler, only: bloom_stats_t, bloom_stats_init_report_loop
        use dmqmc_data, only: dmqmc_estimates_t

        type(bloom_stats_t), intent(inout) :: bloom_stats
        type(dmqmc_estimates_t), intent(inout) :: dmqmc_estimates
        real(p), intent(out) :: rspawn

        call bloom_stats_init_report_loop(bloom_stats)

        rspawn = 0.0_p

        if (allocated(dmqmc_estimates%excit_dist)) dmqmc_estimates%excit_dist = 0.0_p
        dmqmc_estimates%trace = 0.0_p
        dmqmc_estimates%numerators = 0.0_p

    end subroutine init_dmqmc_report_loop

end module dmqmc

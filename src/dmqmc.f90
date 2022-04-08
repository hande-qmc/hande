module dmqmc

! Main loop for performing DMQMC calculations.
! [todo] - Top level comments explaining various dmqmc algorithms.

use qmc_io
use proc_pointers
implicit none

contains

    subroutine do_dmqmc(sys, qmc_in, dmqmc_in, dmqmc_estimates, restart_in, load_bal_in, reference_in, qs, sampling_probs, &
                        qmc_state_restart)

        ! Run DMQMC calculation. We run from a beta=0 to a value of beta
        ! specified by the user and then repeat this main loop beta_loops
        ! times, to accumulate statistics for each value for beta.

        ! In:
        !    sys: system being studied.
        !    qmc_in: input options relating to QMC methods.
        !    dmqmc_in: input options relating to DMQMC.
        !    restart_in: input options for HDF5 restart files.
        !    load_bal_in: input options for load balancing.
        !    reference_in: current reference determinant.  If not set (ie
        !       components allocated) then a best guess is made based upon the
        !       desired spin/symmetry.
        ! In/Out:
        !    qmc_state_restart (optional): if present, restart from a previous fciqmc calculation.
        !       Deallocated on exit.
        !    dmqmc_estimates: type containing all DMQMC estimates.
        ! Out:
        !    qs: qmc_state for use if restarting the calculation
        !    sampling_probs: sampling probabilities found via the output_and_alter_weights algorithm.  Only allocated and set if
        !       dmqmc_in%find_weights is true.

        use parallel
        use json_out
        use checking, only: check_allocate
        use annihilation, only: direct_annihilation
        use bit_utils, only: count_set_bits
        use bloom_handler, only: init_bloom_stats_t, bloom_mode_fixedn, bloom_stats_warning, &
                                 bloom_stats_t, accumulate_bloom_stats, write_bloom_report
        use death, only: stochastic_death
        use determinants, only: alloc_det_info_t, dealloc_det_info_t, decode_det, update_sys_spin_info
        use determinant_data, only: det_info_t
        use dmqmc_estimators
        use dmqmc_procedures
        use dmqmc_initialise_dm, only: create_initial_density_matrix
        use excitations, only: excit_t, connection_exists
        use qmc, only: init_qmc, init_proc_pointers
        use qmc_common
        use restart_hdf5, only: restart_info_t, dump_restart_hdf5, init_restart_info_t
        use system
        use dSFMT_interface, only: dSFMT_t, dSFMT_end
        use qmc_data, only: qmc_in_t, restart_in_t, load_bal_in_t, annihilation_flags_t, qmc_state_t, &
                            qmc_in_t_json, restart_in_t_json, load_bal_in_t_json
        use logging, only: logging_t
        use reference_determinant, only: reference_t, reference_t_json
        use dmqmc_data, only: dmqmc_in_t, dmqmc_estimates_t, dmqmc_weighted_sampling_t, dmqmc_in_t_json, ipdmqmc_in_t_json, &
                              rdm_in_t_json, operators_in_t_json, hartree_fock_dm
        use qmc_io, only: write_dmqmc_report_header, write_dmqmc_report
        use check_input, only: check_qmc_opts, check_dmqmc_opts
        use spawn_data, only: write_memcheck_report, dealloc_spawn_t
        use idmqmc, only: set_parent_flag_dmqmc
        use hash_table, only: free_hash_table
        use chem_pot, only: find_chem_pot
        use errors, only: stop_all
        use calc, only: doing_dmqmc_calc, dmqmc_ref_proj_energy
        use state_histograms

        type(sys_t), intent(inout) :: sys
        type(qmc_in_t), intent(in) :: qmc_in
        type(dmqmc_in_t), intent(inout) :: dmqmc_in
        type(dmqmc_estimates_t), intent(inout) :: dmqmc_estimates
        type(restart_in_t), intent(in) :: restart_in
        type(load_bal_in_t), intent(in) :: load_bal_in
        type(reference_t), intent(in) :: reference_in
        type(qmc_state_t), intent(out), target :: qs
        type(qmc_state_t), intent(inout), optional :: qmc_state_restart
        real(p), intent(out), allocatable :: sampling_probs(:)
        real(p) :: mu(1:2)

        integer :: idet, ireport, icycle, iteration, ireplica, ierr
        integer :: beta_cycle, nreport, prop_switch_report
        integer :: unused_int_1 = -1, unused_int_2 = 0
        integer(int_64) :: init_tot_nparticles
        real(dp), allocatable :: tot_nparticles_old(:)
        real(p), allocatable :: real_population(:)
        integer(int_64) :: nattempts
        integer :: nel_temp
        type(det_info_t) :: cdet1, cdet2
        integer(int_p) :: ndeath
        integer :: nspawn_events
        logical :: imag, calc_ref_proj_energy, piecewise_propagation
        logical :: soft_exit, write_restart_shift, update_tau
        logical :: error, rdm_error, attempt_spawning, restarting
        real :: t1, t2
        real(p) :: energy_shift
        type(dSFMT_t) :: rng
        type(bloom_stats_t) :: bloom_stats
        type(annihilation_flags_t) :: annihilation_flags
        type(dmqmc_weighted_sampling_t) :: weighted_sampling
        type(restart_info_t) :: ri
        type(sys_t) :: sys_copy
        type(json_out_t) :: js
        type(qmc_in_t) :: qmc_in_loc
        type(dmqmc_in_t) :: dmqmc_in_loc
        type(state_histograms_data_t) :: state_hist
        character(36) :: uuid_restart

        type(logging_t) :: logging_info
        integer :: iunit, restart_version_restart

        iunit = 6

        if (parent) then
            write (iunit,'(1X,"DMQMC")')
            write (iunit,'(1X,"-----",/)')
        end if

        if (parent) then
            restarting = present(qmc_state_restart) .or. restart_in%read_restart
            call check_qmc_opts(qmc_in, sys, .not. present(qmc_state_restart), restarting)
            call check_dmqmc_opts(sys, dmqmc_in, qmc_in)
        end if

        ! Initialise data.
        call init_qmc(sys, qmc_in, restart_in, load_bal_in, reference_in, 6, annihilation_flags, qs, uuid_restart, &
                      restart_version_restart, dmqmc_in=dmqmc_in, qmc_state_restart=qmc_state_restart, state_hist=state_hist)

        allocate(tot_nparticles_old(qs%psip_list%nspaces), stat=ierr)
        call check_allocate('tot_nparticles_old', qs%psip_list%nspaces, ierr)
        allocate(real_population(qs%psip_list%nspaces), stat=ierr)
        call check_allocate('real_population', qs%psip_list%nspaces, ierr)

        ! It is convenient to setup calculations using a final beta and
        ! target beta instead of calculating the nreports from tau and mc_cycles.
        ! So we do the nreport calculation here instead.
        ! Additionally we use piecewise_propagation and prop_switch_report
        ! when running piecewise ip-dmqmc, this is the report at which the propagator changes.
        ! Otherwise, if neither beta values are defined we default to
        ! the standard behavior of using nreport from the input.
        piecewise_propagation = .false.
        prop_switch_report = -1
        if (dmqmc_in%final_beta > 0.0_p) then
            nreport = int(ceiling(dmqmc_in%final_beta/(qmc_in%ncycles*qmc_in%tau)))
        else
            nreport = qmc_in%nreport
        end if

        if (dmqmc_in%ipdmqmc) then
            if (dmqmc_in%final_beta < dmqmc_in%target_beta) then
                nreport = int(ceiling(dmqmc_in%target_beta/(qmc_in%ncycles*qmc_in%tau)))
            else
                piecewise_propagation = .true.
                prop_switch_report = int(ceiling(dmqmc_in%target_beta/(qmc_in%ncycles*qmc_in%tau))) + 1
            end if
        end if

        ! When we accumulate data throughout a run, we are actually accumulating
        ! results from the psips distribution from the previous iteration.
        ! For example, in the first iteration, the trace calculated will be that
        ! of the initial distribution, which corresponds to beta=0. Hence, in the
        ! output we subtract one from the iteration number, and run for one more
        ! report loop, asimplemented in the line of code below.
        nreport = nreport + 1

        ! Initialise all the required arrays, ie to store thermal quantities,
        ! and to initalise reduced density matrix quantities if necessary.
        call init_dmqmc(sys, qmc_in, dmqmc_in, qs%psip_list%nspaces, qs, dmqmc_estimates, weighted_sampling)
        ! Determine the chemical potential if doing the ip-dmqmc algorithm with 
        ! grand canonical initialisation enabled.
        if (dmqmc_in%ipdmqmc .and. dmqmc_in%grand_canonical_initialisation) then
            mu = find_chem_pot(sys, qs%target_beta)
            if (dmqmc_in%initial_matrix == hartree_fock_dm) then
                energy_shift = energy_diff_ptr(sys, qs%ref%occ_list0)
            else
                energy_shift = 0.0_p
            end if
        else
            mu = 0.0_p
        end if

        if (parent) then
            call json_object_init(js, tag=.true.)
            call sys_t_json(js, sys)
            ! The default values of pattempt_* are not in qmc_in
            qmc_in_loc = qmc_in
            qmc_in_loc%pattempt_single = qs%excit_gen_data%pattempt_single
            qmc_in_loc%pattempt_double = qs%excit_gen_data%pattempt_double
            qmc_in_loc%shift_damping = qs%shift_damping
            qmc_in_loc%pattempt_parallel = qs%excit_gen_data%pattempt_parallel
            ! Not actually used in DMQMC!
            qmc_in_loc%quasi_newton_threshold = qs%propagator%quasi_newton_threshold
            qmc_in_loc%quasi_newton_value = qs%propagator%quasi_newton_value
            qmc_in_loc%quasi_newton_pop_control = qs%propagator%quasi_newton_pop_control
            ! Initialize the harmonic shift damping algorithm parameters
            qmc_in_loc%shift_harmonic_forcing = qs%shift_harmonic_forcing
            if (dmqmc_in%walker_scale_factor > 0.0_p) then
                ! If we are scaling the initial walker population then
                ! update the reported initial population to be accurate.
                qs%target_particles = qmc_in%target_particles * real(int(dmqmc_in%walker_scale_factor, p), p)
                qmc_in_loc%target_particles = qs%target_particles
                qmc_in_loc%D0_population = nint(qmc_in%D0_population, int_64) * real(int(dmqmc_in%walker_scale_factor, p), p)
            end if
            call qmc_in_t_json(js, qmc_in_loc)
            call dmqmc_in_t_json(js, dmqmc_in)
            dmqmc_in_loc = dmqmc_in
            dmqmc_in_loc%chem_pot = mu
            call ipdmqmc_in_t_json(js, dmqmc_in_loc)
            call rdm_in_t_json(js, dmqmc_in%rdm)
            call operators_in_t_json(js, dmqmc_in)
            call restart_in_t_json(js, restart_in, uuid_restart)
            call load_bal_in_t_json(js, load_bal_in)
            call reference_t_json(js, qs%ref, sys, terminal=.true.)
            call json_object_end(js, terminal=.true., tag=.true.)
            write (js%io, '()')
        end if

        ! Allocate det_info_t components. We need two cdet objects for each 'end'
        ! which may be spawned from in the DMQMC algorithm.
        call alloc_det_info_t(sys, cdet1, .false.)
        call alloc_det_info_t(sys, cdet2, .false.)

        ! Initialise bloom_stats components to the following parameters.
        call init_bloom_stats_t(bloom_stats, mode=bloom_mode_fixedn, encoding_factor=qs%psip_list%pop_real_factor)

        ! Main DMQMC loop.
        if (parent) call write_dmqmc_report_header(qs%psip_list%nspaces, dmqmc_in, sys%max_number_excitations, &
                                                  &dmqmc_estimates, sys%read_in%comp)
        ! Initialise timer.
        call cpu_time(t1)

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
        calc_ref_proj_energy = .false.
        if (doing_dmqmc_calc(dmqmc_ref_proj_energy)) calc_ref_proj_energy = .true.

        outer_loop: do beta_cycle = 1, dmqmc_in%beta_loops

            call init_dmqmc_beta_loop(sys, reference_in, rng, qmc_in, dmqmc_in, dmqmc_estimates, qs, beta_cycle, &
                                      qs%psip_list%nstates, qs%psip_list%nparticles, qs%spawn_store%spawn, &
                                      weighted_sampling%probs, annihilation_flags, piecewise_propagation, state_hist)

            ! Distribute psips uniformly along the diagonal of the density matrix.
            call create_initial_density_matrix(rng, sys, qmc_in, dmqmc_in, qs, annihilation_flags, &
                                               init_tot_nparticles, qs%psip_list, qs%spawn_store%spawn, &
                                               mu, energy_shift)

            ! Allow the shift to vary from the very start of the beta loop, if
            ! shift_harmonic_forcing is not equal to zero. 
            if (qmc_in_loc%shift_harmonic_forcing .ne. 0.00_p) then
                qs%vary_shift = .true.
            else
            ! Allow the shift to vary from the very start of the beta loop, if
            ! this condition is met.
                qs%vary_shift = qs%psip_list%tot_nparticles >= qs%target_particles 
            end if

            ! DMQMC quasi-newton not functional, so we artificially set this value to 0 which should not affect non-QN calcs.
            qs%estimators%proj_energy_old = 0_p

            do ireport = 1, nreport

                call init_dmqmc_report_loop(calc_ref_proj_energy, dmqmc_in%calc_excit_dist, dmqmc_in%calc_mom_dist, &
                                            dmqmc_in%calc_struc_fac, bloom_stats, dmqmc_estimates, qs%spawn_store%rspawn)
                tot_nparticles_old = qs%psip_list%tot_nparticles

                if (ireport .eq. 1) then
                    ! If we have just initalized the density matrix, print out
                    ! the number of diagonal elemenets that were occupied.
                    call dmqmc_estimate_comms(dmqmc_in, error, nspawn_events, sys%max_number_excitations, qmc_in%ncycles, &
                                              qs%psip_list, qs, weighted_sampling%probs_old, dmqmc_estimates)
                    if (parent) write (iunit,'(1X,"# Initial diagonal density matrix element(s):",1X,I0)') qs%estimators%tot_nstates
                end if

                ! If we are propagating past IP-DMQMC target beta.
                ! Change over the appropriate parameters and pointers to continue on with sampling.
                if (piecewise_propagation .and. prop_switch_report == ireport) then
                    call propagator_change(sys, qs, dmqmc_in, annihilation_flags, iunit)
                    call init_proc_pointers(sys, qmc_in, reference_in, iunit, dmqmc_in)
                end if

                do icycle = 1, qmc_in%ncycles

                    call init_mc_cycle(qs%psip_list, qs%spawn_store%spawn, nattempts, ndeath, complx=sys%read_in%comp)

                    iteration = (ireport-1)*qmc_in%ncycles + icycle

                    ! Store (beta-tau)/2 for use in symmetric ip-dmqmc spawning probabilities.
                    if (dmqmc_in%ipdmqmc .and. dmqmc_in%symmetric) &
                            & weighted_sampling%probs(sys%max_number_excitations+1) = &
                            & 0.5*(qs%target_beta-(iteration-1)*qs%dmqmc_factor*qs%tau)

                    do idet = 1, qs%psip_list%nstates ! loop over walkers/dets

                        ! f points to the bitstring that is spawning, f2 to the
                        ! other bit string.
                        cdet1%f => qs%psip_list%states(:sys%basis%tot_string_len,idet)
                        cdet1%f2 => qs%psip_list%states((sys%basis%tot_string_len+1):(2*sys%basis%tot_string_len),idet)
                        cdet1%data => qs%psip_list%dat(:,idet)
                        cdet2%f => qs%psip_list%states((sys%basis%tot_string_len+1):(2*sys%basis%tot_string_len),idet)
                        cdet2%f2 => qs%psip_list%states(:sys%basis%tot_string_len,idet)
                        ! Decode and store the the relevant information for
                        ! both bitstrings. Both of these bitstrings are required
                        ! to refer to the correct element in the density matrix.
                        call decoder_ptr(sys, cdet1%f, cdet1)
                        call decoder_ptr(sys, cdet2%f, cdet2)
                        ! Quasi-Newton not yet implemented -- see also stochastic_death call.
                        cdet1%fock_sum = 1.0_p
                        cdet2%fock_sum = 1.0_p

                        ! If using multiple symmetry sectors then find the
                        ! symmetry labels of this particular det.
                        if (dmqmc_in%all_spin_sectors) call update_sys_spin_info(cdet1, sys)

                        ! Extract the real signs from the encoded signs.
                        real_population = real(qs%psip_list%pops(:,idet),p)/qs%psip_list%pop_real_factor
                        ! Is this density matrix element an initiator?
                        if (qmc_in%initiator_approx) then
                            call set_parent_flag_dmqmc(real_population(1), qmc_in%initiator_pop, cdet1%f, cdet1%f2, &
                                                       dmqmc_in%initiator_level, cdet1%initiator_flag)
                            cdet2%initiator_flag = cdet1%initiator_flag
                        end if

                        ! Call wrapper function which calls routines to update
                        ! all estimators being calculated, and also always
                        ! updates the trace separately.
                        ! Note DMQMC averages over multiple loops over
                        ! temperature/imaginary time so only get data from one
                        ! temperature value per ncycles.
                        if (icycle == 1) then
                            ! If we are accumulating state histogram data, update
                            ! the walker/det count from this state.
                            if (qmc_in%state_histograms) then
                                call update_state_histogram(qs, cdet1%f, cdet1%f2, real_population(1), state_hist, &
                                                           icycle, ireport, final_report=ireport == nreport)
                            end if

                            call update_dmqmc_estimators(sys, dmqmc_in, idet, iteration, cdet1, qs%ref%H00, &
                                                         qs%psip_list, dmqmc_estimates, weighted_sampling, rdm_error, &
                                                         qs%ref%f0)
                        end if

                        ! Only attempt spawning if a valid connection exists.
                        attempt_spawning = connection_exists(sys)

                        do ireplica = 1, qs%psip_list%nspaces

                            ! If the system is complex, ireplica will be:
                            ! { replica1_real, replica1_imag, replica2_real, replica2_imag }
                            imag = sys%read_in%comp .and. mod(ireplica, 2) == 0

                            ! Only attempt spawning if a valid excitation exists.
                            if (attempt_spawning) &
                                call do_dmqmc_spawning_attempt(rng, qs%spawn_store%spawn, bloom_stats, sys, qs, &
                                                              &decide_nattempts(rng, real_population(ireplica)), &
                                                              &cdet1, cdet2, qs%psip_list%pops(ireplica, idet), &
                                                              &imag, dmqmc_in%symmetric, ireplica, &
                                                              &weighted_sampling%probs, logging_info)

                            ! Clone or die.
                            ! We have contributions to the clone/death step from
                            ! both ends of the current walker. We do both of
                            ! these at once by using qs%psip_list%dat(:,idet) which,
                            ! when running a DMQMC algorithm, stores the average
                            ! of the two diagonal elements corresponding to the
                            ! two indicies of the density matrix.
                            call stochastic_death(rng, qs, cdet1%fock_sum, qs%psip_list%dat(1, idet), &
                                                  qs%shift(ireplica), qs%estimators(1)%proj_energy_old, logging_info, &
                                                  qs%psip_list%pops(ireplica,idet), qs%psip_list%nparticles(ireplica), ndeath)
                        end do
                    end do

                    ! Perform the annihilation step where the spawned walker
                    ! list is merged with the main walker list, and walkers of
                    ! opposite sign on the same sites are annihilated.
                    call direct_annihilation(sys, rng, qs%ref, annihilation_flags, qs%psip_list, qs%spawn_store%spawn, &
                                             nspawn_events)

                    call end_mc_cycle(nspawn_events, ndeath, qs%psip_list%pop_real_factor, nattempts, qs%spawn_store%rspawn)

                    ! If doing importance sampling *and* varying the weights of
                    ! the trial function, call a routine to update these weights
                    ! and alter the number of psips on each excitation level
                    ! accordingly.
                    if (dmqmc_in%vary_weights .and. iteration <= dmqmc_in%finish_varying_weights) &
                        call update_sampling_weights(rng, sys%basis, qmc_in, qs%psip_list, &
                                                     sys%max_number_excitations, weighted_sampling)

                end do

                ! Sum all quantities being considered across all MPI processes.
                error = qs%spawn_store%spawn%error .or. qs%psip_list%error .or. rdm_error
                call dmqmc_estimate_comms(dmqmc_in, error, nspawn_events, sys%max_number_excitations, qmc_in%ncycles, &
                                          &qs%psip_list, qs, weighted_sampling%probs_old, dmqmc_estimates)
                if (error) exit outer_loop

                call update_shift_dmqmc(qmc_in, qs, qs%psip_list%tot_nparticles, tot_nparticles_old)

                ! Forcibly disable update_tau as need to average over multiple loops over beta
                ! and hence want to use the same timestep throughout.
                update_tau = .false.
                call end_report_loop(iunit, qmc_in, iteration, update_tau, qs, tot_nparticles_old, &
                                     nspawn_events, unused_int_1, unused_int_2, soft_exit, &
                                     load_bal_in, .false., bloom_stats=bloom_stats)

                call cpu_time(t2)

                if (qmc_in%state_histograms) then
                    call comm_and_report_state_histogram(state_hist, ireport, iunit, &
                                                         final_report=ireport == nreport .or. soft_exit)
                end if

                if (parent) then
                    if (bloom_stats%nblooms_curr > 0) call bloom_stats_warning(bloom_stats)
                    call write_dmqmc_report(sys, qmc_in, qs, ireport, tot_nparticles_old, t2-t1, .false., &
                                            dmqmc_in, dmqmc_estimates, sys%read_in%comp)
                end if

                ! Update the time for the start of the next iteration.
                t1 = t2

                if (soft_exit) exit outer_loop

            end do

            ! Calculate and output all requested estimators based on the reduced
            ! density matrix. This is for ground-state RDMs only.
            if (dmqmc_in%rdm%calc_ground_rdm) call call_ground_rdm_procedures(dmqmc_estimates, beta_cycle, dmqmc_in%rdm)
            ! Calculate and output new weights based on the psip distriubtion in
            ! the previous loop.
            if (dmqmc_in%find_weights) call output_and_alter_weights(sys, dmqmc_in, dmqmc_estimates%excit_dist, weighted_sampling)

        end do outer_loop

        if (parent) write (iunit,'()')
        call write_bloom_report(bloom_stats)
        call load_balancing_report(qs%psip_list%nparticles, qs%psip_list%nstates, qmc_in%use_mpi_barriers, &
                                   qs%spawn_store%spawn%mpi_time)
        call write_memcheck_report(qs%spawn_store%spawn)

        if (soft_exit .or. error) then
            qs%mc_cycles_done = qs%mc_cycles_done + qmc_in%ncycles*ireport
        else
            qs%mc_cycles_done = qs%mc_cycles_done + qmc_in%ncycles*nreport
        end if

        if (restart_in%write_restart) then
            call dump_restart_hdf5(qs, qs%mc_cycles_done, qs%psip_list%tot_nparticles, sys%basis%nbasis, .false., &
                                    sys%basis%info_string_len)
            if (parent) write (iunit,'()')
        end if

        if (dmqmc_in%find_weights) then
            allocate(sampling_probs(size(weighted_sampling%sampling_probs)), stat=ierr)
            call check_allocate('sampling_probs', size(weighted_sampling%sampling_probs), ierr)
            sampling_probs = weighted_sampling%sampling_probs
        end if

        call copy_sys_spin_info(sys_copy, sys)
        call dealloc_det_info_t(cdet1, .false.)
        call dealloc_det_info_t(cdet2, .false.)

        if (qmc_in%state_histograms) call deallocate_histogram_t(state_hist)

        call dSFMT_end(rng)

        if (allocated(dmqmc_estimates%inst_rdm%spawn)) then
            do idet = 1, size(dmqmc_estimates%subsys_info)
                call free_hash_table(dmqmc_estimates%inst_rdm%spawn(idet)%ht)
                call dealloc_spawn_t(dmqmc_estimates%inst_rdm%spawn(idet)%spawn)
            end do
        end if

    end subroutine do_dmqmc

    subroutine do_dmqmc_spawning_attempt(rng, spawn, bloom_stats, sys, qs, nattempts_current_det, cdet1, cdet2, pop, &
                                        &imag_parent, symmetric, ireplica, weighted_probs, logging_info)

        ! Perform spawning from a given matrix element in a given space.

        ! In:
        !   sys: information on system under consideration.
        !   logging_info: information on current logging settings.
        !   weighted_sampling: stores sampling weights information when using
        !       importance sampling.
        !   nattempts_current_det: total number of spawning attempts
        !       to make on this matrix element.
        !   ireplica: replica (and its re/im part) currently under
        !       consideration.
        !   imag_parent: true if spawning from psips within an imaginary
        !       space.
        !   pop: population of given density matrix element in given space.
        !   symmetric: if we should spawn into rows and columns
        !       symmetrically.
        ! In/Out:
        !   rng: random number generator.
        !   qs: qmc_state_t derived type with information on
        !       current calculation. In/Out for compatibility with some
        !       excitation generator features (not currently used in DMQMC).
        !   bloom_stats: information on blooms during calculation.
        !   spawn: stored information on spawning.
        !   cdet[1/2]: determinants corresponding to the row and column of 
        !       density matrix element spawning is originating from.

        use dSFMT_interface, only: dSFMT_t
        use system, only: sys_t
        use logging, only: logging_t
        use qmc_data, only: qmc_state_t
        use determinants, only: det_info_t
        use bloom_handler, only: bloom_stats_t, accumulate_bloom_stats
        use excitations, only: excit_t
        use spawn_data, only: spawn_t

        type(sys_t), intent(in) :: sys
        type(qmc_state_t), intent(inout) :: qs
        type(logging_t), intent(in) :: logging_info
        type(det_info_t), intent(inout) :: cdet1, cdet2
        integer, intent(in) :: nattempts_current_det, ireplica
        logical, intent(in) :: imag_parent, symmetric
        integer(int_p), intent(in) :: pop

        type(dSFMT_t), intent(inout) :: rng
        type(spawn_t), intent(inout) :: spawn
        type(bloom_stats_t), intent(inout) :: bloom_stats
        real(p), intent(in), allocatable :: weighted_probs(:)

        type(excit_t) :: connection
        integer(int_p) :: nspawned, nspawned_im
        integer :: iparticle, replica_real, replica_imag, spawning_end
        integer(int_p) :: scratch

        ! First, determine the particle types possibly created by spawning.
        ! NB this implicitly assumes the replicas are ordered (real, imaginary)

        if (imag_parent) then
            replica_imag = ireplica
            replica_real = ireplica - 1
        else
            replica_real = ireplica
            replica_imag = ireplica + 1
        end if

        do iparticle = 1, nattempts_current_det
            ! When using importance sampling in DMQMC we only spawn from 
            ! one end of a density matrix element.

            ! Spawn from the first end.
            spawning_end = 1
            ! Attempt to spawn.
            call spawner_ptr(rng, sys, qs, qs%spawn_store%spawn%cutoff, qs%psip_list%pop_real_factor, cdet1, pop, &
                            &gen_excit_ptr, weighted_probs, logging_info, nspawned, nspawned_im, connection)
            if (imag_parent) then
                ! If has imaginary parent have to factor into resulting signs/components
                scratch = nspawned_im
                nspawned_im = nspawned
                nspawned = -scratch
            end if
            ! Spawn if attempt was successful.
            if (nspawned /= 0_int_p) then ! real component
                call create_spawned_particle_dm_ptr(sys%basis, qs%ref, cdet1, connection, nspawned, &
                                                   &spawning_end, replica_real, spawn)
                call accumulate_bloom_stats(bloom_stats, nspawned)
            end if
            if (nspawned_im /= 0_int_p) then
                call create_spawned_particle_dm_ptr(sys%basis, qs%ref, cdet1, connection, nspawned_im, &
                                                   &spawning_end, replica_imag, spawn)
                call accumulate_bloom_stats(bloom_stats, nspawned_im)
            end if

            ! Now attempt to spawn from the second end.
            if (symmetric) then
                spawning_end = 2
                call spawner_ptr(rng, sys, qs, qs%spawn_store%spawn%cutoff, qs%psip_list%pop_real_factor, cdet2, pop, &
                                &gen_excit_ptr, weighted_probs, logging_info, nspawned, nspawned_im, connection)
                if (imag_parent) then
                    scratch = nspawned_im
                    nspawned_im = nspawned
                    nspawned = -scratch
                end if
                if (nspawned /= 0_int_p) then
                    call create_spawned_particle_dm_ptr(sys%basis, qs%ref, cdet2, connection, nspawned, &
                                                       &spawning_end, replica_real, spawn)
                    call accumulate_bloom_stats(bloom_stats, nspawned)
                end if
                if (nspawned_im /= 0_int_p) then
                    call create_spawned_particle_dm_ptr(sys%basis, qs%ref, cdet2, connection, nspawned_im, &
                                                       &spawning_end, replica_imag, spawn)
                    call accumulate_bloom_stats(bloom_stats, nspawned_im)
                end if
            end if
        end do

    end subroutine do_dmqmc_spawning_attempt

    subroutine init_dmqmc_beta_loop(sys, reference_in, rng, qmc_in, dmqmc_in, dmqmc_estimates, qs, beta_cycle, &
                                    nstates_active, nparticles, spawn, accumulated_probs, annihilation_flags, &
                                    piecewise_propagation, state_hist)

        ! Initialise/reset DMQMC data for a new run over the temperature range.

        ! In/Out:
        !    rng: random number generator.
        !    spawn: spawn_t object.  Reset on exit.
        !    dmqmc_estimates: type containing dmqmc estimates.
        !    qs: state of QMC calculation. Shift is reset on exit.
        !    state_hist: A type containing information
        !       for state histograms
        ! In:
        !    qmc_in: input options relating to QMC calculations.
        !    beta_cycle: The index of the beta loop about to be started.
        !    dmqmc_in: input options for DMQMC.
        !    reference_in: current reference determinant.  If not set (ie
        !       components allocated) then a best guess is made based upon the
        !       desired spin/symmetry.
        !    piecewise_propagation: A boolean indiciating if peicewise propagation
        !       is being used. In this case we may need to reset some parameters/pointers.
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
        use system
        use qmc_data, only: qmc_in_t, qmc_state_t, annihilation_flags_t
        use dmqmc_data, only: dmqmc_in_t, dmqmc_estimates_t
        use spawn_data, only: spawn_t
        use utils, only: int_fmt
        use dmqmc_estimators, only: dmqmc_energy_and_trace_propagate
        use dmqmc_procedures, only: propagator_restore
        use reference_determinant, only: reference_t, reference_t_json
        use qmc, only: init_proc_pointers
        use state_histograms

        type(sys_t), intent(inout) :: sys
        type(dSFMT_t), intent(inout) :: rng
        type(spawn_t), intent(inout) :: spawn
        type(dmqmc_in_t), intent(inout) :: dmqmc_in
        type(dmqmc_estimates_t), intent(inout) :: dmqmc_estimates
        type(qmc_in_t), intent(in) :: qmc_in
        type(reference_t), intent(in) :: reference_in
        integer, intent(in) :: beta_cycle
        type(qmc_state_t), intent(inout) :: qs
        type(state_histograms_data_t), intent(inout) :: state_hist
        integer, intent(out) :: nstates_active
        real(dp), intent(out) :: nparticles(:)
        real(p), intent(out) :: accumulated_probs(:)
        integer :: new_seed, iunit
        type(annihilation_flags_t) :: annihilation_flags
        logical, intent(in) :: piecewise_propagation

        iunit = 6

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
            write (iunit,'(a32,'//int_fmt(beta_cycle,1)//')') " # Resetting beta... Beta loop =", beta_cycle
            write (iunit,'(a52,'//int_fmt(new_seed,1)//',a1)') " # Resetting random number generator with a seed of:", new_seed, "."
        end if

        ! if we are running piecewise ip-dmqmc restore the relevant values 
        ! and pointers back to ip-dmqmc.
        if (beta_cycle /= 1 .and. piecewise_propagation) then
            call propagator_restore(qs, dmqmc_in, annihilation_flags, iunit)
            call init_proc_pointers(sys, qmc_in, reference_in, iunit, dmqmc_in)
        end if

        if (beta_cycle /= 1 .and. qmc_in%state_histograms) then
            state_hist%current_seed = new_seed
        end if

        ! Reset the random number generator with new_seed = old_seed +
        ! nprocs (each beta loop)
        call dSFMT_init(new_seed, 50000, rng)

    end subroutine init_dmqmc_beta_loop

    subroutine init_dmqmc_report_loop(calc_ref_proj_energy, calc_excit_dist, calc_mom_dist, &
                                      calc_struc_fac, bloom_stats, dmqmc_estimates, rspawn)

        ! Initialise a report loop (basically zero quantities accumulated over
        ! a report loop).

        ! In:
        !    calc_ref_proj_energy: true if generating the projected energy
        !        from the reference row (or column) of the density matrix.
        !    calc_excit_dist: true if the excitation distribution is being
        !        calculated at each report loop.
        !    calc_mom_dist: true if the momentum distribution is being calculated.
        ! In/Out:
        !    bloom_stats: type containing information regarding blooming events.
        !    dmqmc_estimates: type containing estimates of observables.
        ! Out:
        !    rspawn: spawning rate.

        use bloom_handler, only: bloom_stats_t, bloom_stats_init_report_loop
        use dmqmc_data, only: dmqmc_estimates_t

        logical, intent(in) :: calc_excit_dist
        logical, intent(in) :: calc_mom_dist
        logical, intent(in) :: calc_struc_fac
        logical, intent(in) :: calc_ref_proj_energy
        type(bloom_stats_t), intent(inout) :: bloom_stats
        type(dmqmc_estimates_t), intent(inout) :: dmqmc_estimates
        real(p), intent(out) :: rspawn

        call bloom_stats_init_report_loop(bloom_stats)

        rspawn = 0.0_p
        if (calc_excit_dist) dmqmc_estimates%excit_dist = 0.0_p
        dmqmc_estimates%trace = 0.0_p
        dmqmc_estimates%numerators = 0.0_p
        if (calc_mom_dist) dmqmc_estimates%mom_dist%f_k = 0.0_p
        if (calc_struc_fac) dmqmc_estimates%struc_fac%f_k = 0.0_p
        if (calc_ref_proj_energy) dmqmc_estimates%ref_trace = 0.0_p
        if (calc_ref_proj_energy) dmqmc_estimates%ref_D0j_particles = 0.0_p

    end subroutine init_dmqmc_report_loop

end module dmqmc

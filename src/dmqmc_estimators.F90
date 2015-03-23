module dmqmc_estimators

use const

implicit none

! indices for corresponding data items in estimator/population/etc buffers (see
! comments below).
enum, bind(c)
    enumerator :: rspawn_ind = 1
    enumerator :: nocc_states_ind
    enumerator :: nspawned_ind
    enumerator :: nparticles_ind
    enumerator :: trace_ind
    enumerator :: operators_ind
    enumerator :: excit_dist_ind
    enumerator :: rdm_trace_ind
    enumerator :: rdm_r2_ind
    enumerator :: final_ind ! ensure this remains the last index.
end enum

contains

    subroutine dmqmc_estimate_comms(dmqmc_in, nspawn_events, max_num_excits, ncycles, psip_list, qs, dmqmc_estimates)

        ! Sum together the contributions to the various DMQMC estimators (and
        ! some other non-physical quantities such as the rate of spawning and
        ! total number of walkers) across all MPI processes.

        ! This is called every report loop in a DMQMC calculation.

        ! In:
        !    dmqmc_in: input options relating to DMQMC.
        !    nspawn_events: The total number of spawning events to this process.
        !    max_num_excits: The maximum excitation level for the system being
        !        studied.
        !    ncycles: the number of monte carlo cycles.
        ! In/Out:
        !    psip_list: particle information.  On output total (ie not
        !        per-processor) quantities are updated.
        !    dmqmc_estimates: type containing dmqmc estimates.

        use checking, only: check_allocate, check_deallocate
        use fciqmc_data, only: num_dmqmc_operators, calc_inst_rdm, nrdms
        use qmc_data, only: particle_t, qmc_state_t
        use parallel
        use dmqmc_data, only: dmqmc_in_t, dmqmc_estimates_t

        type(dmqmc_in_t), intent(in) :: dmqmc_in
        integer, intent(in) :: nspawn_events, max_num_excits, ncycles
        type(particle_t), intent(inout) :: psip_list
        type(qmc_state_t), intent(inout) :: qs
        type(dmqmc_estimates_t), intent(inout) :: dmqmc_estimates

        real(dp), allocatable :: rep_loop_loc(:)
        real(dp), allocatable :: rep_loop_sum(:)
        integer :: nelems(final_ind-1), min_ind(final_ind-1), max_ind(final_ind-1)
        integer :: tot_nelems, i, ierr

        ! If calculating instantaneous RDM estimates then call this routine to
        ! perform annihilation for RDM particles, and calculate RDM estimates
        ! before they are summed over processors. below.
        if (calc_inst_rdm) call communicate_inst_rdms()

        ! How big is each variable to be communicated?
        nelems(rspawn_ind) = 1
        nelems(nocc_states_ind) = 1
        nelems(nspawned_ind) = 1
        nelems(nparticles_ind) = psip_list%nspaces
        nelems(trace_ind) = psip_list%nspaces
        nelems(operators_ind) = num_dmqmc_operators 
        nelems(excit_dist_ind) = max_num_excits + 1
        nelems(rdm_trace_ind) = psip_list%nspaces*nrdms
        nelems(rdm_r2_ind) = nrdms

        ! The total number of elements in the array to be communicated.
        tot_nelems = sum(nelems)

        ! Create min_ind and max_ind, which hold the minimum and maximum
        ! indices that each object takes in rep_loop_loc and rep_loop_sum.
        do i = 1, final_ind-1
            min_ind(i) = sum(nelems(1:i-1)) + 1
            max_ind(i) = sum(nelems(1:i))
        end do

        allocate(rep_loop_loc(1:tot_nelems), stat=ierr)
        call check_allocate('rep_loop_loc', tot_nelems, ierr)
        allocate(rep_loop_sum(1:tot_nelems), stat=ierr)
        call check_allocate('rep_loop_sum', tot_nelems, ierr)

        ! Move the variables to be communicated to rep_loop_loc.
        call local_dmqmc_estimators(dmqmc_in, dmqmc_estimates, rep_loop_loc, min_ind, max_ind, psip_list%nparticles, &
                                    psip_list%nstates, nspawn_events)

#ifdef PARALLEL
        call mpi_allreduce(rep_loop_loc, rep_loop_sum, size(rep_loop_loc), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
        rep_loop_sum = rep_loop_loc
        ierr = 0 ! Prevent warning about unused variable in serial so -Werror can be used.
#endif

        ! Move the communicated quantites to the corresponding variables.
        call communicated_dmqmc_estimators(dmqmc_in, rep_loop_sum, min_ind, max_ind, ncycles, psip_list%tot_nparticles, qs, &
                                           dmqmc_estimates)

        ! Clean up.
        deallocate(rep_loop_loc, stat=ierr)
        call check_deallocate('rep_loop_loc', ierr)
        deallocate(rep_loop_sum, stat=ierr)
        call check_deallocate('rep_loop_sum', ierr)

    end subroutine dmqmc_estimate_comms

    subroutine local_dmqmc_estimators(dmqmc_in, dmqmc_estimates, rep_loop_loc, min_ind, max_ind, nparticles, &
                                      nstates_active, nspawn_events)

        ! Enter processor dependent report loop quantites into array for
        ! efficient sending to other processors.

        ! In:
        !    dmqmc_in: input options for DMQMC.
        !    dmqmc_estimates: type containing dmqmc estimates.
        !    min_ind: Array holding the minimum indices of the various
        !        quantities in rep_loop_sum.
        !    max_ind: Array holding the maximum indices of the various
        !        quantities in rep_loop_sum.
        !    nparticles: number of particles in each space on the current processor.
        !    nstates_active: number of occupied states in the particle lists.
        !    nspawn_events: The total number of spawning events to this process.
        ! Out:
        !    rep_loop_loc: array containing local quantities to be communicated.

        use calc, only: doing_dmqmc_calc, dmqmc_rdm_r2
        use fciqmc_data, only: rspawn
        use fciqmc_data, only: trace
        use fciqmc_data, only: excit_dist
        use fciqmc_data, only: trace, rdm_traces, renyi_2
        use fciqmc_data, only: calc_inst_rdm
        use dmqmc_data, only: dmqmc_in_t, dmqmc_estimates_t

        type(dmqmc_in_t), intent(in) :: dmqmc_in
        type(dmqmc_estimates_t), intent(in) :: dmqmc_estimates
        integer, intent(in) :: min_ind(:), max_ind(:)
        integer, intent(in) :: nstates_active
        real(p), intent(in) :: nparticles(:)
        integer, intent(in) :: nspawn_events
        real(dp), intent(out) :: rep_loop_loc(:)
        integer :: nrdms(1)

        rep_loop_loc = 0.0_dp

        rep_loop_loc(rspawn_ind) = rspawn
        rep_loop_loc(nocc_states_ind) = nstates_active
        rep_loop_loc(nspawned_ind) = nspawn_events
        rep_loop_loc(min_ind(nparticles_ind):max_ind(nparticles_ind)) = nparticles
        rep_loop_loc(min_ind(trace_ind):max_ind(trace_ind)) = dmqmc_estimates%trace
        rep_loop_loc(min_ind(operators_ind):max_ind(operators_ind)) = dmqmc_estimates%numerators
        if (dmqmc_in%calc_excit_dist) then
            rep_loop_loc(min_ind(excit_dist_ind):max_ind(excit_dist_ind)) = dmqmc_estimates%excit_dist
        end if
        if (calc_inst_rdm) then
            ! Reshape this 2d array into a 1d array to add it to rep_loop_loc.
            nrdms = max_ind(rdm_trace_ind) - min_ind(rdm_trace_ind) + 1
            rep_loop_loc(min_ind(rdm_trace_ind):max_ind(rdm_trace_ind)) = reshape(rdm_traces, nrdms)
        end if
        if (doing_dmqmc_calc(dmqmc_rdm_r2)) then
            rep_loop_loc(min_ind(rdm_r2_ind):max_ind(rdm_r2_ind)) = renyi_2
        end if

    end subroutine local_dmqmc_estimators

    subroutine communicated_dmqmc_estimators(dmqmc_in, rep_loop_sum, min_ind, max_ind, ncycles, tot_nparticles, qs, &
                                             dmqmc_estimates)

        ! Update report loop quantites with information received from other
        ! processors.

        ! In:
        !    dmqmc_in: input options relating to DMQMC.
        !    rep_loop_sum: array containing quantites which have been
        !        summed over all processors, and are to be moved to their
        !         corresponding report loop quantities.
        !    min_ind: Array holding the minimum indices of the various
        !        quantities in rep_loop_sum.
        !    max_ind: Array holding the maximum indices of the various
        !        quantities in rep_loop_sum.
        !    ncycles: number of monte carlo cycles.
        ! Out:
        !    dmqmc_estimates: type containing dmqmc_estimates.
        !    tot_nparticles: total number of particles of each type across all
        !       processors.

        use calc, only: doing_dmqmc_calc, dmqmc_rdm_r2
        use fciqmc_data, only: rspawn, nrdms
        use fciqmc_data, only: rdm_traces, renyi_2, calc_inst_rdm
        use dmqmc_data, only: dmqmc_in_t, dmqmc_estimates_t
        use parallel, only: nprocs
        use qmc_data, only: qmc_state_t

        type(dmqmc_in_t), intent(in) :: dmqmc_in
        real(dp), intent(in) :: rep_loop_sum(:)
        integer, intent(in) :: min_ind(:), max_ind(:), ncycles
        real(p), intent(out) :: tot_nparticles(:)
        type(qmc_state_t), intent(inout) :: qs
        type(dmqmc_estimates_t), intent(out) :: dmqmc_estimates

        rspawn = rep_loop_sum(rspawn_ind)
        qs%estimators%tot_nstates = rep_loop_sum(nocc_states_ind)
        qs%estimators%tot_nspawn_events = rep_loop_sum(nspawned_ind)
        tot_nparticles = rep_loop_sum(min_ind(nparticles_ind):max_ind(nparticles_ind))
        dmqmc_estimates%trace = rep_loop_sum(min_ind(trace_ind):max_ind(trace_ind))
        dmqmc_estimates%numerators = rep_loop_sum(min_ind(operators_ind):max_ind(operators_ind))
        if (dmqmc_in%calc_excit_dist) then
            dmqmc_estimates%excit_dist = rep_loop_sum(min_ind(excit_dist_ind):max_ind(excit_dist_ind))
        end if
        if (calc_inst_rdm) then
            rdm_traces = reshape(rep_loop_sum(min_ind(rdm_trace_ind):max_ind(rdm_trace_ind)), shape(rdm_traces))
        end if
        if (doing_dmqmc_calc(dmqmc_rdm_r2)) then
            renyi_2 = rep_loop_sum(min_ind(rdm_r2_ind):max_ind(rdm_r2_ind))
        end if

        ! Average the spawning rate.
        rspawn = rspawn/(ncycles*nprocs)

    end subroutine communicated_dmqmc_estimators

    subroutine communicate_inst_rdms()

        ! Perform annihilation between 'RDM particles'. This is done exactly
        ! by calling the normal annihilation routine for a spawned object,
        ! done for each RDM being calculated, which communicates the RDM
        ! particles and annihilates them.

        ! Also, once annihilation has occurred, calculate all instantaneous
        ! RDM estimates.

        use calc, only: doing_dmqmc_calc, dmqmc_rdm_r2
        use fciqmc_data, only: rdm_spawn, rdms, nrdms, rdm_traces, renyi_2
        use hash_table, only: reset_hash_table
        use spawn_data, only: annihilate_wrapper_spawn_t

        integer :: irdm

        ! WARNING: cannot pass rdm_spawn%spawn to procedures expecting an
        ! array of type spawn_t due to a bug in gfortran which results in
        ! memory deallocations!
        ! See https://groups.google.com/forum/#!topic/comp.lang.fortran/VuFvOsLs6hE
        ! and http://gcc.gnu.org/bugzilla/show_bug.cgi?id=58310.
        ! The explicit loop is also meant to be more efficient anyway, as it
        ! prevents any chance of copy-in/copy-out...
        do irdm = 1, nrdms
            call annihilate_wrapper_spawn_t(rdm_spawn(irdm)%spawn, .false.)
            ! Now is also a good time to reset the hash table (otherwise we
            ! attempt to lookup non-existent data in the next cycle!).
            call reset_hash_table(rdm_spawn(irdm)%ht)
            ! spawn_t comms changes the memory used by spawn%sdata.  Make
            ! sure the hash table always uses the currently 'active'
            ! spawning memory.
            rdm_spawn(irdm)%ht%data_label => rdm_spawn(irdm)%spawn%sdata
        end do

        call calculate_rdm_traces(rdms, rdm_spawn%spawn, rdm_traces)
        if (doing_dmqmc_calc(dmqmc_rdm_r2)) call calculate_rdm_renyi_2(rdms, rdm_spawn%spawn, renyi_2)

        do irdm = 1, nrdms
            rdm_spawn(irdm)%spawn%head = rdm_spawn(irdm)%spawn%head_start
        end do

    end subroutine communicate_inst_rdms

    subroutine update_shift_dmqmc(qmc_in, qs, loc_totnparticles, loc_totnparticles_old)

        ! In:
        !    qmc_in: input options relating to QMC methods.
        !    loc_totnparticles: total number (across all processors) of
        !        particles in the simulation at end of the previous report loop.
        !    loc_totnparticles_old: total number (across all processors) of
        !        particles in the simulation currently.

        use energy_evaluation, only: update_shift
        use qmc_data, only: qmc_in_t, qmc_state_t

        type(qmc_in_t), intent(in) :: qmc_in
        type(qmc_state_t), intent(inout) :: qs
        real(p), intent(in) :: loc_totnparticles(:)
        real(p), intent(in) :: loc_totnparticles_old(:)
        integer :: ireplica

        do ireplica = 1, size(loc_totnparticles)
            if (qs%vary_shift(ireplica)) then
                call update_shift(qmc_in, qs, qs%shift(ireplica), loc_totnparticles_old(ireplica), &
                    loc_totnparticles(ireplica), qmc_in%ncycles)
            end if
            if (loc_totnparticles(ireplica) > qmc_in%target_particles .and. (.not. qs%vary_shift(ireplica))) &
                qs%vary_shift(ireplica) = .true.
        end do

    end subroutine update_shift_dmqmc

    subroutine update_dmqmc_estimators(sys, dmqmc_in, idet, iteration, cdet, H00, nload_slots, psip_list, dmqmc_estimates)

        ! This function calls the processes to update the estimators which have
        ! been requested by the user to be calculated. First, calculate the
        ! excitation level between the two bitstrings corresponding to the the
        ! two ends. Then add the contribution from the current density matrix
        ! element to the trace, which is always calculated. Then call other
        ! estimators, as required.

        ! In:
        !    sys: system being studied.
        !    dmqmc_in: input options for DMQMC.
        !    idet: Current position in the main particle list.
        !    iteration: current Monte Carlo cycle.
        !    H00: diagonal Hamiltonian element for the reference.
        !    nload_slots: number of load balancing slots (per processor).
        !    psip_list: particle information/lists.
        ! In/Out:
        !    cdet: det_info_t object containing information of current density
        !        matrix element.
        !    dmqmc_estimates: type containing dmqmc estimates.

        use calc, only: doing_dmqmc_calc, dmqmc_energy, dmqmc_staggered_magnetisation
        use calc, only: dmqmc_energy_squared, dmqmc_correlation, dmqmc_full_r2
        use excitations, only: get_excitation, excit_t
        use fciqmc_data, only: doing_reduced_dm
        use fciqmc_data, only: accumulated_probs
        use fciqmc_data, only: accumulated_probs_old, real_factor
        use proc_pointers, only:  update_dmqmc_energy_and_trace_ptr, update_dmqmc_stag_mag_ptr
        use proc_pointers, only: update_dmqmc_energy_squared_ptr, update_dmqmc_correlation_ptr
        use determinants, only: det_info_t
        use system, only: sys_t
        use qmc_data, only: reference_t, particle_t
        use dmqmc_data, only: dmqmc_in_t, dmqmc_estimates_t, energy_ind, energy_squared_ind, &
                              correlation_fn_ind, staggered_mag_ind, full_r2_ind

        type(sys_t), intent(in) :: sys
        type(dmqmc_in_t), intent(in) :: dmqmc_in
        integer, intent(in) :: idet, iteration
        type(det_info_t), intent(inout) :: cdet
        real(p), intent(in) :: H00
        integer, intent(in) :: nload_slots
        type(particle_t), intent(in) :: psip_list
        type(dmqmc_estimates_t), intent(inout) :: dmqmc_estimates

        type(excit_t) :: excitation
        real(p) :: unweighted_walker_pop(psip_list%nspaces)

        ! Get excitation.
        excitation = get_excitation(sys%nel, sys%basis, psip_list%states(:sys%basis%string_len,idet), &
                         psip_list%states((1+sys%basis%string_len):sys%basis%tensor_label_len,idet))

        ! When performing importance sampling the result is that certain
        ! excitation levels have smaller psips populations than the true density
        ! matrix by some factor. In these cases, we want to multiply the psip
        ! population by this factor to calculate the contribution from these
        ! excitation levels correctly.

        ! In the case of no importance sampling, unweighted_walker_pop = particle_t%pops(1,idet).
        unweighted_walker_pop = real(psip_list%pops(:,idet),p)*accumulated_probs(excitation%nexcit)/&
            real_factor

        ! The following only use the populations with ireplica = 1, so only call
        ! them if the determinant is occupied in the first replica.
        associate (est => dmqmc_estimates)
            if (abs(unweighted_walker_pop(1)) > 0) then
                ! See which estimators are to be calculated, and call the
                ! corresponding procedures.
                ! Energy
                if (doing_dmqmc_calc(dmqmc_energy)) call update_dmqmc_energy_and_trace_ptr&
                        &(sys, excitation, cdet, H00, unweighted_walker_pop(1), psip_list%dat(1, idet), &
                          est%trace, est%numerators(energy_ind))
                ! Energy squared.
                if (doing_dmqmc_calc(dmqmc_energy_squared)) call update_dmqmc_energy_squared_ptr&
                    &(sys, cdet, excitation, H00, unweighted_walker_pop(1), est%numerators(energy_squared_ind))
                ! Spin-spin correlation function.
                if (doing_dmqmc_calc(dmqmc_correlation)) call update_dmqmc_correlation_ptr&
                    &(sys, cdet, excitation, H00, unweighted_walker_pop(1), dmqmc_in%correlation_mask, &
                      est%numerators(correlation_fn_ind))
                ! Staggered magnetisation.
                if (doing_dmqmc_calc(dmqmc_staggered_magnetisation)) call update_dmqmc_stag_mag_ptr&
                    &(sys, cdet, excitation, H00, unweighted_walker_pop(1), est%numerators(staggered_mag_ind))
                ! Excitation distribution.
                if (dmqmc_in%calc_excit_dist) est%excit_dist(excitation%nexcit) = &
                    est%excit_dist(excitation%nexcit) + real(abs(psip_list%pops(1,idet)),p)/real_factor
                ! Excitation distribtuion for calculating importance sampling weights.
                if (dmqmc_in%find_weights .and. iteration > dmqmc_in%start_av_excit_dist) est%excit_dist(excitation%nexcit) = &
                    est%excit_dist(excitation%nexcit) + real(abs(psip_list%pops(1,idet)),p)/real_factor
            end if

            ! Full Renyi entropy (S_2).
            if (doing_dmqmc_calc(dmqmc_full_r2)) call update_full_renyi_2(unweighted_walker_pop, excitation%nexcit, &
                                                                          dmqmc_in%half_density_matrix, est%numerators(full_r2_ind))

            ! Update the contribution to the trace from other replicas
            if (dmqmc_in%replica_tricks .and. excitation%nexcit == 0) then
                est%trace(2) = est%trace(2) + unweighted_walker_pop(2)
            end if
        end associate

        ! Reduced density matrices.
        if (doing_reduced_dm) call update_reduced_density_matrix_heisenberg&
            &(sys%basis, cdet, excitation, psip_list%pops(:,idet), iteration, nload_slots, dmqmc_in%start_av_rdm)

        accumulated_probs_old = accumulated_probs

    end subroutine update_dmqmc_estimators

    subroutine dmqmc_energy_and_trace(sys, excitation, cdet, H00, pop, diagonal_contribution, trace, energy)

        ! Add the contribution for the current density matrix element to the thermal
        ! energy estimate.

        ! In:
        !    sys: system being studied.
        !    excitation: excit_t type variable which stores information on
        !        the excitation between the two bitstring ends, corresponding
        !        to the two labels for the density matrix element.
        !    H00: diagonal Hamiltonian element for the reference.
        !    pop: number of particles on the current density matrix
        !        element.
        !    cdet: det_info_t object containing bit strings of densitry matrix
        !    element under consideration.
        !    diagonal_contribution: <D_i|H|D_i>-<D0|H|D0>
        ! In/Out:
        !    trace: total population on diagonal elements of density matrix
        !    energy: current thermal energy estimate.

        use determinants, only: det_info_t
        use system, only: sys_t
        use excitations, only: excit_t
        use proc_pointers, only: update_proj_energy_ptr

        type(sys_t), intent(in) :: sys
        type(excit_t), intent(inout) :: excitation
        type(det_info_t), intent(in) :: cdet
        real(p), intent(in) :: H00, pop
        real(p), intent(in) :: diagonal_contribution
        real(p), intent(inout) :: trace(:)
        real(p), intent(inout) :: energy

        real(p) :: hmatel

        ! Update trace and off-diagonal contributions to the total enegy
        call update_proj_energy_ptr(sys, cdet%f2, cdet, pop, trace(1), energy, excitation, hmatel)

        ! Update diagaonal contribution to the total energy
        if (excitation%nexcit == 0) energy = energy + (diagonal_contribution+H00)*pop

    end subroutine dmqmc_energy_and_trace

    subroutine dmqmc_energy_and_trace_propagate(sys, excitation, cdet, H00, pop, diagonal_contribution, trace, energy)

        ! Add the contribution for the current density matrix element to the thermal
        ! energy estimate. Routine is specific to when using propagate_to_beta
        ! option.

        ! In:
        !    sys: system being studied.
        !    excitation: excit_t type variable which stores information on
        !        the excitation between the two bitstring ends, corresponding
        !        to the two labels for the density matrix element.
        !    H00: diagonal Hamiltonian element for the reference.
        !    pop: number of particles on the current density matrix
        !        element.
        !    cdet: det_info_t object containing bit strings of densitry matrix
        !    element under consideration.
        !    diagonal_contribution: <D_i|H|D_i>-<D0|H|D0>
        ! In/Out:
        !    trace: total population on diagonal elements of density matrix
        !    energy: current thermal energy estimate.

        use determinants, only: det_info_t
        use system, only: sys_t
        use excitations, only: excit_t
        use proc_pointers, only: update_proj_energy_ptr, sc0_ptr

        type(sys_t), intent(in) :: sys
        type(excit_t), intent(inout) :: excitation
        type(det_info_t), intent(in) :: cdet
        real(p), intent(in) :: H00, pop
        real(p), intent(in) :: diagonal_contribution
        real(p), intent(inout) :: trace(:)
        real(p), intent(inout) :: energy

        real(p) :: hmatel

        ! Update trace and off-diagonal contributions to the total enegy
        call update_proj_energy_ptr(sys, cdet%f2, cdet, pop, trace(1), energy, excitation, hmatel)

        ! Update diagaonal contribution to the total energy
        if (excitation%nexcit == 0) energy = energy + sc0_ptr(sys, cdet%f)*pop

    end subroutine dmqmc_energy_and_trace_propagate

    subroutine dmqmc_energy_squared_heisenberg(sys, cdet, excitation, H00, walker_pop, energy_sq)

        ! For the Heisenberg model only.
        ! Add the contribution from the current density matrix element to the
        ! thermal energy squared estimate.

        ! In:
        !    sys: system being studied.
        !    excitation: excit_t type variable which stores information on the
        !        excitation between the two bitstring ends, corresponding to the
        !        two labels for the density matrix element.
        !    cdet: information on the density matrix element, in particular
        !        cdet%data contains (<D_i|H|D_i>+<D_j|H|D_j>)/2 - H00, where
        !        |D_i> and |D_j> are the states labelling the density matrix
        !        element.
        !    H00: diagonal Hamiltonian element for the reference.
        !    walker_pop: number of particles on the current density matrix
        !        element.
        ! In/Out:
        !    energy_sq: running total of thermal energy squared estimate.

        use excitations, only: excit_t
        use determinants, only: det_info_t
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        type(det_info_t), intent(in) :: cdet
        type(excit_t), intent(in) :: excitation
        real(p), intent(in) :: H00, walker_pop
        real(p), intent(inout) :: energy_sq

        integer :: bit_element1, bit_position1, bit_element2, bit_position2
        real(p) :: sum_H1_H2, J_coupling_squared

        sum_H1_H2 = 0.0_p
        J_coupling_squared = sys%heisenberg%J**2/16

        if (excitation%nexcit == 0) then
            ! If there are 0 excitations then either nothing happens twice, or we
            ! flip the same pair of spins twice. The Hamiltonian element for doing
            ! nothing is just the diagonal element. For each possible pairs of
            ! spins which can be flipped, there is a mtarix element of -J/2, so
            ! we just need to count the number of such pairs, which can be found 
            ! simply from the diagonal element.

            sum_H1_H2 = (cdet%data(1)+H00)**2
            associate(sh=>sys%heisenberg)
                sum_H1_H2 = sum_H1_H2 + 2*J_coupling_squared*sh%nbonds + sh%J*(cdet%data(1)+H00)/2
            end associate

        else if (excitation%nexcit == 1) then
            ! If there is only one excitation (2 spins flipped) then the
            ! contribution to H^2 depend on the positions of the spins relative
            ! to one another. If the the spins are nearest neighbors then we
            ! could either do nothing and then flip the pair, or flip the pair
            ! and then do nothing. If next nearest neighbors and there is only
            ! one two-bond path to get from one spin to the other, we first flip
            ! the pair on the first bond, then flip the pair on the second bond.
            ! This flipping can only be done in exactly one order, not both - the
            ! two spins which change are opposite, so the middle spin will
            ! initially only be the same as one or the other spin. This is nice, 
            ! because we don't have check which way up the intermediate spin is -
            ! there will always be one order which contributes. If there are two
            ! such paths, then this could happen by either paths, but again, the
            ! two intermediate spins will only allow one order of spin flipping
            ! for each path, no matter which way up they are, so we only need to
            ! check if there are two possible paths.

            if (sys%real_lattice%next_nearest_orbs(excitation%from_orb(1),excitation%to_orb(1)) /= 0_i0) then
                ! Contribution for next-nearest neighbors.
                sum_H1_H2 = 4.0_p*J_coupling_squared*sys%real_lattice%next_nearest_orbs(excitation%from_orb(1),excitation%to_orb(1))
            end if
            ! Contributions for nearest neighbors.
            ! Note, for certain lattices, such as the triangular lattice, two
            ! spins can be both nearest neighbors *and* next-nearest neighbors.
            ! Therefore, it is necessary in general to check for both situations.
            bit_position1 = sys%basis%bit_lookup(1,excitation%from_orb(1))
            bit_element1 = sys%basis%bit_lookup(2,excitation%from_orb(1))
            if (btest(sys%real_lattice%connected_orbs(bit_element1, excitation%to_orb(1)), bit_position1)) &
                sum_H1_H2 = sum_H1_H2 - sys%heisenberg%J*(cdet%data(1)+H00)

        else if (excitation%nexcit == 2) then
            ! If there are two excitations (4 spins flipped) then, once again,
            ! the contribution to the thermal energy squared will depend on the
            ! positions of the spins. If there are two pairs of spins flipped
            ! which are separated then there is one way for this to happen - by
            ! flipping one pair, and then the other (this also requires that the
            ! two spins within each neighboring pair are opposite, as ever for
            ! the Heisenberg model). These two flips can happen in either order.
            ! In some cases the spins may be such that we may pair the spins in
            ! more than one way. For example, if the four spins are in a square
            ! shape, or for a 4-by-4 Heisenberg model, the spins could be
            ! connected across the whole lattice, forming a ring due to the
            ! periodic boundaries. In these cases it may be possible to perform
            ! the spin flips by pairing them in either of two ways. To account
            ! for this possibility we have to try and pair the spins in both
            ! ways, so we always check both if statements below. Again, once
            ! these pairings have been chosen, the flips can be performed in
            ! either order.

            bit_position1 = sys%basis%bit_lookup(1,excitation%from_orb(1))
            bit_element1 = sys%basis%bit_lookup(2,excitation%from_orb(1))
            bit_position2 = sys%basis%bit_lookup(1,excitation%from_orb(2))
            bit_element2 = sys%basis%bit_lookup(2,excitation%from_orb(2))
            if (btest(sys%real_lattice%connected_orbs(bit_element1, excitation%to_orb(1)), bit_position1) .and. &
                btest(sys%real_lattice%connected_orbs(bit_element2, excitation%to_orb(2)), bit_position2)) &
                sum_H1_H2 = 8.0*J_coupling_squared
            if (btest(sys%real_lattice%connected_orbs(bit_element1, excitation%to_orb(2)), bit_position1) .and. &
                btest(sys%real_lattice%connected_orbs(bit_element2, excitation%to_orb(1)), bit_position2)) &
                sum_H1_H2 = sum_H1_H2 + 8.0*J_coupling_squared

        end if

        energy_sq = energy_sq + sum_H1_H2*walker_pop

    end subroutine dmqmc_energy_squared_heisenberg

    subroutine dmqmc_correlation_function_heisenberg(sys, cdet, excitation, H00, walker_pop, correlation_mask, &
                                                     correlation_fn)

        ! For the Heisenberg model only.
        ! Add the contribution from the current density matrix element to the
        ! thermal spin correlation function estimator.

        ! In:
        !    sys: system being studied.
        !    cdet: information on density matrix element, in particular cdet%f,
        !       the bit strings of each label.
        !    excitation: excit_t type variable which stores information on
        !       the excitation between the two bitstring ends, corresponding to
        !       the two labels for the density matrix element.
        !    H00: diagonal Hamiltonian element for the reference.
        !    walker_pop: number of particles on the current density matrix
        !       element.
        !    correlation_mask: masks for calculating spin-spin correlation
        !       function.
        ! In/Out:
        !    correlation_fn: running estimate for the correlation function estimator.

        use bit_utils, only: count_set_bits
        use determinants, only: det_info_t
        use excitations, only: excit_t
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        type(det_info_t), intent(in) :: cdet
        type(excit_t), intent(in) :: excitation
        real(p), intent(in) :: H00, walker_pop
        integer(i0), allocatable, intent(in) :: correlation_mask(:)
        real(p), intent(inout) :: correlation_fn

        integer(i0) :: f(sys%basis%string_len)
        integer :: bit_element1, bit_position1, bit_element2, bit_position2
        integer :: sign_factor

        ! If no excitation, we want the diagonal element of the correlation
        ! function operator.
        if (excitation%nexcit == 0) then
            ! If the two spins i and j are the same, the matrix element is +1/4.
            ! If they are different, the matrix element is -1/4. So we want
            ! sign_factor to be +1 in the former case, and -1 in the latter case.
            ! f as calculated below will have 0's at sites other than i and j,
            ! and the same values as cdet%f at i and j. Hence, if f has
            ! two 1's or no 1's, we want sign_factor = +1. Else if we have one 1,
            ! we want sign_factor = -1.
            f = iand(cdet%f(:sys%basis%string_len), correlation_mask)
            ! Count if we have zero, one or two 1's.
            sign_factor = sum(count_set_bits(f))
            ! The operation below will map 0 and 2 to +1, and will map 1 to -1,
            ! as is easily checked.
            sign_factor = (mod(sign_factor+1,2)*2)-1
            ! Hence sign_factor can be used to find the matrix element, as used
            ! below.
            correlation_fn = correlation_fn + (sign_factor*(walker_pop/4))
        else if (excitation%nexcit == 1) then
            ! If not a diagonal element, but only a single excitation, then the
            ! corresponding matrix element will be 1/2 if and only if the two
            ! sites which are flipped are sites i and j, else it will be 0. We
            ! assume that excitations will only be set if i and j are opposite
            ! (else they could not be flipped, for ms=0).
            bit_position1 = sys%basis%bit_lookup(1,excitation%from_orb(1))
            bit_element1 = sys%basis%bit_lookup(2,excitation%from_orb(1))
            bit_position2 = sys%basis%bit_lookup(1,excitation%to_orb(1))
            bit_element2 = sys%basis%bit_lookup(2,excitation%to_orb(1))
            if (btest(correlation_mask(bit_element1), bit_position1) .and. btest(correlation_mask(bit_element2), bit_position2)) &
                correlation_fn = correlation_fn + walker_pop/2
        end if

    end subroutine dmqmc_correlation_function_heisenberg

    subroutine dmqmc_stag_mag_heisenberg(sys, cdet, excitation, H00, walker_pop, staggered_mag)

        ! For the Heisenberg model only.
        ! Add the contribution from the current density matrix element to the
        ! thermal staggered magnetisation estimate.

        ! In:
        !    sys: system being studied.
        !    cdet: information on density matrix element, in particular cdet%f,
        !       the bit strings of each label.
        !    excitation: excit_t type variable which stores information on
        !       the excitation between the two bitstring ends, corresponding
        !       to the two labels for the density matrix element.
        !    H00: diagonal Hamiltonian element for the reference.
        !    walker_pop: number of particles on the current density matrix
        !       element.
        ! In/Out:
        !    staggered_mag: running estimate for the staggered magnitisation.

        use bit_utils, only: count_set_bits
        use determinants, only: det_info_t
        use excitations, only: excit_t
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        type(det_info_t), intent(in) :: cdet
        type(excit_t), intent(in) :: excitation
        real(p), intent(in) :: H00, walker_pop
        real(p), intent(inout) :: staggered_mag

        integer :: bit_element1, bit_position1, bit_element2, bit_position2
        integer(i0) :: f(sys%basis%string_len)
        integer :: n_up_plus
        integer :: total_sum

        total_sum = 0

        ! This is for the staggered magnetisation squared. This is given by
        ! M^2 = M_x^2 + M_y^2 + M_z^2
        ! where M_x = \sum{i}(-1)^i S_i^x, etc... and S_i^x is the spin operator
        ! at site i, in the x direction.
        if (excitation%nexcit == 0) then
            ! Need to calculate the number of spins up on sublattice 1:
            ! N_u(+) - Number up on + sublattice (where (-1)^i above is +1)
            ! Note the number of down spins on a sublattice is easily obtained
            ! from N_u(+) since there are N/2 spins on each - and the number of
            ! spins up on a different sublattice is easily obtained since there
            ! are nel spins up in total. Hence the matrix element will be written
            ! only in terms of the number of up spins on sublattice 1, to save
            ! computation.
            f = iand(cdet%f(:sys%basis%string_len), sys%heisenberg%lattice_mask)
            n_up_plus = sum(count_set_bits(f))
            ! Below, the term in brackets and middle term come from the z
            ! component (the z operator is diagonal) and one nsites/4 factor
            ! comes from the x operator, the other nsites/4 factor from the y
            ! operator.
            total_sum = (2*n_up_plus-sys%nel)**2 + (sys%lattice%nsites/2)
        else if (excitation%nexcit == 1) then
            ! Off-diagonal elements from the y and z operators. For the pair of
            ! spins that are flipped, if they are on the same sublattice, we get
            ! a factor of 1, or if on different sublattices, a factor of -1.
            bit_position1 = sys%basis%bit_lookup(1,excitation%from_orb(1))
            bit_element1 = sys%basis%bit_lookup(2,excitation%from_orb(1))
            bit_position2 = sys%basis%bit_lookup(1,excitation%to_orb(1))
            bit_element2 = sys%basis%bit_lookup(2,excitation%to_orb(1))
            if (btest(sys%heisenberg%lattice_mask(bit_element1), bit_position1)) total_sum = total_sum+1
            if (btest(sys%heisenberg%lattice_mask(bit_element2), bit_position2)) total_sum = total_sum+1
            ! The operation below will map 0 and 2 to +1, and will map 1 to -1,
            ! as is easily checked. We want this - if both or no spins on this
            ! sublattice, then both on same sublattice either way, so plus one.
            ! Else they are on different sublattices, so we want a factor of -1,
            ! as we get.
            total_sum = (mod(total_sum+1,2)*2)-1
        end if

        staggered_mag = staggered_mag + (real(total_sum)/real(sys%lattice%nsites**2))*walker_pop

    end subroutine dmqmc_stag_mag_heisenberg

    subroutine update_full_renyi_2(walker_pop, excit_level, half_density_matrix, full_r2)

        ! Add the contribution from the current density matrix element to the
        ! Renyi entropy (S_2) of the full density matrix.

        ! In:
        !    walker_pop: number of particles on the current density matrix
        !        element, for both replicas.
        !    excit_level: The excitation level between the two bitstrings
        !        contributing to the full density matrix bitstring.
        !    half_density_matrix: reflect psips spawned from lower triangle into
        !        the upper one?
        ! In/Out:
        !    full_r2: running estimate for Renyi entropy of the full density
        !       matrix.

        real(p), intent(in) :: walker_pop(:)
        integer, intent(in) :: excit_level
        logical, intent(in) :: half_density_matrix
        real(p), intent(inout) :: full_r2

        if (half_density_matrix .and. excit_level /= 0) then
            ! With the half-density matrix option, only the upper-half of the
            ! density matrix is stored, but the off-diagonal elements are twice
            ! as large instead. Thus, a product of off-diagonal elements is four
            ! times as large. We want two 'correct size' contributions, so we
            ! need to divide these contirbutions by two to get this.
            full_r2 = full_r2 + walker_pop(1)*walker_pop(2)/2.0_p
        else
            full_r2 = full_r2 + walker_pop(1)*walker_pop(2)
        end if

    end subroutine update_full_renyi_2

    subroutine update_reduced_density_matrix_heisenberg(basis, cdet, excitation, walker_pop, &
                                                        iteration, nload_slots, start_av_rdm)

        ! Add the contribution from the current walker to the reduced density
        ! matrices being sampled. This is performed by 'tracing out' the
        ! subsystem B spins.

        ! Applicable only to the Heisenberg model.

        ! This procedure takes the two determinants bitstrings for the current
        ! psips and, if the two bitstrings of the B subsystem are identical,
        ! adds the walker population to the corresponding reduced density matrix
        ! element.

        ! In:
        !    basis: information about the single-particle basis.
        !    cdet: information on density matrix element, in particular cdet%f,
        !       the bit strings of each label.
        !    excitation: excit_t type variable which stores information on
        !        the excitation between the two bitstring ends, corresponding to
        !        the two labels for the density matrix element.
        !    walker_pop: number of particles on the current density matrix
        !        element. Note that this walker population is still weighted
        !        by the importance sampling factors. These factors must be
        !        removed before any estimates can be calculated.
        !    iteration: interation number.  No accumulation of the RDM is
        !        performed if iteration <= start_av_rdm.
        !    nload_slots: number of load balancing slots (per processor).
        !    start_av_rdm: iteration we start averaging the rdm on.

        use basis_types, only: basis_t
        use determinants, only: det_info_t
        use dmqmc_procedures, only: decode_dm_bitstring
        use excitations, only: excit_t
        use fciqmc_data, only: reduced_density_matrix
        use fciqmc_data, only: calc_inst_rdm, calc_ground_rdm, rdms, nrdms
        use fciqmc_data, only: rdm_spawn, accumulated_probs
        use fciqmc_data, only: nsym_vec, real_factor
        use spawning, only: create_spawned_particle_rdm

        type(basis_t), intent(in) :: basis
        type(det_info_t), intent(in) :: cdet
        integer, intent(in) :: iteration
        integer(int_p), intent(in) :: walker_pop(:)
        type(excit_t), intent(in) :: excitation
        integer, intent(in) :: nload_slots
        integer, intent(in) :: start_av_rdm

        real(p) :: unweighted_walker_pop(size(walker_pop))
        integer :: irdm, isym, ireplica
        integer(i0) :: f1(basis%string_len), f2(basis%string_len)

        if (.not. (iteration > start_av_rdm .or. calc_inst_rdm)) return

        ! Loop over all RDMs to be calculated.
        do irdm = 1, nrdms
        ! Loop over every symmetry-equivalent subsystem for this RDM.
        do isym = 1, nsym_vec

        ! Apply the mask for the B subsystem to set all sites in the A
        ! subsystem to 0.
        f1 = iand(rdms(irdm)%B_masks(:,isym),cdet%f(:basis%string_len))
        f2 = iand(rdms(irdm)%B_masks(:,isym),cdet%f(basis%string_len+1:basis%tensor_label_len))

        ! Once this is done, check if the resulting bitstrings (which can
        ! only possibly have 1's in the B subsystem) are identical. If
        ! they are, then this psip contributes to the reduced density
        ! matrix for subsystem A. This is because we get the reduced
        ! density matrix for A by 'tracing out' over B, which in practice
        ! means only keeping matrix elements that are on the diagonal for
        ! subsystem B.
        if (sum(abs(f1-f2)) == 0_i0) then
            ! Call a function which maps the subsystem A state to two RDM
            ! bitstrings.
            call decode_dm_bitstring(basis, cdet%f,irdm,isym)

            if (calc_ground_rdm) then
                ! The above routine actually maps to numbers between 0
                ! and 2^rdms(1)%A_nsites-1, but the smallest and largest
                ! reduced density matrix indices are one more than these,
                ! so add one.
                rdms(irdm)%end1 = rdms(irdm)%end1 + 1
                rdms(irdm)%end2 = rdms(irdm)%end2 + 1
                unweighted_walker_pop = real(walker_pop,p)*&
                    accumulated_probs(excitation%nexcit)/real_factor
                ! Note, when storing the entire RDM (as done here), the
                ! maximum value of rdms(i)%string_len is 1, so we
                ! only consider this one element here.
                reduced_density_matrix(rdms(irdm)%end1(1),rdms(irdm)%end2(1)) = &
                    reduced_density_matrix(rdms(irdm)%end1(1),rdms(irdm)%end2(1)) + unweighted_walker_pop(1)
            end if

            if (calc_inst_rdm) then
                do ireplica = 1, size(walker_pop)
                if (abs(walker_pop(ireplica)) > 0) then
                    call create_spawned_particle_rdm(rdms(irdm), walker_pop(ireplica), &
                        ireplica, rdm_spawn(irdm), nload_slots)
                end if
                end do
            end if

        end if
        end do

        end do

    end subroutine update_reduced_density_matrix_heisenberg

    subroutine call_ground_rdm_procedures(beta_cycle)

        ! Wrapper for calling ground-state RDM procedures (*not*
        ! beta-dependent RDMs).

        ! In:
        !    beta_cycle: index of the beta loop being performed.

        use checking, only: check_allocate, check_deallocate
        use fciqmc_data, only: rdms, reduced_density_matrix
        use fciqmc_data, only: doing_vn_entropy, doing_concurrence
        use fciqmc_data, only: output_rdm, rdm_unit, rdm_traces
        use parallel
        use utils, only: get_free_unit, append_ext, int_fmt

        integer, intent(in) :: beta_cycle
        real(p), allocatable :: old_rdm_elements(:)
        integer :: i, j, k, ierr, new_unit
        character(10) :: rdm_filename
        ! If in parallel then merge the reduced density matrix onto one
        ! processor.
#ifdef PARALLEL

        real(dp), allocatable :: dm(:,:)
        real(dp), allocatable :: dm_sum(:,:)
        integer :: num_eigv

        num_eigv = 2**rdms(1)%A_nsites

        allocate(dm(num_eigv,num_eigv), stat=ierr)
        call check_allocate('dm',num_eigv**2,ierr)
        allocate(dm_sum(num_eigv,num_eigv), stat=ierr)
        call check_allocate('dm_sum',num_eigv**2,ierr)

        dm = reduced_density_matrix
        call mpi_allreduce(dm, dm_sum, size(dm), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
        reduced_density_matrix = dm_sum

        deallocate(dm)
        call check_deallocate('dm',ierr)
        deallocate(dm_sum)
        call check_deallocate('dm_sum',ierr)

#endif

        rdm_traces = 0.0_p

        if (parent) then
            ! Force the reduced density matrix to be symmetric by averaging the
            ! upper and lower triangles.
            do i = 1, ubound(reduced_density_matrix,1)
                do j = 1, i-1
                    reduced_density_matrix(i,j) = 0.5_p*(reduced_density_matrix(i,j) + reduced_density_matrix(j,i))
                    reduced_density_matrix(j,i) = reduced_density_matrix(i,j)
                end do
                ! Add current contribution to the trace.
                rdm_traces(1,1) = rdm_traces(1,1) + reduced_density_matrix(i,i)
            end do

            ! Call the routines to calculate the desired quantities.
            if (doing_vn_entropy) call calculate_vn_entropy(rdm_traces(1,1))
            if (doing_concurrence) call calculate_concurrence()

            write (6,'(1x,"# RDM trace =",1X,es17.10)') rdm_traces(1,1)

            if (output_rdm) then
                new_unit = get_free_unit()
                call append_ext('rdm', beta_cycle, rdm_filename)
                open(new_unit, file=trim(rdm_filename), status='replace')
                write(new_unit,'(a5,1x,es15.8)') "Trace", rdm_traces(1,1)
                do i = 1, ubound(reduced_density_matrix,1)
                    do j = i, ubound(reduced_density_matrix,1)
                        write(new_unit,'(a1,'//int_fmt(i,0)//',a1,'//int_fmt(j,0)//',a1,1x,es15.8)') &
                            "(",i,",",j,")", reduced_density_matrix(i,j)
                    end do
                end do
                close(new_unit)
            end if

        end if

    end subroutine call_ground_rdm_procedures

    subroutine calculate_vn_entropy(trace_rdm)

        ! Calculate the Von Neumann Entropy. Use lapack to calculate the
        ! eigenvalues {\lambda_j} of the reduced density matrix.
        ! Then VN Entropy S = -\sum_j\lambda_j\log_2{\lambda_j}.

        ! Need to paralellise for large subsystems and introduce test to check
        ! whether diagonalisation should be performed in serial or paralell.

        ! In:
        !    trace_rdm: The trace of the RDM being considered.

        use checking, only: check_allocate, check_deallocate
        use fciqmc_data, only: reduced_density_matrix, rdms

        real(p), intent(in) :: trace_rdm
        integer :: i, rdm_size
        integer :: info, ierr, lwork
        real(p), allocatable :: work(:)
        real(p), allocatable :: dm_tmp(:,:)
        real(p) :: eigv(2**rdms(1)%A_nsites)
        real(p) :: vn_entropy
        logical :: thrown_away

        rdm_size = 2**rdms(1)%A_nsites
        vn_entropy = 0.0_p

        ! Find the optimal size of the workspace.
        allocate(work(1), stat=ierr)
        call check_allocate('work',1,ierr)
#ifdef SINGLE_PRECISION
        call ssyev('N', 'U', rdm_size, reduced_density_matrix, rdm_size, eigv, work, -1, info)
#else
        call dsyev('N', 'U', rdm_size, reduced_density_matrix, rdm_size, eigv, work, -1, info)
#endif

        lwork = nint(work(1))
        deallocate(work)
        call check_deallocate('work',ierr)

        ! Now perform the diagonalisation.
        allocate(work(lwork), stat=ierr)
        call check_allocate('work',lwork,ierr)

        ! The matrix input into the following diagonalisation routines will have
        ! their upper half (including the diagonal) destroyed. We might want
        ! reduced_desntiy_matrix later, so use some temporary space:
        allocate(dm_tmp(rdm_size,rdm_size), stat=ierr)
        call check_allocate('dm_tmp',rdm_size**2,ierr)
        dm_tmp = reduced_density_matrix

#ifdef SINGLE_PRECISION
        call ssyev('N', 'U', rdm_size, dm_tmp, rdm_size, eigv, work, lwork, info)
#else
        call dsyev('N', 'U', rdm_size, dm_tmp, rdm_size, eigv, work, lwork, info)
#endif
        thrown_away = .false.
        write(6,'(1X,"# Eigenvalues thrown away:",1X)',advance='no')
        do i = 1, ubound(eigv,1)
            if (eigv(i) < depsilon) then
                write(6,'(es15.8,2x)',advance='no') eigv(i)
                thrown_away = .true.
                cycle
            end if
            vn_entropy = vn_entropy - eigv(i)*(log(eigv(i))/log(2.0_p))
        end do
        if (thrown_away) then
            write(6,'()',advance='yes')
        else
            write(6,'(1X,"none")',advance='yes')
        end if
        write (6,'(1x,"# Unnormalised von Neumann entropy =",1X,es17.10)') vn_entropy

        deallocate(dm_tmp)
        call check_deallocate('dm_tmp',ierr)

    end subroutine calculate_vn_entropy

    subroutine calculate_concurrence()

        ! Calculate the concurrence of a qubit. For a reduced density matrix
        ! \rho, the concurrence,
        ! C =  max(0, \lamda_1 - \lambda_2 - \lambda_3 -\lambda_4)
        ! where \lambda_i are the eigenvalues of the matrix,
        ! R = \sqrt{\sqrt{\rho}\~{\rho}\sqrt{\rho}},
        ! and
        ! \~\rho = {\sigma_y \otimes \sigma_y} \rho^{\ast} {\sigma_y \otimes \sigma_y}.
        ! \lambda_1 > ... > \lambda_4.

        ! This can be simplified to finding the square root of the eigenvalues
        ! \{\lambda_i\} of \rho\~{\rho} and in the case where \rho is a real,
        ! symmetric matrix then we can further simplify the problem to finding
        ! the eigenvalues of R = \rho \sigma_y \otimes \sigma_y.

        ! Below we have named {\sigma_y \otimes \sigma_y} flip_spin_matrix as in
        ! the literature.

        use checking, only: check_allocate, check_deallocate
        use fciqmc_data, only: reduced_density_matrix
        integer :: info, ierr, lwork
        real(p), allocatable :: work(:)
        real(p) :: reigv(4), ieigv(4)
        real(p) :: concurrence
        real(p) :: rdm_spin_flip(4,4), rdm_spin_flip_tmp(4,4), VL(4,4), VR(4,4)

        ! This will store the 4x4 flip spin matrix \sigma_y \otimes \sigma_y if
        ! concurrence is to be calculated.
        real(p), parameter :: flip_spin_matrix(4,4) = reshape( &
                                    [  0.0_p,  0.0_p, 0.0_p, -1.0_p,  &
                                       0.0_p,  0.0_p, 1.0_p,  0.0_p,  &
                                       0.0_p,  1.0_p, 0.0_p,  0.0_p,  &
                                       -1.0_p,  0.0_p, 0.0_p,  0.0_p  ], shape(flip_spin_matrix))

        ! Make rdm_spin_flip_tmp because sgeev and dgeev delete input matrix.
        rdm_spin_flip_tmp = matmul(reduced_density_matrix, flip_spin_matrix)
        rdm_spin_flip = rdm_spin_flip_tmp

        ! Find the optimal size of the workspace.
        allocate(work(1), stat=ierr)
        call check_allocate('work',1,ierr)
#ifdef SINGLE_PRECISION
        call sgeev('N', 'N', 4, rdm_spin_flip_tmp, 4, reigv, ieigv, VL, 1, VR, 1, work, -1, info)
#else
        call dgeev('N', 'N', 4, rdm_spin_flip_tmp, 4, reigv, ieigv, VL, 1, VR, 1,  work, -1, info)
#endif
        lwork = nint(work(1))
        deallocate(work)
        call check_deallocate('work',ierr)

        ! Now perform the diagonalisation.
        allocate(work(lwork), stat=ierr)
        call check_allocate('work',lwork,ierr)
#ifdef SINGLE_PRECISION
        call sgeev('N', 'N', 4, rdm_spin_flip, 4, reigv, ieigv, VL, 1, VR, 1, work, lwork, info)
#else
        call dgeev('N', 'N', 4, rdm_spin_flip, 4, reigv, ieigv, VL, 1, VR, 1, work, lwork, info)
#endif
        ! Calculate the concurrence. Take abs of eigenvalues so that this is
        ! equivalant to sqauring and then square-rooting.
        concurrence = 2.0_p*maxval(abs(reigv)) - sum(abs(reigv)) 
        concurrence = max(0.0_p, concurrence)
        write (6,'(1x,"# Unnormalised concurrence =",1X,es17.10)') concurrence

    end subroutine calculate_concurrence

    subroutine calculate_rdm_traces(rdm_data, rdm_lists, traces)

        ! In:
        !    rdm_data: Array of rdm_t derived types, holding information about
        !        the various subsystems for which RDMs are being estimated.
        !    rdm_lists: Array of rdm_spawn_t derived types, which hold all of
        !        the RDM psips which belong to this processor.
        ! Out:
        !    r2: The calculated RDM traces.

        use fciqmc_data, only: rdm_t
        use excitations, only: get_excitation_level
        use spawn_data, only: spawn_t

        type(rdm_t), intent(in) :: rdm_data(:)
        type(spawn_t), intent(in) :: rdm_lists(:)
        real(p), intent(out) :: traces(:,:)
        integer :: irdm, i, rdm_bl
        integer, parameter :: thread_id = 0

        traces = 0.0_p

        ! Loop over all RDMs being calculated.
        do irdm = 1, size(rdm_data)
            rdm_bl = rdm_data(irdm)%string_len
            ! Loop over the total population of RDM psips on this processor.
            do i = 1, rdm_lists(irdm)%head(thread_id,0)
                ! If on the diagonal of the RDM...
                if (all( rdm_lists(irdm)%sdata(1:rdm_bl,i) == rdm_lists(irdm)%sdata(rdm_bl+1:2*rdm_bl,i))) then
                    traces(:,irdm) = traces(:,irdm) + &
                        real(rdm_lists(irdm)%sdata(rdm_lists(irdm)%bit_str_len+1:rdm_lists(irdm)%element_len,i),p)
                end if
            end do
        end do

    end subroutine calculate_rdm_traces

    subroutine calculate_rdm_renyi_2(rdm_data, rdm_lists, r2)

        ! Calculate the Renyi entropy (S_2) for all instantaneous RDMs being
        ! calculated.

        ! In:
        !    rdm_data: Array of rdm_t derived types, holding information about the
        !        various subsystems for which RDMs are being estimated.
        !    rdm_lists: Array of rdm_spawn_t derived types, which hold all of
        !        the RDM psips which belong to this processor.
        ! Out:
        !    r2: The calculated Renyi entropies (S_2).

        use fciqmc_data, only: rdm_t
        use excitations, only: get_excitation_level
        use fciqmc_data, only: accumulated_probs_old
        use spawn_data, only: spawn_t

        type(rdm_t), intent(in) :: rdm_data(:)
        type(spawn_t), intent(in) :: rdm_lists(:)
        real(p), intent(out) :: r2(:)
        integer :: i, irdm, excit_level, rdm_bl
        real(p) :: unweighted_pop_1, unweighted_pop_2
        integer, parameter :: thread_id = 0

        r2 = 0.0_p

        ! Loop over all RDMs being calculated.
        do irdm = 1, size(rdm_data)
            rdm_bl = rdm_data(irdm)%string_len
            ! Loop over the total population of RDM psips on this processor.
            do i = 1, rdm_lists(irdm)%head(thread_id,0)

                excit_level = get_excitation_level(int(rdm_lists(irdm)%sdata(1:rdm_bl,i), i0), &
                    int(rdm_lists(irdm)%sdata(rdm_bl+1:2*rdm_bl,i), i0) )

                ! Renormalise the psip populations to correct for the importance
                ! sampling procedure used.
                unweighted_pop_1 = rdm_lists(irdm)%sdata(rdm_lists(irdm)%bit_str_len+1,i)*accumulated_probs_old(excit_level)
                unweighted_pop_2 = rdm_lists(irdm)%sdata(rdm_lists(irdm)%bit_str_len+2,i)*accumulated_probs_old(excit_level)

                ! As we only hold RDM elements above the diagonal, off-diagonal
                ! elements must be counted twice.
                if (excit_level == 0) then
                    r2(irdm) = r2(irdm) + unweighted_pop_1*unweighted_pop_2
                else
                    r2(irdm) = r2(irdm) + 2*unweighted_pop_1*unweighted_pop_2
                end if

            end do
        end do

    end subroutine calculate_rdm_renyi_2

end module dmqmc_estimators

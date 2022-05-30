module dmqmc_estimators

use const

implicit none

! indices for corresponding data items in estimator/population/etc buffers (see
! comments below).
enum, bind(c)
    enumerator :: rspawn_ind = 1
    enumerator :: error_ind
    enumerator :: nocc_states_ind
    enumerator :: nspawned_ind
    enumerator :: nparticles_ind
    enumerator :: trace_ind
    enumerator :: operators_ind
    enumerator :: excit_dist_ind
    enumerator :: ground_rdm_trace_ind
    enumerator :: inst_rdm_trace_ind
    enumerator :: rdm_r2_ind
    enumerator :: mom_dist_ind
    enumerator :: struc_fac_ind
    enumerator :: ref_trace_ind
    enumerator :: ref_D0j_ind
    enumerator :: final_ind ! ensure this remains the last index.
end enum

contains

    subroutine dmqmc_estimate_comms(dmqmc_in, error, nspawn_events, max_num_excits, ncycles, &
                                    psip_list, qs, accumulated_probs_old, dmqmc_estimates)

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
        !    accumulated_probs_old: the value of accumulated_probs on the last report cycle.
        ! In/Out:
        !    error: true if an error has occured on this processor. On output, true if an
        !        error has occured on any processor.
        !    psip_list: particle information.  On output total (ie not
        !        per-processor) quantities are updated.
        !    qs: QMC state containing quantities not specific to DMQMC.
        !    dmqmc_estimates: type containing dmqmc estimates.

        use checking, only: check_allocate, check_deallocate
        use qmc_data, only: particle_t, qmc_state_t
        use parallel
        use dmqmc_data, only: dmqmc_in_t, dmqmc_estimates_t, num_dmqmc_operators
        use calc, only: doing_dmqmc_calc, dmqmc_ref_proj_energy

        type(dmqmc_in_t), intent(in) :: dmqmc_in
        logical, intent(inout) :: error
        integer, intent(in) :: nspawn_events, max_num_excits, ncycles
        real(p), intent(in) :: accumulated_probs_old(0:)
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
        if (dmqmc_in%rdm%calc_inst_rdm) call communicate_inst_rdms(accumulated_probs_old, dmqmc_estimates%subsys_info, &
                                                                    dmqmc_estimates%inst_rdm)

        ! How big is each variable to be communicated?
        nelems(rspawn_ind) = 1
        nelems(error_ind) = 1
        nelems(nocc_states_ind) = 1
        nelems(nspawned_ind) = 1
        nelems(nparticles_ind) = psip_list%nspaces
        nelems(trace_ind) = psip_list%nspaces
        nelems(operators_ind) = num_dmqmc_operators
        nelems(excit_dist_ind) = max_num_excits + 1
        nelems(ground_rdm_trace_ind) = 1
        nelems(inst_rdm_trace_ind) = psip_list%nspaces*dmqmc_estimates%inst_rdm%nrdms
        nelems(rdm_r2_ind) = dmqmc_estimates%inst_rdm%nrdms
        if (dmqmc_in%calc_mom_dist) then
            nelems(mom_dist_ind) = size(dmqmc_estimates%mom_dist%f_k)
        else
            nelems(mom_dist_ind) = 0
        end if
        if (dmqmc_in%calc_struc_fac) then
            nelems(struc_fac_ind) = size(dmqmc_estimates%struc_fac%f_k)
        else
            nelems(struc_fac_ind) = 0
        end if

        if (doing_dmqmc_calc(dmqmc_ref_proj_energy)) then
            nelems(ref_trace_ind) = psip_list%nspaces
            nelems(ref_D0j_ind) = psip_list%nspaces
        else
            nelems(ref_trace_ind) = 0
            nelems(ref_D0j_ind) = 0
        end if

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
                                    psip_list%nstates, nspawn_events, qs%spawn_store%rspawn, error)

#ifdef PARALLEL
        call mpi_allreduce(rep_loop_loc, rep_loop_sum, size(rep_loop_loc), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
        rep_loop_sum = rep_loop_loc
        ierr = 0 ! Prevent warning about unused variable in serial so -Werror can be used.
#endif

        ! Move the communicated quantites to the corresponding variables.
        call communicated_dmqmc_estimators(dmqmc_in, rep_loop_sum, min_ind, max_ind, ncycles, psip_list%tot_nparticles, qs, &
                                           dmqmc_estimates, error)

        ! Clean up.
        deallocate(rep_loop_loc, stat=ierr)
        call check_deallocate('rep_loop_loc', ierr)
        deallocate(rep_loop_sum, stat=ierr)
        call check_deallocate('rep_loop_sum', ierr)

    end subroutine dmqmc_estimate_comms

    subroutine local_dmqmc_estimators(dmqmc_in, dmqmc_estimates, rep_loop_loc, min_ind, max_ind, nparticles, &
                                      nstates_active, nspawn_events, rspawn, error)

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
        !    rspawn: rate of spawning on this processor.
        !    error: whether an error has occured on this processor.
        ! Out:
        !    rep_loop_loc: array containing local quantities to be communicated.

        use calc, only: doing_dmqmc_calc, dmqmc_rdm_r2, dmqmc_ref_proj_energy
        use dmqmc_data, only: dmqmc_in_t, dmqmc_estimates_t

        type(dmqmc_in_t), intent(in) :: dmqmc_in
        type(dmqmc_estimates_t), intent(in) :: dmqmc_estimates
        integer, intent(in) :: min_ind(:), max_ind(:)
        integer, intent(in) :: nstates_active
        real(dp), intent(in) :: nparticles(:)
        integer, intent(in) :: nspawn_events
        real(p), intent(in) :: rspawn
        logical, intent(in) :: error
        real(dp), intent(out) :: rep_loop_loc(:)

        integer :: nrdms(1)

        rep_loop_loc = 0.0_dp

        rep_loop_loc(rspawn_ind) = rspawn
        if (error) rep_loop_loc(error_ind) = 1.0_p
        rep_loop_loc(nocc_states_ind) = nstates_active
        rep_loop_loc(nspawned_ind) = nspawn_events
        rep_loop_loc(min_ind(nparticles_ind):max_ind(nparticles_ind)) = nparticles
        rep_loop_loc(min_ind(trace_ind):max_ind(trace_ind)) = dmqmc_estimates%trace
        rep_loop_loc(min_ind(operators_ind):max_ind(operators_ind)) = dmqmc_estimates%numerators
        if (dmqmc_in%calc_excit_dist) then
            rep_loop_loc(min_ind(excit_dist_ind):max_ind(excit_dist_ind)) = dmqmc_estimates%excit_dist
        end if
        if (dmqmc_in%rdm%calc_ground_rdm) then
            rep_loop_loc(min_ind(ground_rdm_trace_ind)) = dmqmc_estimates%ground_rdm%trace
        end if
        if (dmqmc_in%rdm%calc_inst_rdm) then
            ! Reshape this 2d array into a 1d array to add it to rep_loop_loc.
            nrdms = max_ind(inst_rdm_trace_ind) - min_ind(inst_rdm_trace_ind) + 1
            rep_loop_loc(min_ind(inst_rdm_trace_ind):max_ind(inst_rdm_trace_ind)) = reshape(dmqmc_estimates%inst_rdm%traces, nrdms)
        end if
        if (doing_dmqmc_calc(dmqmc_rdm_r2)) then
            rep_loop_loc(min_ind(rdm_r2_ind):max_ind(rdm_r2_ind)) = dmqmc_estimates%inst_rdm%renyi_2
        end if
        if (dmqmc_in%calc_mom_dist) then
            rep_loop_loc(min_ind(mom_dist_ind):max_ind(mom_dist_ind)) = dmqmc_estimates%mom_dist%f_k
        end if
        if (dmqmc_in%calc_struc_fac) then
            rep_loop_loc(min_ind(struc_fac_ind):max_ind(struc_fac_ind)) = dmqmc_estimates%struc_fac%f_k
        end if
        if (doing_dmqmc_calc(dmqmc_ref_proj_energy)) then
            rep_loop_loc(min_ind(ref_trace_ind):max_ind(ref_trace_ind)) = dmqmc_estimates%ref_trace
            rep_loop_loc(min_ind(ref_D0j_ind):max_ind(ref_D0j_ind)) = dmqmc_estimates%ref_D0j_particles
        end if

    end subroutine local_dmqmc_estimators

    subroutine communicated_dmqmc_estimators(dmqmc_in, rep_loop_sum, min_ind, max_ind, ncycles, tot_nparticles, qs, &
                                             dmqmc_estimates, error)

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
        ! In/Out:
        !    qs: QMC state containing quantities not specific to DMQMC.
        !    dmqmc_estimates: type containing dmqmc_estimates.
        ! Out:
        !    tot_nparticles: total number of particles of each type across all
        !       processors.
        !    error: whether an error occured on any processor.

        use calc, only: doing_dmqmc_calc, dmqmc_rdm_r2, dmqmc_ref_proj_energy
        use dmqmc_data, only: dmqmc_in_t, dmqmc_estimates_t
        use parallel, only: nprocs
        use qmc_data, only: qmc_state_t

        type(dmqmc_in_t), intent(in) :: dmqmc_in
        real(dp), intent(in) :: rep_loop_sum(:)
        integer, intent(in) :: min_ind(:), max_ind(:), ncycles
        real(dp), intent(out) :: tot_nparticles(:)
        type(qmc_state_t), intent(inout) :: qs
        type(dmqmc_estimates_t), intent(inout) :: dmqmc_estimates
        logical, intent(out) :: error

        qs%spawn_store%rspawn = real(rep_loop_sum(rspawn_ind), p)
        error = abs(rep_loop_sum(error_ind)) > depsilon
        qs%estimators%tot_nstates = nint(rep_loop_sum(nocc_states_ind))
        qs%estimators%tot_nspawn_events = nint(rep_loop_sum(nspawned_ind))
        tot_nparticles = rep_loop_sum(min_ind(nparticles_ind):max_ind(nparticles_ind))
        dmqmc_estimates%trace = real(rep_loop_sum(min_ind(trace_ind):max_ind(trace_ind)), p)
        dmqmc_estimates%numerators = real(rep_loop_sum(min_ind(operators_ind):max_ind(operators_ind)), p)
        if (dmqmc_in%calc_excit_dist) then
            dmqmc_estimates%excit_dist = real(rep_loop_sum(min_ind(excit_dist_ind):max_ind(excit_dist_ind)), p)
        end if
        if (dmqmc_in%rdm%calc_ground_rdm) then
            dmqmc_estimates%ground_rdm%trace = real(rep_loop_sum(min_ind(ground_rdm_trace_ind)), p)
        end if
        if (dmqmc_in%rdm%calc_inst_rdm) then
            dmqmc_estimates%inst_rdm%traces = real(reshape(rep_loop_sum(min_ind(inst_rdm_trace_ind):max_ind(inst_rdm_trace_ind)), &
                                                      shape(dmqmc_estimates%inst_rdm%traces)), p)
        end if
        if (doing_dmqmc_calc(dmqmc_rdm_r2)) then
            dmqmc_estimates%inst_rdm%renyi_2 = real(rep_loop_sum(min_ind(rdm_r2_ind):max_ind(rdm_r2_ind)), p)
        end if
        if (dmqmc_in%calc_mom_dist) then
            dmqmc_estimates%mom_dist%f_k = real(rep_loop_sum(min_ind(mom_dist_ind):max_ind(mom_dist_ind)), p)
        end if
        if (dmqmc_in%calc_struc_fac) then
            dmqmc_estimates%struc_fac%f_k = real(rep_loop_sum(min_ind(struc_fac_ind):max_ind(struc_fac_ind)), p)
        end if
        if (doing_dmqmc_calc(dmqmc_ref_proj_energy)) then
            dmqmc_estimates%ref_trace = real(rep_loop_sum(min_ind(ref_trace_ind):max_ind(ref_trace_ind)), p)
            dmqmc_estimates%ref_D0j_particles = real(rep_loop_sum(min_ind(ref_D0j_ind):max_ind(ref_D0j_ind)), p)
        end if
        ! Average the spawning rate.
        qs%spawn_store%rspawn = qs%spawn_store%rspawn/(ncycles*nprocs)

    end subroutine communicated_dmqmc_estimators

    subroutine communicate_inst_rdms(accumulated_probs_old, subsys_info, inst_rdms)

        ! Perform annihilation between 'RDM particles'. This is done exactly
        ! by calling the normal annihilation routine for a spawned object,
        ! done for each RDM being calculated, which communicates the RDM
        ! particles and annihilates them.

        ! Also, once annihilation has occurred, calculate all instantaneous
        ! RDM estimates.

        ! In:
        !    accumulated_probs_old: The value of accumulated_probs on the last
        !        report cycle.
        !    subsys_info: information relating to the subsystems being studied
        !        through the RDMs.
        ! In/Out:
        !    inst_rdms: The instantaneous RDM estimates which are to be
        !        communicated.

        use calc, only: doing_dmqmc_calc, dmqmc_rdm_r2
        use dmqmc_data, only: subsys_t, dmqmc_inst_rdms_t
        use hash_table, only: reset_hash_table
        use spawn_data, only: annihilate_wrapper_spawn_t

        real(p), intent(in) :: accumulated_probs_old(0:)
        type(subsys_t), intent(in) :: subsys_info(:)
        type(dmqmc_inst_rdms_t), intent(inout) :: inst_rdms

        integer :: irdm, nrdms

        nrdms = size(subsys_info)

        ! WARNING: cannot pass inst_rdm%spawn%spawn to procedures expecting an
        ! array of type spawn_t due to a bug in gfortran which results in
        ! memory deallocations!
        ! See https://groups.google.com/forum/#!topic/comp.lang.fortran/VuFvOsLs6hE
        ! and http://gcc.gnu.org/bugzilla/show_bug.cgi?id=58310.
        ! The explicit loop is also meant to be more efficient anyway, as it
        ! prevents any chance of copy-in/copy-out...
        do irdm = 1, nrdms
            call annihilate_wrapper_spawn_t(inst_rdms%spawn(irdm)%spawn, .false., .false.)
            ! Now is also a good time to reset the hash table (otherwise we
            ! attempt to lookup non-existent data in the next cycle!).
            call reset_hash_table(inst_rdms%spawn(irdm)%ht)
            ! spawn_t comms changes the memory used by spawn%sdata.  Make
            ! sure the hash table always uses the currently 'active'
            ! spawning memory.
            inst_rdms%spawn(irdm)%ht%data_label => inst_rdms%spawn(irdm)%spawn%sdata
        end do

        call calculate_rdm_traces(subsys_info, inst_rdms%spawn%spawn, inst_rdms%traces)
        if (doing_dmqmc_calc(dmqmc_rdm_r2)) call calculate_rdm_renyi_2(subsys_info, inst_rdms%spawn%spawn, &
                                                                       accumulated_probs_old, inst_rdms%renyi_2)

        do irdm = 1, nrdms
            inst_rdms%spawn(irdm)%spawn%head = inst_rdms%spawn(irdm)%spawn%head_start
        end do

    end subroutine communicate_inst_rdms

    subroutine update_shift_dmqmc(qmc_in, qs, loc_totnparticles, loc_totnparticles_old, comp)

        ! In:
        !    qmc_in: input options relating to QMC methods.
        !    loc_totnparticles: total number (across all processors) of
        !        particles in the simulation at end of the previous report loop.
        !    loc_totnparticles_old: total number (across all processors) of
        !        particles in the simulation currently.
        !    comp: if what we're doing a complex calculation.
        ! In/Out:
        !    qs: QMC state. Shift is updated.

        use energy_evaluation, only: update_shift
        use qmc_data, only: qmc_in_t, qmc_state_t

        type(qmc_in_t), intent(in) :: qmc_in
        type(qmc_state_t), intent(inout) :: qs
        real(dp), intent(in) :: loc_totnparticles(:)
        real(dp), intent(in) :: loc_totnparticles_old(:)
        logical, intent(in), optional :: comp
        logical :: comp_set
        integer :: ireplica

        comp_set = .false.
        if (present(comp)) comp_set = comp

        do ireplica = 1, size(loc_totnparticles)
            if (.not.comp_set .or. mod(ireplica, 2) /= 0) then
                if (qs%vary_shift(ireplica)) then
                    call update_shift(qs, qs%shift(ireplica), loc_totnparticles_old(ireplica), &
                                     &loc_totnparticles(ireplica), qmc_in%ncycles)
                end if
                if (loc_totnparticles(ireplica) > qs%target_particles .and. (.not. qs%vary_shift(ireplica))) &
                    qs%vary_shift(ireplica) = .true.
            else
                ! imaginary parts only needs to follow its real counterpart
                qs%shift(ireplica) = qs%shift(ireplica-1)
                qs%vary_shift(ireplica) = qs%vary_shift(ireplica-1)
            end if
        end do

    end subroutine update_shift_dmqmc

    subroutine update_dmqmc_estimators(sys, dmqmc_in, idet, iteration, cdet, H00, psip_list, &
                                       dmqmc_estimates, weighted_sampling, rdm_error, ref_f0)

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
        !    ref_f0: the reference bitstring of the system.
        ! In/Out:
        !    cdet: det_info_t object containing information of current density
        !        matrix element.
        !    dmqmc_estimates: type containing dmqmc estimates.
        !    weighted_sampling: type containing weighted sampling parameters.
        !    rdm_error: true if an error has occured and we need to quit at the
        !        end of the report loop.

        use calc, only: doing_dmqmc_calc, dmqmc_energy, dmqmc_staggered_magnetisation
        use calc, only: dmqmc_energy_squared, dmqmc_correlation, dmqmc_full_r2, dmqmc_kinetic_energy
        use calc, only: dmqmc_H0_energy, dmqmc_potential_energy, dmqmc_HI_energy
        use calc, only: dmqmc_ref_proj_energy
        use excitations, only: get_excitation, excit_t
        use proc_pointers, only:  update_dmqmc_energy_and_trace_ptr, update_dmqmc_stag_mag_ptr
        use proc_pointers, only: update_dmqmc_energy_squared_ptr, update_dmqmc_correlation_ptr
        use proc_pointers, only: update_dmqmc_kinetic_energy_ptr
        use determinant_data, only: det_info_t
        use system, only: sys_t
        use qmc_data, only: particle_t
        use dmqmc_data, only: dmqmc_in_t, dmqmc_estimates_t, energy_ind, energy_imag_ind, energy_squared_ind, &
                              correlation_fn_ind, staggered_mag_ind, full_r2_ind, full_r2_imag_ind, dmqmc_weighted_sampling_t, &
                              kinetic_ind, potential_ind, H0_ind, H0_imag_ind, HI_ind, HI_imag_ind, ref_proj_ind, ref_proj_imag_ind

        type(sys_t), intent(in) :: sys
        type(dmqmc_in_t), intent(in) :: dmqmc_in
        integer, intent(in) :: idet, iteration
        integer(i0), intent(in) :: ref_f0(sys%basis%tot_string_len)
        type(det_info_t), intent(inout) :: cdet
        real(p), intent(in) :: H00
        type(particle_t), intent(in) :: psip_list
        type(dmqmc_estimates_t), intent(inout) :: dmqmc_estimates
        type(dmqmc_weighted_sampling_t), intent(inout) :: weighted_sampling
        logical, intent(inout) :: rdm_error

        logical :: nonzero_pop
        type(excit_t) :: excitation
        real(p) :: unweighted_walker_pop(psip_list%nspaces)

        ! Get excitation.
        excitation = get_excitation(sys%nel, sys%basis, psip_list%states(:sys%basis%tot_string_len,idet), &
                         psip_list%states((1+sys%basis%tot_string_len):sys%basis%tensor_label_len,idet))

        ! When performing importance sampling the result is that certain
        ! excitation levels have smaller psips populations than the true density
        ! matrix by some factor. In these cases, we want to multiply the psip
        ! population by this factor to calculate the contribution from these
        ! excitation levels correctly.

        ! In the case of no importance sampling, unweighted_walker_pop = particle_t%pops(1,idet).
        unweighted_walker_pop = real(psip_list%pops(:,idet),p)*weighted_sampling%probs(excitation%nexcit)/psip_list%pop_real_factor

        ! The following only use the populations with ireplica = 1, so only call
        ! them if the determinant is occupied in the first replica.
        associate (est => dmqmc_estimates, comp => sys%read_in%comp)

            ! Check if we have a non-zero population.
            nonzero_pop = abs(unweighted_walker_pop(1)) > 0
            if (comp) nonzero_pop = nonzero_pop .or. abs(unweighted_walker_pop(2)) > 0

            if (nonzero_pop) then
                ! See which estimators are to be calculated, and call the
                ! corresponding procedures.
                ! Energy.
                if (doing_dmqmc_calc(dmqmc_energy)) &
                    call update_dmqmc_energy_and_trace_ptr(sys, excitation, cdet, H00, unweighted_walker_pop, &
                                                           psip_list%dat(1, idet), est%trace, &
                                                           est%numerators(energy_ind:energy_imag_ind), comp)
                ! Energy squared.
                if (doing_dmqmc_calc(dmqmc_energy_squared)) &
                    call update_dmqmc_energy_squared_ptr(sys, cdet, excitation, H00, unweighted_walker_pop(1), &
                                                         est%numerators(energy_squared_ind))
                ! Spin-spin correlation function.
                if (doing_dmqmc_calc(dmqmc_correlation)) &
                    call update_dmqmc_correlation_ptr(sys, cdet, excitation, H00, unweighted_walker_pop(1), &
                                                      est%correlation_mask, est%numerators(correlation_fn_ind))
                ! Staggered magnetisation.
                if (doing_dmqmc_calc(dmqmc_staggered_magnetisation)) &
                    call update_dmqmc_stag_mag_ptr(sys, cdet, excitation, H00, unweighted_walker_pop(1), &
                                                   est%numerators(staggered_mag_ind))
                ! Kinetic energy.
                if (doing_dmqmc_calc(dmqmc_kinetic_energy)) &
                    call update_dmqmc_kinetic_energy_ptr(sys, cdet, excitation, H00, unweighted_walker_pop(1), &
                                                         est%numerators(kinetic_ind))
                ! Potential energy.
                if (doing_dmqmc_calc(dmqmc_potential_energy)) &
                    call update_dmqmc_potential_energy(sys, cdet, excitation, unweighted_walker_pop(1), &
                                                       est%numerators(potential_ind))
                ! Stores the projected energy from only considering the
                ! reference row of the density matrix.
                if (doing_dmqmc_calc(dmqmc_ref_proj_energy)) &
                    call update_dmqmc_ref_proj_energy(sys, excitation, cdet, H00, unweighted_walker_pop, &
                                                      psip_list%dat(1, idet), est%ref_trace, est%ref_D0j_particles, &
                                                      est%numerators(ref_proj_ind:ref_proj_imag_ind), ref_f0, comp)
                ! H^0 energy, where H^0 = H - V. See subroutines interface
                ! comments for description.
                if (doing_dmqmc_calc(dmqmc_H0_energy)) &
                    call update_dmqmc_H0_energy(sys, cdet, excitation, unweighted_walker_pop, &
                                                est%numerators(H0_ind:H0_imag_ind), comp)
                ! HI energy, HI(tau-beta) = e^{-0.5(beta-tau)H^0} H e^{0.5(beta-tau)H^0}
                if (doing_dmqmc_calc(dmqmc_HI_energy)) &
                    call update_dmqmc_HI_energy(sys, cdet, excitation, unweighted_walker_pop, &
                                                weighted_sampling%probs(sys%max_number_excitations+1), &
                                                est%numerators(HI_ind:HI_imag_ind), comp)
                ! Momentum distribution.
                if (dmqmc_in%calc_mom_dist) &
                    call update_dmqmc_momentum_distribution(sys, cdet, excitation, H00, unweighted_walker_pop(1), &
                                                            est%mom_dist%f_k)
                ! UEG structure factor.
                if (dmqmc_in%calc_struc_fac) &
                    call update_dmqmc_structure_factor_ueg(sys, cdet, excitation, H00, unweighted_walker_pop(1), &
                                                           est%struc_fac%f_k)
                ! Excitation distribution. (Should also be evaluated if we're 
                ! calculating importance sampling weights.)
                if (dmqmc_in%calc_excit_dist .or. (dmqmc_in%find_weights .and. iteration > dmqmc_in%find_weights_start)) then
                    if (.not.comp) then
                        est%excit_dist(excitation%nexcit) = est%excit_dist(excitation%nexcit) + &
                            &real(abs(psip_list%pops(1,idet)),p)/psip_list%pop_real_factor
                    else 
                        ! Now we just sum the particles from the real and 
                        ! imaginary components for complex systems.
                        est%excit_dist(excitation%nexcit) = est%excit_dist(excitation%nexcit) + &
                            &sum(real(abs(psip_list%pops(1:2,idet)),p))/psip_list%pop_real_factor
                    end if
                end if
            end if

            ! Full Renyi entropy (S_2).
            if (doing_dmqmc_calc(dmqmc_full_r2)) &
                call update_full_renyi_2(unweighted_walker_pop, excitation%nexcit, dmqmc_in%half_density_matrix, &
                                         est%numerators(full_r2_ind:full_r2_imag_ind), complx=comp)

            ! Update the contribution to the trace from other replicas
            if (dmqmc_in%replica_tricks .and. excitation%nexcit == 0) then
                if (.not.comp) then
                    est%trace(2) = est%trace(2) + unweighted_walker_pop(2)
                else
                    est%trace(3:4) = est%trace(3:4) + unweighted_walker_pop(3:4)
                end if
            end if

            ! Reduced density matrices.
            if (dmqmc_in%rdm%doing_rdm) then
                call update_reduced_density_matrix_heisenberg(sys%basis, est, dmqmc_in%rdm, cdet, excitation, &
                                                              psip_list%pops(:,idet), psip_list%pop_real_factor, iteration, &
                                                              dmqmc_in%start_av_rdm, weighted_sampling%probs, rdm_error)
            end if

            weighted_sampling%probs_old = weighted_sampling%probs

        end associate

    end subroutine update_dmqmc_estimators

    subroutine dmqmc_energy_and_trace(sys, excitation, cdet, H00, pop, diagonal_contribution, trace, energy, &
                                     &complx)

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
        !    complex: if the system is complex-valued.
        ! In/Out:
        !    trace: total population on diagonal elements of density matrix
        !    energy: current thermal energy estimate.

        use determinant_data, only: det_info_t
        use system, only: sys_t
        use excitations, only: excit_t
        use hamiltonian_data

        type(sys_t), intent(in) :: sys
        type(excit_t), intent(inout) :: excitation
        type(det_info_t), intent(in) :: cdet
        real(p), intent(in) :: H00, diagonal_contribution
        real(p), intent(in) :: pop(:)
        logical, intent(in), optional :: complx
        real(p), intent(inout) :: trace(:)
        real(p), intent(inout) :: energy(:)

        logical :: complx_set
        type(hmatel_t) :: hmatel
        ! Importance sampling (in the FCIQMC-sense) isn't used in DMQMC...
        real(p) :: trial_wfn_dat(0)

        complx_set = .false.
        if (present(complx)) complx_set = complx

        ! Update trace and off-diagonal contributions to the total enegy
        call update_proj_energy_dmqmc(sys, cdet%f2, trial_wfn_dat, cdet, pop, trace, energy, excitation, hmatel, complx)

        ! Update diagaonal contribution to the total energy. 
        ! Only the first replica is considered.
        ! NB The diagonal elements of Hermitian operators (namely, all observables)
        !    should be all real (cf. read_in.f90 and hamiltonian_periodic_complex.f90)
        if (excitation%nexcit == 0) then 
            if (.not.complx_set) then
                energy(1) = energy(1) + (diagonal_contribution+H00)*pop(1)
            else
                energy(1:2) = energy(1:2) + (diagonal_contribution+H00)*pop(1:2)
            end if
        end if

    end subroutine dmqmc_energy_and_trace

    subroutine dmqmc_energy_and_trace_propagate(sys, excitation, cdet, H00, pop, diagonal_contribution, trace, energy, &
                                               &complx)

        ! Add the contribution for the current density matrix element to the thermal
        ! energy estimate. Routine is specific to when using ipdmqmc
        ! option.

        ! In:
        !    sys: system being studied.
        !    excitation: excit_t type variable which stores information on
        !        the excitation between the two bitstring ends, corresponding
        !        to the two labels for the density matrix element.
        !    H00: diagonal Hamiltonian element for the reference.
        !    pop: number of particles on the current density matrix
        !        element.
        !    pop_im: imaginary population, not used here.
        !    cdet: det_info_t object containing bit strings of densitry matrix
        !        element under consideration.
        !    diagonal_contribution: <D_i|H|D_i>-<D0|H|D0>
        ! In/Out:
        !    trace: total population on diagonal elements of density matrix
        !    energy: current thermal energy estimate.

        use determinant_data, only: det_info_t
        use system, only: sys_t
        use excitations, only: excit_t
        use proc_pointers, only: sc0_ptr
        use hamiltonian_data

        type(sys_t), intent(in) :: sys
        type(excit_t), intent(inout) :: excitation
        type(det_info_t), intent(in) :: cdet
        real(p), intent(in) :: H00, diagonal_contribution
        real(p), intent(in) :: pop(:)
        logical, intent(in), optional :: complx
        real(p), intent(inout) :: trace(:)
        real(p), intent(inout) :: energy(:)

        logical :: complx_set
        type(hmatel_t) :: hmatel
        ! Importance sampling (in the FCIQMC-sense) isn't used in DMQMC...
        real(p) :: trial_wfn_dat(0)

        complx_set = .false.
        if (present(complx)) complx_set = complx

        ! Update trace and off-diagonal contributions to the total enegy.
        call update_proj_energy_dmqmc(sys, cdet%f2, trial_wfn_dat, cdet, pop, trace, energy, excitation, hmatel, complx)

        ! Update diagaonal contribution to the total energy.
        if (excitation%nexcit == 0) then 
            if (.not.complx_set) then
                energy(1) = energy(1) + sc0_ptr(sys, cdet%f)*pop(1)
            else
                energy(1:2) = energy(1:2) + sc0_ptr(sys, cdet%f)*pop(1:2)
            end if
        end if

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
        use determinant_data, only: det_info_t
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
        use determinant_data, only: det_info_t
        use excitations, only: excit_t
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        type(det_info_t), intent(in) :: cdet
        type(excit_t), intent(in) :: excitation
        real(p), intent(in) :: H00, walker_pop
        integer(i0), allocatable, intent(in) :: correlation_mask(:)
        real(p), intent(inout) :: correlation_fn

        integer(i0) :: f(sys%basis%tot_string_len)
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
            f = iand(cdet%f(:sys%basis%tot_string_len), correlation_mask)
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
        use determinant_data, only: det_info_t
        use excitations, only: excit_t
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        type(det_info_t), intent(in) :: cdet
        type(excit_t), intent(in) :: excitation
        real(p), intent(in) :: H00, walker_pop
        real(p), intent(inout) :: staggered_mag

        integer :: bit_element1, bit_position1, bit_element2, bit_position2
        integer(i0) :: f(sys%basis%tot_string_len)
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
            f = iand(cdet%f(:sys%basis%tot_string_len), sys%heisenberg%lattice_mask)
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

        staggered_mag = staggered_mag + (real(total_sum,p)/real(sys%lattice%nsites**2,p))*walker_pop

    end subroutine dmqmc_stag_mag_heisenberg

    subroutine update_full_renyi_2(walker_pop, excit_level, half_density_matrix, full_r2, complx)

        ! Add the contribution from the current density matrix element to the
        ! Renyi entropy (S_2) of the full density matrix.

        ! In:
        !    walker_pop: number of particles on the current density matrix
        !        element, for both replicas.
        !    excit_level: The excitation level between the two bitstrings
        !        contributing to the full density matrix bitstring.
        !    half_density_matrix: reflect psips spawned from lower triangle into
        !        the upper one?
        !    complx: if we are sampling a complex read_in system.
        ! In/Out:
        !    full_r2: running estimate for Renyi entropy of the full density
        !       matrix.

        real(p), intent(in) :: walker_pop(:)
        integer, intent(in) :: excit_level
        logical, intent(in) :: half_density_matrix
        logical, intent(in), optional :: complx
        real(p), intent(inout) :: full_r2(:)
        real(p) :: delta_r2(size(full_r2))

        delta_r2 = 0
        if (present(complx) .and. complx) then
            delta_r2(1) = walker_pop(1)*walker_pop(3) - walker_pop(2)*walker_pop(4)
            delta_r2(2) = walker_pop(1)*walker_pop(4) + walker_pop(2)*walker_pop(3)
        else
            delta_r2(1) = walker_pop(1)*walker_pop(2)
        end if

        if (half_density_matrix .and. excit_level /= 0) then
            ! With the half-density matrix option, only the upper-half of the
            ! density matrix is stored, but the off-diagonal elements are twice
            ! as large instead. Thus, a product of off-diagonal elements is four
            ! times as large. We want two 'correct size' contributions, so we
            ! need to divide these contirbutions by two to get this.
            full_r2 = full_r2 + delta_r2/2.0_p
        else
            full_r2 = full_r2 + delta_r2
        end if

    end subroutine update_full_renyi_2

    subroutine update_reduced_density_matrix_heisenberg(basis, dmqmc_estimates, rdm_in, cdet, excitation, walker_pop, &
                                                        pop_real_factor, iteration, start_av_rdm, accumulated_probs, rdm_error)

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
        !    rdm_in: input options relating to reduced density matrices.
        !    cdet: information on density matrix element, in particular cdet%f,
        !       the bit strings of each label.
        !    excitation: excit_t type variable which stores information on
        !        the excitation between the two bitstring ends, corresponding to
        !        the two labels for the density matrix element.
        !    walker_pop: number of particles on the current density matrix
        !        element. Note that this walker population is still weighted
        !        by the importance sampling factors. These factors must be
        !        removed before any estimates can be calculated.
        !    pop_real_factor: The factor by which populations are multiplied to
        !        enable non-integer populations.
        !    iteration: interation number.  No accumulation of the RDM is
        !        performed if iteration <= start_av_rdm.
        !    start_av_rdm: iteration we start averaging the rdm on.
        !    accumulated_probs: factors by which the population on each
        !        excitation level are reduced.
        ! In/Out:
        !    dmqmc_estimates: type containing dmqmc estimates.
        !    rdm_error: true on output if we run of memory in an RDM spawning
        !        array. True on input if an this already happened before this
        !        call.

        use basis_types, only: basis_t
        use determinant_data, only: det_info_t
        use dmqmc_data, only: dmqmc_rdm_in_t, dmqmc_estimates_t
        use dmqmc_procedures, only: decode_dm_bitstring
        use excitations, only: excit_t
        use spawning, only: create_spawned_particle_rdm

        type(basis_t), intent(in) :: basis
        type(dmqmc_estimates_t), intent(inout) :: dmqmc_estimates
        type(dmqmc_rdm_in_t), intent(in) :: rdm_in
        type(det_info_t), intent(in) :: cdet
        integer, intent(in) :: iteration
        integer(int_p), intent(in) :: walker_pop(:)
        integer(int_p), intent(in) :: pop_real_factor
        type(excit_t), intent(in) :: excitation
        integer, intent(in) :: start_av_rdm
        real(p), intent(in) :: accumulated_probs(0:)
        logical, intent(inout) :: rdm_error

        real(p) :: unweighted_walker_pop(size(walker_pop))
        integer :: irdm, isym, ireplica, nrdms, nsym_vecs
        integer(i0) :: f1(basis%tot_string_len), f2(basis%tot_string_len)
        integer(i0) :: f3(basis%tensor_label_len)
        integer(i0) :: rdm_f1(basis%tot_string_len), rdm_f2(basis%tot_string_len)

        if (.not. (iteration > start_av_rdm .or. rdm_in%calc_inst_rdm)) return

        ! Combined bitstring.
        f3(1:basis%tot_string_len) = cdet%f(:basis%tot_string_len)
        f3(basis%tot_string_len+1:) = cdet%f2(:basis%tot_string_len)

        nrdms = size(dmqmc_estimates%subsys_info)
        nsym_vecs = size(dmqmc_estimates%subsys_info(1)%B_masks, 2)

        ! Loop over all RDMs to be calculated.
        do irdm = 1, nrdms
        ! Loop over every symmetry-equivalent subsystem for this RDM.
        do isym = 1, nsym_vecs

        rdm_f1 = 0_i0
        rdm_f2 = 0_i0

        associate(rdm=>dmqmc_estimates%ground_rdm%rdm, end1=>rdm_f1(1), end2=>rdm_f2(1), &
                  subsys_info=>dmqmc_estimates%subsys_info)

        ! Apply the mask for the B subsystem to set all sites in the A
        ! subsystem to 0.
        f1 = iand(subsys_info(irdm)%B_masks(:,isym),cdet%f(:basis%tot_string_len))
        f2 = iand(subsys_info(irdm)%B_masks(:,isym),cdet%f2(:basis%tot_string_len))

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
            call decode_dm_bitstring(basis, f3, isym, subsys_info(irdm), &
                                     rdm_f1(1:subsys_info(irdm)%string_len), rdm_f2(1:subsys_info(irdm)%string_len))

            if (rdm_in%calc_ground_rdm) then
                ! The above routine actually maps to numbers between 0
                ! and 2^subsys_info(1)%A_nsites-1, but the smallest and largest
                ! reduced density matrix indices are one more than these,
                ! so add one.
                rdm_f1 = rdm_f1 + 1
                rdm_f2 = rdm_f2 + 1
                unweighted_walker_pop = real(walker_pop,p)*&
                    accumulated_probs(excitation%nexcit)/pop_real_factor
                ! Note, when storing the entire RDM (as done here), the
                ! maximum value of subsys_info(i)%tot_string_len is 1, so we
                ! only consider this one element here.
                rdm(end1, end2) = rdm(end1, end2) + unweighted_walker_pop(1)
            end if

            if (rdm_in%calc_inst_rdm) then
                do ireplica = 1, size(walker_pop)
                    if (abs(walker_pop(ireplica)) > 0) then
                        call create_spawned_particle_rdm(rdm_f1, rdm_f2, walker_pop(ireplica), ireplica, &
                                                         dmqmc_estimates%inst_rdm%spawn(irdm))
                        ! Did we run out of the space in the RDM spawning array.
                        rdm_error = rdm_error .or. dmqmc_estimates%inst_rdm%spawn(irdm)%spawn%error
                    end if
                end do
            end if

        end if
        end associate
        end do

        end do

    end subroutine update_reduced_density_matrix_heisenberg

    subroutine call_ground_rdm_procedures(dmqmc_estimates, beta_cycle, rdm_in)

        ! Wrapper for calling ground-state RDM procedures (*not*
        ! beta-dependent RDMs).

        ! In:
        !    beta_cycle: index of the beta loop being performed.
        !    rdm_in: input options relating to reduced density matrices.
        ! In/Out:
        !    dmqmc_estimates: type containing dmqmc estimates.

        use checking, only: check_allocate, check_deallocate
        use dmqmc_data, only: dmqmc_rdm_in_t, dmqmc_estimates_t
        use parallel
        use utils, only: append_ext, int_fmt

        type(dmqmc_estimates_t), intent(inout) :: dmqmc_estimates
        integer, intent(in) :: beta_cycle
        type(dmqmc_rdm_in_t), intent(in) :: rdm_in

        integer :: i, j, new_unit, iunit
        character(10) :: rdm_filename
        ! If in parallel then merge the reduced density matrix onto one
        ! processor.
#ifdef PARALLEL

        real(p), allocatable :: dm(:,:)
        real(p), allocatable :: dm_sum(:,:)
        integer :: num_eigv, ierr

        num_eigv = 2**dmqmc_estimates%subsys_info(1)%A_nsites

        allocate(dm(num_eigv,num_eigv), stat=ierr)
        call check_allocate('dm',num_eigv**2,ierr)
        allocate(dm_sum(num_eigv,num_eigv), stat=ierr)
        call check_allocate('dm_sum',num_eigv**2,ierr)

        dm = dmqmc_estimates%ground_rdm%rdm
        call mpi_allreduce(dm, dm_sum, size(dm), MPI_PREAL, MPI_SUM, MPI_COMM_WORLD, ierr)
        dmqmc_estimates%ground_rdm%rdm = dm_sum

        deallocate(dm)
        call check_deallocate('dm',ierr)
        deallocate(dm_sum)
        call check_deallocate('dm_sum',ierr)

#endif
        iunit = 6

        associate(rdm=>dmqmc_estimates%ground_rdm%rdm, trace=>dmqmc_estimates%ground_rdm%trace)

            trace = 0.0_p

            if (parent) then
                ! Force the reduced density matrix to be symmetric by averaging the
                ! upper and lower triangles.
                do i = 1, ubound(rdm,1)
                    do j = 1, i-1
                        rdm(i,j) = 0.5_p*(rdm(i,j) + rdm(j,i))
                        rdm(j,i) = rdm(i,j)
                    end do
                    ! Add current contribution to the trace.
                    trace = trace + rdm(i,i)
                end do

                ! Call the routines to calculate the desired quantities.
                if (rdm_in%doing_vn_entropy) call calculate_vn_entropy(rdm, dmqmc_estimates%subsys_info)
                if (rdm_in%doing_concurrence) call calculate_concurrence(rdm)

                write (iunit,'(1x,"# RDM trace =",1X,es17.10)') trace

                if (rdm_in%output_rdm) then
                    call append_ext('rdm', beta_cycle, rdm_filename)
                    open(newunit=new_unit, file=trim(rdm_filename), status='replace')
                    write(new_unit,'(a5,1x,es15.8)') "Trace", trace
                    do i = 1, ubound(rdm,1)
                        do j = i, ubound(rdm,1)
                            write(new_unit,'(a1,'//int_fmt(i,0)//',a1,'//int_fmt(j,0)//',a1,1x,es15.8)') &
                                "(",i,",",j,")", rdm(i,j)
                        end do
                    end do
                    close(new_unit)
                end if

            end if

        end associate

    end subroutine call_ground_rdm_procedures

    subroutine calculate_vn_entropy(rdm, subsys_info)

        ! Calculate the Von Neumann Entropy. Use lapack to calculate the
        ! eigenvalues {\lambda_j} of the reduced density matrix.
        ! Then VN Entropy S = -\sum_j\lambda_j\log_2{\lambda_j}.

        ! Need to paralellise for large subsystems and introduce test to check
        ! whether diagonalisation should be performed in serial or paralell.

        ! In:
        !    subsys_info: information relating to the subsystems being studied
        !        through the RDMs.
        ! In/Out:
        !     rdm: Reduced density matrix estimate.

        use checking, only: check_allocate, check_deallocate
        use linalg, only: syev_wrapper
        use dmqmc_data, only: subsys_t

        real(p), intent(inout) :: rdm(:,:)
        type(subsys_t) :: subsys_info(:)

        integer :: i, rdm_size
        integer :: info, ierr
        real(p), allocatable :: dm_tmp(:,:)
        real(p) :: eigv(2**subsys_info(1)%A_nsites)
        real(p) :: vn_entropy
        logical :: thrown_away
        integer :: iunit

        iunit = 6

        rdm_size = 2**subsys_info(1)%A_nsites
        vn_entropy = 0.0_p

        ! The matrix input into the following diagonalisation routines will have
        ! their upper half (including the diagonal) destroyed. We might want
        ! reduced_desntiy_matrix later, so use some temporary space:
        allocate(dm_tmp(rdm_size,rdm_size), stat=ierr)
        call check_allocate('dm_tmp',rdm_size**2,ierr)
        dm_tmp = rdm

        call syev_wrapper('N', 'U', rdm_size, dm_tmp, rdm_size, eigv, info)
        thrown_away = .false.
        write(iunit,'(1X,"# Eigenvalues thrown away:",1X)',advance='no')
        do i = 1, ubound(eigv,1)
            if (eigv(i) < depsilon) then
                write(iunit,'(es15.8,2x)',advance='no') eigv(i)
                thrown_away = .true.
                cycle
            end if
            vn_entropy = vn_entropy - eigv(i)*(log(eigv(i))/log(2.0_p))
        end do
        if (thrown_away) then
            write(iunit,'()',advance='yes')
        else
            write(iunit,'(1X,"none")',advance='yes')
        end if
        write (iunit,'(1x,"# Unnormalised von Neumann entropy =",1X,es17.10)') vn_entropy

        deallocate(dm_tmp)
        call check_deallocate('dm_tmp',ierr)

    end subroutine calculate_vn_entropy

    subroutine calculate_concurrence(rdm)

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

        ! In:
        !     rdm: Reduced density matrix estimate.

        use linalg, only: geev_wrapper

        real(p), intent(in) :: rdm(:,:)

        integer :: info
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
        integer :: iunit
        iunit = 6

        ! Make rdm_spin_flip_tmp because sgeev and dgeev delete input matrix.
        rdm_spin_flip_tmp = matmul(rdm, flip_spin_matrix)
        rdm_spin_flip = rdm_spin_flip_tmp

        call geev_wrapper('N', 'N', 4, rdm_spin_flip, 4, reigv, ieigv, VL, 1, VR, 1, info)
        ! Calculate the concurrence. Take abs of eigenvalues so that this is
        ! equivalant to sqauring and then square-rooting.
        concurrence = 2.0_p*maxval(abs(reigv)) - sum(abs(reigv))
        concurrence = max(0.0_p, concurrence)
        write (iunit,'(1x,"# Unnormalised concurrence =",1X,es17.10)') concurrence

    end subroutine calculate_concurrence

    subroutine calculate_rdm_traces(subsys_info, rdm_lists, traces)

        ! In:
        !    subsys_info: Array of subsys_t derived types, holding information
        !        about the various subsystems for which RDMs are being
        !        estimated.
        !    rdm_lists: Array of rdm_spawn_t derived types, which hold all of
        !        the RDM psips which belong to this processor.
        ! Out:
        !    traces: The calculated RDM traces.

        use dmqmc_data, only: subsys_t
        use excitations, only: get_excitation_level
        use spawn_data, only: spawn_t

        type(subsys_t), intent(in) :: subsys_info(:)
        type(spawn_t), intent(in) :: rdm_lists(:)
        real(p), intent(out) :: traces(:,:)
        integer :: irdm, i, rdm_bl
        integer, parameter :: thread_id = 0

        traces = 0.0_p

        ! Loop over all RDMs being calculated.
        do irdm = 1, size(subsys_info)
            rdm_bl = subsys_info(irdm)%string_len
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

    subroutine calculate_rdm_renyi_2(subsys_info, rdm_lists, accumulated_probs_old, r2)

        ! Calculate the Renyi entropy (S_2) for all instantaneous RDMs being
        ! calculated.

        ! In:
        !    subsys_info: Array of subsys_t derived types, holding information
        !        about the various subsystems for which RDMs are being
        !        estimated.
        !    rdm_lists: Array of rdm_spawn_t derived types, which hold all of
        !        the RDM psips which belong to this processor.
        !    accumulated_probs_old: value of accumulated probs on the last
        !        report cycle.
        ! Out:
        !    r2: The calculated Renyi entropies (S_2).

        use dmqmc_data, only: subsys_t
        use excitations, only: get_excitation_level
        use spawn_data, only: spawn_t

        type(subsys_t), intent(in) :: subsys_info(:)
        type(spawn_t), intent(in) :: rdm_lists(:)
        real(p), intent(in) :: accumulated_probs_old(0:)
        real(p), intent(out) :: r2(:)

        integer :: i, irdm, excit_level, rdm_bl
        real(p) :: unweighted_pop_1, unweighted_pop_2
        integer, parameter :: thread_id = 0

        r2 = 0.0_p

        ! Loop over all RDMs being calculated.
        do irdm = 1, size(subsys_info)
            rdm_bl = subsys_info(irdm)%string_len
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

    subroutine dmqmc_kinetic_energy_diag(sys, cdet, excitation, H00, pop, kinetic_energy)

        ! Add the contribution for the current density matrix element to the thermal
        ! kinetic energy estimate.

        ! In:
        !    sys: system being studied.
        !    cdet: det_info_t object containing bit strings of densitry matrix
        !       element under consideration.
        !    excitation: excit_t type variable which stores information on
        !        the excitation between the two bitstring ends, corresponding
        !        to the two labels for the density matrix element.
        !    H00: diagonal hamiltonian element for the reference. only for
        !       interface consistency, not used.
        !    pop: number of particles on the current density matrix
        !        element.
        ! In/Out:
        !    kinetic_energy: current thermal kinetic energy estimate.

        use determinant_data, only: det_info_t
        use system, only: sys_t
        use excitations, only: excit_t
        use proc_pointers, only: kinetic_diag_ptr

        type(sys_t), intent(in) :: sys
        type(det_info_t), intent(in) :: cdet
        type(excit_t), intent(in) :: excitation
        real(p), intent(in) :: H00, pop
        real(p), intent(inout) :: kinetic_energy

        if (excitation%nexcit == 0) kinetic_energy = kinetic_energy + pop*kinetic_diag_ptr(sys, cdet%f)

    end subroutine dmqmc_kinetic_energy_diag

    subroutine update_dmqmc_momentum_distribution(sys, cdet, excitation, H00, pop, momentum_dist)

        ! Add the contribution for the current density matrix element to the thermal
        ! momentum distribution estimate.

        ! In:
        !    sys: system being studied.
        !    cdet: det_info_t object containing bit strings of densitry matrix
        !       element under consideration.
        !    excitation: excit_t type variable which stores information on
        !        the excitation between the two bitstring ends, corresponding
        !        to the two labels for the density matrix element.
        !    H00: diagonal hamiltonian element for the reference. only for
        !       interface consistency, not used.
        !    pop: number of particles on the current density matrix
        !        element.
        ! In/Out:
        !    momentum_dist: array expressing the number operator numerator
        !        for each k point


        use determinant_data, only: det_info_t
        use system, only: sys_t
        use excitations, only: excit_t

        type(sys_t), intent(in) :: sys
        type(det_info_t), intent(in) :: cdet
        type(excit_t), intent(in) :: excitation
        real(p), intent(in) :: H00, pop
        real(p), intent(inout) :: momentum_dist(:)
        integer :: iocc, oc_orb, pol

        pol = (sys%nalpha-sys%nbeta) / sys%nel

        if (excitation%nexcit == 0) then
            do iocc = 1, sys%nel
                ! Really want spin averaged momentum distribution, so only consider up index.
                oc_orb = cdet%occ_list(iocc)/2 + mod(cdet%occ_list(iocc),2)
                if (oc_orb <= size(momentum_dist)) then
                    ! Divide by two for unpolarised system.
                    momentum_dist(oc_orb) = momentum_dist(oc_orb) + pop / (2-pol)
                end if
            end do
        endif

    end subroutine update_dmqmc_momentum_distribution

    subroutine update_dmqmc_structure_factor_ueg(sys, cdet, excitation, H00, pop, structure_factor)

        ! Add the contribution for the current density matrix element to the thermal
        ! structure factor estimate.

        ! The structure factor is given as:
        !   S_{sigma,sigma'}(q) = N_{sigma}delta_{sigma,sigma'}/N +
        !   1/N sum_{kp} <c^{dagger}_{k+q,sigma}c^{dagger}_{p-q,sigma'}c_{p,sigma'}c_{k,sigma}>
        ! This is closely related to the the evaluation of the potential energy as
        !     V = -N/(2V) \sum_{q} v_q (S(q)-1)
        ! Note the constant term and the factors of 1/N are accounted for during output and
        ! in post-processing.

        ! In:
        !    sys: system being studied.
        !    cdet: det_info_t object containing bit strings of densitry matrix
        !       element under consideration.
        !    excitation: excit_t type variable which stores information on
        !        the excitation between the two bitstring ends, corresponding
        !        to the two labels for the density matrix element.
        !    H00: diagonal hamiltonian element for the reference. only for
        !       interface consistency, not used.
        !    pop: number of particles on the current density matrix
        !        element.
        ! In/Out:
        !    structure_factor : array expressing the numerator for the estimator of the
        !        structure factor (up to a few constants (see above)) for each k point.

        use determinant_data, only: det_info_t
        use system, only: sys_t
        use excitations, only: excit_t

        type(sys_t), intent(in) :: sys
        type(det_info_t), intent(in) :: cdet
        type(excit_t), intent(in) :: excitation
        real(p), intent(in) :: H00, pop
        real(p), intent(inout) :: structure_factor(:)

        integer :: i, j, a, b, indx, sign_factor
        integer :: q(sys%lattice%ndim)
        logical :: in_basis

        ! Due to momentum conservation only need to consider diagonal contributions
        ! and those from double excitations.
        ! We store S_{up,up}(q) + S_{down,down} in odd indices, and
        ! S_{up,down}(q) + S_{down,up} in even indices i.e.,
        ! S = (S_{uu}(k1)+S_{dd}(k1), S_{ud}(k1) + S_{du}(k1),...)
        ! The index of S_{uu}+S_{dd} is ueg_basis%lookup(q) + 0,
        ! while for S_{ud} + S_{du} is ueg_basis%lookup(q) + 1
        ! Note ueg_basis%lookup only returns the up spin index.
        associate(ueg_basis => sys%ueg%basis, bs => sys%basis)
            select case(excitation%nexcit)
            case(0)
                ! Hartree--Fock-like contribution to the structure factor.
                do i = 1, sys%nel
                    do j = i+1, sys%nel
                        if (mod(cdet%occ_list(i),2) == mod(cdet%occ_list(j),2)) then
                            q = bs%basis_fns(cdet%occ_list(i))%l - bs%basis_fns(cdet%occ_list(j))%l
                            call q_in_basis(ueg_basis, q, in_basis, indx)
                            if (in_basis) then
                                ! The factor of two accounts for only summing j > i.
                                structure_factor(indx) = structure_factor(indx) - 2.0_p*pop
                            end if
                        end if
                    end do
                end do
            case(2)
                ! This is an allowed excitation by virtue of it exisiting as a density
                ! matrix element so don't need to check for momentum conservation.
                i = excitation%from_orb(1)
                j = excitation%from_orb(2)
                a = excitation%to_orb(1)
                b = excitation%to_orb(2)
                if (excitation%perm) then
                    sign_factor = -1
                else
                    sign_factor = 1
                end if
                ! The factors of two account for the possible ways the creation an
                ! annihilation operators appearing in the expectation value can be
                ! ordered yielding a non-zero result for a given excitation ij->ab.
                if (bs%basis_fns(i)%ms == bs%basis_fns(j)%ms) then
                    ! There are two contributions for matrix elements involving like spins:
                    ! q = k_i - k_a, and q = k_i - k_b
                    q = bs%basis_fns(i)%l - bs%basis_fns(a)%l
                    call q_in_basis(ueg_basis, q, in_basis, indx)
                    if (in_basis) then
                        ! Contribution from like spins goes in the odd indices.
                        structure_factor(indx) = structure_factor(indx) + 2.0_p*sign_factor*pop
                    end if
                    q = bs%basis_fns(i)%l - bs%basis_fns(b)%l
                    call q_in_basis(ueg_basis, q, in_basis, indx)
                    if (in_basis) then
                        structure_factor(indx) = structure_factor(indx) - 2.0_p*sign_factor*pop
                    end if
                else
                    ! Only one contribution for unlike spins.
                    if (bs%basis_fns(i)%ms == bs%basis_fns(a)%ms) then
                        q = bs%basis_fns(i)%l - bs%basis_fns(a)%l
                        call q_in_basis(ueg_basis, q, in_basis, indx)
                        if (in_basis) then
                            ! Contribution from un-like spins goes in the even indices.
                            structure_factor(indx+1) = structure_factor(indx+1) + 2.0_p*sign_factor*pop
                        end if
                    else if (bs%basis_fns(i)%ms == bs%basis_fns(b)%ms) then
                        q = bs%basis_fns(i)%l - bs%basis_fns(b)%l
                        call q_in_basis(ueg_basis, q, in_basis, indx)
                        if (in_basis) then
                            structure_factor(indx+1) = structure_factor(indx+1) - 2.0_p*sign_factor*pop
                        end if
                    end if
                end if
            end select
        end associate

        contains

            subroutine q_in_basis(ueg_basis, q, in_basis, indx)

                ! Check if momentum transfer is in our basis set, and less than
                ! the size of the structure factor.

                ! In:
                !    ueg_basis: ueg_basis_t
                !    q: momentum transfer vector
                ! Out:
                !    in_basis: true if transfer is in the basis set
                !    indx: index of momentum transfer

                use ueg_types, only: ueg_basis_t
                type(ueg_basis_t), intent(in) :: ueg_basis
                integer, intent(in) :: q(:)
                logical, intent(out) :: in_basis
                integer, intent(out) :: indx

                in_basis = .false.

                if (all(abs(q) <= ueg_basis%kmax)) then
                    indx = ueg_basis%lookup(dot_product(q,ueg_basis%offset_inds) + ueg_basis%offset)
                    if (indx > 0 .and. indx <= size(structure_factor)) then
                        in_basis = .true.
                    end if
                end if

            end subroutine q_in_basis

    end subroutine update_dmqmc_structure_factor_ueg

    subroutine update_dmqmc_H0_energy(sys, cdet, excitation, pop, H0_energy, complx)

        ! Add the contribution for the current density matrix element to the thermal
        ! zeroth-order Hamiltonian (H^0) energy estimate used.

        ! Usually one writes H = T + U, where T is the kinetic energy and U is the potential
        ! energy. However we can in principle partition H in many different
        ! ways. When working in the interaction pitcture using DMQMC it is
        ! useful to split H = H^0 + V where H^0 is the zeroth order Hamiltonian
        ! and V is some perturbation. Currently two partitions are implemented
        ! so that if dmqmc_in%initial_matrix = 'free_electron' then the
        ! partitioning H^0 = T, V = U is used, and if dmqmc_in%initial_matrix =
        ! 'hartree_fock' H^0 = \sum_i |D_i> <D_i|H|D_i> <D_i| and V = H - H^0.

        ! In:
        !    sys: system being studied.
        !    cdet: det_info_t object containing bit strings of densitry matrix
        !       element under consideration.
        !    excitation: excit_t type variable which stores information on
        !        the excitation between the two bitstring ends, corresponding
        !        to the two labels for the density matrix element.
        !    pop: number of particles on the current density matrix
        !        element.
        ! In/Out:
        !    H0_energy: current thermal zeroth order Hamiltonian energy estimate.

        use determinant_data, only: det_info_t
        use system, only: sys_t
        use excitations, only: excit_t
        use proc_pointers, only: h0_ptr

        type(sys_t), intent(in) :: sys
        type(det_info_t), intent(in) :: cdet
        type(excit_t), intent(in) :: excitation
        real(p), intent(in) :: pop(:)
        logical, intent(in), optional :: complx
        real(p), intent(inout) :: H0_energy(:)

        logical :: complx_set

        complx_set = .false.
        if (present(complx)) complx_set = complx

        if (excitation%nexcit == 0) then 
            if (.not.complx_set) then
                H0_energy = H0_energy + pop(1)*h0_ptr(sys, cdet%f)
            else
                ! Although H0 should be real in all cases, in complex 
                ! simulations, we cannot guarantee the that the trace doesn't
                ! have a phase shift.
                H0_energy = H0_energy + pop(1:2)*h0_ptr(sys, cdet%f)
            end if
        end if

    end subroutine update_dmqmc_H0_energy

    subroutine update_dmqmc_ref_proj_energy(sys, excitation, cdet, H00, pop, diagonal_contribution, trace, D0j, &
                                            energy, ref_f0, complx)

        ! If the current density matrix element is on the reference row
        ! then store the componenant of this element for the projector.
        ! Otherwise we do not include an element that is not on the referencei row.
        ! This is really a wrapper to check the matrix element is on the
        ! reference row, then uses the normal energy updating procedures
        ! to update the "trace" and numerator of the energy esitmate on this row.

        ! In:
        !    sys: system being studied.
        !    excitation: excit_t type variable which stores information on
        !        the excitation between the two bitstring ends, corresponding
        !        to the two labels for the density matrix element.
        !    H00: diagonal Hamiltonian element for the reference.
        !    pop: number of particles on the current density matrix
        !        element.
        !    pop_im: imaginary population, not used here.
        !    cdet: det_info_t object containing bit strings of densitry matrix
        !        element under consideration.
        !    diagonal_contribution: <D_i|H|D_i>-<D0|H|D0>
        !    ref_f0: the bitstring of the reference determinant used to check
        !        if we are on the reference row of the density matrix.
        !    complx: if we are sampling a complex read_in system.
        ! In/Out:
        !    energy: current reference projected energy numerator.
        !    trace: total population on diagonal elements of density matrix
        !       reference row.

        use determinant_data, only: det_info_t
        use system, only: sys_t
        use excitations, only: get_excitation_level, excit_t
        use proc_pointers, only: update_dmqmc_energy_and_trace_ptr

        type(sys_t), intent(in) :: sys
        type(excit_t), intent(inout) :: excitation
        type(det_info_t), intent(in) :: cdet
        real(p), intent(in) :: H00, diagonal_contribution
        real(p), intent(in) :: pop(:)
        integer(i0), intent(in) :: ref_f0(sys%basis%tot_string_len)
        logical, intent(in), optional :: complx
        real(p), intent(inout) :: trace(:)
        real(p), intent(inout) :: D0j(:)
        real(p), intent(inout) :: energy(:)
        integer :: row_col_excit_level

        logical :: complx_set
        ! Importance sampling (in the FCIQMC-sense) isn't used in DMQMC...
        real(p) :: trial_wfn_dat(0)

        if (present(complx)) complx_set = complx
        row_col_excit_level = get_excitation_level(ref_f0,cdet%f2)

        if (row_col_excit_level == 0) then
            D0j = D0j + abs(pop(1))
            call update_dmqmc_energy_and_trace_ptr(sys, excitation, cdet, H00, pop, &
                                                   diagonal_contribution, trace, &
                                                   energy, complx)
        end if

    end subroutine update_dmqmc_ref_proj_energy

    subroutine update_dmqmc_HI_energy(sys, cdet, excitation, pop, tdiff, HI_energy, complx)

        ! Add the contribution for the current density matrix element to the thermal
        ! interaction picture Hamiltonian (H^I) energy estimate used.

        ! H^I(tau-beta) = e^{-(beta-tau)/2 H^0} H e^{(beta-tau)/2 H^0}.

        ! Usually one writes H = T + U, where T is the kinetic energy and U is the potential
        ! energy. However we can in principle partition H in many different
        ! ways. When working in the interaction picture using DMQMC it is
        ! useful to split H = H^0 + V where H^0 is the zeroth order Hamiltonian
        ! and V is some perturbation. Currently two partitions are implemented
        ! so that if dmqmc_in%initial_matrix = 'free_electron' then the
        ! partitioning H^0 = T, V = U is used, and if dmqmc_in%initial_matrix =
        ! 'hartree_fock' H^0 = \sum_i |D_i> <D_i|H|D_i> <D_i| and V = H - H^0.

        ! In:
        !    sys: system being studied.
        !    cdet: det_info_t object containing bit strings of densitry matrix
        !       element under consideration.
        !    excitation: excit_t type variable which stores information on
        !        the excitation between the two bitstring ends, corresponding
        !        to the two labels for the density matrix element.
        !    pop: number of particles on the current density matrix
        !        element.
        !    tdiff: 0.5*(beta-tau).
        ! In/Out:
        !    HI_energy: current thermal interaction picture Hamiltonian energy estimate.

        use determinant_data, only: det_info_t
        use system, only: sys_t
        use excitations, only: excit_t
        use proc_pointers, only: sc0_ptr, h0_ptr
        use hamiltonian_data

        type(sys_t), intent(in) :: sys
        type(det_info_t), intent(in) :: cdet
        type(excit_t), intent(inout) :: excitation
        real(p), intent(in) :: pop(:)
        real(p), intent(in) :: tdiff
        logical, intent(in), optional :: complx
        real(p), intent(inout) :: HI_energy(:)

        real(p) :: diff_ijab, trace(2), energy(2)
        type(hmatel_t) :: hmatel
        logical :: complx_set
        ! Importance sampling (in the FCIQMC-sense) isn't used in DMQMC...
        real(p) :: trial_wfn_dat(0)

        diff_ijab = 0.0_p
        trace = 0.0_p
        energy = 0.0_p
        complx_set = .false.
        if (present(complx)) complx_set = complx

        ! Hamiltonian matrix element.
        call update_proj_energy_dmqmc(sys, cdet%f2, trial_wfn_dat, cdet, pop, trace, energy, excitation, hmatel, complx=complx)

        if (excitation%nexcit == 0) then
            ! Diagonal Hamiltonian elements.
            hmatel%r = sc0_ptr(sys, cdet%f)
            ! No additional test for complex system required here as the value 
            ! is always real and only one component actually used below.
            hmatel%c = hmatel%r
        else
            diff_ijab = h0_ptr(sys, cdet%f) - h0_ptr(sys, cdet%f2)
        end if

        if (.not.complx_set) then
            HI_energy(1) = HI_energy(1) + exp(-tdiff*diff_ijab)*hmatel%r*pop(1)
        else
            HI_energy(1) = HI_energy(1) + real(exp(-tdiff*diff_ijab)*hmatel%c, p)*pop(1) - &
                                        &aimag(exp(-tdiff*diff_ijab)*hmatel%c)*pop(2)
            HI_energy(2) = HI_energy(2) + real(exp(-tdiff*diff_ijab)*hmatel%c, p)*pop(2) + &
                                        &aimag(exp(-tdiff*diff_ijab)*hmatel%c)*pop(1)
        end if

    end subroutine update_dmqmc_HI_energy

    subroutine update_dmqmc_potential_energy(sys, cdet, excitation, pop, potential_energy)

        ! Add the contribution for the current density matrix element to the thermal
        ! estimate for the (electronic) potential energy.

        ! In:
        !    sys: system being studied.
        !    cdet: det_info_t object containing bit strings of densitry matrix
        !       element under consideration.
        !    excitation: excit_t type variable which stores information on
        !        the excitation between the two bitstring ends, corresponding
        !        to the two labels for the density matrix element.
        !    pop: number of particles on the current density matrix
        !        element.
        ! In/Out:
        !    potential_energy: current thermal potential energy estimate.

        use determinant_data, only: det_info_t
        use system, only: sys_t
        use excitations, only: excit_t
        use proc_pointers, only: potential_energy_ptr

        type(sys_t), intent(in) :: sys
        type(det_info_t), intent(in) :: cdet
        type(excit_t), intent(in) :: excitation
        real(p), intent(in) :: pop
        real(p), intent(inout) :: potential_energy

        potential_energy = potential_energy + pop*potential_energy_ptr(sys, cdet%f, excitation)

    end subroutine update_dmqmc_potential_energy

    subroutine update_proj_energy_dmqmc(sys, f, trial_wfn_dat, cdet, pop, trace, energy, excitation, hmatel, complx)

        ! Update projected energy contribution from current density matrix element.

        ! In:
        !    sys: system being studied.
        !    f: one bit string of density matrix element under consideration.
        !    trial_wfn_dat: info on trial wavefunction in use.
        !    cdet: det_info_t object containing bit strings of densitry matrix
        !    pop: number of particles on the current density matrix
        !        element.
        !    complx(optional): true if the system is complex
        ! In/Out:
        !    excitation: excit_t type variable which stores information on
        !        the excitation between the two bitstring ends, corresponding
        !    trace: total population on diagonal elements of density matrix
        !    energy: current thermal energy estimate.
        ! Out:
        !        to the two labels for the density matrix element.
        !    hmatel: Hamiltonian matrix element between f and the reference.

        use determinant_data, only: det_info_t
        use system, only: sys_t
        use excitations, only: excit_t
        use proc_pointers, only: update_proj_energy_ptr
        use hamiltonian, only: hmatel_t
        use qmc_data, only: estimators_t
        use const, only: i0

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f(:)
        type(det_info_t), intent(in) :: cdet
        type(excit_t), intent(inout) :: excitation
        real(p), intent(in) :: pop(:)
        real(p), intent(inout) :: trace(:), energy(:)
        type(hmatel_t), intent(out) :: hmatel
        real(p), intent(in) :: trial_wfn_dat(0)
        logical, intent(in), optional :: complx

        type(estimators_t) :: dummy_estimators
        logical :: complx_set

        complx_set = .false.
        if (present(complx)) complx_set = complx

        if (complx_set) then
            dummy_estimators%proj_energy_comp = cmplx(energy(1), energy(2), p)
            dummy_estimators%D0_population_comp = cmplx(trace(1), trace(2), p)
        else
            dummy_estimators%proj_energy = energy(1)
            dummy_estimators%D0_population = trace(1)
        end if

        call update_proj_energy_ptr(sys, f, trial_wfn_dat, cdet, pop, dummy_estimators, excitation, hmatel)

        if (complx_set) then
            energy(1) = real(dummy_estimators%proj_energy_comp, p)
            energy(2) = aimag(dummy_estimators%proj_energy_comp)
            trace(1) = real(dummy_estimators%D0_population_comp, p)
            trace(2) = aimag(dummy_estimators%D0_population_comp)
        else
            energy(1) = dummy_estimators%proj_energy
            trace(1) = dummy_estimators%D0_population
        end if

    end subroutine update_proj_energy_dmqmc

end module dmqmc_estimators

module energy_evaluation

! This module contains procedure for evaluating and estimating the energy of
! a system based upon the population dynamics of an FCIQMC calculation.

use const

implicit none

! indices for corresponding data items in estimator/population/etc buffers (see
! comments below).
enum, bind(c)
    enumerator :: proj_energy_ind = 1
    enumerator :: D0_pop_ind
    enumerator :: rspawn_ind
    enumerator :: update_tau_ind
    enumerator :: bloom_tot_ind
    enumerator :: bloom_num_ind
    enumerator :: hf_signed_pop_ind
    enumerator :: hf_proj_O_ind
    enumerator :: hf_proj_H_ind
    enumerator :: hf_D0_pop_ind
    enumerator :: proj_energy_imag_ind
    enumerator :: D0_pop_imag_ind
    enumerator :: nocc_states_ind
    enumerator :: nspawned_ind
    enumerator :: comms_found_ind
    enumerator :: error_ind
    enumerator :: nparticles_start_ind ! ensure this is always the last enumerator
end enum

contains

! --- Estimator/population/etc communication ---

    ! In order to avoid doing repeated MPI calls, we place information that must be
    ! summed over processors into one buffer and sum that buffer in one call.
    ! The buffer (typically with a name begging with rep_loop) contains:

    ! buffer(1:nparticles_start_ind-1)
    !   single-valued data given by the (descriptive) name in the enumerator
    !   above.
    ! buffer(nparticles_start_ind+iproc*particle_t%nspaces:nparticles_start_ind-1+(iproc+1)*particle_t%nspaces)
    !   the total population of each particle type on processor iproc.  Prior to
    !   communication, each processor only sets its only value.

    ! To add a data item to the buffer:
    ! 1. add it (before nparticles_start_ind) to the enum above.
    ! 2. set the buffer (with appropriate index) in local_energy_estimators.
    ! 3. set the variable to the summed version after the comm call (i.e. in
    !    communicated_energy_estimators).

    ! All other elements are set to zero.

    subroutine update_energy_estimators(qmc_in, qs, nspawn_events, ntot_particles_old, load_bal_in, doing_lb, &
                                        comms_found, error, update_tau, bloom_stats, comp)

        ! Update the shift and average the shift and projected energy
        ! estimators.

        ! Should be called every report loop in an FCIQMC/iFCIQMC calculation.

        ! In:
        !    qmc_in: input options relating to QMC methods.
        !    nspawn_events: The total number of spawning events to this process.
        !    doing_lb (optional): true if performing load balancing.
        !    load_bal_in: input options for load balancing.
        ! In/Out:
        !    qs: qmc state. Estimators are updated on output.
        !    ntot_particles_old: total number (across all processors) of
        !        particles in the simulation at end of the previous report loop.
        !        Returns the current total number of particles for use in the
        !        next report loop.
        !    comms_found: On input, true if HANDE.COMM file present on this processor. On output
        !        true if present on any processor.
        !    error (optional): true if an error has occured on this processor. On output, true if an
        !        error has occured on any processor.
        !    update_tau (optional): On input, true if the processor wants to automatically rescale
        !       tau.  On output the logical or of this across all processors (i.e. true if
        !       any one processor wants to rescale tau).
        !    bloom_stats (optional): Bloom stats.  The report loop quantities are accumulated into
        !       their respective components.
        !    comp (optional): true if running calculation with real and complex walkers. This means that:
        !       - pairwise combine total populations when calculating shift variation and generate a
        !         combined shift for both.
        !       - have complex proj energy estimator to pass real and imaginary components of via mpi.

        use bloom_handler, only: bloom_stats_t
        use qmc_data, only: qmc_in_t, load_bal_in_t, qmc_state_t

        use parallel

        type(qmc_in_t), intent(in) :: qmc_in
        integer, intent(in) :: nspawn_events
        type(qmc_state_t), intent(inout) :: qs
        real(dp), intent(inout) :: ntot_particles_old(qs%psip_list%nspaces)
        logical, optional, intent(in) :: doing_lb, comp
        logical, intent(inout) :: comms_found
        logical, intent(inout), optional :: error
        type(bloom_stats_t), intent(inout), optional :: bloom_stats
        logical, intent(inout), optional :: update_tau
        type(load_bal_in_t), intent(in) :: load_bal_in

        real(dp) :: rep_loop_loc(qs%psip_list%nspaces*nprocs+nparticles_start_ind-1)
        real(dp) :: rep_loop_sum(qs%psip_list%nspaces*nprocs+nparticles_start_ind-1)
        integer :: ierr

        call local_energy_estimators(qs, rep_loop_loc, nspawn_events, comms_found, error, update_tau, &
                                    bloom_stats, comp = comp)
        ! Don't bother to optimise for running in serial.  This is a fast
        ! routine and is run only once per report loop anyway!

#ifdef PARALLEL
        call mpi_allreduce(rep_loop_loc, rep_loop_sum, size(rep_loop_loc), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
        rep_loop_sum = rep_loop_loc
        ierr = 0 ! Prevent warning about unused variable in serial so -Werror can be used.
#endif

        call communicated_energy_estimators(qmc_in, qs, rep_loop_sum, ntot_particles_old, load_bal_in, doing_lb, &
                                            qs%psip_list%nparticles_proc, comms_found, error, update_tau, &
                                            bloom_stats, comp)

    end subroutine update_energy_estimators

    subroutine update_energy_estimators_send(rep_comm)

        ! Send report loop quantities to other processors.
        ! These won't be received until the next report loop.

        ! In/Out:
        !     rep_comm: nb_rep_t object containing report loop information.

        use qmc_data, only: nb_rep_t
        use parallel

        type(nb_rep_t), intent(inout) :: rep_comm

#ifdef PARALLEL
        integer :: i, ierr
        do i = 0, nprocs-1
           call MPI_ISend(rep_comm%rep_info, size(rep_comm%rep_info), MPI_REAL8, &
                          i, 789, MPI_COMM_WORLD, rep_comm%request(i), ierr)
        end do
#endif

    end subroutine update_energy_estimators_send

    subroutine update_energy_estimators_recv(qmc_in, qs, ntypes, rep_request_s, ntot_particles_old, nparticles_proc, &
                                             load_bal_in, doing_lb, comms_found, error, update_tau, bloom_stats, comp)

        ! Receive report loop quantities from all other processors and reduce.

        ! In:
        !    qmc_in: input options relating to QMC methods.
        !    ntypes: number of particle types/spaces being sampled.
        !    load_bal_in: input options for load balancing.
        ! In/Out:
        !    qs: qmc_state_t object containing estimators that are upadted.
        !    rep_request_s: array of requests initialised during non-blocking send of
        !        information.
        !    ntot_particles_old: total number (across all processors) of
        !        particles in the simulation at end of the previous report loop.
        !        Returns the current total number of particles for use in the
        !        next report loop.
        ! Out:
        !    nparticles_proc: number of particles in each space on each processor.
        ! In (optional):
        !    doing_lb: true if performing load balancing.
        !    comp: if true running a calculation with complex walkers and so complex
        !        projected energy to pass for reporting.
        ! In/Out (optional):
        !    bloom_stats: Bloom stats.  The report loop quantities are accumulated into
        !       their respective components.
        ! Out (optional):
        !    comms_found: whether HANDE.COMM exists
        !    error: whether an error occured on any processor.
        !    update_tau: if true, tau should be automatically rescaled.

        use bloom_handler, only: bloom_stats_t
        use qmc_data, only: qmc_in_t, load_bal_in_t, qmc_state_t
        use parallel

        type(qmc_in_t), intent(in) :: qmc_in
        type(qmc_state_t), intent(inout) :: qs
        integer :: ntypes
        integer, intent(inout) :: rep_request_s(:)
        real(dp), intent(inout) :: ntot_particles_old(:)
        type(load_bal_in_t), intent(in) :: load_bal_in
        real(dp), intent(out) :: nparticles_proc(:,:)
        logical, optional, intent(in) :: doing_lb, comp
        logical, intent(out), optional :: comms_found, error
        logical, intent(out), optional :: update_tau
        type(bloom_stats_t), intent(inout), optional :: bloom_stats

        real(dp) :: rep_info_sum(nprocs*ntypes+nparticles_start_ind-1)
        real(dp) :: rep_loop_reduce(nprocs*(nprocs*ntypes+nparticles_start_ind-1))
#ifdef PARALLEL
        integer :: rep_request_r(0:nprocs-1)
        integer :: stat_ir_s(MPI_STATUS_SIZE, nprocs), stat_ir_r(MPI_STATUS_SIZE, nprocs), ierr
#endif
        integer :: i, j, data_size
        logical :: comp_loc

        comp_loc = .false.
        if (present(comp)) comp_loc = comp

        data_size = nprocs*ntypes+ nparticles_start_ind-1
#ifdef PARALLEL
        do i = 0, nprocs-1
            call MPI_IRecv(rep_loop_reduce(i*data_size+1:(i+1)*data_size), data_size, MPI_REAL8, &
                           i, 789, MPI_COMM_WORLD, rep_request_r(i), ierr)
        end do
        call MPI_Waitall(nprocs, rep_request_r, stat_ir_r, ierr)
        call MPI_Waitall(nprocs, rep_request_s, stat_ir_s, ierr)
#endif
        ! Reduce quantities across processors.
        ! Each processor sent its buffer (see comments at top of this section of
        ! procedures for format), which are now concatenated together in
        ! rep_loop_reduce.
        do i = 1, nparticles_start_ind-1
            rep_info_sum(i) = sum(rep_loop_reduce(i::data_size))
        end do
        forall (i=nparticles_start_ind:nparticles_start_ind+ntypes-1,j=0:nprocs-1) &
                rep_info_sum(i+j) = sum(rep_loop_reduce(i+j::data_size))

        call communicated_energy_estimators(qmc_in, qs, rep_info_sum, ntot_particles_old, load_bal_in, doing_lb, &
                                    nparticles_proc, comms_found, error, update_tau, bloom_stats, comp_loc)

    end subroutine update_energy_estimators_recv

    subroutine local_energy_estimators(qs, rep_loop_loc, nspawn_events, comms_found, error, update_tau, &
                                       bloom_stats, spawn_elsewhere, comp)

        ! Enter processor dependent report loop quantites into array for
        ! efficient sending to other processors.

        ! In:
        !    qs: qmc_state_t object containing estimators.
        ! In (optional):
        !    nspawn_events: The total number of spawning events to this process.
        !    comms_found: whether HANDE.COMM exists on this processor
        !    error: whether an error has occured on this processor.
        !    update_tau: if true, then the current processor wants to automatically rescale tau.
        !    bloom_stats: Bloom stats.  The report loop quantities must be set.
        !    spawn_elsewhere: number of walkers spawned from current processor
        !        not including those spawned to current processor.
        !    comp: if true running a calculation with complex walkers and so complex
        !        projected energy to pass for reporting.
        ! Out:
        !    rep_loop_loc: array containing local quantities required for energy
        !       evaluation.

        use bloom_handler, only: bloom_stats_t
        use calc, only: doing_calc, hfs_fciqmc_calc
        use parallel, only: iproc
        use qmc_data, only: qmc_state_t

        type(qmc_state_t), intent(in) :: qs
        real(dp), intent(out) :: rep_loop_loc(:)
        type(bloom_stats_t), intent(in), optional :: bloom_stats
        integer, intent(in), optional :: nspawn_events
        logical, intent(in), optional :: update_tau
        integer, intent(in) , optional :: spawn_elsewhere
        logical, intent(in), optional :: comms_found, error, comp

        integer :: offset
        logical :: comp_param

        rep_loop_loc = 0.0_dp

        comp_param = .false.
        if (present(comp)) comp_param = comp

        if (comp_param) then
            ! [review] - JSS: I think it's fine to do this branch in all cases
            ! [review] - JSS: just for convenience.  The additional work is
            ! [review] - JSS: minimal and the imaginary quantities will just be zero anyway.
            ! [reply] - CJCS: I'm unsure; if %proj_energy and %D0_population were both completely
            ! [reply] - CJCS: set to 0 in the complex case we could just do addition, but since we're
            ! [reply] - CJCS: setting the real values to the real component of the complex for
            ! [reply] - CJCS: convenience in communicated_energy_estimators we can't guarantee
            ! [reply] - CJCS: the one we don't want will be unset.
            rep_loop_loc(proj_energy_ind) = real(qs%estimators%proj_energy_comp, p)
            rep_loop_loc(D0_pop_ind) = real(qs%estimators%D0_population_comp, p)
            rep_loop_loc(proj_energy_imag_ind) = aimag(qs%estimators%proj_energy_comp)
            rep_loop_loc(D0_pop_imag_ind) = aimag(qs%estimators%D0_population_comp)
        else
            rep_loop_loc(proj_energy_ind) = qs%estimators%proj_energy
            rep_loop_loc(D0_pop_ind) = qs%estimators%D0_population
        end if
        rep_loop_loc(rspawn_ind) = qs%spawn_store%rspawn
        if (present(update_tau)) then
            if (update_tau) rep_loop_loc(update_tau_ind) = 1.0_p
        end if
        if (present(bloom_stats)) then
            rep_loop_loc(bloom_tot_ind) = bloom_stats%tot_bloom_curr
            rep_loop_loc(bloom_num_ind) = bloom_stats%nblooms_curr
        end if
        if (doing_calc(hfs_fciqmc_calc)) &
            rep_loop_loc(hf_signed_pop_ind) = calculate_hf_signed_pop(qs%psip_list)
        rep_loop_loc(hf_proj_O_ind) = qs%estimators%proj_hf_O_hpsip
        rep_loop_loc(hf_proj_H_ind) = qs%estimators%proj_hf_H_hfpsip
        rep_loop_loc(hf_D0_pop_ind) = qs%estimators%D0_hf_population
        rep_loop_loc(nocc_states_ind) = qs%psip_list%nstates

        if (present(nspawn_events)) rep_loop_loc(nspawned_ind) = nspawn_events

        offset = nparticles_start_ind-1 + iproc*qs%psip_list%nspaces
        if (present(spawn_elsewhere)) then
            rep_loop_loc(offset+1:offset+qs%psip_list%nspaces) = qs%psip_list%nparticles + spawn_elsewhere
        else
            rep_loop_loc(offset+1:offset+qs%psip_list%nspaces) = qs%psip_list%nparticles
        end if
        if (present(comms_found)) then
            if (comms_found) rep_loop_loc(comms_found_ind) = 1.0_p
        end if
        if (present(error)) then
            if (error) rep_loop_loc(error_ind) = 1.0_p
        end if

    end subroutine local_energy_estimators

    subroutine communicated_energy_estimators(qmc_in, qs, rep_loop_sum, ntot_particles_old, load_bal_in, &
                                              doing_lb, nparticles_proc, comms_found, error, update_tau, &
                                              bloom_stats, comp)

        ! Update report loop quantites with information received from other
        ! processors.

        ! In:
        !    qmc_in: input options relating to QMC methods.
        !    rep_loop_sum: array containing quantites required for energy
        !        evaluation.
        !    load_bal_in: input options for load balancing.
        ! In/Out:
        !    qs: QMC state containing estimators.
        !    ntot_particles_old: total number (across all processors) of
        !        particles in the simulation at end of the previous report loop.
        !        Returns the current total number of particles for use in the
        !        next report loop.
        ! Out:
        !    nparticles_proc: number of particles in each space on each processor.
        ! In (optional):
        !    doing_lb: true if performing load balancing.
        !    comp: true if qmc calculation with real and imaginary walkers.
        ! In/Out (optional):
        !    bloom_stats: blooming stats.
        ! Out (optional):
        !    update_tau: if true, tau should be automatically rescaled.
        !    comms_found: whether a HANDE.COMM file exists
        !    error: whether an error occured on any processor.

        use load_balancing, only: check_imbalance
        use bloom_handler, only: bloom_stats_t
        use calc, only: doing_calc, hfs_fciqmc_calc
        use parallel, only: nprocs
        use qmc_data, only: qmc_in_t, load_bal_in_t, qmc_state_t

        type(qmc_in_t), intent(in) :: qmc_in
        type(qmc_state_t), intent(inout) :: qs
        real(dp), intent(in) :: rep_loop_sum(:)
        real(dp), intent(inout) :: ntot_particles_old(:)
        real(dp), intent(out) :: nparticles_proc(:,:)
        type(load_bal_in_t), intent(in) :: load_bal_in
        logical, optional, intent(in) :: doing_lb, comp
        type(bloom_stats_t), intent(inout), optional :: bloom_stats
        logical, intent(out), optional :: comms_found, error
        logical, intent(out), optional :: update_tau

        real(dp) :: ntot_particles(size(ntot_particles_old)), new_hf_signed_pop
        real(p) :: pop_av
        real(dp) :: nparticles_wfn
        integer :: i, ntypes
        logical :: comp_param

        comp_param = .false.
        if (present(comp)) comp_param = comp

        ntypes = size(ntot_particles_old) ! Just to save passing in another parameter...

        qs%estimators%proj_energy_comp = cmplx(rep_loop_sum(proj_energy_ind), &
                                                rep_loop_sum(proj_energy_imag_ind), p)
        qs%estimators%D0_population_comp = cmplx(rep_loop_sum(D0_pop_ind), &
                                                rep_loop_sum(D0_pop_imag_ind), p)
        qs%estimators%proj_energy = rep_loop_sum(proj_energy_ind)
        qs%estimators%D0_population = rep_loop_sum(D0_pop_ind)

        qs%spawn_store%rspawn = real(rep_loop_sum(rspawn_ind), p)
        if (present(update_tau)) then
            update_tau = abs(rep_loop_sum(update_tau_ind)) > depsilon
        end if
        if (present(bloom_stats)) then
            bloom_stats%tot_bloom_curr = real(rep_loop_sum(bloom_tot_ind), p)
            bloom_stats%nblooms_curr = nint(rep_loop_sum(bloom_num_ind))
            ! Also add to running totals.
            bloom_stats%tot_bloom = bloom_stats%tot_bloom + bloom_stats%tot_bloom_curr 
            bloom_stats%nblooms = bloom_stats%nblooms + bloom_stats%nblooms_curr
        end if
        new_hf_signed_pop = rep_loop_sum(hf_signed_pop_ind)
        qs%estimators%proj_hf_O_hpsip = real(rep_loop_sum(hf_proj_O_ind), p)
        qs%estimators%proj_hf_H_hfpsip = real(rep_loop_sum(hf_proj_H_ind), p)
        qs%estimators%D0_hf_population = real(rep_loop_sum(hf_D0_pop_ind), p)
        qs%estimators%tot_nstates = nint(rep_loop_sum(nocc_states_ind))
        qs%estimators%tot_nspawn_events = nint(rep_loop_sum(nspawned_ind))
        if (present(comms_found)) then
            comms_found = abs(rep_loop_sum(comms_found_ind)) > depsilon
        end if
        if (present(error)) then
            error = abs(rep_loop_sum(error_ind)) > depsilon
        end if

        do i = 1, ntypes
            nparticles_proc(i,:nprocs) = rep_loop_sum(nparticles_start_ind-1+i::ntypes)
            ntot_particles(i) = sum(nparticles_proc(i,:nprocs))
        end do

        associate(lb=>qs%par_info%load)
            if (present(doing_lb)) then
                if (doing_lb .and. ntot_particles(1) > load_bal_in%pop .and. lb%nattempts < load_bal_in%max_attempts) then
                    pop_av = real(sum(nparticles_proc(1,:nprocs))/nprocs, p)
                    ! Check if there is at least one processor with load imbalance.
                    call check_imbalance(real(nparticles_proc,p), pop_av, load_bal_in%percent, lb%needed)
                end if
            end if
        end associate

        if (qs%vary_shift(1)) then
            if (comp_param) then
                call update_shift(qmc_in, qs, qs%shift(1), ntot_particles_old(1) + ntot_particles_old(2), &
                                    ntot_particles(1) + ntot_particles(2), qmc_in%ncycles)
                ! [review] - JSS: shift(2) is not used in the complex code.  What is it meant to represent?
                ! [reply] - CJCS: Shift in imaginary space, which is the same as in real space.
                ! [reply] - CJCS: Mainly just definined for compatibility with general multiple space routines
                ! [reply] -  eg. if we were doing replica tricks + complex could just call:
                ! do idet = 1, ndets
                !   do i = 1, nspaces
                !       call death( pop(i, idet), qs%shift(i)...)
                !   etc
                ! [reply] - CJCS: rather than having a more complicated structure with shift have size nspaces/2.
                ! [reply] - CJCS: The gain from removing the duplicate information seems negligible compared to
                ! [reply] - CJCS: the gain in compatiblity.
                qs%shift(2) = qs%shift(1)
            else
                call update_shift(qmc_in, qs, qs%shift(1), ntot_particles_old(1), ntot_particles(1), &
                                    qmc_in%ncycles)
                if (doing_calc(hfs_fciqmc_calc)) then
                    call update_hf_shift(qmc_in, qs, qs%shift(2), ntot_particles_old(1), ntot_particles(1), &
                                         qs%estimators%hf_signed_pop, new_hf_signed_pop, qmc_in%ncycles)
                end if
            end if
        end if
        ntot_particles_old = ntot_particles
        qs%estimators%hf_signed_pop = new_hf_signed_pop
        !      2. I am not sure that setting the target population to be the
        !         total of the real and complex space is ideal.
        ! [reply] - CJCS: What's the alternative? And the benefit to be gained? So long as
        ! [reply] - CJCS: we have an effective measure to base population control from and
        ! [reply] - CJCS: don't bias the wavefunction one way or the other it shouldn't have
        ! [reply] - CJCS: much impact. If we were to use only the real or imaginary population
        ! [reply] - CJCS: we could change our sampling of the space depending upon the phase
        ! [reply] - CJCS: the wavefunction settles on (which may not be consistent), making
        ! [reply] - CJCS: array sizing much more difficult.
        ! [reply] - CJCS: We could use the sum of complex magnitudes as a measure, which would
        ! [reply] - CJCS: be robust to potential phase rotations from the standpoint of the
        ! [reply] - CJCS: derivation of shift control (exponential growth of the groundstate
        ! [reply] - CJCS: wavefunction, so ignoring wavefunction phase). However, this hasn't
        ! [reply] - CJCS: shown itself to be a major problem thus far.
        !      3. Why set the shift from the real(proj/N_0) - this should be
        !         equivalent to just the ratio of the real parts (though only on
        !         average).  The key is to set the shift to something closer to
        !         the true energy rather so this shouldn't matter too much.
        ! [reply] - CJCS: Depending upon the system, we can end up with the phase of the
        ! [reply] - CJCS: reference being fairly far away from 0 (although stable in a reasonable
        ! [reply] - CJCS: calculation). In this case there's a fairly big difference between the
        ! [reply] - CJCS: two values.

        if (comp_param) then
            nparticles_wfn = ntot_particles(1) + ntot_particles(2)
        else
            nparticles_wfn = ntot_particles(1)
        end if

        if (nparticles_wfn > qs%target_particles .and. .not.qs%vary_shift(1)) then
            qs%vary_shift = .true.
            if (qmc_in%vary_shift_from_proje) then
                ! Set shift to be instantaneous projected energy.
                associate(est=>qs%estimators)
                    if (comp_param) then
                        qs%shift = real(est%proj_energy_comp/est%D0_population_comp, p)
                    else
                        qs%shift = est%proj_energy/est%D0_population
                    end if
                    if (doing_calc(hfs_fciqmc_calc)) then
                        qs%shift(2) = est%proj_hf_O_hpsip/est%D0_population + est%proj_hf_H_hfpsip/est%D0_population &
                                                             - (est%proj_energy*est%D0_hf_population)/est%D0_population**2
                    end if
                end associate
            else
                qs%shift = qmc_in%vary_shift_from
            end if
        end if

        ! average energy quantities over report loop.
        qs%estimators%proj_energy = qs%estimators%proj_energy/qmc_in%ncycles
        qs%estimators%D0_population = qs%estimators%D0_population/qmc_in%ncycles
        ! Similarly for the HFS estimator
        qs%estimators%D0_hf_population = qs%estimators%D0_hf_population/qmc_in%ncycles
        qs%estimators%proj_hf_O_hpsip = qs%estimators%proj_hf_O_hpsip/qmc_in%ncycles
        qs%estimators%proj_hf_H_hfpsip = qs%estimators%proj_hf_H_hfpsip/qmc_in%ncycles
        ! Similarly for complex quantities.
        qs%estimators%proj_energy_comp = qs%estimators%proj_energy_comp/qmc_in%ncycles
        qs%estimators%D0_population_comp = qs%estimators%D0_population_comp/qmc_in%ncycles
        ! average spawning rate over report loop and processor.
        qs%spawn_store%rspawn = qs%spawn_store%rspawn/(qmc_in%ncycles*nprocs)

    end subroutine communicated_energy_estimators

!--- Shift updates ---

    subroutine update_shift(qmc_in, qs, loc_shift, nparticles_old, nparticles, nupdate_steps)

        ! Update the shift according to:
        !  shift(beta) = shift(beta-A*tau) - xi*log(N_w(tau)/N_w(beta-A*tau))/(A*tau)
        ! where
        !  * shift(beta) is the shift at imaginary time beta;
        !  * A*tau is the amount of imaginary time between shift-updates (=# of
        !    Monte Carlo cycles between updating the shift);
        !  * xi is a damping factor (0.05-0.10 is appropriate) to prevent large fluctations;
        !  * N_w(beta) is the total number of particles at imaginary time beta.
        ! The running average of the shift is also updated.

        ! In:
        !    qmc_in: Input options relating to QMC methods.
        !    qs: qmc state.
        !    nparticles_old: N_w(beta-A*tau).
        !    nparticles: N_w(beta).
        !    nupdate_steps: number of iterations between shift updates.
        ! In/Out:
        !    loc_shift: energy shift/offset.  Set to S(beta-A*tau) on input and S(beta) on output.

        use qmc_data, only: qmc_in_t, qmc_state_t

        type(qmc_in_t), intent(in) :: qmc_in
        type(qmc_state_t), intent(in) :: qs
        real(p), intent(inout) :: loc_shift
        real(dp), intent(in) :: nparticles_old, nparticles
        integer, intent(in) :: nupdate_steps

        ! dmqmc_factor is included to account for a factor of 1/2 introduced into tau in
        ! DMQMC calculations. In all other calculation types, it is set to 1, and so can be ignored.
        loc_shift = loc_shift - real(log(nparticles/nparticles_old)*qmc_in%shift_damping/(qs%dmqmc_factor*qs%tau*nupdate_steps) ,p)

    end subroutine update_shift

    subroutine update_hf_shift(qmc_in, qs, hf_shift, nparticles_old, nparticles, nhf_particles_old, nhf_particles, nupdate_steps)

        ! Update the Hellmann-Feynman shift, \tilde{S}.

        ! In:
        !    qmc_in: Input options relating to QMC methods.
        !    qs: qmc state.
        !    nparticles_old: N_w(beta-A*tau); total Hamiltonian population at beta-Atau.
        !    nparticles: N_w(beta); total Hamiltonian population at beta.
        !    nhf_particles_old: N_w(beta-A*tau); total Hellmann-Feynman (signed) population at beta-Atau.
        !    nhf_particles: N_w(beta); total Hellmann-Feynman (signed) population at beta.
        !    nupdate_steps: number of iterations between shift updates.
        ! In/Out:
        !    hf_shift: Hellmann--Feynman shift.  Set to \tilde{S}(beta-A*tau) on input and
        !       \tilde{S}(beta) on output.

        ! WARNING:
        ! The Hellmann-Feynman signed population is not simply the sum over
        ! Hellmann-Feynman walkers but also involves the Hamiltonian walkers and
        ! *must* be calculated using calculate_hf_signed_pop.

        use qmc_data, only: qmc_in_t, qmc_state_t

        type(qmc_in_t), intent(in) :: qmc_in
        type(qmc_state_t), intent(in) :: qs
        real(dp), intent(in) :: nparticles_old, nparticles, nhf_particles_old, nhf_particles
        integer, intent(in) :: nupdate_steps
        real(p), intent(inout) :: hf_shift

        ! Given the definition of the shift, S, \tilde{S} \equiv \frac{dS}{d\alpha}|_{\alpha=0}.
        ! Hence \tilde{S}(\beta) =
        !           \tilde{S}(\beta-A\tau)
        !           - \frac{\xi}{A\tau} [ \frac{\tilde{N}_w(\beta)}{N_w(\beta)}
        !                                 - \frac{\tilde{N}_w(\beta-A\tau)}{N_w(\beta-A\tau)} ]
        ! where N_w(\beta) is the total population of (Hamiltonian) walkers at
        ! imaginary time \beta and \tilde{N}_w = \frac{dN_w}{d\alpha}|_{\alpha=0}.
        ! The latter quantity is calculated in calculate_hf_signed fpop.

        hf_shift = hf_shift - &
                 real((qmc_in%shift_damping/(qs%tau*nupdate_steps)) &
                 *(nhf_particles/nparticles - nhf_particles_old/nparticles_old),p)

    end subroutine update_hf_shift

    function calculate_hf_signed_pop(psip_list) result(hf_signed_pop)

        ! Find
        !    \sum_j sign(N_j(\beta)) \tilde{N}_j(\beta)
        ! where N_j(\beta) is the Hamiltonian population on j at imaginary time
        ! \beta and \tilde{N}_j(\beta) is the Hellmann-Feynman population on
        ! j at imaginary time \beta.

        ! In:
        !    psip_list: particle_t object containing occupied states and their populations in both the Hamiltonian and
        !       Hellmann--Feynman spaces.

        use qmc_data, only: particle_t

        type(particle_t), intent(in) :: psip_list

        real(p) :: hf_signed_pop
        real(p) :: real_population(psip_list%nspaces)

        integer :: i

        hf_signed_pop = 0.0_dp
        do i = 1, psip_list%nstates
            real_population = real(abs(psip_list%pops(:,i)),p)/psip_list%pop_real_factor
            if (psip_list%pops(1,i) == 0_int_p) then
                ! An approximation: let alpha->0_+ (using alpha->0_- appears not to make much of a difference).
                ! letting alpha->0_+
                hf_signed_pop = hf_signed_pop + real_population(2)
            else
                hf_signed_pop = hf_signed_pop + sign(1.0_p, real_population(1))*&
                                                 real_population(2)
            end if
        end do

    end function calculate_hf_signed_pop

!--- Projected estimator updates ---

    pure subroutine update_proj_energy_hub_k(sys, f0, wfn_dat, cdet, pop, D0_pop_sum, proj_energy_sum, excitation, hmatel)

        ! Add the contribution of the current determinant to the projected
        ! energy.
        ! The correlation energy given by the projected energy is:
        !   \sum_{i \neq 0} <D_i|H|D_0> N_i/N_0
        ! where N_i is the population on the i-th determinant, D_i,
        ! and 0 refers to the reference determinant.
        ! During a MC cycle we store N_0 and \sum_{i \neq 0} <D_i|H|D_0> N_i.
        ! This procedure is for the Hubbard model in momentum space only.

        ! In:
        !    sys: system being studied.
        !    f0: reference determinant.
        !    wfn_dat: trial wavefunction data (unused, included for interface compatibility).
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.  Only the bit string field needs to be set.
        !    pop: population on current determinant.
        ! In/Out:
        !    D0_pop_sum: running total of N_0, the population on the reference
        !        determinant, |D_0>.  Updated only if cdet is |D_0>.
        !    proj_energy_sum: running total of \sum_{i \neq 0} <D_i|H|D_0> N_i.
        !        Updated only if <D_i|H|D_0> is non-zero.
        !    excitation: excitation connecting the determinant to the reference determinant.
        ! Out:
        !    hmatel: <D_i|H|D_0>, the Hamiltonian matrix element between the
        !       determinant and the reference determinant.

        ! NOTE: it is the programmer's responsibility to ensure D0_pop_sum and
        ! proj_energy_sum are zero before the first call.

        use determinants, only: det_info_t
        use excitations, only: excit_t
        use hamiltonian_hub_k, only: slater_condon2_hub_k
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f0(:)
        real(p), intent(in) :: wfn_dat(:)
        type(det_info_t), intent(in) :: cdet
        real(p), intent(in) :: pop
        real(p), intent(inout) :: D0_pop_sum, proj_energy_sum
        type(excit_t), intent(inout) :: excitation
        real(p), intent(out) :: hmatel

        if (excitation%nexcit == 0) then
            ! Have reference determinant.
            D0_pop_sum = D0_pop_sum + pop
        else if (excitation%nexcit == 2) then
            ! Have a determinant connected to the reference determinant: add to
            ! projected energy.
            hmatel = slater_condon2_hub_k(sys, excitation%from_orb(1), excitation%from_orb(2), &
                                       & excitation%to_orb(1), excitation%to_orb(2),excitation%perm)
            proj_energy_sum = proj_energy_sum + hmatel*pop
        end if

    end subroutine update_proj_energy_hub_k

    pure subroutine update_proj_energy_hub_real(sys, f0, wfn_dat, cdet, pop, D0_pop_sum, proj_energy_sum, excitation, hmatel)

        ! Add the contribution of the current determinant to the projected
        ! energy.
        ! The correlation energy given by the projected energy is:
        !   \sum_{i \neq 0} <D_i|H|D_0> N_i/N_0
        ! where N_i is the population on the i-th determinant, D_i,
        ! and 0 refers to the reference determinant.
        ! During a MC cycle we store N_0 and \sum_{i \neq 0} <D_i|H|D_0> N_i.
        ! This procedure is for the Hubbard model in real space only.

        ! In:
        !    sys: system being studied.
        !    f0: reference determinant.
        !    wfn_dat: trial wavefunction data (unused, included for interface compatibility).
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.  Only the bit string field needs to be set.
        !    pop: population on current determinant.
        ! In/Out:
        !    D0_pop_sum: running total of N_0, the population on the reference
        !        determinant, |D_0>.  Updated only if cdet is |D_0>.
        !    proj_energy_sum: running total of \sum_{i \neq 0} <D_i|H|D_0> N_i.
        !        Updated only if <D_i|H|D_0> is non-zero.
        !    excitation: excitation connecting the determinant to the reference determinant.
        ! Out:
        !    hmatel: <D_i|H|D_0>, the Hamiltonian matrix element between the
        !       determinant and the reference determinant.

        ! NOTE: it is the programmer's responsibility to ensure D0_pop_sum and
        ! proj_energy_sum are zero before the first call.

        use determinants, only: det_info_t
        use excitations, only: excit_t
        use hamiltonian_hub_real, only: slater_condon1_hub_real
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f0(:)
        real(p), intent(in) :: wfn_dat(:)
        type(det_info_t), intent(in) :: cdet
        real(p), intent(in) :: pop
        real(p), intent(inout) :: D0_pop_sum, proj_energy_sum
        type(excit_t), intent(inout) :: excitation
        real(p), intent(out) :: hmatel

        if (excitation%nexcit == 0) then
            ! Have reference determinant.
            D0_pop_sum = D0_pop_sum + pop
        else if (excitation%nexcit == 1) then
            ! Have a determinant connected to the reference determinant: add to
            ! projected energy.
            hmatel = slater_condon1_hub_real(sys, excitation%from_orb(1), excitation%to_orb(1), excitation%perm)
            proj_energy_sum = proj_energy_sum + hmatel*pop
        end if

    end subroutine update_proj_energy_hub_real

    pure subroutine update_proj_energy_mol(sys, f0, wfn_dat, cdet, pop, D0_pop_sum, proj_energy_sum, excitation, hmatel)

        ! Add the contribution of the current determinant to the projected
        ! energy.
        ! The correlation energy given by the projected energy is:
        !   \sum_{i \neq 0} <D_i|H|D_0> N_i/N_0
        ! where N_i is the population on the i-th determinant, D_i,
        ! and 0 refers to the reference determinant.
        ! During a MC cycle we store N_0 and \sum_{i \neq 0} <D_i|H|D_0> N_i.
        ! This procedure is for molecular systems (i.e. those defined by an
        ! FCIDUMP file).

        ! In:
        !    sys: system being studied.
        !    f0: reference determinant.
        !    wfn_dat: trial wavefunction data (unused, included for interface compatibility).
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.  Only the bit string field needs to be set.
        !    pop: population on current determinant.
        ! In/Out:
        !    D0_pop_sum: running total of N_0, the population on the reference
        !        determinant, |D_0>.  Updated only if cdet is |D_0>.
        !    proj_energy_sum: running total of \sum_{i \neq 0} <D_i|H|D_0> N_i.
        !        Updated only if <D_i|H|D_0> is non-zero.
        !    excitation: excitation connecting the determinant to the reference determinant.
        ! Out:
        !    hmatel: <D_i|H|D_0>, the Hamiltonian matrix element between the
        !       determinant and the reference determinant.

        ! NOTE: it is the programmer's responsibility to ensure D0_pop_sum and
        ! proj_energy_sum are zero before the first call.

        use determinants, only: det_info_t
        use excitations, only: excit_t
        use hamiltonian_molecular, only: slater_condon1_mol_excit, slater_condon2_mol_excit
        use point_group_symmetry, only: cross_product_pg_basis
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f0(:)
        real(p), intent(in) :: wfn_dat(:)
        type(det_info_t), intent(in) :: cdet
        real(p), intent(in) :: pop
        real(p), intent(inout) :: D0_pop_sum, proj_energy_sum
        type(excit_t), intent(inout) :: excitation
        real(p), intent(out) :: hmatel

        integer :: ij_sym, ab_sym

        hmatel = 0.0_p

        select case(excitation%nexcit)
        case (0)
            ! Have reference determinant.
            D0_pop_sum = D0_pop_sum + pop
        case(1)
            ! Have a determinant connected to the reference determinant by
            ! a single excitation: add to projected energy.
            ! Is excitation symmetry allowed?
            if (sys%basis%basis_fns(excitation%from_orb(1))%Ms == sys%basis%basis_fns(excitation%to_orb(1))%Ms .and. &
                    sys%basis%basis_fns(excitation%from_orb(1))%sym == sys%basis%basis_fns(excitation%to_orb(1))%sym) then
                hmatel = slater_condon1_mol_excit(sys, cdet%occ_list, excitation%from_orb(1), excitation%to_orb(1), &
                                                  excitation%perm)
                proj_energy_sum = proj_energy_sum + hmatel*pop
            end if
        case(2)
            ! Have a determinant connected to the reference determinant by
            ! a double excitation: add to projected energy.
            ! Is excitation symmetry allowed?
            if (sys%basis%basis_fns(excitation%from_orb(1))%Ms+sys%basis%basis_fns(excitation%from_orb(2))%Ms == &
                    sys%basis%basis_fns(excitation%to_orb(1))%Ms+sys%basis%basis_fns(excitation%to_orb(2))%Ms) then
                ij_sym = cross_product_pg_basis(sys%read_in%pg_sym, excitation%from_orb(1), excitation%from_orb(2), &
                                                sys%basis%basis_fns)
                ab_sym = cross_product_pg_basis(sys%read_in%pg_sym, excitation%to_orb(1), excitation%to_orb(2), sys%basis%basis_fns)
                if (ij_sym == ab_sym) then
                    hmatel = slater_condon2_mol_excit(sys, excitation%from_orb(1), excitation%from_orb(2), &
                                                      excitation%to_orb(1), excitation%to_orb(2),     &
                                                      excitation%perm)
                    proj_energy_sum = proj_energy_sum + hmatel*pop
                end if
            end if
        end select

    end subroutine update_proj_energy_mol

    pure subroutine update_proj_energy_mol_complex(sys, f0, wfn_dat, cdet, pop, D0_pop_sum_comp,&
                             proj_energy_sum_comp, excitation, hmatel)

        ! Add the contribution of the current determinant to the projected
        ! energy. This function is specifically for systems with complex
        ! hamiltonian elements and coefficients.
        ! The correlation energy given by the projected energy is:
        !   \sum_{i \neq 0} <D_i|H|D_0> N_i/N_0
        ! where N_i is the population on the i-th determinant, D_i,
        ! and 0 refers to the reference determinant.
        ! During a MC cycle we store N_0 and \sum_{i \neq 0} <D_i|H|D_0> N_i.
        ! This procedure is for molecular systems (i.e. those defined by an
        ! FCIDUMP file).

        ! In:
        !    sys: system being studied.
        !    f0: reference determinant.
        !    wfn_dat: trial wavefunction data (unused, included for interface compatibility).
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.  Only the bit string field needs to be set.
        !    pop: population on current determinant.
        ! In/Out:
        !    D0_pop_sum_comp: running total of N_0, the population on the reference
        !        determinant, |D_0>.  Updated only if cdet is |D_0>.
        !    proj_energy_sum_comp: running total of \sum_{i \neq 0} <D_i|H|D_0> N_i.
        !        Updated only if <D_i|H|D_0> is non-zero.
        !    excitation: excitation connecting the determinant to the reference determinant.
        ! Out:
        !    hmatel_comp: <D_i|H|D_0>, the Hamiltonian matrix element between the
        !       determinant and the reference determinant.

        ! NOTE: it is the programmer's responsibility to ensure D0_pop_sum and
        ! proj_energy_sum are zero before the first call.

        use determinants, only: det_info_t
        use excitations, only: excit_t
        use hamiltonian_molecular_complex, only: slater_condon1_mol_excit_complex, slater_condon2_mol_excit_complex
        use point_group_symmetry, only: cross_product_pg_basis
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f0(:)
        real(p), intent(in) :: wfn_dat(:)
        type(det_info_t), intent(in) :: cdet
        complex(p), intent(in) :: pop
        complex(p), intent(inout) :: D0_pop_sum_comp, proj_energy_sum_comp
        type(excit_t), intent(inout) :: excitation
        complex(p), intent(out) :: hmatel

        integer :: ij_sym, ab_sym

        hmatel = cmplx(0.0, 0.0, p)

        select case(excitation%nexcit)
        case (0)
            ! Have reference determinant.
            D0_pop_sum_comp = D0_pop_sum_comp + pop
        case(1)
            ! Have a determinant connected to the reference determinant by
            ! a single excitation: add to projected energy.
            ! Is excitation symmetry allowed?
            if (sys%basis%basis_fns(excitation%from_orb(1))%Ms == sys%basis%basis_fns(excitation%to_orb(1))%Ms .and. &
                    sys%basis%basis_fns(excitation%from_orb(1))%sym == sys%basis%basis_fns(excitation%to_orb(1))%sym) then
                hmatel = slater_condon1_mol_excit_complex(sys, cdet%occ_list, excitation%from_orb(1), excitation%to_orb(1), &
                                                  excitation%perm)
                proj_energy_sum_comp = proj_energy_sum_comp + hmatel * pop
            end if
        case(2)
            ! Have a determinant connected to the reference determinant by
            ! a double excitation: add to projected energy.
            ! Is excitation symmetry allowed?
            if (sys%basis%basis_fns(excitation%from_orb(1))%Ms+sys%basis%basis_fns(excitation%from_orb(2))%Ms == &
                    sys%basis%basis_fns(excitation%to_orb(1))%Ms+sys%basis%basis_fns(excitation%to_orb(2))%Ms) then
                ij_sym = cross_product_pg_basis(sys%read_in%pg_sym, excitation%from_orb(1), excitation%from_orb(2), &
                                                sys%basis%basis_fns)
                ab_sym = cross_product_pg_basis(sys%read_in%pg_sym, excitation%to_orb(1), excitation%to_orb(2), sys%basis%basis_fns)
                if (ij_sym == ab_sym) then
                    hmatel = slater_condon2_mol_excit_complex(sys, excitation%from_orb(1), excitation%from_orb(2), &
                                                      excitation%to_orb(1), excitation%to_orb(2),     &
                                                      excitation%perm)
                    proj_energy_sum_comp = proj_energy_sum_comp + hmatel * pop
                end if
            end if
        end select

    end subroutine update_proj_energy_mol_complex

    pure subroutine update_proj_energy_ueg(sys, f0, wfn_dat, cdet, pop, D0_pop_sum, proj_energy_sum, excitation, hmatel)

        ! Add the contribution of the current determinant to the projected
        ! energy.
        ! The correlation energy given by the projected energy is:
        !   \sum_{i \neq 0} <D_i|H|D_0> N_i/N_0
        ! where N_i is the population on the i-th determinant, D_i,
        ! and 0 refers to the reference determinant.
        ! During a MC cycle we store N_0 and \sum_{i \neq 0} <D_i|H|D_0> N_i.
        ! This procedure is for the electron gas only.
        ! In:
        !    sys: system being studied.
        !    f0: reference determinant.
        !    wfn_dat: trial wavefunction data (unused, included for interface compatibility).
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.  Only the bit string field needs to be set.
        !    pop: population on current determinant.
        ! In/Out:
        !    D0_pop_sum: running total of N_0, the population on the reference
        !        determinant, |D_0>.  Updated only if cdet is |D_0>.
        !    proj_energy_sum: running total of \sum_{i \neq 0} <D_i|H|D_0> N_i.
        !        Updated only if <D_i|H|D_0> is non-zero.
        !    excitation: excitation connecting the determinant to the reference determinant.
        ! Out:
        !    hmatel: <D_i|H|D_0>, the Hamiltonian matrix element between the
        !       determinant and the reference determinant.

        ! NOTE: it is the programmer's responsibility to ensure D0_pop_sum and
        ! proj_energy_sum are zero before the first call.

        use determinants, only: det_info_t
        use excitations, only: excit_t
        use hamiltonian_ueg, only: slater_condon2_ueg
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f0(:)
        real(p), intent(in) :: wfn_dat(:)
        type(det_info_t), intent(in) :: cdet
        real(p), intent(in) :: pop
        real(p), intent(inout) :: D0_pop_sum, proj_energy_sum
        type(excit_t), intent(inout) :: excitation
        real(p), intent(out) :: hmatel

        hmatel = 0.0_p

        if (excitation%nexcit == 0) then
            ! Have reference determinant.
            D0_pop_sum = D0_pop_sum + pop
        else if (excitation%nexcit == 2) then
            ! Have a determinant connected to the reference determinant: add to
            ! projected energy.
            hmatel = slater_condon2_ueg(sys, excitation%from_orb(1), excitation%from_orb(2), &
                                       & excitation%to_orb(1), excitation%to_orb(2),excitation%perm)
            proj_energy_sum = proj_energy_sum + hmatel*pop
        end if

    end subroutine update_proj_energy_ueg

    pure subroutine update_proj_energy_ringium(sys, f0, wfn_dat, cdet, pop, D0_pop_sum, proj_energy_sum, excitation, hmatel)

        ! Add the contribution of the current determinant to the projected
        ! energy.
        ! The correlation energy given by the projected energy is:
        !   \sum_{i \neq 0} <D_i|H|D_0> N_i/N_0
        ! where N_i is the population on the i-th determinant, D_i,
        ! and 0 refers to the reference determinant.
        ! During a MC cycle we store N_0 and \sum_{i \neq 0} <D_i|H|D_0> N_i.
        ! This procedure is for ringium.
        ! In:
        !    sys: system being studied.
        !    f0: reference determinant.
        !    wfn_dat: trial wavefunction data (unused, included for interface compatibility).
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.  Only the bit string field needs to be set.
        !    pop: population on current determinant.
        ! In/Out:
        !    D0_pop_sum: running total of N_0, the population on the reference
        !        determinant, |D_0>.  Updated only if cdet is |D_0>.
        !    proj_energy_sum: running total of \sum_{i \neq 0} <D_i|H|D_0> N_i.
        !        Updated only if <D_i|H|D_0> is non-zero.
        ! Out:
        !    excitation: excitation connecting the determinant to the reference determinant.
        !    hmatel: <D_i|H|D_0>, the Hamiltonian matrix element between the
        !       determinant and the reference determinant.

        ! NOTE: it is the programmer's responsibility to ensure D0_pop_sum and
        ! proj_energy_sum are zero before the first call.

        use determinants, only: det_info_t
        use excitations, only: excit_t, get_excitation
        use hamiltonian_ringium, only: slater_condon2_ringium
        use system, only: sys_t
        use bit_utils, only: count_set_bits

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f0(:)
        real(p), intent(in) :: wfn_dat(:)
        type(det_info_t), intent(in) :: cdet
        real(p), intent(in) :: pop
        real(p), intent(inout) :: D0_pop_sum, proj_energy_sum
        type(excit_t), intent(inout) :: excitation
        real(p), intent(out) :: hmatel

        excitation = get_excitation(sys%nel, sys%basis, cdet%f, f0)
        hmatel = 0.0_p

        if (excitation%nexcit == 0) then
            ! Have reference determinant.
            D0_pop_sum = D0_pop_sum + pop
        else if (excitation%nexcit == 2) then
            ! Have a determinant connected to the reference determinant: add to
            ! projected energy.
            hmatel = slater_condon2_ringium(sys, excitation%from_orb(1), excitation%from_orb(2), &
                                       & excitation%to_orb(1), excitation%to_orb(2),excitation%perm)
            proj_energy_sum = proj_energy_sum + hmatel*pop
        end if

    end subroutine update_proj_energy_ringium


    subroutine update_proj_hfs_hamiltonian(sys, f, fpop, f_hfpop, fdata, excitation, hmatel, &
                                           D0_hf_pop,proj_hf_O_hpsip, proj_hf_H_hfpsip)

        ! Add the contribution of the current determinant to the projected
        ! energy in an identical way to update_proj_energy_hub_k.

        ! Also add the contribution of the current determinant to the running
        ! total of the projected Hellmann--Feynman estimators.

        ! For debugging purposes, this procedure is for when we are sampling
        ! O=H.

        ! This procedure is for the Hubbard model in momentum space only.

        ! In:
        !    sys: system being studied.  Unused.
        !    f(string_len): bit string representation of the Slater determinant, D_i.
        !    fpop: Hamiltonian population on the determinant.
        !    f_hfpop: Hellmann-Feynman population on the determinant.
        !    fdata(:): additional information about the determinant (unused, for
        !       interface compatibility only).
        !    excitation: excitation connecting the determinant to the reference determinant.
        !    hmatel: <D_i|H|D_0>, the Hamiltonian matrix element between the
        !       determinant and the reference determinant.
        ! In/Out:
        !    D0_hf_population: running total of the Hellmann-Feynman population
        !       on the reference.  Only updated if D_i *is* the reference determinant.
        !    proj_hf_O_hpsip: running total of the numerator, \sum_{i \neq 0}
        !       <D_i|O|D_0> N_i, where N_i is the Hamiltonian population on D_i.
        !    proj_hf_H_fhpsip: running total of the numerator, \sum_{i \neq 0}
        !       <D_i|H|D_0> \tilde{N}_i, where \tilde{N}_i is the
        !       Hellmann-Feynman population on D_i.

        use excitations, only: excit_t
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f(:)
        integer, intent(in) :: fpop, f_hfpop
        real(p), intent(in) :: fdata(:), hmatel
        type(excit_t), intent(in) :: excitation
        real(p), intent(inout) :: D0_hf_pop, proj_hf_O_hpsip, proj_hf_H_hfpsip

        if (excitation%nexcit == 0) then
            ! Have reference determinant.
            D0_hf_pop = D0_hf_pop + f_hfpop
        else if (excitation%nexcit <= 2) then
            ! Have a determinant connected to the reference determinant: add to
            ! projected estimators.
            ! DEBUG/TESTING: For now, just using O=H
            ! In this case, \sum_j O_0j c_j = proj_energy
            proj_hf_O_hpsip = proj_hf_O_hpsip + hmatel*fpop
            ! \sum_j H_0j \tilde{c}_j is similarly easy to evaluate
            proj_hf_H_hfpsip = proj_hf_H_hfpsip + hmatel*f_hfpop
        end if

    end subroutine update_proj_hfs_hamiltonian

    subroutine update_proj_hfs_diagonal(sys, f, fpop, f_hfpop, fdata, excitation, hmatel, &
                                              D0_hf_pop,proj_hf_O_hpsip, proj_hf_H_hfpsip)

        ! Add the contribution of the current determinant to the running
        ! total of the projected Hellmann--Feynman estimator.

        ! This procedure is for when we are sampling an operator, O, which is
        ! diagonal in the Slater determinant space.

        ! In:
        !    sys: system being studied.  Unused.
        !    f(string_len): bit string representation of the Slater determinant, D_i
        !       (unused, for interface compatibility only).
        !    fpop: Hamiltonian population on the determinant (unused, for interface
        !       compatibility only).
        !    f_hfpop: Hellmann-Feynman population on the determinant.
        !    fdata(:): additional information about the determinant (unused, for
        !       interface compatibility only).
        !    excitation: excitation connecting the determinant to the reference determinant.
        !    hmatel: <D_i|H|D_0>, the Hamiltonian matrix element between the
        !       determinant and the reference determinant.
        ! In/Out:
        !    D0_hf_population: running total of the Hellmann-Feynman population
        !       on the reference.  Only updated if D_i *is* the reference determinant.
        !    proj_hf_O_hpsip: running total of the numerator, \sum_{i \neq 0}
        !       <D_i|O|D_0> N_i, where N_i is the Hamiltonian population on D_i.
        !    proj_hf_H_fhpsip: running total of the numerator, \sum_{i \neq 0}
        !       <D_i|H|D_0> \tilde{N}_i, where \tilde{N}_i is the
        !       Hellmann-Feynman population on D_i.

        use excitations, only: excit_t
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f(:)
        integer, intent(in) :: fpop, f_hfpop
        real(p), intent(in) :: fdata(:), hmatel
        type(excit_t), intent(in) :: excitation
        real(p), intent(inout) :: D0_hf_pop, proj_hf_O_hpsip, proj_hf_H_hfpsip

        if (excitation%nexcit == 0) then
            ! Have reference determinant.
            D0_hf_pop = D0_hf_pop + f_hfpop
        else if (excitation%nexcit <= 2) then
            ! Have a determinant connected to the reference determinant: add to
            ! projected energy.

            ! O is diagonal in the determinant basis.  As we are actually
            ! sampling O - <D0|O|D0>, this means that \sum_j O_j0 c_j = 0.

            ! \sum_j H_0j \tilde{c}_j is similarly easy to evaluate
            proj_hf_H_hfpsip = proj_hf_H_hfpsip + hmatel*f_hfpop
        end if

    end subroutine update_proj_hfs_diagonal

    subroutine update_proj_hfs_double_occ_hub_k(sys, f, fpop, f_hfpop, fdata, excitation, hmatel, &
                                              D0_hf_pop,proj_hf_O_hpsip, proj_hf_H_hfpsip)

        ! Add the contribution of the current determinant to the running
        ! total of the projected Hellmann--Feynman estimator.

        ! This procedure is for when we are sampling D, the double occupancy
        ! operator, in the Bloch (momentum) basis set for the Hubbard model.

        ! In:
        !    sys: system being studied.  Requires hubbard%u and lattice%nsites.
        !    f(string_len): bit string representation of the Slater determinant, D_i
        !       (unused, for interface compatibility only).
        !    fpop: Hamiltonian population on the determinant.
        !    f_hfpop: Hellmann-Feynman population on the determinant.
        !    fdata(:): additional information about the determinant (unused, for
        !       interface compatibility only).
        !    excitation: excitation connecting the determinant to the reference determinant.
        !    hmatel: <D_i|H|D_0>, the Hamiltonian matrix element between the
        !       determinant and the reference determinant.
        ! In/Out:
        !    D0_hf_population: running total of the Hellmann-Feynman population
        !       on the reference.  Only updated if D_i *is* the reference determinant.
        !    proj_hf_O_hpsip: running total of the numerator, \sum_{i \neq 0}
        !       <D_i|O|D_0> N_i, where N_i is the Hamiltonian population on D_i.
        !    proj_hf_H_fhpsip: running total of the numerator, \sum_{i \neq 0}
        !       <D_i|H|D_0> \tilde{N}_i, where \tilde{N}_i is the
        !       Hellmann-Feynman population on D_i.

        use excitations, only: excit_t
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f(:)
        integer, intent(in) :: fpop, f_hfpop
        real(p), intent(in) :: fdata(:), hmatel
        type(excit_t), intent(in) :: excitation
        real(p), intent(inout) :: D0_hf_pop, proj_hf_O_hpsip, proj_hf_H_hfpsip

        ! Note: two-electron operator.

        select case(excitation%nexcit)
        case(0)
            ! Have reference determinant.
            D0_hf_pop = D0_hf_pop + f_hfpop
        case(2)
            ! Have a determinant connected to the reference determinant: add to
            ! projected energy.

            !\hat{O}_0j = H_0j / (U L), where L is the number of sites.
            ! sampling \hat{O} - <D0|O|D0>, this means that \sum_j O_j0 c_j = 0.
            proj_hf_O_hpsip = proj_hf_O_hpsip + (hmatel/(sys%hubbard%u*sys%lattice%nsites))*fpop

            ! \sum_j H_0j \tilde{c}_j is similarly easy to evaluate
            proj_hf_H_hfpsip = proj_hf_H_hfpsip + hmatel*f_hfpop
        end select

    end subroutine update_proj_hfs_double_occ_hub_k

    subroutine update_proj_hfs_one_body_mol(sys, f, fpop, f_hfpop, fdata, excitation, hmatel, &
                                              D0_hf_pop,proj_hf_O_hpsip, proj_hf_H_hfpsip)

        ! Add the contribution of the current determinant to the running
        ! total of the projected Hellmann--Feynman estimator.

        ! This procedure is for when we are sampling O_1, a one-body operator in
        ! a molecular system (i.e. where the integrals have been read in).

        ! In:
        !    sys: system being studied.  Unused.
        !    f(string_len): bit string representation of the Slater determinant, D_i
        !       (unused, for interface compatibility only).
        !    fpop: Hamiltonian population on the determinant.
        !    f_hfpop: Hellmann-Feynman population on the determinant.
        !    fdata(:): additional information about the determinant (unused, for
        !       interface compatibility only).
        !    excitation: excitation connecting the determinant to the reference determinant.
        !    hmatel: <D_i|H|D_0>, the Hamiltonian matrix element between the
        !       determinant and the reference determinant.
        ! In/Out:
        !    D0_hf_population: running total of the Hellmann-Feynman population
        !       on the reference.  Only updated if D_i *is* the reference determinant.
        !    proj_hf_O_hpsip: running total of the numerator, \sum_{i \neq 0}
        !       <D_i|O_1|D_0> N_i, where N_i is the Hamiltonian population on D_i.
        !    proj_hf_H_fhpsip: running total of the numerator, \sum_{i \neq 0}
        !       <D_i|H|D_0> \tilde{N}_i, where \tilde{N}_i is the
        !       Hellmann-Feynman population on D_i.

        use excitations, only: excit_t
        use operators, only: one_body1_mol
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f(:)
        integer, intent(in) :: fpop, f_hfpop
        real(p), intent(in) :: fdata(:), hmatel
        type(excit_t), intent(in) :: excitation
        real(p), intent(inout) :: D0_hf_pop, proj_hf_O_hpsip, proj_hf_H_hfpsip

        real(p) :: matel

        ! Note: one-electron operator.

        select case(excitation%nexcit)
        case(0)
            ! Have reference determinant.
            D0_hf_pop = D0_hf_pop + f_hfpop
        case(1)
            ! Have a determinant connected to the reference determinant: add to
            ! projected energy.

            ! \sum_j O_0j c_j
            matel = one_body1_mol(sys, excitation%from_orb(1), excitation%to_orb(1), excitation%perm)
            proj_hf_O_hpsip = proj_hf_O_hpsip + matel*fpop

            ! \sum_j H_0j \tilde{c}_j
            proj_hf_H_hfpsip = proj_hf_H_hfpsip + hmatel*f_hfpop
        case(2)
            ! O is a one-body operator => no contributions from double
            ! excitations to \sum_j O_0j c_j.

            ! \sum_j H_0j \tilde{c}_j
            proj_hf_H_hfpsip = proj_hf_H_hfpsip + hmatel*f_hfpop
        end select

    end subroutine update_proj_hfs_one_body_mol

end module energy_evaluation

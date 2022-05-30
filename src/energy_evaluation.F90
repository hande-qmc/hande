module energy_evaluation

! This module contains procedure for evaluating and estimating the energy of
! a system based upon the population dynamics of an FCIQMC calculation.

use const
use hamiltonian_data
use qmc_data, only: estimators_t

implicit none

! indices for corresponding data items in estimator/population/etc buffers (see
! comments below).
enum, bind(c)
    enumerator :: proj_energy_ind = 1
    enumerator :: D0_pop_ind
    enumerator :: D0_noncomp_pop_ind
    ! [todo] - having a separate index for each space is not very general.
    enumerator :: proj_energy_replica_ind
    enumerator :: D0_pop_replica_ind
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
    enumerator :: proj_energy_imag_replica_ind
    enumerator :: D0_pop_imag_replica_ind
    enumerator :: nocc_states_ind
    enumerator :: nspawned_ind
    enumerator :: comms_found_ind
    enumerator :: rdm_energy_ind
    enumerator :: rdm_trace_ind
    enumerator :: nattempts_ind
    enumerator :: reblock_done_ind
    enumerator :: error_ind
end enum
! Index of last single-data entry in estimator buffer.
! WARNING: ensure this is the last entry in the enum above.
integer, parameter :: est_buf_data_size = error_ind

! Per-processor data - leave order unchanged and keep at end of enum.
! This data is sent all-to-all via a MPI_SUM and setting slots correspond to other processors to 0 initially.
! Note this means that (bar the next entry) the real start of subsequent entries are shifted.
enum, bind(c)
    enumerator :: bloom_max_indx = 1 ! Keep this first.
    enumerator :: nparticles_indx
end enum

! Index of first data
! WARNING: ensure this is the first entry in the enum above.
integer, parameter :: est_buf_start_per_proc = bloom_max_indx 

! Index of last data
! WARNING: ensure this is the last entry in the enum above.
integer, parameter :: est_buf_n_per_proc = nparticles_indx

contains

! --- Estimator/population/etc communication ---

    ! In order to avoid doing repeated MPI calls, we place information that must be summed over processors into one
    ! buffer and sum that buffer in one call.  The buffer (typically with a name beginning with rep_loop) contains:

    ! buffer(1:est_buf_data_size)
    !   single-valued data given by the (descriptive) name in the enumerator
    !   above.
    ! buffer(est_buf_data_size+1:)
    !   per-processor quantities which are broadcast to all other processors.  This is done by setting values for all
    !   bar the current processor to 0 and then performing MPI_SUM.  The (small) additional extra amount of data is more
    !   than compensated for by avoiding the latency of multiple MPI calls.  The start and length of each entry can be found
    !   from proc_data_info(:,indx), where indx is an entry in the per-processor enum above and proc_data_info is produced by
    !   get_comm_processor_indx.

    ! To add a data item to the buffer:
    ! 1. add it (preferably before error_ind and, if not, update est_buf_data_size) to the enum above.
    ! 2. set the buffer (with appropriate index) in local_energy_estimators.
    ! 3. set the variable to the summed version after the comm call (i.e. in
    !    communicated_energy_estimators).

    ! To add a per-processor quantity:
    ! 1. add it to per-processor enum after bloom_max_indx and preferably before nparticles_indx (if not, update est_buf_n_per_proc).
    ! 2. set the entry corresponding to the new enum value in get_comm_processor_indx for length.
    ! 3. follow steps 2 and 3 as above, possibly with custom handling to further sum over all processors.

    ! All other elements are set to zero.

    pure subroutine get_comm_processor_indx(nspaces, proc_data_info, ntot_proc_data)

        ! Work out the start position and length for each per-processor quantity in the comms buffer.
        ! This is sufficiently cheap that it can just be called whenever the information is required.

        ! In:
        !    nspaces: number of spaces being sampled.
        ! Out:
        !    proc_data_info(2,est_buf_n_per_proc): one entry for each per-processor quantity containing the
        !       start position of that entry in the comms buffer (proc_data_info(1,indx)) and length of
        !       that entry (proc_data_info(2,indx).
        !    ntot_proc_data: total amount of per-processor information in the comms buffer.

        use parallel, only: nprocs

        integer, intent(in) :: nspaces
        integer, intent(out) :: proc_data_info(:,:), ntot_proc_data
        integer :: i

        ! bloom_max requires 1 entry per processor.
        proc_data_info(2,bloom_max_indx) = nprocs
        ! nparticles requires nspaces entries per processor.
        proc_data_info(2,nparticles_indx) = nspaces*nprocs

        ! Start positions.  est_buf_start_per_proc is the first per-processor item.
        proc_data_info(1,est_buf_start_per_proc) = est_buf_data_size+1
        do i = est_buf_start_per_proc+1, est_buf_n_per_proc
            proc_data_info(1,i) = sum(proc_data_info(:,i-1))
        end do

        ntot_proc_data = sum(proc_data_info(2,:est_buf_n_per_proc))

    end subroutine get_comm_processor_indx

    subroutine update_energy_estimators(qmc_in, qs, nspawn_events, ntot_particles_old, load_bal_in, doing_lb, &
                                        comms_found, error, update_tau, bloom_stats, vary_shift_reference, comp)

        ! Update the shift and average the shift and projected energy
        ! estimators.

        ! Should be called every report loop in an FCIQMC/iFCIQMC calculation.

        ! In:
        !    qmc_in: input options relating to QMC methods.
        !    nspawn_events: The total number of spawning events to this process.
        !    doing_lb (optional): true if performing load balancing.
        !    load_bal_in: input options for load balancing.
        !    vary_shift_reference (optional): if true, vary shift to control reference, not total, population
        !    comp (optional): true if running calculation with real and complex walkers. This means that:
        !       - pairwise combine total populations when calculating shift variation and generate a
        !         combined shift for both.
        !       - have complex proj energy estimator to pass real and imaginary components of via mpi.
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
        logical, intent(in), optional :: vary_shift_reference

        real(dp), allocatable :: rep_loop_loc(:)
        real(dp), allocatable :: rep_loop_sum(:)
        integer :: proc_data_info(2,est_buf_n_per_proc), ntot_proc_data, ierr

        call get_comm_processor_indx(qs%psip_list%nspaces, proc_data_info, ntot_proc_data)
        allocate(rep_loop_loc(ntot_proc_data+est_buf_data_size))
        allocate(rep_loop_sum(ntot_proc_data+est_buf_data_size))

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
                                            qs%psip_list%nparticles_proc, comms_found, error, update_tau, bloom_stats, &
                                            vary_shift_reference, comp)

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

    subroutine update_energy_estimators_recv(qmc_in, qs, rep_request_s, ntot_particles_old, nparticles_proc, &
                                             load_bal_in, doing_lb, comms_found, error, update_tau, bloom_stats, comp)

        ! Receive report loop quantities from all other processors and reduce.

        ! In:
        !    qmc_in: input options relating to QMC methods.
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
        integer, intent(inout) :: rep_request_s(:)
        real(dp), intent(inout) :: ntot_particles_old(:)
        type(load_bal_in_t), intent(in) :: load_bal_in
        real(dp), intent(out) :: nparticles_proc(:,:)
        logical, optional, intent(in) :: doing_lb, comp
        logical, intent(out), optional :: comms_found, error
        logical, intent(out), optional :: update_tau
        type(bloom_stats_t), intent(inout), optional :: bloom_stats

        real(dp), allocatable :: rep_info_sum(:)
        real(dp), allocatable :: rep_loop_reduce(:)
#ifdef PARALLEL
        integer :: rep_request_r(0:nprocs-1)
        integer :: stat_ir_s(MPI_STATUS_SIZE, nprocs), stat_ir_r(MPI_STATUS_SIZE, nprocs), ierr
#endif
        integer :: i, j, data_size, proc_data_info(2,est_buf_n_per_proc), ntot_proc_data
        logical :: comp_loc

        comp_loc = .false.
        if (present(comp)) comp_loc = comp

        call get_comm_processor_indx(qs%psip_list%nspaces, proc_data_info, ntot_proc_data)
        allocate(rep_info_sum(ntot_proc_data+est_buf_data_size))
        data_size = size(rep_info_sum)

        allocate(rep_loop_reduce(nprocs*data_size))

#ifdef PARALLEL
        do i = 0, nprocs-1
            ! 789 is just a tag for this call to label it.
            call MPI_IRecv(rep_loop_reduce(i*data_size+1:(i+1)*data_size), data_size, MPI_REAL8, &
                           i, 789, MPI_COMM_WORLD, rep_request_r(i), ierr)
        end do
        call MPI_Waitall(nprocs, rep_request_r, stat_ir_r, ierr)
        call MPI_Waitall(nprocs, rep_request_s, stat_ir_s, ierr)
#endif
        ! Reduce quantities across processors.
        ! Each processor sent its buffer (length data_size) (see comments at top of this section of
        ! procedures for format), which are now concatenated together in
        ! rep_loop_reduce.
        !   NB proc_data_info(1,indx)) is start position of that entry in the comms buffer 

        ! First reduce all the non-particle-number data
        do i = 1, proc_data_info(1, nparticles_indx) - 1
            rep_info_sum(i) = sum(rep_loop_reduce(i::data_size))
        end do
        ! Each particle type's population (there are nspaces of these) is listed in each proc's buffer,
        !  starting at proc_data_info(1, nparticles_indx)
        associate(real_start_ind => proc_data_info(1, nparticles_indx))
            forall (i=real_start_ind:real_start_ind+qs%psip_list%nspaces-1,j=0:nprocs-1) &
                    rep_info_sum(i+j) = sum(rep_loop_reduce(i+j::data_size))
        end associate

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
        integer :: ntot_proc_data, proc_data_info(2,est_buf_n_per_proc)
        logical :: comp_param

        rep_loop_loc = 0.0_dp

        comp_param = .false.
        if (present(comp)) comp_param = comp

        if (comp_param) then
            rep_loop_loc(proj_energy_ind) = real(qs%estimators(1)%proj_energy_comp, p)
            rep_loop_loc(D0_pop_ind) = real(qs%estimators(1)%D0_population_comp, p)
            rep_loop_loc(proj_energy_imag_ind) = aimag(qs%estimators(1)%proj_energy_comp)
            rep_loop_loc(D0_pop_imag_ind) = aimag(qs%estimators(1)%D0_population_comp)
            if (qs%psip_list%nspaces > 2) then
                rep_loop_loc(proj_energy_replica_ind) = real(qs%estimators(3)%proj_energy_comp, p)
                rep_loop_loc(D0_pop_replica_ind) = real(qs%estimators(3)%D0_population_comp, p)
                rep_loop_loc(proj_energy_imag_replica_ind) = aimag(qs%estimators(3)%proj_energy_comp)
                rep_loop_loc(D0_pop_imag_replica_ind) = aimag(qs%estimators(3)%D0_population_comp)
            end if
        else
            rep_loop_loc(proj_energy_ind) = qs%estimators(1)%proj_energy
            rep_loop_loc(D0_pop_ind) = qs%estimators(1)%D0_population
            rep_loop_loc(D0_noncomp_pop_ind) = qs%estimators(1)%D0_noncomposite_population
            if (qs%psip_list%nspaces > 1) then
                rep_loop_loc(proj_energy_replica_ind) = qs%estimators(2)%proj_energy
                rep_loop_loc(D0_pop_replica_ind) = qs%estimators(2)%D0_population
            end if
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
        rep_loop_loc(hf_proj_O_ind) = qs%estimators(1)%proj_hf_O_hpsip
        rep_loop_loc(hf_proj_H_ind) = qs%estimators(1)%proj_hf_H_hfpsip
        rep_loop_loc(hf_D0_pop_ind) = qs%estimators(1)%D0_hf_population
        rep_loop_loc(nocc_states_ind) = qs%psip_list%nstates

        if (present(nspawn_events)) rep_loop_loc(nspawned_ind) = nspawn_events
        rep_loop_loc(rdm_energy_ind) = qs%estimators(1)%rdm_energy
        rep_loop_loc(rdm_trace_ind) = qs%estimators(1)%rdm_trace

        rep_loop_loc(nattempts_ind) = real(qs%estimators(1)%nattempts, dp)

        if (present(comms_found)) then
            if (comms_found) rep_loop_loc(comms_found_ind) = 1.0_p
        end if
        if (present(error)) then
            if (error) rep_loop_loc(error_ind) = 1.0_p
        end if
        if (qs%reblock_done) rep_loop_loc(reblock_done_ind) = 1.0_p

        ! [todo]  Include bloom_max_indx data

        ! Multi-dimensional and per-processor data.
        call get_comm_processor_indx(qs%psip_list%nspaces, proc_data_info, ntot_proc_data)
        associate(offset => proc_data_info(1,nparticles_indx) + iproc*qs%psip_list%nspaces)
            if (present(spawn_elsewhere)) then
                rep_loop_loc(offset:offset+qs%psip_list%nspaces-1) = qs%psip_list%nparticles + spawn_elsewhere
            else
                rep_loop_loc(offset:offset+qs%psip_list%nspaces-1) = qs%psip_list%nparticles
            end if
        end associate

    end subroutine local_energy_estimators

    subroutine communicated_energy_estimators(qmc_in, qs, rep_loop_sum, ntot_particles_old, load_bal_in, doing_lb, &
                                              nparticles_proc, comms_found, error, update_tau, bloom_stats, &
                                              vary_shift_reference, comp)

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
        !    vary_shift_reference: if true, vary the shift to control the reference, rather than the
        !    total, population.
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
        logical, optional, intent(in) :: doing_lb, vary_shift_reference, comp
        type(bloom_stats_t), intent(inout), optional :: bloom_stats
        logical, intent(out), optional :: comms_found, error
        logical, intent(out), optional :: update_tau

        real(dp) :: ntot_particles(size(ntot_particles_old)), new_hf_signed_pop
        real(p) :: pop_av
        integer :: i, proc_data_info(2,est_buf_n_per_proc), ntot_proc_data
        logical :: comp_param, vary_shift_reference_loc

        comp_param = .false.
        if (present(comp)) comp_param = comp

        vary_shift_reference_loc = .false.
        if (present(vary_shift_reference)) vary_shift_reference_loc = vary_shift_reference

        call get_comm_processor_indx(qs%psip_list%nspaces, proc_data_info, ntot_proc_data)

        qs%estimators(1)%proj_energy_comp = cmplx(rep_loop_sum(proj_energy_ind), &
                                                rep_loop_sum(proj_energy_imag_ind), p)
        qs%estimators(1)%D0_population_comp = cmplx(rep_loop_sum(D0_pop_ind), &
                                                rep_loop_sum(D0_pop_imag_ind), p)
        if (size(qs%estimators) > 1) then
            qs%estimators(2)%proj_energy_comp = cmplx(rep_loop_sum(proj_energy_ind), &
                                                    rep_loop_sum(proj_energy_imag_ind), p)
            qs%estimators(2)%D0_population_comp = cmplx(rep_loop_sum(D0_pop_ind), &
                                                    rep_loop_sum(D0_pop_imag_ind), p)
        end if
        if (size(qs%estimators) > 2) then
            qs%estimators(3)%proj_energy_comp = cmplx(rep_loop_sum(proj_energy_replica_ind), &
                                                    rep_loop_sum(proj_energy_imag_replica_ind), p)
            qs%estimators(4)%proj_energy_comp = cmplx(rep_loop_sum(proj_energy_replica_ind), &
                                                    rep_loop_sum(proj_energy_imag_replica_ind), p)
            qs%estimators(3)%D0_population_comp = cmplx(rep_loop_sum(D0_pop_replica_ind), &
                                                    rep_loop_sum(D0_pop_imag_replica_ind), p)
            qs%estimators(4)%D0_population_comp = cmplx(rep_loop_sum(D0_pop_replica_ind), &
                                                    rep_loop_sum(D0_pop_imag_replica_ind), p)
        end if

        qs%estimators(1)%proj_energy = real(rep_loop_sum(proj_energy_ind), p)
        qs%estimators(1)%D0_population = real(rep_loop_sum(D0_pop_ind), p)
        qs%estimators(1)%D0_noncomposite_population = real(rep_loop_sum(D0_noncomp_pop_ind), p)
        if (size(qs%estimators) > 1) then
            qs%estimators(2)%proj_energy = real(rep_loop_sum(proj_energy_replica_ind), p)
            qs%estimators(2)%D0_population = real(rep_loop_sum(D0_pop_replica_ind), p)
        end if

        qs%spawn_store%rspawn = real(rep_loop_sum(rspawn_ind), p)
        if (present(update_tau)) then
            update_tau = abs(rep_loop_sum(update_tau_ind)) > depsilon
        end if
        if (present(bloom_stats)) then
            bloom_stats%tot_bloom_curr = real(rep_loop_sum(bloom_tot_ind), p)
            bloom_stats%nblooms_curr = nint(rep_loop_sum(bloom_num_ind), int_64)
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
        qs%estimators%rdm_energy = real(rep_loop_sum(rdm_energy_ind), p)
        qs%estimators%rdm_trace = real(rep_loop_sum(rdm_trace_ind), p)
        qs%estimators%nattempts = int(rep_loop_sum(nattempts_ind), int_64)

        if (present(comms_found)) then
            comms_found = abs(rep_loop_sum(comms_found_ind)) > depsilon
        end if
        if (present(error)) then
            error = abs(rep_loop_sum(error_ind)) > depsilon
        end if

        qs%reblock_done = abs(rep_loop_sum(reblock_done_ind)) > depsilon


        ! [todo]  Extra bloom_max_indx data

        associate(nparticle_start=>proc_data_info(1,nparticles_indx))
            do i = 1, qs%psip_list%nspaces
                nparticles_proc(i,:nprocs) = rep_loop_sum(nparticle_start+i-1::qs%psip_list%nspaces)
                ntot_particles(i) = sum(nparticles_proc(i,:nprocs))
            end do
        end associate

        associate(lb=>qs%par_info%load)
            if (present(doing_lb)) then
                if (doing_lb .and. ntot_particles(1) > load_bal_in%pop .and. lb%nattempts < load_bal_in%max_attempts) then
                    pop_av = real(sum(nparticles_proc(1,:nprocs))/nprocs, p)
                    ! Check if there is at least one processor with load imbalance.
                    call check_imbalance(real(nparticles_proc,p), pop_av, load_bal_in%percent, lb%needed)
                end if
            end if
        end associate

        ! average energy quantities over report loop.
        qs%estimators%proj_energy = qs%estimators%proj_energy/qmc_in%ncycles
        qs%estimators%D0_population = qs%estimators%D0_population/qmc_in%ncycles
        qs%estimators%D0_noncomposite_population = qs%estimators%D0_noncomposite_population/qmc_in%ncycles
        ! Similarly for the HFS estimator
        qs%estimators%D0_hf_population = qs%estimators%D0_hf_population/qmc_in%ncycles
        qs%estimators%proj_hf_O_hpsip = qs%estimators%proj_hf_O_hpsip/qmc_in%ncycles
        qs%estimators%proj_hf_H_hfpsip = qs%estimators%proj_hf_H_hfpsip/qmc_in%ncycles
        ! Similarly for complex quantities.
        qs%estimators%proj_energy_comp = qs%estimators%proj_energy_comp/qmc_in%ncycles
        qs%estimators%D0_population_comp = qs%estimators%D0_population_comp/qmc_in%ncycles
        ! average spawning rate over report loop and processor.
        qs%spawn_store%rspawn = qs%spawn_store%rspawn/(qmc_in%ncycles*nprocs)

        if (doing_calc(hfs_fciqmc_calc)) then
            if (qs%vary_shift(1)) then
                call update_shift(qs, qs%shift(1), ntot_particles_old(1), ntot_particles(1), qmc_in%ncycles)
                call update_hf_shift(qmc_in, qs, qs%shift(2), ntot_particles_old(1), ntot_particles(1), &
                                     qs%estimators(1)%hf_signed_pop, new_hf_signed_pop, qmc_in%ncycles)
            end if
        else if (comp_param) then
            do i = 1, qs%psip_list%nspaces, 2
                if (qs%vary_shift(i)) then
                    call update_shift(qs, qs%shift(i), ntot_particles_old(i) + ntot_particles_old(i+1), &
                                        ntot_particles(i) + ntot_particles(i+1), qmc_in%ncycles)
                    qs%shift(i+1) = qs%shift(i)
                end if
            end do
        else
            do i = 1, qs%psip_list%nspaces
                if (qs%vary_shift(i)) then
                    if (vary_shift_reference_loc) then
                        call update_shift(qs, qs%shift(i), real(qs%estimators(i)%D0_population_old, dp), &
                                          real(qs%estimators(i)%D0_population, dp), qmc_in%ncycles)
                    else
                        call update_shift(qs, qs%shift(i), ntot_particles_old(i), ntot_particles(i), qmc_in%ncycles)
                    end if
                end if
            end do
        end if

        qs%estimators%D0_population_old = qs%estimators%D0_population
        ntot_particles_old = ntot_particles
        qs%estimators(1)%hf_signed_pop = new_hf_signed_pop

        associate (est=>qs%estimators)
            if (comp_param) then
                do i = 1, qs%psip_list%nspaces, 2
                    if (.not. qs%vary_shift(i) .and. sum(ntot_particles(i:i+1)) > qs%target_particles) then
                        qs%vary_shift(i) = .true.
                        if (qmc_in%vary_shift_from_proje) then
                            ! Set shift to be instantaneous projected energy.
                            qs%shift(i) = real(est(i)%proj_energy_comp/est(i)%D0_population_comp, p)
                        else
                            qs%shift(i) = qmc_in%vary_shift_from
                        end if
                        qs%shift(i+1) = qs%shift(i)
                    end if
                end do
            else if (doing_calc(hfs_fciqmc_calc)) then
                if (.not. qs%vary_shift(1) .and. ntot_particles(1) > qs%target_particles) then
                    qs%vary_shift = .true.
                    if (qmc_in%vary_shift_from_proje) then
                        qs%shift(1) = est(1)%proj_energy/est(1)%D0_population
                        qs%shift(2) = est(1)%proj_hf_O_hpsip/est(1)%D0_population + est(1)%proj_hf_H_hfpsip/est(1)%D0_population &
                                                             - (est(1)%proj_energy*est(1)%D0_hf_population)/est(1)%D0_population**2
                    else
                        qs%shift = qmc_in%vary_shift_from
                    end if
                end if
            else
                do i = 1, qs%psip_list%nspaces 
                    if (.not. qs%vary_shift(i)) then
                        if ((ntot_particles(i) > qs%target_particles .and. .not. qmc_in%target_reference) .or. &
                            (est(i)%D0_population > qs%target_particles .and. qmc_in%target_reference)) then
                            qs%vary_shift(i) = .true.
                            if (qmc_in%vary_shift_from_proje) then
                                ! Set shift to be instantaneous projected energy.
                                qs%shift(i) = est(i)%proj_energy/est(i)%D0_population
                            else
                                qs%shift(i) = qmc_in%vary_shift_from
                            end if
                        end if
                    end if
                end do
            end if
        end associate

    end subroutine communicated_energy_estimators

!--- Shift updates ---

    subroutine update_shift(qs, loc_shift, nparticles_old, nparticles, nupdate_steps)

        ! Update the shift according to Yang, Pahl and Brand's harmonic damping
        ! population control algorithm (J. Chem. Phys. 153, 174103 (2020)
        ! (DOI:10.1063/5.0023088)):
        !  shift(beta) = shift(beta-A*tau) - xi*log(N_w(beta)/N_w(beta-A*tau))/(A*tau) 
        ! - zeta*(log(N_w(beta)/N_t))/(A*tau) 
        ! where
        !  * shift(beta) is the shift at imaginary time beta;
        !  * A*tau is the amount of imaginary time between shift-updates (=# of
        !    Monte Carlo cycles between updating the shift);
        !  * xi is a damping factor (0.05-0.10 is appropriate) to prevent large fluctations;
        !  * N_w(beta) is the total number of particles at imaginary time beta.
        !  * zeta is a dimensionless forcing strength parameter. 
        !  * N_t is the target population of the simulation. 
        ! The running average of the shift is also updated.
        ! Note that when zeta is equal to zero, the original algorithm is
        ! obtained. 

        ! In:
        !    qs: qmc state.
        !    nparticles_old: N_w(beta-A*tau).
        !    nparticles: N_w(beta).
        !    nupdate_steps: number of iterations between shift updates.
        ! In/Out:
        !    loc_shift: energy shift/offset.  Set to S(beta-A*tau) on input and S(beta) on output.

        use qmc_data, only: qmc_in_t, qmc_state_t
 
        type(qmc_state_t), intent(in) :: qs
        real(p), intent(inout) :: loc_shift
        real(dp), intent(in) :: nparticles_old, nparticles
        integer, intent(in) :: nupdate_steps

        ! dmqmc_factor is included to account for a factor of 1/2 introduced into tau in
        ! DMQMC calculations. In all other calculation types, it is set to 1, and so can be ignored.
        if (qs%target_particles .le. 0.00_p) then
            loc_shift = loc_shift - real(log(nparticles/nparticles_old)*qs%shift_damping/(qs%dmqmc_factor*qs%tau*nupdate_steps) ,p)
        else 
            loc_shift = loc_shift - real(log(nparticles/nparticles_old)*qs%shift_damping/(qs%dmqmc_factor*qs%tau*nupdate_steps),p) & 
                  - real(log(nparticles/qs%target_particles)*(qs%shift_harmonic_forcing)/ &
                  (qs%dmqmc_factor*qs%tau*nupdate_steps),p)
        end if 
    
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

    pure subroutine update_proj_energy_hub_k(sys, f0, wfn_dat, cdet, pop, estimators, excitation, hmatel)

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
        !    estimators: estimators_t object containing running totals of N_0
        !        and proj energy contribution.
        !    excitation: excitation connecting the determinant to the reference determinant.
        ! Out:
        !    hmatel: <D_i|H|D_0>, the Hamiltonian matrix element between the
        !       determinant and the reference determinant.

        ! NOTE: it is the programmer's responsibility to ensure D0_pop_sum and
        ! proj_energy_sum are zero before the first call.

        use determinant_data, only: det_info_t
        use excitations, only: excit_t
        use hamiltonian_hub_k, only: slater_condon2_hub_k
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f0(:)
        real(p), intent(in) :: wfn_dat(:)
        type(det_info_t), intent(in) :: cdet
        real(p), intent(in) :: pop(:)
        type(estimators_t), intent(inout) :: estimators
        type(excit_t), intent(inout) :: excitation
        type(hmatel_t), intent(out) :: hmatel

        if (excitation%nexcit == 0) then
            ! Have reference determinant.
            estimators%D0_population = estimators%D0_population + pop(1)
        else if (excitation%nexcit == 2) then
            ! Have a determinant connected to the reference determinant: add to
            ! projected energy.
            hmatel%r = slater_condon2_hub_k(sys, excitation%from_orb(1), excitation%from_orb(2), &
                                       & excitation%to_orb(1), excitation%to_orb(2),excitation%perm)
            estimators%proj_energy = estimators%proj_energy + hmatel%r*pop(1)
        end if

    end subroutine update_proj_energy_hub_k

    pure subroutine update_proj_energy_hub_real(sys, f0, wfn_dat, cdet, pop, estimators, excitation, hmatel)

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
        !    estimators: estimators_t object containing running totals of N_0
        !        and proj energy contribution.
        !    excitation: excitation connecting the determinant to the reference determinant.
        ! Out:
        !    hmatel: <D_i|H|D_0>, the Hamiltonian matrix element between the
        !       determinant and the reference determinant.

        ! NOTE: it is the programmer's responsibility to ensure D0_pop_sum and
        ! proj_energy_sum are zero before the first call.

        use determinant_data, only: det_info_t
        use excitations, only: excit_t
        use hamiltonian_hub_real, only: slater_condon1_hub_real
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f0(:)
        real(p), intent(in) :: wfn_dat(:)
        type(det_info_t), intent(in) :: cdet
        real(p), intent(in) :: pop(:)
        type(estimators_t), intent(inout) :: estimators
        type(excit_t), intent(inout) :: excitation
        type(hmatel_t), intent(out) :: hmatel

        if (excitation%nexcit == 0) then
            ! Have reference determinant.
            estimators%D0_population = estimators%D0_population + pop(1)
        else if (excitation%nexcit == 1) then
            ! Have a determinant connected to the reference determinant: add to
            ! projected energy.
            hmatel%r = slater_condon1_hub_real(sys, excitation%from_orb(1), excitation%to_orb(1), excitation%perm)
            estimators%proj_energy = estimators%proj_energy + hmatel%r*pop(1)
        end if

    end subroutine update_proj_energy_hub_real

    pure subroutine update_proj_energy_mol(sys, f0, wfn_dat, cdet, pop, estimators, excitation, hmatel)

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
        !    estimators: estimators_t object containing running totals of N_0
        !        and proj energy contribution.
        !    excitation: excitation connecting the determinant to the reference determinant.
        ! Out:
        !    hmatel: <D_i|H|D_0>, the Hamiltonian matrix element between the
        !       determinant and the reference determinant.

        ! NOTE: it is the programmer's responsibility to ensure D0_pop_sum and
        ! proj_energy_sum are zero before the first call.

        use determinant_data, only: det_info_t
        use excitations, only: excit_t
        use hamiltonian_molecular, only: slater_condon1_mol_excit, slater_condon2_mol_excit
        use read_in_symmetry, only: cross_product_basis_read_in
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f0(:)
        real(p), intent(in) :: wfn_dat(:)
        type(det_info_t), intent(in) :: cdet
        real(p), intent(in) :: pop(:)
        type(estimators_t), intent(inout) :: estimators
        type(excit_t), intent(inout) :: excitation
        type(hmatel_t), intent(out) :: hmatel

        integer :: ij_sym, ab_sym

        hmatel%r = 0.0_p

        select case(excitation%nexcit)
        case (0)
            ! Have reference determinant.
            estimators%D0_population = estimators%D0_population + pop(1)
        case(1)
            ! Have a determinant connected to the reference determinant by
            ! a single excitation: add to projected energy.
            ! Is excitation symmetry allowed?
            if (sys%basis%basis_fns(excitation%from_orb(1))%Ms == sys%basis%basis_fns(excitation%to_orb(1))%Ms .and. &
                    sys%basis%basis_fns(excitation%from_orb(1))%sym == sys%basis%basis_fns(excitation%to_orb(1))%sym) then
                hmatel = slater_condon1_mol_excit(sys, cdet%occ_list, excitation%from_orb(1), excitation%to_orb(1), &
                                                  excitation%perm)
                estimators%proj_energy = estimators%proj_energy + hmatel%r*pop(1)
            end if
        case(2)
            ! Have a determinant connected to the reference determinant by
            ! a double excitation: add to projected energy.
            ! Is excitation symmetry allowed?
            if (sys%basis%basis_fns(excitation%from_orb(1))%Ms+sys%basis%basis_fns(excitation%from_orb(2))%Ms == &
                    sys%basis%basis_fns(excitation%to_orb(1))%Ms+sys%basis%basis_fns(excitation%to_orb(2))%Ms) then
                ij_sym = cross_product_basis_read_in(sys, excitation%from_orb(1), excitation%from_orb(2))
                ab_sym = cross_product_basis_read_in(sys, excitation%to_orb(1), excitation%to_orb(2))
                if (ij_sym == ab_sym) then
                    hmatel = slater_condon2_mol_excit(sys, excitation%from_orb(1), excitation%from_orb(2), &
                                                      excitation%to_orb(1), excitation%to_orb(2),     &
                                                      excitation%perm)
                    estimators%proj_energy = estimators%proj_energy + hmatel%r*pop(1)
                end if
            end if
        end select

    end subroutine update_proj_energy_mol

    pure subroutine update_proj_energy_periodic_complex(sys, f0, wfn_dat, cdet, pop, estimators, excitation, hmatel)

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
        !    estimators: estimators_t object containing running totals of N_0
        !        and proj energy contribution.
        !    excitation: excitation connecting the determinant to the reference determinant.
        ! Out:
        !    hmatel_comp: <D_i|H|D_0>, the Hamiltonian matrix element between the
        !       determinant and the reference determinant.

        ! NOTE: it is the programmer's responsibility to ensure D0_pop_sum and
        ! proj_energy_sum are zero before the first call.

        use determinant_data, only: det_info_t
        use excitations, only: excit_t
        use hamiltonian_periodic_complex, only: slater_condon1_periodic_excit_complex, &
                                                slater_condon2_periodic_excit_complex
        use read_in_symmetry, only: cross_product_basis_read_in
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f0(:)
        real(p), intent(in) :: wfn_dat(:)
        type(det_info_t), intent(in) :: cdet
        real(p), intent(in) :: pop(:)
        type(estimators_t), intent(inout) :: estimators
        type(excit_t), intent(inout) :: excitation
        type(hmatel_t), intent(out) :: hmatel

        integer :: ij_sym, ab_sym

        hmatel%c = cmplx(0.0, 0.0, p)

        select case(excitation%nexcit)
        case (0)
            ! Have reference determinant.
            estimators%D0_population_comp = estimators%D0_population_comp + cmplx(pop(1),pop(2),p)
        case(1)
            ! Have a determinant connected to the reference determinant by
            ! a single excitation: add to projected energy.
            ! Is excitation symmetry allowed?
            if (sys%basis%basis_fns(excitation%from_orb(1))%Ms == sys%basis%basis_fns(excitation%to_orb(1))%Ms .and. &
                    sys%basis%basis_fns(excitation%from_orb(1))%sym == sys%basis%basis_fns(excitation%to_orb(1))%sym) then
                hmatel = slater_condon1_periodic_excit_complex(sys, cdet%occ_list, excitation%from_orb(1), excitation%to_orb(1), &
                                                  excitation%perm)
                estimators%proj_energy_comp = estimators%proj_energy_comp + hmatel%c*cmplx(pop(1),pop(2),p)
            end if
        case(2)
            ! Have a determinant connected to the reference determinant by
            ! a double excitation: add to projected energy.
            ! Is excitation symmetry allowed?
            if (sys%basis%basis_fns(excitation%from_orb(1))%Ms+sys%basis%basis_fns(excitation%from_orb(2))%Ms == &
                    sys%basis%basis_fns(excitation%to_orb(1))%Ms+sys%basis%basis_fns(excitation%to_orb(2))%Ms) then
                ij_sym = cross_product_basis_read_in(sys, excitation%from_orb(1), excitation%from_orb(2))
                ab_sym = cross_product_basis_read_in(sys, excitation%to_orb(1), excitation%to_orb(2))
                if (ij_sym == ab_sym) then
                    hmatel = slater_condon2_periodic_excit_complex(sys, excitation%from_orb(1), excitation%from_orb(2), &
                                                      excitation%to_orb(1), excitation%to_orb(2),     &
                                                      excitation%perm)
                    estimators%proj_energy_comp = estimators%proj_energy_comp + hmatel%c*cmplx(pop(1),pop(2),p)
                end if
            end if
        end select

    end subroutine update_proj_energy_periodic_complex

    pure subroutine update_proj_energy_ueg(sys, f0, wfn_dat, cdet, pop, estimators, excitation, hmatel)

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
        !    estimators: estimators_t object containing running totals of N_0
        !        and proj energy contribution.
        !    excitation: excitation connecting the determinant to the reference determinant.
        ! Out:
        !    hmatel: <D_i|H|D_0>, the Hamiltonian matrix element between the
        !       determinant and the reference determinant.

        ! NOTE: it is the programmer's responsibility to ensure D0_pop_sum and
        ! proj_energy_sum are zero before the first call.

        use determinant_data, only: det_info_t
        use excitations, only: excit_t
        use hamiltonian_ueg, only: slater_condon2_ueg
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f0(:)
        real(p), intent(in) :: wfn_dat(:)
        type(det_info_t), intent(in) :: cdet
        real(p), intent(in) :: pop(:)
        type(estimators_t), intent(inout) :: estimators
        type(excit_t), intent(inout) :: excitation
        type(hmatel_t), intent(out) :: hmatel

        hmatel%r = 0.0_p

        if (excitation%nexcit == 0) then
            ! Have reference determinant.
            estimators%D0_population = estimators%D0_population + pop(1)
        else if (excitation%nexcit == 2) then
            ! Have a determinant connected to the reference determinant: add to
            ! projected energy.
            hmatel%r = slater_condon2_ueg(sys, excitation%from_orb(1), excitation%from_orb(2), &
                                       & excitation%to_orb(1), excitation%to_orb(2),excitation%perm)
            estimators%proj_energy = estimators%proj_energy + hmatel%r*pop(1)
        end if

    end subroutine update_proj_energy_ueg

    pure subroutine update_proj_energy_ringium(sys, f0, wfn_dat, cdet, pop, estimators, excitation, hmatel)

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
        !    estimators: estimators_t object containing running totals of N_0
        !        and proj energy contribution.
        ! Out:
        !    excitation: excitation connecting the determinant to the reference determinant.
        !    hmatel: <D_i|H|D_0>, the Hamiltonian matrix element between the
        !       determinant and the reference determinant.

        ! NOTE: it is the programmer's responsibility to ensure D0_pop_sum and
        ! proj_energy_sum are zero before the first call.

        use determinant_data, only: det_info_t
        use excitations, only: excit_t, get_excitation
        use hamiltonian_ringium, only: slater_condon2_ringium
        use system, only: sys_t
        use bit_utils, only: count_set_bits

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f0(:)
        real(p), intent(in) :: wfn_dat(:)
        type(det_info_t), intent(in) :: cdet
        real(p), intent(in) :: pop(:)
        type(estimators_t), intent(inout) :: estimators
        type(excit_t), intent(inout) :: excitation
        type(hmatel_t), intent(out) :: hmatel

        excitation = get_excitation(sys%nel, sys%basis, cdet%f, f0)
        hmatel%r = 0.0_p

        if (excitation%nexcit == 0) then
            ! Have reference determinant.
            estimators%D0_population = estimators%D0_population + pop(1)
        else if (excitation%nexcit == 2) then
            ! Have a determinant connected to the reference determinant: add to
            ! projected energy.
            hmatel%r = slater_condon2_ringium(sys, excitation%from_orb(1), excitation%from_orb(2), &
                                       & excitation%to_orb(1), excitation%to_orb(2),excitation%perm)
            estimators%proj_energy = estimators%proj_energy + hmatel%r*pop(1)
        end if

    end subroutine update_proj_energy_ringium


    subroutine update_proj_hfs_hamiltonian(sys, f, fpop, f_hfpop, fdata, excitation, hmatel, &
                                           D0_hf_pop, proj_hf_O_hpsip, proj_hf_H_hfpsip)

        ! Add the contribution of the current determinant to the projected
        ! energy in an identical way to update_proj_energy_hub_k.

        ! Also add the contribution of the current determinant to the running
        ! total of the projected Hellmann--Feynman estimators.

        ! For debugging purposes, this procedure is for when we are sampling
        ! O=H.

        ! This procedure is for the Hubbard model in momentum space only.

        ! In:
        !    sys: system being studied.  Unused.
        !    f(tot_string_len): bit string representation of the Slater determinant, D_i.
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
        real(p), intent(in) :: fdata(:)
        type(hmatel_t), intent(in) :: hmatel
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
            proj_hf_O_hpsip = proj_hf_O_hpsip + hmatel%r*fpop
            ! \sum_j H_0j \tilde{c}_j is similarly easy to evaluate
            proj_hf_H_hfpsip = proj_hf_H_hfpsip + hmatel%r*f_hfpop
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
        !    f(tot_string_len): bit string representation of the Slater determinant, D_i
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
        real(p), intent(in) :: fdata(:)
        type(hmatel_t), intent(in) :: hmatel
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
            proj_hf_H_hfpsip = proj_hf_H_hfpsip + hmatel%r*f_hfpop
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
        !    f(tot_string_len): bit string representation of the Slater determinant, D_i
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
        real(p), intent(in) :: fdata(:)
        type(hmatel_t), intent(in) :: hmatel
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
            proj_hf_O_hpsip = proj_hf_O_hpsip + (hmatel%r/(sys%hubbard%u*sys%lattice%nsites))*fpop

            ! \sum_j H_0j \tilde{c}_j is similarly easy to evaluate
            proj_hf_H_hfpsip = proj_hf_H_hfpsip + hmatel%r*f_hfpop
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
        !    f(tot_string_len): bit string representation of the Slater determinant, D_i
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
        real(p), intent(in) :: fdata(:)
        type(hmatel_t), intent(in) :: hmatel
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
            proj_hf_H_hfpsip = proj_hf_H_hfpsip + hmatel%r*f_hfpop
        case(2)
            ! O is a one-body operator => no contributions from double
            ! excitations to \sum_j O_0j c_j.

            ! \sum_j H_0j \tilde{c}_j
            proj_hf_H_hfpsip = proj_hf_H_hfpsip + hmatel%r*f_hfpop
        end select

    end subroutine update_proj_hfs_one_body_mol

    pure function get_sanitized_projected_energy(qs) result(proje)
        
        ! From a qmc_state, qs, return either the value of the projected energy,
        ! or 0 if this is undefined.

        ! In:
        !    qs: qmc state containing estimators.

        ! Returns:
        !   real containing (instantaneous) projected energy.

        use qmc_data, only: qmc_state_t
        type(qmc_state_t), intent(in) :: qs
        real(p) ::  proje(size(qs%estimators))
 
        where (abs(qs%estimators%D0_population) < tiny(0.0_p))
            proje = 0.0_p
        elsewhere
            proje = qs%estimators%proj_energy/qs%estimators%D0_population
        end where

    end function get_sanitized_projected_energy 

    pure function get_sanitized_projected_energy_cmplx(qs) result(proje)

        ! From a qmc_state, qs, return either the value of the projected energy,
        ! or 0 if this is undefined.
        ! Returns the real component of the projected energy. Since we have assumed
        ! elsewhere the Hamiltonian is hermitian, we can safely assume a real energy
        ! estimator.

        ! In:
        !    qs: qmc state containing estimators.

        ! Returns:
        !   real containing real component of (complex, instantaneous) projected energy.

        use qmc_data, only: qmc_state_t
        type(qmc_state_t), intent(in) :: qs
        real(p) ::  proje(size(qs%estimators))

        where (abs(qs%estimators%D0_population_comp)<tiny(0.0_p))
            proje = 0.0_p
        elsewhere
            proje = real(qs%estimators%proj_energy_comp/qs%estimators%D0_population_comp, p)
        end where

    end function get_sanitized_projected_energy_cmplx

end module energy_evaluation

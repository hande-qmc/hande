module trotterized_uccmc

! Module for performing trotterized (disentangled) unitary coupled cluster Monte Carlo (UCCMC) calculations.

! The approach is very similar to conventional CCMC and many routines are shared between the
! two algorithms. For details of general computational considerations, see comments in ccmc.f90.

! The main difference between unitary and conventional CC stems from the form of the exponential
! ansatz. In CC, the wavefunction is given by
!
!     \Psi_CC = exp(T)|D_0>
!
! and the cluster operator, T, is::
!
!     T = \sum_{ia} t_i^a(t) a_i^a + 1/2!^2 \sum_{ijab} t_{ij}^{ab}(t) a_{ij}^{ab} + ... 
!
! In UCC, the wavefunction becomes:: 
!
!     \Psi_CC = exp(T - T^dagger)|D_0>
!
! For tUCC, we take the further approximation that::
!
!    exp(\sum t) = \prod exp(t)
!
! This introduces an additional degree of freedom in the ordering of operators in the tUCC expansion which
! can affect the quality of the results.

use const, only: i0, int_p, int_64, p, dp, debug

implicit none

contains

    subroutine do_trot_uccmc(sys, qmc_in, ccmc_in, uccmc_in, restart_in, load_bal_in, reference_in, &
                        logging_in, io_unit, qs, qmc_state_restart)

        ! This subroutine is derived from do_ccmc in ccmc.f90. 
        ! Some of the features of CCMC have not yet been implemented for UCCMC: even selection, blocking_in etc.

        ! Run the Trotterized UCCMC algorithm starting from the initial walker distribution
        ! using the timestep algorithm.

        ! See notes about the implementation of this using function pointers
        ! in fciqmc_main.

        ! In:
        !    sys: system being studied.
        !    qmc_in: input options relating to QMC methods.
        !    ccmc_in: input options relating to CCMC.
        !    uccmc_in: input options relating to UCCMC.
        !    restart_in: input options for HDF5 restart files.
        !    load_bal_in: input options for load balancing.
        !    reference_in: current reference determinant.  If not set (ie
        !       components allocated) then a best guess is made based upon the
        !       desired spin/symmetry.
        !    logging_in: input options for debug logs.
        !    io_unit: input unit to write all output to.
        ! In/Out:
        !    qmc_state_restart (optional): if present, restart from a previous calculation.
        !       Deallocated on exit.
        ! Out:
        !    qs: qmc_state for use if restarting the calculation

        use checking, only: check_allocate, check_deallocate
        use dSFMT_interface, only: dSFMT_t, dSFMT_init, dSFMT_end, dSFMT_state_t_to_dSFMT_t, dSFMT_t_to_dSFMT_state_t, &
                                   free_dSFMT_state_t
        use errors, only: stop_all, warning
        use parallel
        use restart_hdf5, only: dump_restart_hdf5, restart_info_t, init_restart_info_t, dump_restart_file_wrapper

        use annihilation, only: insert_new_walker, direct_annihilation
        use bloom_handler, only: init_bloom_stats_t, bloom_stats_t, bloom_mode_fractionn, bloom_mode_fixedn, &
                                 write_bloom_report, bloom_stats_warning, update_bloom_threshold_prop
        use ccmc, only: do_ccmc_accumulation, perform_ccmc_spawning_attempt, do_stochastic_ccmc_propagation, &
                        do_nc_ccmc_propagation
        use ccmc_data
        use ccmc_death_spawning, only: stochastic_ccmc_death_nc
        use ccmc_selection, only: create_null_cluster
        use ccmc_selection, only: init_selection_data, update_selection_probabilities, set_cluster_selections, &
                                  init_amp_psel_accumulation, select_nc_cluster
        use ccmc_utils, only: get_D0_info, init_contrib, dealloc_contrib, cumulative_population, & 
                              regenerate_ex_levels_psip_list
        use determinants, only: alloc_det_info_t, dealloc_det_info_t, sum_sp_eigenvalues_occ_list, &
                                sum_fock_values_bit_string, decode_det
        use determinant_data, only: det_info_t
        use excitations, only: excit_t, get_excitation_level, get_excitation
        use qmc_io, only: write_qmc_report, write_qmc_report_header, write_qmc_var
        use qmc, only: init_qmc, init_secondary_references
        use qmc_common, only: initial_qmc_status, initial_cc_projected_energy, load_balancing_report, init_report_loop, &
                              init_mc_cycle, end_report_loop, end_mc_cycle, redistribute_particles, rescale_tau
        use proc_pointers
        use spawning, only: assign_particle_processor
        use system, only: sys_t, sys_t_json
        use spawn_data, only: calc_events_spawn_t, write_memcheck_report, spawn_t, alloc_spawn_t
        use replica_rdm, only: update_rdm, calc_rdm_energy, write_final_rdm

        use qmc_data, only: qmc_in_t, ccmc_in_t, uccmc_in_t, restart_in_t

        use qmc_data, only: load_bal_in_t, qmc_state_t, annihilation_flags_t, estimators_t, particle_t
        use qmc_data, only: qmc_in_t_json, ccmc_in_t_json, uccmc_in_t_json, restart_in_t_json
        use qmc_data, only: excit_gen_power_pitzer_orderN, excit_gen_heat_bath
        use reference_determinant, only: reference_t, reference_t_json
        use check_input, only: check_qmc_opts, check_uccmc_opts
        use json_out, only: json_out_t, json_object_init, json_object_end
        use hamiltonian_data
        use energy_evaluation, only: get_sanitized_projected_energy, get_sanitized_projected_energy_cmplx

        use logging, only: init_logging, end_logging, prep_logging_mc_cycle, write_logging_calc_ccmc
        use logging, only: logging_in_t, logging_t, logging_in_t_json, logging_t_json, write_logging_select_ccmc
        use report, only: write_date_time_close
        use excit_gens, only: p_single_double_coll_t
        use particle_t_utils, only: init_particle_t
        use search, only: binary_search
        use uccmc_utils, only: add_info_str_trot, add_ci_contribution, allocate_time_average_lists, &
                               var_energy_uccmc, initialise_average_wfn, write_average_wfn_trot, add_t_contributions
        type(sys_t), intent(in) :: sys
        type(qmc_in_t), intent(in) :: qmc_in
        type(ccmc_in_t), intent(in) :: ccmc_in
        type(uccmc_in_t), intent(in) :: uccmc_in
        type(restart_in_t), intent(in) :: restart_in
        type(load_bal_in_t), intent(in) :: load_bal_in
        type(reference_t), intent(in) :: reference_in
        type(logging_in_t), intent(in) :: logging_in
        type(qmc_state_t), target, intent(out) :: qs
        type(qmc_state_t), intent(inout), optional :: qmc_state_restart
        integer, intent(in) :: io_unit

        integer :: i, ireport, icycle, iter, it
        integer(int_64) :: iattempt
        integer(int_64) :: nattempts_spawn
        real(dp), allocatable :: nparticles_old(:), nparticles_change(:)
        type(det_info_t) :: ref_det

        integer(int_p) :: ndeath, ndeath_nc
        integer :: nspawn_events, ierr
        type(wfn_contrib_t), allocatable :: contrib(:)
        type(multispawn_stats_t), allocatable :: ms_stats(:)
        type(p_single_double_coll_t), allocatable :: ps_stats(:)
        type(dSFMT_t), allocatable :: rng(:)
        type(json_out_t) :: js
        type(qmc_in_t) :: qmc_in_loc
        type(logging_t) :: logging_info
        type(selection_data_t) :: selection_data

        logical :: soft_exit, dump_restart_shift, restarting, restart_proj_est

        real(p), allocatable :: cumulative_abs_real_pops(:)
        integer :: D0_proc, D0_pos, nD0_proc, min_cluster_size, max_cluster_size, iexcip_pos
        real(p) :: tot_abs_real_pop
        complex(p) :: D0_normalisation
        type(bloom_stats_t) :: bloom_stats
        type(annihilation_flags_t) :: annihilation_flags
        type(restart_info_t) :: ri, ri_shift
        character(36) :: uuid_restart

        real :: t1, t2

        logical :: update_tau, error

        logical :: seen_D0, regenerate_info
        real(p) :: dfock
        complex(p) :: D0_population_cycle, proj_energy_cycle
        real(p) :: D0_population_noncomp_cycle

        real(p), allocatable :: rdm(:,:)

        integer :: iunit, restart_version_restart
        integer :: date_values(8)
        character(:), allocatable :: err_msg
 
        real(p), allocatable :: time_avg_psip_list_sq(:,:), time_avg_psip_list_pops(:), time_avg_psip_list_ci_pops(:)
        integer(i0), allocatable :: time_avg_psip_list_states(:,:), time_avg_psip_list_ci_states(:,:)
        integer :: nstates_sq, nstates_ci
        integer :: semi_stoch_it, pos, j, k
        logical :: hit
        integer(i0), allocatable :: state(:)
        real(p) :: population
        real(p) :: real_population, var_energy
        logical :: old_vary
        integer :: avg_start
        real(p) :: p_ref, pcumul
        real(p) :: p_avg, p_complement_avg
        integer :: count_select
        logical :: variational
        real(p) :: cluster_pop

        variational = .false.
        if (uccmc_in%variational_energy) variational = .true.
        ! old_vary encodes whether the shift started varying before the current iteration. Used to determine when
        ! to start taking averages of the wavefunction.
        old_vary = .false.
        avg_start = 0
        count_select = 0

        if (parent) then
            write (io_unit,'(1X,"Trotterized UCCMC")')
            write (io_unit,'(1X,"----",/)')
        end if

        restarting = present(qmc_state_restart) .or. restart_in%read_restart
        ! Check input options.
        if (parent) then
            if (present(qmc_state_restart)) then
                call check_qmc_opts(qmc_in, sys, .not.present(qmc_state_restart), restarting, &
                    qmc_state_restart=qmc_state_restart)
            else
                call check_qmc_opts(qmc_in, sys, .not.present(qmc_state_restart), restarting)
            end if
            call check_uccmc_opts(sys, ccmc_in, uccmc_in, qmc_in)
        end if
        ! Initialise data.
        call init_qmc(sys, qmc_in, restart_in, load_bal_in, reference_in, io_unit, annihilation_flags, qs, &
                      uuid_restart, restart_version_restart, qmc_state_restart=qmc_state_restart, &
                      regenerate_info=regenerate_info)

        ! Add information strings to the psip_list and the reference determinant.
        call regenerate_trot_info_psip_list(sys%basis, sys%nel, qs)
        qs%psip_list%descending=.true.
        call add_trot_info_reference(qs%ref%f0, sys)

        !Allocate memory for time averaged populations and variational energy computation.
        allocate(state(sys%basis%tot_string_len))

        if(uccmc_in%variational_energy) then
             population = 0
             call allocate_time_average_lists(qs%psip_list, time_avg_psip_list_ci_states, time_avg_psip_list_ci_pops, nstates_ci)
             time_avg_psip_list_ci_states(:,1) = qs%psip_list%states(:,1)
             time_avg_psip_list_ci_pops(1) = (real(qs%psip_list%pops(1,1))/qs%psip_list%pop_real_factor)
             nstates_ci = 1
             ![todo] deal with restarting
        end if

        if (uccmc_in%average_wfn) then
            call allocate_time_average_lists(qs%psip_list, time_avg_psip_list_states, time_avg_psip_list_pops, nstates_sq)
            allocate(time_avg_psip_list_sq(sys%basis%tot_string_len + 1 ,size(qs%psip_list%states(1,:))))
            time_avg_psip_list_sq(:sys%basis%tot_string_len,:qs%psip_list%nstates) = qs%psip_list%states(:,:qs%psip_list%nstates)
            time_avg_psip_list_sq(sys%basis%tot_string_len+1,:qs%psip_list%nstates) = &
                (real(qs%psip_list%pops(1,:qs%psip_list%nstates))/qs%psip_list%pop_real_factor)**2
            nstates_sq = 1
        end if
        
        if (ccmc_in%multiref) then
            ! Initialise multireference CCMC specific data.
            qs%multiref = .true.
            qs%mr_acceptance_search = ccmc_in%mr_acceptance_search
            qs%n_secondary_ref = ccmc_in%n_secondary_ref
            if(ccmc_in%mr_read_in) then
                qs%mr_read_in = ccmc_in%mr_read_in
                qs%mr_secref_file = ccmc_in%mr_secref_file
                qs%mr_n_frozen = ccmc_in%mr_n_frozen
                qs%mr_excit_lvl = ccmc_in%mr_excit_lvl
            endif

            allocate (qs%secondary_refs(qs%n_secondary_ref))
            call init_secondary_references(sys, ccmc_in%secondary_refs, io_unit, qs)
        else 
            qs%ref%max_ex_level = qs%ref%ex_level
        end if

        if (debug) call init_logging(logging_in, logging_info, qs%ref%ex_level)

        if (parent) then
            call json_object_init(js, tag=.true., io=io_unit)
            call sys_t_json(js, sys)
            ! The default values of pattempt_* are not in qmc_in
            qmc_in_loc = qmc_in
            ! [todo] -  This is repeated in DMQMC (and FCIQMC?).  Should it be a subroutine?
            qmc_in_loc%pattempt_single = qs%excit_gen_data%pattempt_single
            qmc_in_loc%pattempt_double = qs%excit_gen_data%pattempt_double
            qmc_in_loc%shift_damping = qs%shift_damping
            qmc_in_loc%pattempt_parallel = qs%excit_gen_data%pattempt_parallel
            call qmc_in_t_json(js, qmc_in_loc)
            call ccmc_in_t_json(js, ccmc_in)
            call uccmc_in_t_json(js, uccmc_in)
            call restart_in_t_json(js, restart_in, uuid_restart)
            call reference_t_json(js, qs%ref, sys)
            call logging_in_t_json(js, logging_in)
            call logging_t_json(js, logging_info, terminal=.true.)
            call json_object_end(js, terminal=.true., tag=.true.)
            write (js%io, '()')
        end if

        allocate(nparticles_old(qs%psip_list%nspaces), stat=ierr)
        call check_allocate('nparticles_old', qs%psip_list%nspaces, ierr)
        allocate(nparticles_change(qs%psip_list%nspaces), stat=ierr)
        call check_allocate('nparticles_change', qs%psip_list%nspaces, ierr)

        ! Initialise bloom_stats components to the following parameters.
        call init_bloom_stats_t(bloom_stats, mode=bloom_mode_fractionn, encoding_factor=qs%psip_list%pop_real_factor)

        ! Allocate and initialise per thread...
        allocate(rng(0:nthreads-1), stat=ierr)
        call check_allocate('rng', nthreads, ierr)
        allocate(ms_stats(0:nthreads-1), stat=ierr)
        call check_allocate('ms_stats', nthreads, ierr)
        allocate(ps_stats(0:nthreads-1), stat=ierr)
        call check_allocate('ps_stats', nthreads, ierr)

        call init_contrib(sys, uccmc_in%pow_trunc, ccmc_in%linked, contrib)

        do i = 0, nthreads-1
            ! Initialise and allocate RNG store.
            call dSFMT_init(qmc_in%seed+iproc+i*nprocs, 50000, rng(i))
        end do

        if (restart_in%restart_rng .and. allocated(qs%rng_state%dsfmt_state)) then
            call dSFMT_state_t_to_dSFMT_t(rng(0), qs%rng_state, err_msg=err_msg)
            if (allocated(err_msg)) call stop_all('do_trot_uccmc', 'Failed to reset RNG state: '//err_msg)
            call free_dSFMT_state_t(qs%rng_state)
        end if

        ! ...and scratch space for calculative cumulative probabilities.

        allocate(cumulative_abs_real_pops(size(qs%psip_list%states,dim=2)), stat=ierr)
        call check_allocate('cumulative_abs_real_pops', size(qs%psip_list%states, dim=2), ierr)

        if (debug) call init_amp_psel_accumulation(qs%ref%max_ex_level+2, logging_info, ccmc_in%linked, selection_data)

        nparticles_old = qs%psip_list%tot_nparticles

        ! Initialise D0_pos to be somewhere (anywhere) in the list.
        D0_pos = 1

        associate(pl=>qs%psip_list, spawn=>qs%spawn_store%spawn)
            ! Initialise hash shift if restarting...
            spawn%hash_shift = qs%mc_cycles_done
            ! NOTE: currently hash_seed is not exposed and so cannot change unless the hard-coded value changes. Therefore, as we
            ! have not evolved the particles since the were written out (i.e. hash_shift hasn't changed) the only parameter
            ! which can be altered which can change an excitors location since the restart files were written is move_freq.
            if (ccmc_in%move_freq /= spawn%move_freq .and. nprocs > 1) then
                spawn%move_freq = ccmc_in%move_freq
                if (restarting) then
                    if (parent) call warning('do_trot_uccmc', 'move_freq is different from that in the restart file. &
                                            &Reassigning processors. Please check for equilibration effects.')
                    ! Cannot rely on spawn%head being initialised to indicate an empty spawn_t object so reset it before
                    ! redistributing,.
                    spawn%head = spawn%head_start
                    call redistribute_particles(pl%states, pl%pop_real_factor, pl%pops, pl%nstates, pl%nparticles, spawn)
                    call direct_annihilation(sys, rng(0), qs%ref, annihilation_flags, pl, spawn)
                end if
            end if
        end associate

        if (parent) then
            call write_qmc_report_header(qs%psip_list%nspaces, cmplx_est=sys%read_in%comp, rdm_energy=ccmc_in%density_matrices, &
                                         nattempts=.true., io_unit=io_unit)
        end if

        restart_proj_est = present(qmc_state_restart) .or. (restart_in%read_restart .and. restart_version_restart >= 2)
        if (.not.restart_proj_est) then
            call initial_cc_projected_energy(sys, qs, qmc_in%seed+iproc, logging_info, cumulative_abs_real_pops, nparticles_old)
        end if

        call initial_qmc_status(sys, qmc_in, qs, nparticles_old, doing_ccmc=.true., io_unit=io_unit)

        ! Initialise timer.
        call cpu_time(t1)

        ! Should we dump a restart file just before the shift is turned on?
        dump_restart_shift = restart_in%write_restart_shift
        call init_restart_info_t(ri, write_id=restart_in%write_id)
        call init_restart_info_t(ri_shift, write_id=restart_in%write_shift_id)

        do ireport = 1, qmc_in%nreport

            ! Projected energy from last report loop to correct death
            qs%estimators%proj_energy_old = get_sanitized_projected_energy(qs)
                
            call init_report_loop(qs, bloom_stats)

            do icycle = 1, qmc_in%ncycles

                iter = qs%mc_cycles_done + (ireport-1)*qmc_in%ncycles + icycle

                if (all(qs%vary_shift) .and. (.not. old_vary) .and. (uccmc_in%average_wfn .or. uccmc_in%variational_energy)) then
                    ! On the first iteration after the shift has started varying, begin
                    ! storing average wavefunction.
                    old_vary = all(qs%vary_shift) 
                    avg_start = iter
                end if
                if (all(qs%vary_shift) .and. (.not. old_vary) .and. uccmc_in%average_wfn) then
                    time_avg_psip_list_pops(:qs%psip_list%nstates) = &
                        real(qs%psip_list%pops(1,:qs%psip_list%nstates))/qs%psip_list%pop_real_factor
                    time_avg_psip_list_states(:,:qs%psip_list%nstates) = qs%psip_list%states(:,:qs%psip_list%nstates)
                    time_avg_psip_list_sq(:sys%basis%tot_string_len,:qs%psip_list%nstates) &
                        = qs%psip_list%states(:,:qs%psip_list%nstates)
                    time_avg_psip_list_sq(sys%basis%tot_string_len+1,:qs%psip_list%nstates) &
                        = (real(qs%psip_list%pops(1,:qs%psip_list%nstates))/qs%psip_list%pop_real_factor)**2
                    nstates_sq = qs%psip_list%nstates
                end if
                  
                ! Recover total number of particles added to time-averages 
                if(uccmc_in%variational_energy .and. all(qs%vary_shift)) then
                          time_avg_psip_list_ci_pops(:nstates_ci) =  time_avg_psip_list_ci_pops(:nstates_ci)*(iter-avg_start)
                end if


                if (debug) call prep_logging_mc_cycle(iter, logging_in, logging_info, sys%read_in%comp, &
                                                        min(sys%nel, qs%ref%ex_level+2))

                call get_D0_info_trot(qs, sys%read_in%comp, D0_proc, D0_pos, nD0_proc, D0_normalisation)

                ! Update the shift of the excitor locations to be the end of this
                ! current iteration.
                qs%spawn_store%spawn%hash_shift = qs%spawn_store%spawn%hash_shift + 1

                ! Maximum possible cluster size that we can generate.
                ! If only the reference is populated, we can only generate clusters of size 0.
                ! Otherwise we can generate clusters up to our chosen truncation level.

                if(qs%psip_list%nstates-nD0_proc == 0) then
                    max_cluster_size = 0
                else
                    max_cluster_size = qs%psip_list%nstates-1
                end if

                call init_mc_cycle(qs%psip_list, qs%spawn_store%spawn, qs%estimators(1)%nattempts, ndeath, &
                                   min_attempts=nint(abs(D0_normalisation), kind=int_64), &
                                   complx=sys%read_in%comp)

                nparticles_change = 0.0_p
                ! We need to count spawning attempts differently as there may be multiple spawns
                ! per cluster

                nattempts_spawn=0

                call update_bloom_threshold_prop(bloom_stats, nparticles_old(1))

                ! Three options for evolution:

                ! * Original CCMC algorithm
                !       + The number of excips on this processor determines the number
                !         of cluster generations, each of which can spawn and die.
                !         non-composite clusters therefore are seldom selected.
                ! * 'full non-composite' algorithm, where spawning and death are split into two tranches.
                !       + non-composite clusters (i.e. consisting of a single excitor):
                !         enumerate explicitly (this is just the list of excitors)
                !       + composite clusters, which must be selected stochastically (as in
                !         the original algorithm for all clusters).  We sample the space
                !         of composite clusters, choosing nattempts samples.  For convenience
                !         nattempts = # excitors not on the reference (i.e. the number of
                !         excitors which can actually be involved in a composite cluster).
                ! * 'even selection' algorithm, where all clusters are selected with probability
                !         proportional to their contribution to the final wavefunction.
                !       + non-composite cluster enumerated as in full non-composite algorithm.
                !       + composite clusters more complicated selection probability required.

                !Initially for UCC we will simply use a modification of the original algorithm.
                call cumulative_population(qs%psip_list%pops, qs%psip_list%states(sys%basis%tot_string_len,:), &
                                           qs%psip_list%nstates, D0_proc, D0_pos, qs%psip_list%pop_real_factor, &
                                           ccmc_in%even_selection, sys%read_in%comp, cumulative_abs_real_pops, &
                                           tot_abs_real_pop)

                call set_cluster_selections(selection_data, qs%estimators(1)%nattempts, min_cluster_size, &
                                            max_cluster_size, D0_normalisation, tot_abs_real_pop, qs%psip_list%nstates, &
                                            ccmc_in%full_nc, .false.)
                call zero_ps_stats(ps_stats, qs%excit_gen_data%p_single_double%rep_accum%overflow_loc)

                ! OpenMP chunk size determined completely empirically from a single
                ! test.  Please feel free to improve...
                ! NOTE: we can't refer to procedure pointers in shared blocks so
                ! can't use default(none).  I *strongly* recommend turning
                ! default(none) on when making changes and ensure that the only
                ! errors relate to the procedure pointers...

                proj_energy_cycle = cmplx(0.0, 0.0, p)
                D0_population_cycle = cmplx(0.0, 0.0, p)
                D0_population_noncomp_cycle = 0.0_p
                !$omp parallel default(none) &
                !$omp private(it, iexcip_pos, i, seen_D0, hit, pos, population, real_population,k, cluster_pop, &
                !$omp state,annihilation_flags) &
                !$omp shared(rng, cumulative_abs_real_pops, tot_abs_real_pop,  &
                !$omp        max_cluster_size, contrib, D0_normalisation, D0_pos, rdm,    &
                !$omp        qs, sys, bloom_stats, min_cluster_size, ref_det,             &
                !$omp         selection_data,      &
                !$omp        uccmc_in, ccmc_in, nprocs, ms_stats, ps_stats, qmc_in, load_bal_in, &
                !$omp        ndeath_nc,   pcumul, nstates_ci, &
                !$omp        nparticles_change,logging_info, &
                !$omp        time_avg_psip_list_ci_states, time_avg_psip_list_ci_pops, &
                !$omp        time_avg_psip_list_states, time_avg_psip_list_pops) &
                !$omp reduction(+:D0_population_cycle,proj_energy_cycle, D0_population_noncomp_cycle, nattempts_spawn,ndeath)

                it = get_thread_id()
                iexcip_pos = 0
                seen_D0 = .false.
                pcumul = (tot_abs_real_pop**(max_cluster_size+1)-1) /(tot_abs_real_pop-1)
                cluster_pop = 1
                do i = 1, qs%psip_list%nstates
                if (i /= D0_pos) &
                    cluster_pop = cluster_pop * cos((real(qs%psip_list%pops(1, i))/real(qs%psip_list%pop_real_factor))&
                                                    /real(D0_normalisation,p))                     
                end do

                !$omp do schedule(dynamic,200) 
                do iattempt = 1, selection_data%nsingle_excitors
                    ! For OpenMP scalability, have this test inside a single loop rather
                    ! than attempt to parallelise over three separate loops.
                   
                        ! As noncomposite clusters can't be above truncation level or linked-only all can accumulate +
                        ! propagate. Only need to check not selecting the reference as we treat it separately.
                        if (iattempt /= D0_pos) then
                            ! Deterministically select each excip as a non-composite cluster.
                            call select_nc_cluster_trot(sys, qs%psip_list, qs%ref%f0, &
                                        iattempt, qmc_in%initiator_pop, .false., &
                                        contrib(it)%cdet, contrib(it)%cluster, qs%excit_gen_data, D0_normalisation, D0_pos)

                            if (qs%propagator%quasi_newton) contrib(it)%cdet%fock_sum = &
                                            sum_sp_eigenvalues_occ_list(sys, contrib(it)%cdet%occ_list) - qs%ref%fock_sum
                            ! [VAN]: This is quite dangerous when using OpenMP as selection_data is shared but updated here if
                            ! [VAN]: in debug mode. However, this updated selection_data will only be used if selection logging
                            ! [VAN]: according to comments. And logging cannot be used with openmp. Dangerous though.
                            call do_ccmc_accumulation(sys, qs, contrib(it)%cdet, contrib(it)%cluster, logging_info, &
                                                    D0_population_cycle, proj_energy_cycle, ccmc_in, ref_det, rdm, &
                                                    selection_data, D0_population_noncomp_cycle)
                            if (uccmc_in%variational_energy .and. all(qs%vary_shift) .and. & 
                                contrib(it)%cluster%excitation_level <= qs%ref%ex_level)  then
                                call add_ci_contribution(contrib(it)%cluster, contrib(it)%cdet, &
                                time_avg_psip_list_ci_states, time_avg_psip_list_ci_pops, nstates_ci)
                            end if
                            call do_nc_ccmc_propagation(rng(it), sys, qs, ccmc_in, logging_info, bloom_stats, &
                                                                contrib(it), nattempts_spawn, ps_stats(it), uccmc_in)
                        end if
                end do
                !$omp end do
                !$omp do schedule(dynamic,200) 
                do iattempt = 1, selection_data%nstochastic_clusters

                        call select_ucc_trot_cluster(rng(it), sys, qs%psip_list, qs%ref%f0, qs%ref%max_ex_level, &
                                            selection_data%nstochastic_clusters, &
                                            D0_normalisation, qmc_in%initiator_pop, D0_pos, cumulative_abs_real_pops,&
                                            tot_abs_real_pop, min_cluster_size, qs%psip_list%nstates-1, &
                                            logging_info, contrib(it)%cdet, contrib(it)%cluster, qs%excit_gen_data, cluster_pop)

                        call add_info_str_trot(sys%basis, qs%ref%f0, sys%nel, contrib(it)%cdet%f)

                            !Add selected cluster contribution to CI wavefunction estimator.
                            if (uccmc_in%variational_energy .and. (.not. all(contrib(it)%cdet%f==0)) .and. &
                                all(qs%vary_shift).and. &
                                contrib(it)%cluster%excitation_level <= qs%ref%ex_level)  then
                                call add_ci_contribution(contrib(it)%cluster, contrib(it)%cdet, &
                                time_avg_psip_list_ci_states, time_avg_psip_list_ci_pops, nstates_ci)
                            end if
                    

                        if (contrib(it)%cluster%excitation_level <= qs%ref%max_ex_level+2) then
                            ! cluster%excitation_level == huge(0) indicates a cluster
                            ! where two excitors share an elementary operator
                            if (qs%propagator%quasi_newton) contrib(it)%cdet%fock_sum = &
                                            sum_sp_eigenvalues_occ_list(sys, contrib(it)%cdet%occ_list) - qs%ref%fock_sum

                            call do_ccmc_accumulation(sys, qs, contrib(it)%cdet, contrib(it)%cluster, logging_info, &
                                                    D0_population_cycle, proj_energy_cycle, ccmc_in, ref_det, rdm, selection_data,&
                                                    D0_population_noncomp_cycle)
                            call do_stochastic_ccmc_propagation(rng(it), sys, qs, &
                                                                ccmc_in, logging_info, ms_stats(it), bloom_stats, &
                                                                contrib(it), nattempts_spawn, ndeath, ps_stats(it), uccmc_in)
                        end if
                end do
                !$omp end do

                !$omp do schedule(dynamic,200) 
                do iattempt = 1, selection_data%nD0_select
                        if (.not. seen_D0) then
                            ! This is the first time this thread is spawning from D0, so it
                            ! needs to be converted into a det_info_t object for the excitation
                            ! generators. On subsequent calls, cdet does not need to change.
                            seen_D0 = .true.
                            call create_null_cluster(sys, qs%ref%f0, nprocs*real(selection_data%nD0_select,p), &
                                                     D0_normalisation*cluster_pop, qmc_in%initiator_pop, contrib(it)%cdet, &
                                                     contrib(it)%cluster, qs%excit_gen_data)
                        end if
                        if (uccmc_in%variational_energy .and. all(qs%vary_shift) .and. &
                            contrib(it)%cluster%excitation_level <= qs%ref%ex_level)  then
                            call add_ci_contribution(contrib(it)%cluster, contrib(it)%cdet, &
                            time_avg_psip_list_ci_states, time_avg_psip_list_ci_pops, nstates_ci)
                        end if
                        if (qs%propagator%quasi_newton) contrib(it)%cdet%fock_sum = &
                                            sum_sp_eigenvalues_occ_list(sys, contrib(it)%cdet%occ_list) - qs%ref%fock_sum

                        call do_ccmc_accumulation(sys, qs, contrib(it)%cdet, contrib(it)%cluster, logging_info, &
                                                    D0_population_cycle, proj_energy_cycle, ccmc_in, ref_det, rdm, selection_data,&
                                                    D0_population_noncomp_cycle)
                        nattempts_spawn = nattempts_spawn + 1
                        call perform_ccmc_spawning_attempt(rng(it), sys, qs, ccmc_in, logging_info, bloom_stats, contrib(it), 1, &
                                                        ps_stats(it), uccmc_in)

                end do
                !$omp end do

                ndeath_nc=0
                if (ccmc_in%full_nc .and. qs%psip_list%nstates > 0) then
                    ! Do death exactly and directly for non-composite clusters
                    !$omp do schedule(dynamic,200) private(dfock) reduction(+:ndeath_nc,nparticles_change)
                    do iattempt = 1, qs%psip_list%nstates
                        ! Note we use the (encoded) population directly in stochastic_ccmc_death_nc
                        ! (unlike the stochastic_ccmc_death) to avoid unnecessary decoding/encoding
                        ! steps (cf comments in stochastic_death for FCIQMC).
                        if (qs%propagator%quasi_newton) then
                            dfock = sum_fock_values_bit_string(sys, qs%propagator%sp_fock, qs%psip_list%states(:,iattempt)) &
                                - qs%ref%fock_sum
                        end if
                        if (iattempt == D0_pos) then
                            call stochastic_trot_uccmc_death_nc(rng(it), ccmc_in%linked, sys, qs, iattempt==D0_pos, dfock, &
                                            qs%psip_list%dat(1,iattempt), qs%estimators(1)%proj_energy_old, &
                                            qs%psip_list%pops(1, iattempt), qs%psip_list%pops(1, iattempt) * cluster_pop, & 
                                            nparticles_change(1), ndeath_nc, logging_info)
                        else
                            call stochastic_trot_uccmc_death_nc(rng(it), ccmc_in%linked, sys, qs, iattempt==D0_pos, dfock, &
                                            qs%psip_list%dat(1,iattempt), qs%estimators(1)%proj_energy_old, &
                                            qs%psip_list%pops(1, iattempt), get_cluster_population(sys, qs%psip_list, D0_pos, &
                                            iattempt, real(D0_normalisation, p), qs%ref%f0) * qs%psip_list%pops(1, D0_pos), &
                                            nparticles_change(1), ndeath_nc, logging_info)
                        end if
                    end do
                    !$omp end do
                end if
                !$omp end parallel
                count_select = 0
                ! Add the accumulated ps_stats data to qs%excit_gen_data%p_single_double.

                if (qs%excit_gen_data%p_single_double%vary_psingles) then
                    call ps_stats_reduction_update(qs%excit_gen_data%p_single_double%rep_accum, ps_stats)
                end if

                qs%psip_list%nparticles = qs%psip_list%nparticles + nparticles_change
                qs%estimators%D0_population_comp = qs%estimators%D0_population_comp + D0_population_cycle
                qs%estimators%proj_energy_comp = qs%estimators%proj_energy_comp + proj_energy_cycle
                do j = 1, qs%psip_list%nstates
                    if (j/=D0_pos) then
                        D0_population_noncomp_cycle = &
                            D0_population_noncomp_cycle/cos((qs%psip_list%pops(1,&
                                j)/real(qs%psip_list%pop_real_factor))/D0_normalisation)
                    end if
                end do
                qs%estimators%D0_noncomposite_population = qs%estimators%D0_noncomposite_population + D0_population_noncomp_cycle

                ! Calculate the number of spawning events before the particles are redistributed,
                ! otherwise sending particles to other processors is counted as a spawning event.
                nspawn_events = calc_events_spawn_t(qs%spawn_store%spawn)
                ! Redistribute excips to new processors.
                ! The spawned excips were sent to the correct processors with
                ! the current hash shift, so it's just those in the main list
                ! that we need to deal with.
                associate(pl=>qs%psip_list, spawn=>qs%spawn_store%spawn)

                    if (nprocs > 1) call redistribute_particles(pl%states, pl%pop_real_factor, pl%pops, pl%nstates, &
                                                                pl%nparticles, spawn)
                    call direct_annihilation(sys, rng(0), qs%ref, annihilation_flags, pl, spawn)
                end associate

                if (debug) call write_logging_calc_ccmc(logging_info, iter, nspawn_events, ndeath + ndeath_nc, &
                                                        selection_data%nD0_select, &
                                                        selection_data%nclusters, selection_data%nstochastic_clusters, &
                                                       selection_data%nsingle_excitors)

                !Take updated time average of CI wavefunction.
                if(uccmc_in%variational_energy .and. all(qs%vary_shift)) then
                          time_avg_psip_list_ci_pops(:nstates_ci) =  time_avg_psip_list_ci_pops(:nstates_ci)/(iter-avg_start+1)
                end if

                ! Add new contributions to time-average cluster populations.
                if (all(qs%vary_shift) .and. old_vary .and. uccmc_in%average_wfn) then
                    call add_t_contributions(qs%psip_list, time_avg_psip_list_states, time_avg_psip_list_pops, &
                                             time_avg_psip_list_sq, nstates_sq, .true.)
                end if
                call end_mc_cycle(nspawn_events, ndeath_nc, qs%psip_list%pop_real_factor, nattempts_spawn, qs%spawn_store%rspawn)
            end do

            update_tau = bloom_stats%nblooms_curr > 0

            error = qs%spawn_store%spawn%error .or. qs%psip_list%error

            qs%estimators%D0_population = real(qs%estimators%D0_population_comp,p)
            qs%estimators%proj_energy = real(qs%estimators%proj_energy_comp,p)
            if (uccmc_in%variational_energy) then
                call var_energy_uccmc(sys, time_avg_psip_list_ci_states, time_avg_psip_list_ci_pops, nstates_ci, var_energy, &
                                      real(D0_normalisation,p))
                !qs%estimators%var_energy = var_energy
            end if 
            if (debug) call write_logging_select_ccmc(logging_info, iter, selection_data)
          
            call end_report_loop(io_unit, qmc_in, iter, update_tau, qs, nparticles_old, nspawn_events, &
                                 -1, semi_stoch_it, soft_exit=soft_exit, &
                                 load_bal_in=load_bal_in, bloom_stats=bloom_stats, comp=sys%read_in%comp, &
                                 error=error, vary_shift_reference=ccmc_in%vary_shift_reference)
            if (error) exit

            call cpu_time(t2)
            
            if (parent) then
                if (bloom_stats%nblooms_curr > 0) call bloom_stats_warning(bloom_stats, io_unit=io_unit)
                call write_qmc_report(qmc_in, qs, ireport, nparticles_old, t2-t1, .false., .false., &
                                        io_unit=io_unit, cmplx_est=sys%read_in%comp, rdm_energy=ccmc_in%density_matrices, &
                                        nattempts=.true.)
            end if


            ! Update the time for the start of the next iteration.
            t1 = t2

            call dump_restart_file_wrapper(qs, dump_restart_shift, restart_in%write_freq, nparticles_old, ireport, &
                                           qmc_in%ncycles, sys%basis%nbasis, ri, ri_shift, .false., sys%basis%info_string_len, &
                                           rng(0))

            qs%psip_list%tot_nparticles = nparticles_old

            if (soft_exit) exit

            if (update_tau) call rescale_tau(qs%tau)


        end do

        if (parent) write (io_unit,'()')

        if (parent .and. uccmc_in%average_wfn) then
            ! Take average of wavefunction.
            time_avg_psip_list_pops(:nstates_sq) =  time_avg_psip_list_pops(:nstates_sq)/(iter-avg_start+1)
            time_avg_psip_list_sq(sys%basis%tot_string_len+1,:nstates_sq) = &
                time_avg_psip_list_sq(sys%basis%tot_string_len+1,:nstates_sq)/(iter-avg_start+1)
            call write_average_wfn_trot(sys, time_avg_psip_list_pops, time_avg_psip_list_sq, io_unit, &
                    time_avg_psip_list_states, nstates_sq)
        end if

        call dSFMT_t_to_dSFMT_state_t(rng(0), qs%rng_state)

        if (parent) write (io_unit,'()')
        call write_bloom_report(bloom_stats, io_unit=io_unit)
        call multispawn_stats_report(ms_stats, io_unit=io_unit)

        call load_balancing_report(qs%psip_list%nparticles, qs%psip_list%nstates, qmc_in%use_mpi_barriers,&
                                   qs%spawn_store%spawn%mpi_time, io_unit=io_unit)
        call write_memcheck_report(qs%spawn_store%spawn, io_unit)

        if (soft_exit .or. error) then
            qs%mc_cycles_done = qs%mc_cycles_done + qmc_in%ncycles*ireport
        else
            qs%mc_cycles_done = qs%mc_cycles_done + qmc_in%ncycles*qmc_in%nreport
        end if

        if (restart_in%write_restart) then
            call dump_restart_hdf5(qs, qs%mc_cycles_done, nparticles_old, sys%basis%nbasis, .false., sys%basis%info_string_len)
            if (parent) write (io_unit,'()')
        end if

        if(uccmc_in%variational_energy) then
            call var_energy_uccmc(sys, time_avg_psip_list_ci_states, time_avg_psip_list_ci_pops, nstates_ci, var_energy, &
                                  real(D0_normalisation,p))
            print*, 'Variational energy: ', var_energy
        end if 
        
        if (debug) call end_logging(logging_info)
        if (debug) call end_selection_data(selection_data)

        call dealloc_contrib(contrib, ccmc_in%linked)
        do i = 0, nthreads-1
            call dSFMT_end(rng(i))
        end do

        ! Deallocate arrays for OpenMP threads.
        deallocate(rng, stat=ierr)
        call check_deallocate('rng', ierr)
        deallocate(ms_stats, stat=ierr)
        call check_deallocate('ms_stats', ierr)
        deallocate(ps_stats, stat=ierr)
        call check_deallocate('ps_stats', ierr)
        deallocate(state, stat=ierr)
        call check_deallocate('state', ierr)
    
        if (uccmc_in%variational_energy) then
            deallocate(time_avg_psip_list_ci_states, stat = ierr)
            call check_deallocate('time_avg_psip_list_ci_states', ierr)
            deallocate(time_avg_psip_list_ci_pops, stat = ierr)
            call check_deallocate('time_avg_psip_list_ci_pops', ierr)
        end if

        if (uccmc_in%average_wfn) then
            deallocate(time_avg_psip_list_states, stat = ierr)
            call check_deallocate('time_avg_psip_list_states', ierr)
            deallocate(time_avg_psip_list_pops, stat = ierr)
            call check_deallocate('time_avg_psip_list_pops', ierr)
            deallocate(time_avg_psip_list_sq, stat = ierr)
            call check_deallocate('time_avg_psip_list_sq', ierr)
        end if

    end subroutine do_trot_uccmc

    subroutine select_ucc_trot_cluster(rng, sys, psip_list, f0, ex_level, nattempts, normalisation, &
                              initiator_pop, D0_pos, cumulative_excip_pop, tot_excip_pop, min_size, max_size, &
                              logging_info, cdet, cluster, excit_gen_data, cluster_pop)

        ! Based on select_cluster (without the linked_cluster parts) and with information about de-excitors.
        ! Select a random cluster of excitors from the excitors on the
        ! processor.  A cluster of excitors is itself an excitor.  For clarity
        ! (if not technical accuracy) in comments we shall distinguish between
        ! the cluster of excitors and a single excitor, from a set of which the
        ! cluster is formed.

        ! In:
        !    sys: system being studied
        !    psip_list: particle_t object containing current excip distribution on
        !       this processor.
        !    f0: bit string of the reference.
        !    ex_level: max number of excitations from the reference to include in
        !        the Hilbert space.
        !    nattempts: the number of times (on this processor) a random cluster
        !        of excitors is generated in the current timestep.
        !    normalisation: intermediate normalisation factor, N_0, where we use the
        !       wavefunction ansatz |\Psi_{CC}> = N_0 e^{T/N_0} | D_0 >.
        !    initiator_pop: the population above which a determinant is an initiator.
        !    D0_pos: position in the excip list of the reference.  Must be negative
        !       if the reference is not on the processor.
        !    cumulative_excip_pop: running cumulative excip population on
        !        all excitors; i.e. cumulative_excip_population(i) = sum(particle_t%pops(1:i)).
        !    tot_excip_pop: total excip population.
        !    min_size: the minimum size cluster to allow.
        !    max_size: the maximum size cluster to allow.
        !    logging_info: derived type containing information on currently logging status
        !    excit_gen_data: information about excitation generators

        ! NOTE: cumulative_excip_pop and tot_excip_pop ignore the population on the
        ! reference as excips on the reference cannot form a cluster and the rounds the
        ! population on all other excitors to the nearest integer (for convenience--see
        ! comments in do_ccmc).  Both these quantities should be generated by
        ! cumulative_population (or be in the same format).

        ! In/Out:
        !    rng: random number generator.
        !    cdet: information about the cluster of excitors applied to the
        !        reference determinant.  This is a bare det_info_t variable on input
        !        with only the relevant fields allocated.  On output the
        !        appropriate (system-specific) fields have been filled by
        !        decoding the bit string of the determinant formed from applying
        !        the cluster to the reference determinant.
        !    cluster:
        !        Additional information about the cluster of excitors.  On
        !        input this is a bare cluster_t variable with the excitors array
        !        allocated to the maximum number of excitors in a cluster.  On
        !        output all fields in cluster have been set.

        use checking, only: check_deallocate
        use determinant_data, only: det_info_t
        use ccmc_data, only: cluster_t
        use ccmc_utils, only: convert_excitor_to_determinant
        use excitations, only: get_excitation_level
        use dSFMT_interface, only: dSFMT_t, get_rand_close_open
        use qmc_data, only: particle_t
        use proc_pointers, only: decoder_ptr
        use utils, only: factorial
        use search, only: binary_search
        use sort, only: insert_sort
        use parallel, only: nprocs
        use system, only: sys_t
        use logging, only: write_logging_stoch_selection, logging_t
        use excit_gens, only: excit_gen_data_t
        use const, only: depsilon
        use ccmc_selection, only: create_null_cluster
        use uccmc, only: ucc_collapse_cluster


        type(sys_t), intent(in) :: sys
        type(particle_t), intent(in), target :: psip_list
        integer(i0), intent(in) :: f0(sys%basis%tot_string_len)
        integer, intent(in) :: ex_level
        integer(int_64), intent(in) :: nattempts
        integer, intent(in) :: D0_pos
        complex(p), intent(in) :: normalisation
        real(p), intent(in) :: initiator_pop
        real(p), intent(in) :: cumulative_excip_pop(:), tot_excip_pop, cluster_pop
        integer :: min_size, max_size
        type(dSFMT_t), intent(inout) :: rng
        type(det_info_t), intent(inout) :: cdet
        type(cluster_t), intent(inout) :: cluster
        type(logging_t), intent(in) :: logging_info
        type(excit_gen_data_t), intent(in) :: excit_gen_data

        real(dp) :: rand
        real(p) :: psize, tot_excip_local
        complex(p) :: cluster_population, excitor_pop
        integer :: i, pos, prev_pos, ierr, current_excit
        real(p) :: pop(max_size), ref_real, pop_real
        integer :: poses(max_size)
        logical :: hit, allowed
        logical :: conjugate
       
        ! We shall accumulate the factors which comprise cluster%pselect as we go.
        !   cluster%pselect = n_sel p_size p_clust
        ! where
        !   n_sel   is the number of cluster selections made;
        !   p_size  is the probability of choosing a cluster of that size;
        !   p_clust is the probability of choosing a specific cluster given
        !           the choice of size.
        ! Each processor does nattempts.
        ! However:
        ! * if min_size=0, then each processor is allowed to select the reference (on
        !   average) nattempts/2 times.  Hence in order to have selection probabilities
        !   consistent and independent of the number of processors being used (which
        !   amounts to a processor-dependent timestep scaling), we need to multiply the
        !   probability the reference is selected by nprocs.
        ! * assuming each excitor spends (on average) the same amount of time on each
        !   processor, the probability that X different excitors are on the same processor
        !   at a given timestep is 1/nprocs^{X-1).
        ! The easiest way to handle both of these is to multiply the number of attempts by
        ! the number of processors here and then deal with additional factors of 1/nprocs
        ! when creating composite clusters.
        ! NB within a processor those nattempts can be split amongst OpenMP
        ! threads though that doesn't affect this probability.
        cluster%pselect = real(nattempts*nprocs, p)
        ! Select the cluster size, i.e. the number of excitors in a cluster.
        ! For a given truncation level, only clusters containing at most
        ! ex_level+2 excitors.
        ! Following the process described by Thom in 'Initiator Stochastic
        ! Coupled Cluster Theory' (unpublished), each size, n_s, has probability
        ! p(n_s) = 1/2^(n_s+1), n_s=0,ex_level and p(ex_level+2)
        ! is such that \sum_{n_s=0}^{ex_level+2} p(n_s) = 1.

        ! This procedure is modified so that clusters of size min_size+n_s
        ! has probability 1/2^(n_s+1), and the max_size picks up the remaining
        ! probability from the series.
        rand = get_rand_close_open(rng)
        psize = 0.0_p
        cluster%nexcitors = -1
        cluster_population = 1
        cdet%f = f0
        do i = 0, max_size-min_size-1
            psize = psize + 1.0_p/2_int_64**(i+1)
            if (rand < psize) then
                ! Found size!
                cluster%nexcitors = i+min_size
                cluster%pselect = cluster%pselect/2_int_64**(i+1)
                exit
            end if
        end do
        ! If not set, then must be the largest possible cluster
        if (cluster%nexcitors == -1) then
            cluster%nexcitors = max_size
            cluster%pselect = cluster%pselect*(1.0_p - psize)
        end if

        ! If could be using logging set to easily identifiable nonsense value.
        if (debug) pop = -1_int_p

        ! Initiator approximation.
        ! This is sufficiently quick that we'll just do it in all cases, even
        ! when not using the initiator approximation.  This matches the approach
        ! used by Alex Thom in 'Initiator Stochastic Coupled Cluster Theory'
        ! (unpublished).
        ! Assume all excitors in the cluster are initiators (initiator_flag=0)
        ! until proven otherwise (initiator_flag=1).
        cdet%initiator_flag = 0

        ! Assume cluster is allowed unless collapse_cluster finds out otherwise
        ! when collapsing/combining excitors or if it could never have been
        ! valid
        allowed = min_size <= max_size
        ! For linked coupled cluster we keep building the cluster after a
        ! disallowed excitation so need to know if there has been a disallowed
        ! excitation at all
        ref_real = real(psip_list%pops(1, D0_pos))/real(psip_list%pop_real_factor)

        select case(cluster%nexcitors)
        case(0)
            call create_null_cluster(sys, f0, cluster%pselect, normalisation*cluster_pop, initiator_pop, &
                                    cdet, cluster, excit_gen_data)
        case default
            ! Select cluster from the excitors on the current processor with
            ! probability for choosing an excitor proportional to the excip
            ! population on that excitor.  (For convenience, we use a probability
            ! proportional to the ceiling(pop), as it makes finding the right excitor
            ! much easier, especially for the non-composite algorithm, as well as
            ! selecting excitors with the correct (relative) probability.  The
            ! additional fractional weight is taken into account in the amplitude.)
            !
            ! Rather than selecting one excitor at a time and adding it to the
            ! cluster, select all excitors and then find their locations and
            ! apply them.  This allows us to sort by population first (as the
            ! number of excitors is small) and hence allows for a more efficient
            ! searching of the cumulative population list.

            do i = 1, cluster%nexcitors
                ! Select nexcitors different positions in the excitors list.
                if (i == 1) then 
                    pop(i) = get_rand_close_open(rng)*tot_excip_pop
                    call binary_search(cumulative_excip_pop, pop(i), 1, psip_list%nstates, hit, poses(i))
                    do
                        if (poses(i) == 1) then
                            exit
                        end if
                        if (abs(cumulative_excip_pop(poses(i)) - cumulative_excip_pop(poses(i)-1)) > depsilon) exit
                        poses(i) = poses(i) - 1
                    end do
                else
                     pop(i) = get_rand_close_open(rng)*tot_excip_pop
                     call binary_search(cumulative_excip_pop, pop(i), 1, psip_list%nstates, hit, poses(i))
                     do
                         if (poses(i) == 1) then
                             exit
                         end if
                         if (abs(cumulative_excip_pop(poses(i)) - cumulative_excip_pop(poses(i)-1)) > depsilon) exit
                         poses(i) = poses(i) - 1
                     end do
                     if (any(poses(1:i-1)==poses(i))) then
                        ! Clusters with the same excitor applied multiple times are not allowed.
                        allowed = .false.
                        exit
                     end if
                end if
            end do
            tot_excip_local = tot_excip_pop
            call insert_sort(pop(:cluster%nexcitors))
            prev_pos = 1
            current_excit = 1
            if (allowed) then
                ! For each excitor in psip_list check if it can be applied as an excitation/deexcitation operator and if it is
                ! in one of the positions selected before.
                do i = 1, psip_list%nstates
                   if (i /= D0_pos) then
                   pop_real = real(psip_list%pops(1, i))/real(psip_list%pop_real_factor)
                   conjugate = .false.
                   if (deexcitation_possible(f0(:sys%basis%bit_string_len),psip_list%states(:sys%basis%bit_string_len,i), &
                       cdet%f(:sys%basis%bit_string_len))) then
                           conjugate = .true.
                           if (current_excit <= cluster%nexcitors) then
                              if (pop(current_excit) <= cumulative_excip_pop(i) .and. &
                                  pop(current_excit) > cumulative_excip_pop(i-1)) then
                                ! If the excitor can be applied as a de-excitation operator and is in the selected list, 
                                ! apply it to current cluster and multiply population by -sin(pop_real/ref_real) = -sin(t)
                                excitor_pop = -sin(pop_real/ref_real)
                                cluster%pselect = cluster%pselect*abs(pop_real)/tot_excip_local
                                call ucc_collapse_cluster(sys%basis, f0, psip_list%states(:,i), excitor_pop, cdet%f, &
                                              cluster_population, allowed, conjugate)
                                cluster%excitors(current_excit)%f => psip_list%states(:,i)
                                if (abs(excitor_pop) <= initiator_pop) cdet%initiator_flag = 3
                                current_excit = current_excit + 1
                              else
                                ! If the excitor can be applied as a de-excitation operator and is NOT in the selected list, 
                                ! multiply population by cos(pop_real/ref_real)
                                excitor_pop = cos(pop_real/ref_real)
                                cluster_population = cluster_population*excitor_pop
                              end if
                           else
                              excitor_pop = cos(pop_real/ref_real)
                              cluster_population = cluster_population*excitor_pop
                           end if
                        else if (excitation_possible(f0(:sys%basis%bit_string_len),psip_list%states(:sys%basis%bit_string_len,i), &
                                 cdet%f(:sys%basis%bit_string_len))) then
                           if (current_excit <= cluster%nexcitors) then
                              if (pop(current_excit) <= cumulative_excip_pop(i) .and. &
                                  pop(current_excit) > cumulative_excip_pop(i-1)) then
                                ! If the excitor can be applied as an excitation operator and is in the selected list, 
                                ! apply it to current cluster and multiply population by sin(pop_real/ref_real) = sin(t)
                                if (current_excit == 1) cdet%data => psip_list%dat(:,i)
                                excitor_pop = sin(pop_real/ref_real)
                                cluster%pselect = cluster%pselect*abs(pop_real)/tot_excip_local
                                call ucc_collapse_cluster(sys%basis, f0, psip_list%states(:,i), excitor_pop, cdet%f, &
                                              cluster_population, allowed, conjugate)
                                cluster%excitors(current_excit)%f => psip_list%states(:,i)
                                if (abs(excitor_pop) <= initiator_pop) cdet%initiator_flag = 3
                                current_excit = current_excit + 1
                              else
                                ! If the excitor can be applied as an excitation operator and is NOT in the selected list, multiply
                                ! population by cos(pop_real/ref_real)
                                excitor_pop = cos(pop_real/ref_real)
                                cluster_population = cluster_population*excitor_pop
                              end if
                           else
                              excitor_pop = cos(pop_real/ref_real)
                              cluster_population = cluster_population*excitor_pop
                           end if
                        else 
                           ! If excitor CANNOT be applied and is in the list, cluster is not allowed so exit. 
                           ! Otherwise, just move onto next excitor.
                           if (current_excit <= cluster%nexcitors) then
                               if (pop(current_excit) <= cumulative_excip_pop(i)) then
                                allowed = .false.
                                exit
                               end if
                           end if
                        end if
                    end if
                end do
            end if

            if (allowed) then
                cluster%excitation_level = get_excitation_level(f0, cdet%f)
                ! To contribute the cluster must be within a double excitation of
                ! the maximum excitation included in the CC wavefunction.
                allowed = cluster%excitation_level <= ex_level+2
            end if

            if (allowed) then
                cluster%pselect = cluster%pselect*factorial(cluster%nexcitors)
                ! Sign change due to difference between determinant
                ! representation and excitors and excitation level.
                call convert_excitor_to_determinant(cdet%f, cluster%excitation_level, cluster%cluster_to_det_sign, f0)
                call decoder_ptr(sys, cdet%f, cdet, excit_gen_data)

                ! Normalisation factor for cluster%amplitudes...
                cluster%amplitude = cluster_population*normalisation
            else
                ! Simply set excitation level to a too high (fake) level to avoid
                ! this cluster being used.
                cluster%excitation_level = huge(0)
            end if
        end select

        if (debug) call write_logging_stoch_selection(logging_info, cluster%nexcitors, cluster%excitation_level, pop, &
                max_size, cluster%pselect, cluster%amplitude, allowed)

    end subroutine select_ucc_trot_cluster

    subroutine regenerate_trot_info_psip_list(basis, nel, qs)

        ! To obtain an ordering consisten with https://doi.org/10.1063/1.5133059
        ! the psip_list is ordered in descending order, 
        ! first by highest orbital excited from(relative to f0), then by (nelec - excitation level).
        ! e.g. D0, t_2^4, t_12^34, t_1^4 
        ! Regenerates excitation level and lowest unoccupied orbital information
        ! stored at start of bit string within states in psip list.
        ! For use when restarting from a restart file
        ! not containing this information.
        ! Also sorts the list, as ordering will change.
        ! Hashing only uses nbasis bits, so should be unaffected by the additional
        ! information at the start of the bit string.
        ! [todo] figure out a way to double check this is the case.

        ! In:
        !   basis: information on single-particle basis in use.
        !   nel: number of electrons in the system
        ! In/Out:
        !   qs: qmc_state_t object with information on current state of calculation. We update and
        !       reorder the bit strings within the psip list, using the reference
        !       determinant bit string stored within qs%ref%f0.

        use basis_types, only: basis_t
        use qmc_data, only: qmc_state_t
        use sort, only: qsort
        use uccmc_utils, only: add_info_str_trot

        type(basis_t), intent(in) :: basis
        type(qmc_state_t), intent(inout) :: qs
        integer, intent(in) :: nel

        integer :: istate

        do istate = 1, qs%psip_list%nstates
            call add_info_str_trot(basis, qs%ref%f0, nel, qs%psip_list%states(:,istate))
        end do

        associate(pl=>qs%psip_list)
             call qsort(pl%nstates, pl%states, pl%pops, pl%dat, .true.)
        end associate

    end subroutine regenerate_trot_info_psip_list

    subroutine find_D0_trot(psip_list, f0, D0_pos)

        ! Find the reference determinant in the list of walkers labelled for trotterized UCC.

        ! In:
        !    psip_list: particle_t object containing current excip distribution on
        !       this processor.
        !    f0: bit string representing the reference.
        ! In/Out:
        !    D0_pos: on input, the position of the reference in
        !       particle_t%states in the previous iteration (or -1 if it was
        !       not on this processor).  On output, the current position.

        use bit_utils, only: bit_str_cmp
        use search, only: binary_search
        use qmc_data, only: particle_t
        use errors, only: stop_all

        type(particle_t), intent(in) :: psip_list
        integer(i0), intent(in) :: f0(:)
        integer, intent(inout) :: D0_pos

        integer :: D0_pos_old
        logical :: hit

        if (D0_pos == -1) then
            ! D0 was just moved to this processor.  No idea where it might be...
            call binary_search(psip_list%states, f0, 1, psip_list%nstates, hit, D0_pos, .true.)
        else
            D0_pos_old = D0_pos
            select case(-bit_str_cmp(f0, psip_list%states(:,D0_pos)))
            case(0)
                ! D0 hasn't moved.
                hit = .true.
            case(1)
                ! D0 < psip_list%states(:,D0_pos) -- it has moved to earlier in
                ! the list and the old D0_pos is an upper bound.
                call binary_search(psip_list%states, f0, 1, D0_pos_old, hit, D0_pos, .true.)
            case(-1)
                ! D0 > psip_list%states(:,D0_pos) -- it has moved to later in
                ! the list and the old D0_pos is a lower bound.
                call binary_search(psip_list%states, f0, D0_pos_old, psip_list%nstates, hit, D0_pos, .true.)
            end select
        end if
        if (.not.hit) call stop_all('find_D0', 'Cannot find reference!')

    end subroutine find_D0_trot

    subroutine get_D0_info_trot(qs, complx, D0_proc, D0_pos, nD0_proc, D0_normalisation)

        ! In:
        !    qs: qmc_state_t object describing the current CCMC state.
        !    complx: true if system has a complex wavefunction (i.e. sys_t%sys_read_in_t%comp).
        ! Out:
        !    D0_proc: the processor index on which the reference currently resides.
        !    D0_pos: the position within the excip list of the reference. Set to -1 if iproc != D0_proc.
        !    nD0_proc: 1 if iproc == D0_proc and 0 otherwise.
        !    D0_normalisation: population of the reference.

        use parallel
        use qmc_data, only: qmc_state_t
        use spawning, only: assign_particle_processor

        type(qmc_state_t), intent(in) :: qs
        logical, intent(in) :: complx
        integer, intent(out) :: D0_proc, D0_pos, nD0_proc
        complex(p), intent(out) :: D0_normalisation
        integer :: slot
#ifdef PARALLEL
        integer :: ierr
#endif

        associate(spawn=>qs%spawn_store%spawn, pm=>qs%spawn_store%spawn%proc_map)
            call assign_particle_processor(qs%ref%f0, spawn%bit_str_nbits, spawn%hash_seed, spawn%hash_shift, &
                                           spawn%move_freq, nprocs, D0_proc, slot, pm%map, pm%nslots)

        end associate

        if (iproc == D0_proc) then

            ! Population on reference determinant.
            ! As we might select the reference determinant multiple times in
            ! a cycle, the running total of D0_population is incorrect (by
            ! a factor of the number of times it was selected).
            call find_D0_trot(qs%psip_list, qs%ref%f0, D0_pos)
            if (complx) then
                D0_normalisation = cmplx(qs%psip_list%pops(1,D0_pos), qs%psip_list%pops(2,D0_pos), p)/qs%psip_list%pop_real_factor
            else
                D0_normalisation = real(qs%psip_list%pops(1,D0_pos),p)/qs%psip_list%pop_real_factor
            end if
            nD0_proc = 1
        else

            ! Can't find D0 on this processor.  (See how D0_pos is used
            ! in select_cluster.)
            D0_pos = -1
            nD0_proc = 0 ! No reference excitor on the processor.

        end if

#ifdef PARALLEL
        call mpi_bcast(D0_normalisation, 1, mpi_pcomplex, D0_proc, MPI_COMM_WORLD, ierr)
#endif

    end subroutine get_D0_info_trot

    pure function deexcitation_possible(f0, excit, cdet_f) result (allowed)
        ! Function to check whether it is possible to apply a particular
        ! deexcitation operator to the current cluster.

        ! In:
        ! f0: bit string corresponding to the reference determinant.
        ! excit: bit string corresponding to the effect of applying a given excitor
        ! to the reference.
        ! cdet_f: bit string corresponding to the current cluter.

        integer(i0), intent(in) :: f0(:), excit(:), cdet_f(:)
        logical :: allowed

        allowed = all(iand(ieor(f0(:),excit(:)), &
                           ieor(f0(:),cdet_f(:))) == &
                           ieor(f0(:),excit(:)))
    end function deexcitation_possible

    pure function excitation_possible(f0, excit, cdet_f) result (allowed)
        ! Function to check whether it is possible to apply a particular
        ! excitation operator to the current cluster.

        ! In:
        ! f0: bit string corresponding to the reference determinant.
        ! excit: bit string corresponding to the effect of applying a given excitor
        ! to the reference.
        ! cdet_f: bit string corresponding to the current cluter.

        integer(i0), intent(in) :: f0(:), excit(:), cdet_f(:)
        logical :: allowed

        allowed = all(iand(ieor(f0(:),excit(:)), &
                           ieor(f0(:),cdet_f(:))) == 0)        
    end function excitation_possible

    subroutine select_nc_cluster_trot(sys, psip_list, f0, iexcitor, initiator_pop, ex_lvl_sort, &
                                            cdet, cluster, excit_gen_data, D0_normalisation, D0_pos)

        ! Select (deterministically) the non-composite cluster containing only
        ! the single excitor iexcitor and set the same information as select_ucc_trot_cluster.

        ! In:
        !    sys: system being studied
        !    psip_list: particle_t object containing current excip distribution on
        !       this processor.
        !    f0: bit string of the reference
        !    iexcitor: the index (in range [1,nstates]) of the excitor to select.
        !    initiator_pop: the population above which a determinant is an initiator.
        !    ex_lvl_sort: if true, excitors are sorted by excitation level and so have
        !       information about the excitation level also encoded in their bit string
        !       that should be removed before being passed on.
        !    excit_gen_data: data for excitation generators.

        ! In/Out:
        !    cdet: information about the cluster of excitors applied to the
        !        reference determinant.  This is a bare det_info variable on input
        !        with only the relevant fields allocated.  On output the
        !        appropriate (system-specific) fields have been filled by
        !        decoding the bit string of the determinant formed from applying
        !        the cluster to the reference determinant.
        !    cluster:
        !        Additional information about the cluster of excitors.  On
        !        input this is a bare cluster_t variable with the excitors array
        !        allocated to the maximum number of excitors in a cluster.  On
        !        output all fields in cluster have been set.

        use system, only: sys_t
        use determinant_data, only: det_info_t
        use basis_types, only: reset_extra_info_bit_string
        use ccmc_data, only: cluster_t
        use ccmc_utils, only: convert_excitor_to_determinant
        use excitations, only: get_excitation_level
        use qmc_data, only: particle_t
        use proc_pointers, only: decoder_ptr
        use excit_gens, only: excit_gen_data_t

        type(sys_t), intent(in) :: sys
        type(particle_t), intent(in), target :: psip_list
        integer(i0), intent(in) :: f0(sys%basis%tot_string_len)
        integer(int_64), intent(in) :: iexcitor
        real(p), intent(in) :: initiator_pop
        logical, intent(in) :: ex_lvl_sort
        complex(p), intent(in) :: D0_normalisation
        type(excit_gen_data_t), intent(in) :: excit_gen_data
        integer, intent(in) :: D0_pos
        type(det_info_t), intent(inout) :: cdet
        type(cluster_t), intent(inout) :: cluster
        complex(p) :: excitor_pop

        ! Rather than looping over individual excips we loop over different sites. This is
        ! because we want to stochastically set the number of spawning attempts such that
        ! we on average select each excitor a number of times proportional to the absolute
        ! population. As such, we can't specify the total number of selections beforehand
        ! within do_ccmc.

        ! As iterating deterministically through all noncomposite clusters, pselect = 1
        ! exactly.
        ! This excitor can only be selected on this processor and only one excitor is
        ! selected in the cluster, so unlike selecting the reference or composite
        ! clusters, there are no additional factors of nprocs or 1/nprocs to include.

        cluster%pselect = 1.0_p

        cluster%nexcitors = 1

        ! Initiator approximation.
        ! This is sufficiently quick that we'll just do it in all cases, even
        ! when not using the initiator approximation.  This matches the approach
        ! used by Alex Thom in 'Initiator Stochastic Coupled Cluster Theory'
        ! (unpublished).
        ! Assume all excitors in the cluster are initiators (initiator_flag=0)
        ! until proven otherwise (initiator_flag=1).
        cdet%initiator_flag = 0

        cdet%f = psip_list%states(:,iexcitor)
        cdet%data => psip_list%dat(:,iexcitor)
        cluster%excitors(1)%f => psip_list%states(:,iexcitor)
        if (sys%read_in%comp) then
            excitor_pop = cmplx(psip_list%pops(1,iexcitor), psip_list%pops(2,iexcitor),p)/psip_list%pop_real_factor
        else
            excitor_pop = real(psip_list%pops(1,iexcitor),p)/psip_list%pop_real_factor
        end if

        if (abs(excitor_pop) <= initiator_pop) cdet%initiator_flag = 3

        if (ex_lvl_sort) call reset_extra_info_bit_string(sys%basis, cdet%f)

        cluster%excitation_level = get_excitation_level(f0(:sys%basis%bit_string_len), cdet%f(:sys%basis%bit_string_len))
        cluster%amplitude = get_cluster_population(sys, psip_list, D0_pos, iexcitor, real(D0_normalisation, p), f0)*D0_normalisation

        ! Sign change due to difference between determinant
        ! representation and excitors and excitation level.
        call convert_excitor_to_determinant(cdet%f, cluster%excitation_level, cluster%cluster_to_det_sign, f0)
        call decoder_ptr(sys, cdet%f, cdet, excit_gen_data)

    end subroutine select_nc_cluster_trot

! [review] - Brian: document this
    function get_cluster_population(sys, psip_list, D0_pos, iattempt, ref_real, f0) result(cluster_population)
        ! Function to obtain the effective cluster population of a non-composite cluster in tUCCMC.
        ! The excitor the cluster corresponds to contributes sin(Ni/N0). Every other excitor that
        ! could be applied (but is not) corresponds cos(Ni/N0).

        ! In:
        ! sys: sys_t object for the system studied.
        ! psip_list: particle_t object encoding the current wavefunction.
        ! D0_pos: position of D0 in psip_list.
        ! iattempt: current excitor considered.
        ! ref_real: real population on the reference.
        ! f0: bit string of the reference determinant.

        use qmc_data, only: particle_t
        use system, only: sys_t

        type(particle_t), intent(in) :: psip_list
        type(sys_t), intent(in) :: sys
        integer, intent(in) :: D0_pos
        ! [review] - Brian: should iattempt be int_64? because it seems that all calls are filled by int_64 integers (iexcitor etc)
        integer(int_64), intent(in) :: iattempt
        integer(i0), intent(in) :: f0(:)
        real(p), intent(in) :: ref_real

        integer(i0) :: f(size(f0))
        real(p) :: cluster_population, pop_real, excitor_pop
        integer(int_64) :: i
        cluster_population = 1.0_p

        f = f0
        
        do i = 1, psip_list%nstates
            if (i /= D0_pos) then
                pop_real = real(psip_list%pops(1, i))/real(psip_list%pop_real_factor)
                if (deexcitation_possible(f0(:sys%basis%bit_string_len),psip_list%states(:sys%basis%bit_string_len,i), &
                                          f(:sys%basis%bit_string_len))) then
                    excitor_pop = cos(pop_real/ref_real)
                    cluster_population = cluster_population*excitor_pop
                else if (excitation_possible(f0(:sys%basis%bit_string_len),psip_list%states(:sys%basis%bit_string_len,i), &
                                             f(:sys%basis%bit_string_len))) then
                    if (i == iattempt) then
                        excitor_pop = sin(pop_real/ref_real)
                        f = psip_list%states(:, iattempt)
                    else
                        excitor_pop = cos(pop_real/ref_real)
                    end if
                    cluster_population = cluster_population*excitor_pop
                end if
            end if
        end do
    end function get_cluster_population

    subroutine stochastic_trot_uccmc_death_nc(rng, linked_ccmc,  sys, qs, isD0, dfock, Hii, proj_energy, population, &
                                        trot_population, tot_population, ndeath, logging_info)

        ! Based on stochastic_ccmc_death_nc, accounting for different cluster population definition in tUCCMC.
        ! Attempt to 'die' (ie create an excip on the current excitor, cdet%f)
        ! with probability
        !    \tau |<D_s|H|D_s> A_s|
        !    ----------------------
        !       n_sel p_s p_clust
        ! where |D_s> is the determinant formed by applying the excitor to the
        ! reference determinant and A_s is the amplitude.  See comments in
        ! select_cluster about the probabilities.

        ! When doing linked CCMC, the matrix elements
        !   <D_s|H|D_s> = <D_s|H T^n|D_0>
        ! are replaced by
        !   <D_s|[..[H,T],..T]|D_0>
        ! which changes the death probabilities, and also means the shift only
        ! applies on the reference determinant.

        ! This procedure kills the excips directly, rather than by creating anti-particles
        ! in the spawned list, so only works for non-composite excips (single excitors or
        ! particles on the reference).

        ! In:
        !    linked_ccmc: if true then only sample linked clusters.
        !    qs: qmc_state_t object. The shift and timestep are used.
        !    isD0: true if the current excip is the null (reference) excitor
        !    dfock: difference in the Fock energy of the determinant formed by applying the excitor to
        !       the reference and the Fock energy of the reference.
        !    Hii: the diagonal matrix element of the determinant formed by applying the excitor to the
        !       reference.
        !    proj_energy: projected energy.  This should be the average value from the last
        !        report loop, not the running total in qs%estimators.
        !    trot_population: the effective tUCCMC population on the current excip, as computed from get_cluster_population.
        ! In/Out:
        !    rng: random number generator.
        !    ndeath: running (encoded) total of number of particles killed/cloned.
        !    population: the (encoded) population on the current excip
        !    tot_population: total number of particles.

        use dSFMT_interface, only: dSFMT_t, get_rand_close_open
        use qmc_data, only: qmc_state_t
        use system, only: sys_t
        use spawning, only: calc_qn_weighting
        use logging, only: logging_t, write_logging_death

        type(sys_t), intent(in) :: sys
        logical, intent(in) :: linked_ccmc
        type(qmc_state_t), intent(in) :: qs
        logical, intent(in) :: isD0
        real(p), intent(in) :: Hii, dfock, proj_energy
        type(logging_t), intent(in) :: logging_info
        type(dSFMT_t), intent(inout) :: rng
        integer(int_p), intent(inout) :: population, ndeath
        real(p), intent(in) :: trot_population
        real(dp), intent(inout) :: tot_population

        real(p) :: pdeath, KiiAi
        integer(int_p) :: nkill, old_pop
        real(p) :: invdiagel

        ! Spawning onto the same excitor so no change in sign due to
        ! a difference in the sign of the determinant formed from applying the
        ! parent excitor to the reference and that formed from applying the
        ! child excitor.

        invdiagel = calc_qn_weighting(qs%propagator, dfock)
        if (isD0) then
            KiiAi = ((- proj_energy)*invdiagel + (proj_energy - qs%shift(1)))*trot_population
        else
            KiiAi = ((Hii - proj_energy)*invdiagel + (proj_energy - qs%shift(1)))*trot_population
        end if

        ! Death is attempted exactly once on this cluster regardless of pselect.
        ! Population passed in is in the *encoded* form.
        pdeath = qs%tau*abs(KiiAi)

        ! Number that will definitely die
        nkill = int(pdeath,int_p)
        ! Stochastic death...
        pdeath = pdeath - nkill
        if (pdeath > get_rand_close_open(rng)) nkill = nkill + 1

        if (nkill /= 0) then
            ! Create nkill excips with sign of -K_ii A_i
            if (KiiAi > 0) nkill = -nkill
            ! Kill directly for single excips
            ! This only works in the full non composite algorithm as otherwise the
            ! population on an excip can still be needed if it as selected as (part of)
            ! another cluster. It is also necessary that death is not done until after
            ! all spawning attempts from the excip
            old_pop = population
            population = population + nkill
            ! Also need to update total population
            tot_population = tot_population + real(abs(population)-abs(old_pop),p)/qs%psip_list%pop_real_factor
            ndeath = ndeath + abs(nkill)
        end if

        if (debug) call write_logging_death(logging_info, KiiAi, proj_energy, qs%shift(1), invdiagel, &
                                            nkill, pdeath, real(old_pop,p)/qs%psip_list%pop_real_factor, &
                                            real(population,p)/qs%psip_list%pop_real_factor)

    end subroutine stochastic_trot_uccmc_death_nc

    subroutine add_trot_info_reference(f0, sys)

        use system, only: sys_t
        use uccmc_utils, only: latest_unset

        integer(i0), intent(inout) :: f0(:)
        type(sys_t), intent(in) :: sys
 
        f0(sys%basis%bit_string_len + 2) = latest_unset(f0(:sys%basis%bit_string_len), &
                                                               f0(:sys%basis%bit_string_len), sys%nel, sys%basis) 
        f0(sys%basis%bit_string_len + 1) = sys%nel
    end subroutine add_trot_info_reference
end module

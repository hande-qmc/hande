module uccmc

! Module for performing unitary coupled cluster Monte Carlo (UCCMC) calculations.

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
! This means the ansatz contains both excitation and deexcitation operators. In consequence,
! the size of composite clusters is no longer limited to ex_level + 2 and can in principle be
! infinite. In practice, we truncate at a finite polynomial power in T and find that this
! converges satisfactorily.
!
! Additional care must be taken with the sign of different contributing clusters. See comments
! in ucc_collapse_cluster.
!
! Finally, selection must take into account the different possible orderings of excitation and 
! deexcitation operators in a given composite cluster.

use const, only: i0, int_p, int_64, p, dp, debug

implicit none

contains

    subroutine do_uccmc(sys, qmc_in, ccmc_in, uccmc_in, restart_in, load_bal_in, reference_in, &
                        logging_in, io_unit, qs, qmc_state_restart)

        ! This subroutine is derived from do_ccmc in ccmc.f90. [todo] check for possible shared functions

        ! Run the UCCMC algorithm starting from the initial walker distribution
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
        !    qmc_state_restart (optional): if present, restart from a previous fciqmc calculation.
        !       Deallocated on exit.
        ! Out:
        !    qs: qmc_state for use if restarting the calculation

        use checking, only: check_allocate, check_deallocate
        use dSFMT_interface, only: dSFMT_t, dSFMT_init, dSFMT_end, dSFMT_state_t_to_dSFMT_t, dSFMT_t_to_dSFMT_state_t, &
                                   free_dSFMT_state_t
        use errors, only: stop_all, warning
        use parallel
        use restart_hdf5, only: dump_restart_hdf5, restart_info_t, init_restart_info_t, dump_restart_file_wrapper

        use annihilation, only: direct_annihilation, insert_new_walker
        use bloom_handler, only: init_bloom_stats_t, bloom_stats_t, bloom_mode_fractionn, bloom_mode_fixedn, &
                                 write_bloom_report, bloom_stats_warning, update_bloom_threshold_prop
        use ccmc, only: do_ccmc_accumulation, perform_ccmc_spawning_attempt, do_stochastic_ccmc_propagation, &
                        do_nc_ccmc_propagation
        use ccmc_data
        use ccmc_death_spawning, only: stochastic_ccmc_death_nc
        use ccmc_selection, only: create_null_cluster, select_nc_cluster
        use ccmc_selection, only: init_selection_data, update_selection_probabilities, &
                                  init_amp_psel_accumulation, set_cluster_selections
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
        use spawn_data, only: calc_events_spawn_t, write_memcheck_report
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
        use uccmc_utils, only: allocate_time_average_lists, add_ci_contribution, add_t_contributions, var_energy_uccmc, &
                               initialise_average_wfn, write_average_wfn

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
 
        real(p), allocatable :: time_avg_psip_list_ci_pops(:), time_avg_psip_list_pops(:),time_avg_psip_list_sq(:,:)
        integer(i0), allocatable :: time_avg_psip_list_ci_states(:,:), time_avg_psip_list_states(:,:)
        integer :: semi_stoch_it, pos, j, k, nstates_ci, nstates_sq
        logical :: hit
        real(p) :: population
        real(p) :: real_population, var_energy
        ! old_vary encodes whether the shift started varying before the current iteration. Used to determine when
        ! to start taking averages of the wavefunction.
        logical :: old_vary
        integer :: avg_start
        integer :: count_discard

        count_discard = 0
        old_vary=.false.
        if (parent) then
            write (io_unit,'(1X,"UCCMC")')
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


        if(uccmc_in%variational_energy) then
        ! If calculating a variational estimator (very slow, only for benchmarking purposes),
        ! store average CI wavefunction.
             population = 0
             call allocate_time_average_lists(qs%psip_list, time_avg_psip_list_ci_states, time_avg_psip_list_ci_pops, nstates_ci)
        end if

        ! If computing the average wavefunction, allocate arrays to store it.
        if (uccmc_in%average_wfn) then
            call allocate_time_average_lists(qs%psip_list, time_avg_psip_list_states, time_avg_psip_list_pops, nstates_sq)
            allocate(time_avg_psip_list_sq(sys%basis%tot_string_len+1,size(qs%psip_list%states(1,:))))
            time_avg_psip_list_sq(:sys%basis%tot_string_len,:qs%psip_list%nstates) = qs%psip_list%states(:,:qs%psip_list%nstates)
            time_avg_psip_list_sq(sys%basis%tot_string_len+1,:qs%psip_list%nstates) &
                = (real(qs%psip_list%pops(1,:qs%psip_list%nstates))/qs%psip_list%pop_real_factor)**2
        end if

        if (ccmc_in%multiref) then
            ! Initialise multireference CCMC specific data.
            qs%multiref = .true.
            qs%mr_acceptance_search = ccmc_in%mr_acceptance_search
            call init_secondary_references(sys, ccmc_in, io_unit, qs)
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
            if (allocated(err_msg)) call stop_all('do_uccmc', 'Failed to reset RNG state: '//err_msg)
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
                    if (parent) call warning('do_uccmc', 'move_freq is different from that in the restart file. &
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
            call initial_cc_projected_energy(sys, qs, qmc_in%seed+iproc, logging_info, cumulative_abs_real_pops, nparticles_old, &
                                             ccmc_in)
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

                if(uccmc_in%variational_energy .and. all(qs%vary_shift)) then
                    ! If computing variational energy, store average CI expansion.
                          time_avg_psip_list_ci_pops(:nstates_ci) =  time_avg_psip_list_ci_pops(:nstates_ci)*(iter-avg_start)
                end if

                if (debug) call prep_logging_mc_cycle(iter, logging_in, logging_info, sys%read_in%comp, &
                                                        min(sys%nel, qs%ref%ex_level+2))

                call get_D0_info(qs, sys%read_in%comp, D0_proc, D0_pos, nD0_proc, D0_normalisation)

                ! Update the shift of the excitor locations to be the end of this
                ! current iteration.
                qs%spawn_store%spawn%hash_shift = qs%spawn_store%spawn%hash_shift + 1

                ! Maximum possible cluster size that we can generate.
                ! If only the reference is populated, we can only generate clusters of size 0.
                ! Otherwise we can generate clusters up to our chosen truncation level.

                if(qs%psip_list%nstates-nD0_proc == 0) then
                    max_cluster_size = 0
                else
                    max_cluster_size = uccmc_in%pow_trunc
                end if

                call init_mc_cycle(qs%psip_list, qs%spawn_store%spawn, qs%estimators(1)%nattempts, ndeath, &
                                   min_attempts=nint(abs(D0_normalisation), kind=int_64), &
                                   complx=sys%read_in%comp)

                nparticles_change = 0.0_p

                ! We need to count spawning attempts differently as there may be multiple spawns
                ! per cluster

                nattempts_spawn=0

                ! Find cumulative population...
                ! NOTE: for simplicity we only consider the integer part of the population on each excitor.
                ! (Populations under 1 are stochastically rounded in the annihilation process, so each excitor in the list has
                ! a non-zero integer population.)
                ! Unlike in FCIQMC, where we loop over each determinant and hence can individually decide whether or not to
                ! stochastically attempt another attempt for a fractional population, in CCMC we select excitors based upon their
                ! population and the total number of attempts based upon the total population.  In the non-composite algorith, we
                ! also need to find the determinant of a given excip.  This is painful to do if we use fractional populations
                ! (as we'd need to keep track of how many fractional populations had been rounded up in order to search the
                ! cumulative list correctly).  Instead, we base the number of attempts and the probably of selecting a given excitor
                ! solely upon the nearest integer of the population.  This decouples (slightly) the selection probability and the
                ! amplitude, which uses the exact population (including fractional part) but is fine as we can choose any
                ! (normalised) selection scheme we want...
                ! Given the contribution to the projected energy is divided by the cluster generation probability and
                ! multiplied by the actual weight, doing this has absolutely no effect on the projected energy.

                call cumulative_population(qs%psip_list%pops, qs%psip_list%states(sys%basis%tot_string_len,:), &
                                           qs%psip_list%nstates, D0_proc, D0_pos, qs%psip_list%pop_real_factor, &
                                           ccmc_in%even_selection, sys%read_in%comp, cumulative_abs_real_pops, &
                                           tot_abs_real_pop)

                call update_bloom_threshold_prop(bloom_stats, nparticles_old(1))


                ! Evolution is done via:

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
                !$omp private(it, iexcip_pos, i, seen_D0, hit, pos, population, real_population,k, &
                !$omp annihilation_flags) &
                !$omp shared(rng, cumulative_abs_real_pops, tot_abs_real_pop,  &
                !$omp        max_cluster_size, contrib, D0_normalisation, D0_pos, rdm,    &
                !$omp        qs, sys, bloom_stats, min_cluster_size, ref_det,             &
                !$omp        selection_data,      &
                !$omp        uccmc_in, ccmc_in, nprocs, ms_stats, ps_stats, qmc_in, load_bal_in, &
                !$omp        count_discard, &  
                !$omp        logging_info, nstates_ci, & 
                !$omp        time_avg_psip_list_ci_states, time_avg_psip_list_ci_pops, &
                !$omp        time_avg_psip_list_states, time_avg_psip_list_pops) &
                !$omp reduction(+:D0_population_cycle,proj_energy_cycle,D0_population_noncomp_cycle, &
                !$omp nattempts_spawn,ndeath,nparticles_change,ndeath_nc)
                it = get_thread_id()
                iexcip_pos = 0
                seen_D0 = .false.

                !$omp do schedule(dynamic,200) 
                do iattempt = 1, selection_data%nsingle_excitors + selection_data%nstochastic_clusters
                    if (iattempt <= selection_data%nsingle_excitors) then
                    ! For OpenMP scalability, have this test inside a single loop rather
                    ! than attempt to parallelise over three separate loops.
                        ! As noncomposite clusters can't be above truncation level or linked-only all can accumulate +
                        ! propagate. Only need to check not selecting the reference as we treat it separately.
                        if (iattempt /= D0_pos) then
                            ! Deterministically select each excip as a non-composite cluster.
                            call select_nc_cluster(sys, qs%psip_list, qs%ref%f0, &
                                        iattempt, qmc_in%initiator_pop, .false., &
                                        contrib(it)%cdet, contrib(it)%cluster, qs%excit_gen_data)

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
                        else
                        call select_ucc_cluster(rng(it), sys, qs%psip_list, qs%ref%f0, qs%ref%max_ex_level, &
                                            selection_data%nstochastic_clusters, D0_normalisation, qmc_in%initiator_pop, D0_pos, &
                                            cumulative_abs_real_pops, tot_abs_real_pop, min_cluster_size, max_cluster_size, &
                                            logging_info, contrib(it)%cdet, contrib(it)%cluster, qs%excit_gen_data, &
                                            uccmc_in%threshold, count_discard)

                        ! Add contribution to average CI wfn
                        if (uccmc_in%variational_energy .and. all(qs%vary_shift) .and. &
                            contrib(it)%cluster%excitation_level <= qs%ref%ex_level)  then
                            call add_ci_contribution(contrib(it)%cluster, contrib(it)%cdet, &
                            time_avg_psip_list_ci_states, time_avg_psip_list_ci_pops, nstates_ci)
                        end if

                        if (contrib(it)%cluster%excitation_level <= qs%ref%max_ex_level+2) then
                            if (qs%propagator%quasi_newton) contrib(it)%cdet%fock_sum = &
                                            sum_sp_eigenvalues_occ_list(sys, contrib(it)%cdet%occ_list) - qs%ref%fock_sum

                            call do_ccmc_accumulation(sys, qs, contrib(it)%cdet, contrib(it)%cluster, logging_info, &
                                                    D0_population_cycle, proj_energy_cycle, ccmc_in, ref_det, rdm, &
                                                    selection_data, D0_population_noncomp_cycle)
                            call do_stochastic_ccmc_propagation(rng(it), sys, qs, &
                                                                ccmc_in, logging_info, ms_stats(it), bloom_stats, &
                                                                contrib(it), nattempts_spawn, ndeath, ps_stats(it), uccmc_in)
                        end if
                    end if
                end do
                !$omp end do
                ! See comments below 'if (.not. seen_D0) then' on why this loop needs to be separate from above.
                ! If noncomposite is turned off, this loop will be 'do i = nclusters+1, nclusters', which will be a null 
                ! loop and ignored (as strides at +1 by default)
                !$omp do schedule(dynamic,200) 
                do iattempt = 1, selection_data%nD0_select
                    ! We just select the empty cluster.
                    ! As in the original algorithm, allow this to happen on
                    ! each processor and hence scale the selection
                    ! probability by nprocs.  See comments in select_cluster
                    ! for more details.
                        if (.not. seen_D0) then
                            ! This is the first time this thread is spawning from D0, so it
                            ! needs to be converted into a det_info_t object for the excitation
                            ! generators. On subsequent calls, cdet does not need to change.
                            seen_D0 = .true.
                            call create_null_cluster(sys, qs%ref%f0, nprocs*real(selection_data%nD0_select,p), D0_normalisation, &
                                                     qmc_in%initiator_pop, contrib(it)%cdet, contrib(it)%cluster, qs%excit_gen_data)
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
                    !$omp do schedule(dynamic,200) private(dfock) 
                    do iattempt = 1, qs%psip_list%nstates
                        ! Note we use the (encoded) population directly in stochastic_ccmc_death_nc
                        ! (unlike the stochastic_ccmc_death) to avoid unnecessary decoding/encoding
                        ! steps (cf comments in stochastic_death for FCIQMC).
                        if (qs%propagator%quasi_newton) then
                            dfock = sum_fock_values_bit_string(sys, qs%propagator%sp_fock, qs%psip_list%states(:,iattempt)) &
                                - qs%ref%fock_sum
                        end if
                        call stochastic_ccmc_death_nc(rng(it), ccmc_in%linked, qs, iattempt==D0_pos, dfock, &
                                          qs%psip_list%dat(1,iattempt), qs%estimators(1)%proj_energy_old, &
                                          qs%psip_list%pops(1, iattempt), nparticles_change(1), ndeath_nc, &
                                          logging_info)
                    end do
                    !$omp end do
                end if
                !$omp end parallel

                ! Add the accumulated ps_stats data to qs%excit_gen_data%p_single_double.
                if (qs%excit_gen_data%p_single_double%vary_psingles) then
                    call ps_stats_reduction_update(qs%excit_gen_data%p_single_double%rep_accum, ps_stats)
                end if

                if (ccmc_in%density_matrices .and. qs%vary_shift(1) .and. parent .and. .not. sys%read_in%comp) then
                    ! Add in diagonal contribution to RDM (only once per cycle not each time reference
                    ! is selected as this is O(N^2))
                    call update_rdm(sys, ref_det%f, ref_det%f, ref_det%occ_list, real(D0_normalisation,p), 1.0_p, 1.0_p, rdm)
                end if


                qs%psip_list%nparticles = qs%psip_list%nparticles + nparticles_change
                qs%estimators%D0_population_comp = qs%estimators%D0_population_comp + D0_population_cycle
                qs%estimators%proj_energy_comp = qs%estimators%proj_energy_comp + proj_energy_cycle
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

                if(uccmc_in%variational_energy .and. all(qs%vary_shift)) then
                          time_avg_psip_list_ci_pops(:nstates_ci) =  time_avg_psip_list_ci_pops(:nstates_ci)/(iter-avg_start+1)
                end if

                ! If the shift has started before this iteration, add contributions to the average wavefunction and divide by 
                ! new number of iterations. 
                if(all(qs%vary_shift) .and. old_vary .and. uccmc_in%average_wfn) then
                    ! Add current wfn value average.
                    call add_t_contributions(qs%psip_list, time_avg_psip_list_states, time_avg_psip_list_pops, &
                                             time_avg_psip_list_sq, nstates_sq, .false.)
                end if
                call end_mc_cycle(nspawn_events, ndeath_nc, qs%psip_list%pop_real_factor, nattempts_spawn, qs%spawn_store%rspawn)
            end do

            update_tau = bloom_stats%nblooms_curr > 0

            if (ccmc_in%density_matrices .and. qs%vary_shift(1)) then
                call calc_rdm_energy(sys, qs%ref, rdm, qs%estimators(1)%rdm_energy, qs%estimators(1)%rdm_trace)
            end if

            error = qs%spawn_store%spawn%error .or. qs%psip_list%error

            qs%estimators%D0_population = real(qs%estimators%D0_population_comp,p)
            qs%estimators%proj_energy = real(qs%estimators%proj_energy_comp,p)
 
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

        if (parent .and. uccmc_in%average_wfn) then
            ! Take average of wavefunction.
            time_avg_psip_list_pops(:nstates_sq) =  time_avg_psip_list_pops(:nstates_sq)/(iter-avg_start+1)
            time_avg_psip_list_sq(sys%basis%tot_string_len+1,:nstates_sq) &
                = time_avg_psip_list_sq(sys%basis%tot_string_len+1,:nstates_sq)/(iter-avg_start+1)
        end if

        if (parent) write (io_unit,'()')

        if (parent .and. uccmc_in%average_wfn) then
            call write_average_wfn(sys, time_avg_psip_list_pops, time_avg_psip_list_sq,io_unit,time_avg_psip_list_states,nstates_sq)
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
            call var_energy_uccmc(sys, time_avg_psip_list_ci_states,time_avg_psip_list_ci_pops,nstates_ci, var_energy, &
                                  real(D0_normalisation,p))
            print*, 'Variational energy: ', var_energy
        end if 
        
        if (debug) call end_logging(logging_info)
        if (debug) call end_selection_data(selection_data)

        if (ccmc_in%density_matrices) then
            call write_final_rdm(rdm, sys%nel, sys%basis%nbasis, ccmc_in%density_matrix_file, io_unit)
            call calc_rdm_energy(sys, qs%ref, rdm, qs%estimators(1)%rdm_energy, qs%estimators(1)%rdm_trace)
            if (parent) &
                write (io_unit,'(1x,"# Final energy from RDM",2x,es17.10, /)') qs%estimators%rdm_energy/qs%estimators%rdm_trace
            deallocate(rdm, stat=ierr)
            call check_deallocate('rdm',ierr)
            call dealloc_det_info_t(ref_det)
        end if

        if (uccmc_in%threshold > 0 .and. parent) print*, 'Number of discard events: ', count_discard

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

    end subroutine do_uccmc

    subroutine select_ucc_cluster(rng, sys, psip_list, f0, ex_level, nattempts, normalisation, &
                              initiator_pop, D0_pos, cumulative_excip_pop, tot_excip_pop, min_size, max_size, &
                              logging_info, cdet, cluster, excit_gen_data, threshold, counter)

        ! Based on select_cluster (without the linked_cluster parts) and with information about de-excitors.
        ! Select a random cluster of excitors and deexcitors from the excitors on the
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

        type(sys_t), intent(in) :: sys
        type(particle_t), intent(in), target :: psip_list
        integer(i0), intent(in) :: f0(sys%basis%tot_string_len)
        integer, intent(in) :: ex_level
        integer(int_64), intent(in) :: nattempts
        integer, intent(in) :: D0_pos
        complex(p), intent(in) :: normalisation
        real(p), intent(in) :: initiator_pop
        real(p), intent(in) :: cumulative_excip_pop(:), tot_excip_pop
        integer :: min_size, max_size
        type(dSFMT_t), intent(inout) :: rng
        type(det_info_t), intent(inout) :: cdet
        type(cluster_t), intent(inout) :: cluster
        type(logging_t), intent(in) :: logging_info
        type(excit_gen_data_t), intent(in) :: excit_gen_data
        real(p), intent(in) :: threshold
        integer, intent(inout) :: counter

        real(dp) :: rand
        real(p) :: psize
        complex(p) :: cluster_population, excitor_pop
        integer :: i, pos, prev_pos, deexcit_count, ierr
        real(p) :: pop(max_size)
        logical :: hit, allowed
        logical, allocatable :: deexcitation(:)
        real(p) :: normal
       
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
        deexcit_count = 0

!        normal= 0
!        do i = min_size, max_size
!            normal = normal + ((2.0_p**(i-1))/factorial(i))*(tot_excip_pop/real(normalisation, p))**(i-1)
!        end do
        do i = 0, max_size-min_size-1
            psize = psize + 1.0_p/2_int_64**(i+1)
!        do i = min_size, max_size-1
!            psize = psize + ((2.0_p**(i-1))/factorial(i))*(tot_excip_pop/real(normalisation, p))**(i-1)
            if (rand < psize) then
                ! Found size!
                cluster%nexcitors = min_size + i
                cluster%pselect = cluster%pselect/2_int_64**(i+1)
!                cluster%pselect = cluster%pselect*(((2.0_p**(i-1))/factorial(i))*(tot_excip_pop/real(normalisation, p))**(i-1))/normal
                exit
            end if
        end do
        ! If not set, then must be the largest possible cluster
        if (cluster%nexcitors == -1) then
            cluster%nexcitors = max_size
            cluster%pselect = cluster%pselect*(1.0_p - psize)
        end if

        if(cluster%nexcitors>0) then
              allocate(deexcitation(cluster%nexcitors))
              deexcitation(:) = .false.
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

        select case(cluster%nexcitors)
        case(0)
            call create_null_cluster(sys, f0, cluster%pselect, normalisation, initiator_pop, &
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
            ! For UCCMC, also consider whether each operator acts as an excitation or
            ! de-excitation operator. A string of all-excitation operators gives the
            ! same overall excitor regardless of operator order. This is not true for
            ! mixed excitation-deexcitation strings, so the operators can no longer be
            ! sorted before application. 

            do i = 1, cluster%nexcitors
                ! Select a position in the excitors list.
                pop(i) = get_rand_close_open(rng)*tot_excip_pop
                rand = get_rand_close_open(rng)
                ! Decide if this is to be an excitation or deexcitation operator.
                ! First operator applied must always be an excitor.
                if(i > 1 .and. rand > 0.5) then
                   deexcitation(i) = .true.
                   deexcit_count = deexcit_count + 1
                end if
            end do
            !For the sign of the cluster population we are interested only in the parity of
            !the number of deexcitation operators.
            if (mod(deexcit_count,2)==1) then
               deexcit_count = -1
            else
               deexcit_count = 1
            end if
            !There are 2^(n-1) possible arrangements of excitation/deexcitation operators for
            ! a given cluster, so the probability of selecting one must be divided by this.
            cluster%pselect = cluster%pselect/(2**(cluster%nexcitors-1))

            prev_pos = 1
            do i = 1, cluster%nexcitors
                call binary_search(cumulative_excip_pop, pop(i), 1, psip_list%nstates, hit, pos)
                ! Not allowed to select the reference as it is not an excitor.
                ! Because we treat (for the purposes of the cumulative
                ! population) the reference to have 0 excips, then
                ! cumulative_excip_pop(D0_pos) = cumulative_excip_pop(D0_pos-1).
                ! The binary search algorithm assumes each value in the array
                ! being searched is unique, which is not true, so we can
                ! accidentally find D0_pos.  As we want to find pos such that
                ! cumulative_excip_pop(pos-1) < pop <= cumulative_excip_pop(pos),
                ! then this means we actually need the slot before D0_pos.
                ! Correcting for this accident is much easier than producing an
                ! array explicitly without D0...
                ! If contain multiple spaces we can have this in a more general
                ! case, where an excitor has population in another space but not
                ! that which we're currently concerned with. More general test
                ! should account for this.
                do
                    if (pos == 1) then
                       exit
                    end if
                    if (abs(cumulative_excip_pop(pos) - cumulative_excip_pop(pos-1)) > depsilon) exit
                    pos = pos - 1
                end do

                if (sys%read_in%comp) then
                    excitor_pop = cmplx(psip_list%pops(1,pos), psip_list%pops(2,pos),p)/psip_list%pop_real_factor
                else
                    excitor_pop = real(psip_list%pops(1,pos),p)/psip_list%pop_real_factor
                end if
                if (i == 1) then
                    ! First excitor 'seeds' the cluster:
                    cdet%f = psip_list%states(:,pos)
                    cdet%data => psip_list%dat(:,pos) ! Only use if cluster is non-composite!
                    cluster_population = excitor_pop*deexcit_count
                    ! Counter the additional *nprocs above.
                    cluster%pselect = cluster%pselect/nprocs
                else
                    call ucc_collapse_cluster(sys%basis, f0, psip_list%states(:,pos), excitor_pop, cdet%f, &
                                          cluster_population, allowed, deexcitation(i))
                    if (.not.allowed) then
                        cluster%excitation_level = huge(0)
                        exit
                    end if
                    ! Each excitor spends the same amount of time on each processor on
                    ! average.  If this excitor is different from the previous excitor,
                    ! then the probability this excitor is on the same processor as the
                    ! previous excitor is 1/nprocs.  (Note choosing the same excitor
                    ! multiple times is valid in linked CC.)
                    ![todo] MPI has not been considered for UCCMC. Please implement.
                    if (pos /= prev_pos) cluster%pselect = cluster%pselect/nprocs
                end if
                ! If the excitor's population is below the initiator threshold, we remove the
                ! initiator status for the cluster
                if (abs(excitor_pop) <= initiator_pop) cdet%initiator_flag = 3
                ! Probability of choosing this excitor = pop/tot_pop.
                cluster%pselect = (cluster%pselect*abs(excitor_pop))/tot_excip_pop
                cluster%excitors(i)%f => psip_list%states(:,pos)
                prev_pos = pos
            end do

            if (allowed) then
                cluster%excitation_level = get_excitation_level(f0, cdet%f)
                ! To contribute the cluster must be within a double excitation of
                ! the maximum excitation included in the CC wavefunction.
                allowed = cluster%excitation_level <= ex_level+2
            end if

            if (allowed) then
                ! We chose excitors with a probability proportional to their
                ! occupation.  However, because (for example) the cluster t_X t_Y
                ! and t_Y t_X collapse onto the same excitor (where X and Y each
                ! label an excitor), the probability of selecting a given cluster is
                ! proportional to the number of ways the cluster could have been
                ! formed.  (One can view this factorial contribution as the
                ! factorial prefactors in the series expansion of e^T---see Eq (8)
                ! in the module-level comments.)
                ! If two excitors in the cluster are the same, the factorial
                ! overcounts the number of ways the cluster could have been formed
                ! but the extra factor(s) of 2 are cancelled by a similar
                ! overcounting in the calculation of hmatel.
                cluster%pselect = cluster%pselect*factorial(cluster%nexcitors)

                ! Sign change due to difference between determinant
                ! representation and excitors and excitation level.
                call convert_excitor_to_determinant(cdet%f, cluster%excitation_level, cluster%cluster_to_det_sign, f0)
                call decoder_ptr(sys, cdet%f, cdet, excit_gen_data)

                ! Normalisation factor for cluster%amplitudes...
                cluster%amplitude = cluster_population/(normalisation**(cluster%nexcitors-1))
                
                if (abs(cluster%amplitude)/cluster%pselect > threshold) then
                    allowed = .false.
                    cluster%excitation_level = huge(0)
                    !$omp atomic update
                    counter = counter + 1
                    !$omp end atomic
                end if
            else
                ! Simply set excitation level to a too high (fake) level to avoid
                ! this cluster being used.
                cluster%excitation_level = huge(0)
            end if
        end select

        if (cluster%nexcitors > 0)  then
            deallocate(deexcitation, stat=ierr)
            call check_deallocate('deexcitation',ierr)
        end if 
        if (debug) call write_logging_stoch_selection(logging_info, cluster%nexcitors, cluster%excitation_level, pop, &
                max_size, cluster%pselect, cluster%amplitude, allowed)

    end subroutine select_ucc_cluster

    subroutine ucc_collapse_cluster(basis, f0, excitor, excitor_population, cluster_excitor, cluster_population, &
                    allowed, conjugate)

        ! Collapse two excitors.  The result is returned in-place. Based on collapse_cluster in 
        ! ccmc_selection

        ! In:
        !    basis: information about the single-particle basis.
        !    f0: bit string representation of the reference determinant.
        !    excitor: bit string of the Slater determinant formed by applying
        !        the excitor, e1, to the reference determinant.
        !    excitor_population: number of excips on the excitor e1.
        !    conjugate : true/false, encodes whether the current excitor is to be used as a 
        !        deexcitor or excitor.
        ! In/Out:
        !    cluster_excitor: bit string of the Slater determinant formed by applying
        !        the excitor, e2, to the reference determinant.
        !    cluster_population: number of excips on the 'cluster' excitor, e2.
        ! Out:
        !    allowed: true if excitor e1 can be applied to excitor e2 (i.e. e1
        !       and e2 do not involve exciting from/to the any identical
        !       spin-orbitals).

        ! On input, cluster excitor refers to an existing excitor, e2.  On output,
        ! cluster excitor refers to the excitor formed from applying the excitor
        ! e1 to the cluster e2.
        ! ***WARNING***: if allowed is false then cluster_excitor is *not* updated.

        use basis_types, only: basis_t, reset_extra_info_bit_string

        use bit_utils, only: count_set_bits
        use const, only: i0_end
        use ccmc_utils, only: collapse_excitor_onto_cluster
        use uccmc_utils, only: collapse_deexcitor_onto_cluster

        type(basis_t), intent(in) :: basis
        integer(i0), intent(in) :: f0(basis%tot_string_len)
        integer(i0), intent(in) :: excitor(basis%tot_string_len)
        complex(p), intent(in) :: excitor_population
        integer(i0), intent(inout) :: cluster_excitor(basis%tot_string_len)
        complex(p), intent(inout) :: cluster_population
        logical,  intent(out) :: allowed
        logical, intent(in) :: conjugate

        integer(i0) :: excitor_loc(basis%tot_string_len), cluster_loc(basis%tot_string_len), f0_loc(basis%tot_string_len)

        integer(i0) :: excitor_excitation(basis%tot_string_len)
        integer(i0) :: excitor_annihilation(basis%tot_string_len)
        integer(i0) :: excitor_creation(basis%tot_string_len)
        integer(i0) :: cluster_excitation(basis%tot_string_len)
        integer(i0) :: cluster_annihilation(basis%tot_string_len)
        integer(i0) :: cluster_creation(basis%tot_string_len)
        integer :: ibasis, ibit

        excitor_loc = excitor
        cluster_loc = cluster_excitor
        f0_loc = f0
        call reset_extra_info_bit_string(basis, excitor_loc)
        call reset_extra_info_bit_string(basis, cluster_loc)
        call reset_extra_info_bit_string(basis, f0_loc)

        ! Apply excitor to the cluster of excitors.

        ! orbitals involved in excitation from reference
        excitor_excitation = ieor(f0_loc, excitor_loc)
        cluster_excitation = ieor(f0_loc, cluster_loc)
        ! Annihilation/creation operators are reversed between excitation/deexcitation operators,
        ! so they must be found appropriately
        ! annihilation operators (relative to the reference)
        if (conjugate) then
            excitor_annihilation = excitor_excitation - iand(excitor_excitation, f0_loc)
        else
            excitor_annihilation = iand(excitor_excitation, f0_loc)
        end if
        cluster_annihilation = iand(cluster_excitation, f0_loc)
        ! creation operators (relative to the reference)
        if (conjugate) then
            excitor_creation = iand(excitor_excitation, f0_loc)
        else
            excitor_creation = iand(excitor_excitation, excitor_loc)
        end if
        cluster_creation = iand(cluster_excitation, cluster_loc)

        ! First, let's find out if the excitor is valid...
        if(conjugate) then
           ! Check deexcitation is possible
           allowed = .true.
           do ibasis = 1, basis%bit_string_len
                do ibit = 0, i0_end
                    if (btest(excitor_annihilation(ibasis),ibit)) then
                        if (.not. btest(cluster_creation(ibasis),ibit)) then
                           allowed = .false.
                           exit
                        end if
                    end if
                    if (btest(excitor_creation(ibasis),ibit)) then
                        if (btest(cluster_loc(ibasis),ibit)) then
                           allowed = .false.
                           exit
                        end if
                    end if
                end do
           end do
        
           if (allowed) then
                cluster_population = cluster_population*excitor_population

                ! Now apply the deexcitor to the cluster (which is, in its own right,
                ! an excitor).
                ! Consider a cluster, e.g. t_{ij}^{ab} = a^+_a a^+_b a_j a_i (which 
                ! corresponds to i,j ->a, b, where i<j and a<b).
                ! We wish to collapse e.g. (t_i^a)^+ with it to get a single 
                ! excitor, t_j^b = a^+_b a_j  However (t_i^a)^+ t_{ij}^{ab} = a^+_j a_a a^+_a a^+_b a_j a_i.  
                ! We thus need to permute the creation and annihilation operators to cancel them out.
                !  Each permutation incurs a sign change.
                call collapse_deexcitor_onto_cluster(basis, excitor_excitation, f0_loc, cluster_excitor, &
                                                  cluster_annihilation, cluster_creation, cluster_population, &
                                                  excitor_annihilation, cluster_excitation)
           end if  
        else 
            if (any(iand(excitor_creation,cluster_creation) /= 0) &
                    .or. any(iand(excitor_annihilation,cluster_annihilation) /= 0)) then
                ! excitor attempts to excite from an orbital already excited from by
                ! the cluster or into an orbital already excited into by the
                ! cluster.
                ! => not valid
                allowed = .false.

                ! We still use the cluster in linked ccmc so need its amplitude
                cluster_population = cluster_population*excitor_population
            else
                ! Applying the excitor to the existing cluster of excitors results
                ! in a valid cluster.
                allowed = .true.

                ! Combine amplitudes.
                ! Might need a sign change as well...see below!
                cluster_population = cluster_population*excitor_population

                call collapse_excitor_onto_cluster(basis, excitor_excitation, f0, cluster_excitor, &
                                                  cluster_annihilation, cluster_creation, cluster_population)
            end if
        end if
    end subroutine ucc_collapse_cluster

end module

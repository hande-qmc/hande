module trotterized_uccmc

use const, only: i0, int_p, int_64, p, dp, debug

implicit none

contains

! [review] - AJWT: Put in a brief summary of the differences between this and conventional CC.
    subroutine do_trot_uccmc(sys, qmc_in, uccmc_in, restart_in, load_bal_in, reference_in, &
                        logging_in, io_unit, qs, qmc_state_restart)

! [review] - AJWT: Perhaps note that this is derived from do_ccmc (and similarly in do_ccmc that this derives from it).  Are there any commonalities which can be sensibly shared as a called function?
! [review] - AJWT: What CCMC features are not present? e.g. blocking, even_selection
        ! Run the Trotterized UCCMC algorithm starting from the initial walker distribution
        ! using the timestep algorithm.

        ! See notes about the implementation of this using function pointers
        ! in fciqmc_main.

        ! In:
        !    sys: system being studied.
        !    qmc_in: input options relating to QMC methods.
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

        use annihilation, only: insert_new_walker
        use bloom_handler, only: init_bloom_stats_t, bloom_stats_t, bloom_mode_fractionn, bloom_mode_fixedn, &
                                 write_bloom_report, bloom_stats_warning, update_bloom_threshold_prop
        use ccmc_data
        use ccmc_selection, only: create_null_cluster
        use ccmc_selection, only: init_selection_data, update_selection_probabilities, set_cluster_selections, &
                                  init_amp_psel_accumulation
        use ccmc_utils, only: get_D0_info, init_contrib, dealloc_contrib, cumulative_population, & 
                              regenerate_ex_levels_psip_list
        use determinants, only: alloc_det_info_t, dealloc_det_info_t, sum_sp_eigenvalues_occ_list, &
                                sum_sp_eigenvalues_bit_string, decode_det
        use determinant_data, only: det_info_t
        use excitations, only: excit_t, get_excitation_level, get_excitation
        use qmc_io, only: write_qmc_report, write_qmc_report_header, write_qmc_var
        use qmc, only: init_qmc, init_secondary_reference
        use qmc_common, only: initial_qmc_status, initial_cc_projected_energy, load_balancing_report, init_report_loop, &
                              init_mc_cycle, end_report_loop, end_mc_cycle, redistribute_particles, rescale_tau
        use proc_pointers
        use spawning, only: assign_particle_processor
        use system, only: sys_t, sys_t_json
        use spawn_data, only: calc_events_spawn_t, write_memcheck_report, spawn_t, alloc_spawn_t
        use replica_rdm, only: update_rdm, calc_rdm_energy, write_final_rdm

        use qmc_data, only: qmc_in_t, uccmc_in_t, restart_in_t

        use qmc_data, only: load_bal_in_t, qmc_state_t, annihilation_flags_t, estimators_t, particle_t
        use qmc_data, only: qmc_in_t_json, uccmc_in_t_json, restart_in_t_json
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
        use uccmc, only: var_energy_uccmc, ucc_set_cluster_selections, do_uccmc_accumulation, &
                         do_stochastic_uccmc_propagation
        use uccmc_utils, only: add_info_str_trot, earliest_unset
        type(sys_t), intent(in) :: sys
        type(qmc_in_t), intent(in) :: qmc_in
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
        real(p) :: D0_population_ucc_cycle

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

        variational = .false.
        if (uccmc_in%variational_energy) variational = .true.
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
            call check_uccmc_opts(sys, uccmc_in, qmc_in)
        end if
        ! Initialise data.
        call init_qmc(sys, qmc_in, restart_in, load_bal_in, reference_in, io_unit, annihilation_flags, qs, &
                      uuid_restart, restart_version_restart, qmc_state_restart=qmc_state_restart, &
                      regenerate_info=regenerate_info, uccmc_in=uccmc_in)

!        qs%psip_list%nstates = 8
!        qs%psip_list%states(1,2) = 33
!        qs%psip_list%states(1,3) = 192
!        qs%psip_list%states(1,4) = 132
!        qs%psip_list%states(1,5) = 72
!        qs%psip_list%states(1,6) = 48
!        qs%psip_list%states(1,7) = 12
!        qs%psip_list%states(1,8) = 18
!
!
!        qs%psip_list%pops(1,1) = 837 * qs%psip_list%pop_real_factor 
!        qs%psip_list%pops(1,2) = int(0.91*qs%psip_list%pop_real_factor) 
!        qs%psip_list%pops(1,3) = -33 * qs%psip_list%pop_real_factor 
!        qs%psip_list%pops(1,4) = 24 * qs%psip_list%pop_real_factor
!        qs%psip_list%pops(1,5) = -24 * qs%psip_list%pop_real_factor
!        qs%psip_list%pops(1,6) = -48 * qs%psip_list%pop_real_factor
!        qs%psip_list%pops(1,7) = -30 * qs%psip_list%pop_real_factor
!        qs%psip_list%pops(1,8) = int(0.99 * qs%psip_list%pop_real_factor)

        ! Add information strings to the psip_list and the reference determinant.
        call regenerate_trot_info_psip_list(sys%basis, sys%nel, qs)
        qs%ref%f0(3) = earliest_unset(qs%ref%f0, qs%ref%f0, sys%nel, sys%basis) 
        qs%ref%f0(2) = sys%nel

        !Allocate memory for time averaged populations and variational energy computation.
        allocate(state(sys%basis%tot_string_len))

        if(uccmc_in%variational_energy) then
             population = 0
             allocate(time_avg_psip_list_ci_states(size(qs%psip_list%states(:,1)),size(qs%psip_list%states(1,:))))
             allocate(time_avg_psip_list_ci_pops(size(qs%psip_list%states(1,:))))
             time_avg_psip_list_ci_states(:,1) = qs%psip_list%states(:,1)
             time_avg_psip_list_ci_pops(1) = (real(qs%psip_list%pops(1,1))/qs%psip_list%pop_real_factor)
             nstates_ci = 1
             ![todo] deal with restarting
        end if

        allocate(time_avg_psip_list_states(size(qs%psip_list%states(:,1)),size(qs%psip_list%states(1,:))))
        allocate(time_avg_psip_list_pops(size(qs%psip_list%states(1,:))))
        allocate(time_avg_psip_list_sq(2,size(qs%psip_list%states(1,:))))
        time_avg_psip_list_states(:,1) = qs%psip_list%states(:,1)
        time_avg_psip_list_pops(1) = (real(qs%psip_list%pops(1,1))/qs%psip_list%pop_real_factor)
        time_avg_psip_list_sq(1,1) = qs%psip_list%states(1,1)
        time_avg_psip_list_sq(2,1) = (real(qs%psip_list%pops(1,1))/qs%psip_list%pop_real_factor)**2
        nstates_sq = 1
        
        qs%ref%max_ex_level = qs%ref%ex_level

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

        call init_contrib(sys, uccmc_in%pow_trunc, uccmc_in%linked, contrib)

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

        if (debug) call init_amp_psel_accumulation(qs%ref%max_ex_level+2, logging_info, uccmc_in%linked, selection_data)

        nparticles_old = qs%psip_list%tot_nparticles

        ! Initialise D0_pos to be somewhere (anywhere) in the list.
        D0_pos = 1

        associate(pl=>qs%psip_list, spawn=>qs%spawn_store%spawn)
            ! Initialise hash shift if restarting...
            spawn%hash_shift = qs%mc_cycles_done
            ! NOTE: currently hash_seed is not exposed and so cannot change unless the hard-coded value changes. Therefore, as we
            ! have not evolved the particles since the were written out (i.e. hash_shift hasn't changed) the only parameter
            ! which can be altered which can change an excitors location since the restart files were written is move_freq.
            if (uccmc_in%move_freq /= spawn%move_freq .and. nprocs > 1) then
                spawn%move_freq = uccmc_in%move_freq
                if (restarting) then
                    if (parent) call warning('do_trot_uccmc', 'move_freq is different from that in the restart file. &
                                            &Reassigning processors. Please check for equilibration effects.')
                    ! Cannot rely on spawn%head being initialised to indicate an empty spawn_t object so reset it before
                    ! redistributing,.
                    spawn%head = spawn%head_start
                    call redistribute_particles(pl%states, pl%pop_real_factor, pl%pops, pl%nstates, pl%nparticles, spawn)
                    call direct_annihilation_trot(sys, rng(0), qs%ref, annihilation_flags, pl, spawn)
                end if
            end if
        end associate

        if (parent) then
            call write_qmc_report_header(qs%psip_list%nspaces, cmplx_est=sys%read_in%comp, rdm_energy=uccmc_in%density_matrices, &
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
            ! [review] Verena: gfortran complained about the "not" so I changed it. Please double check that I have not changed it
            ! [review] Verena: in a wrong way.
                
            call init_report_loop(qs, bloom_stats)

            do icycle = 1, qmc_in%ncycles
                iter = qs%mc_cycles_done + (ireport-1)*qmc_in%ncycles + icycle
                if (all(qs%vary_shift) .and. (.not. old_vary)) then
                    avg_start = iter - 1
                    old_vary = .true.
                    time_avg_psip_list_pops(:qs%psip_list%nstates) = &
                        real(qs%psip_list%pops(1,:qs%psip_list%nstates))/qs%psip_list%pop_real_factor
                    time_avg_psip_list_states(:,:qs%psip_list%nstates) = qs%psip_list%states(:,:qs%psip_list%nstates)
                    time_avg_psip_list_sq(1,:qs%psip_list%nstates) = qs%psip_list%states(1,:qs%psip_list%nstates)
                    time_avg_psip_list_sq(2,:qs%psip_list%nstates) = (real(qs%psip_list%pops(1,:qs%psip_list%nstates))/qs%psip_list%pop_real_factor)**2
                    nstates_sq = qs%psip_list%nstates
                end if
                
                
                ! Recover total number of particles added to time-averages 
                if(uccmc_in%variational_energy) then
                          time_avg_psip_list_ci_pops(:nstates_ci) =  time_avg_psip_list_ci_pops(:nstates_ci)*(iter-1)
                end if

                if(all(qs%vary_shift)) then
                    !time_avg_psip_list%nparticles = time_avg_psip_list%nparticles*(iter-avg_start)
                    time_avg_psip_list_pops(:nstates_sq) =  time_avg_psip_list_pops(:nstates_sq)*(iter - avg_start)

                    time_avg_psip_list_sq(2,:nstates_sq) =  time_avg_psip_list_sq(2,:nstates_sq)*(iter-avg_start)
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

                call ucc_set_cluster_selections(selection_data, qs%estimators(1)%nattempts, min_cluster_size)
                call zero_ps_stats(ps_stats, qs%excit_gen_data%p_single_double%rep_accum%overflow_loc)

                ! OpenMP chunk size determined completely empirically from a single
                ! test.  Please feel free to improve...
                ! NOTE: we can't refer to procedure pointers in shared blocks so
                ! can't use default(none).  I *strongly* recommend turning
                ! default(none) on when making changes and ensure that the only
                ! errors relate to the procedure pointers...

                !$omp parallel default(none) &
                !$omp private(it, iexcip_pos, i, seen_D0, hit, pos, population, real_population,k, &
                !$omp state,annihilation_flags) &
                !$omp shared(rng, cumulative_abs_real_pops, tot_abs_real_pop,  &
                !$omp        max_cluster_size, contrib, D0_normalisation, D0_pos, rdm,    &
                !$omp        qs, sys, bloom_stats, min_cluster_size, ref_det,             &
                !$omp        proj_energy_cycle, D0_population_cycle, selection_data,      &
                !$omp        nattempts_spawn, D0_population_ucc_cycle, &
                !$omp        uccmc_in, nprocs, ms_stats, ps_stats, qmc_in, load_bal_in, &
                !$omp        ndeath_nc,   pcumul, nstates_ci, &
                !$omp        nparticles_change, ndeath, logging_info, &
                !$omp        time_avg_psip_list_ci_states, time_avg_psip_list_ci_pops, &
                !$omp        time_avg_psip_list_states, time_avg_psip_list_pops)

                it = get_thread_id()
                iexcip_pos = 0
                seen_D0 = .false.
                proj_energy_cycle = cmplx(0.0, 0.0, p)
                D0_population_cycle = cmplx(0.0, 0.0, p)
                D0_population_ucc_cycle = 0.0_p
                pcumul = (tot_abs_real_pop**(max_cluster_size+1)-1) /(tot_abs_real_pop-1)
                !$omp do schedule(dynamic,200) reduction(+:D0_population_cycle,proj_energy_cycle, D0_population_ucc_cycle, nattempts_spawn,ndeath)
                do iattempt = 1, selection_data%nclusters
                    ! For OpenMP scalability, have this test inside a single loop rather
                    ! than attempt to parallelise over three separate loops.
                   
                    call select_trot_ucc_cluster(rng(it), sys, qs%psip_list, qs%ref%f0, qs%ref%max_ex_level, &
                                            selection_data%nstochastic_clusters, D0_normalisation, qmc_in%initiator_pop, D0_pos, &
                                            logging_info, contrib(it)%cdet, contrib(it)%cluster, qs%excit_gen_data)
                    if (uccmc_in%trot) call add_info_str_trot(sys%basis, qs%ref%f0, sys%nel, contrib(it)%cdet%f)
                    !Add selected cluster contribution to CI wavefunction estimator.
                    if (uccmc_in%variational_energy .and. (.not. all(contrib(it)%cdet%f==0)) .and. &
                        contrib(it)%cluster%excitation_level <= qs%ref%ex_level)  then
                       state = contrib(it)%cdet%f 
                       call binary_search_i0_list_trot(time_avg_psip_list_ci_states, state, 1, nstates_ci, hit, pos)
                       population = &
                           contrib(it)%cluster%amplitude*contrib(it)%cluster%cluster_to_det_sign/contrib(it)%cluster%pselect
                       if (hit) then
                          time_avg_psip_list_ci_pops(pos) = time_avg_psip_list_ci_pops(pos) + population 
                       else
                           do j = nstates_ci, pos, -1
                               ! i is the number of determinants that will be inserted below j.
                               k = j + 1 
                               time_avg_psip_list_ci_states(:,k) = time_avg_psip_list_ci_states(:,j)
                               time_avg_psip_list_ci_pops(k) = time_avg_psip_list_ci_pops(j)
                           end do

                           time_avg_psip_list_ci_states(:,pos) = state
                           time_avg_psip_list_ci_pops(pos) = population
                           nstates_ci = nstates_ci + 1
                           ! Extract the real sign from the encoded sign.

                       end if
                    end if
                    

                    if (contrib(it)%cluster%excitation_level <= qs%ref%max_ex_level+2) then
                            ! cluster%excitation_level == huge(0) indicates a cluster
                            ! where two excitors share an elementary operator
                        if (qs%propagator%quasi_newton) contrib(it)%cdet%fock_sum = &
                                            sum_sp_eigenvalues_occ_list(sys, contrib(it)%cdet%occ_list) - qs%ref%fock_sum

                        call do_uccmc_accumulation(sys, qs, contrib(it)%cdet, contrib(it)%cluster, logging_info, &
                                                    D0_population_cycle, proj_energy_cycle, uccmc_in, ref_det, rdm, selection_data,&
                                                    D0_population_ucc_cycle)
                        ! [review] Verena: make sure qs is not altered unexpectedly here (it is shared).
                        call do_stochastic_uccmc_propagation(rng(it), sys, qs, &
                                                                uccmc_in, logging_info, ms_stats(it), bloom_stats, &
                                                                contrib(it), nattempts_spawn, ndeath, ps_stats(it))
                    end if

                end do
                !$omp end do

                ndeath_nc=0
                !$omp end parallel
                !print*, count_select
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
                        ! [review] Verena: Is this a fraction of the form a/b/c in the cos? Is this correct?
                        D0_population_ucc_cycle = &
                            D0_population_ucc_cycle/cos(qs%psip_list%pops(1,j)/real(qs%psip_list%pop_real_factor)/D0_normalisation)
                    end if
                end do
                qs%estimators%D0_population_ucc = qs%estimators%D0_population_ucc + D0_population_ucc_cycle

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
                    call direct_annihilation_trot(sys, rng(0), qs%ref, annihilation_flags, pl, spawn)
                end associate

                if (debug) call write_logging_calc_ccmc(logging_info, iter, nspawn_events, ndeath + ndeath_nc, &
                                                        selection_data%nD0_select, &
                                                        selection_data%nclusters, selection_data%nstochastic_clusters, &
                                                       selection_data%nsingle_excitors)

                !Take updated time average of CI wavefunction.
                if(uccmc_in%variational_energy) then
                          time_avg_psip_list_ci_pops(:nstates_ci) =  time_avg_psip_list_ci_pops(:nstates_ci)/(iter)
                end if

                ! Add new contributions to time-average cluster populations.
                if (all(qs%vary_shift)) then
                do i = 1, qs%psip_list%nstates
                    state = qs%psip_list%states(:,i) 
                    call binary_search_i0_list_trot(time_avg_psip_list_states, state, 1, nstates_sq, hit, pos)
                    !call binary_search(time_avg_psip_list%states, state, 1, time_avg_psip_list%nstates, hit, pos)
                    if (hit) then
                          time_avg_psip_list_pops(pos) = &
                              time_avg_psip_list_pops(pos) + (real(qs%psip_list%pops(1,i))/qs%psip_list%pop_real_factor)
                          time_avg_psip_list_sq(2,pos) = &
                              time_avg_psip_list_sq(2,pos) + (real(qs%psip_list%pops(1,i))/qs%psip_list%pop_real_factor)**2
                       else
                           do j = nstates_sq, pos, -1
                               ! i is the number of determinants that will be inserted below j.
                               k = j + 1 
                               time_avg_psip_list_states(:,k) = time_avg_psip_list_states(:,j)
                               time_avg_psip_list_pops(k) = time_avg_psip_list_pops(j)
                               time_avg_psip_list_sq(1,k) = time_avg_psip_list_sq(1,j)
                               time_avg_psip_list_sq(2,k) = time_avg_psip_list_sq(2,j)
                           end do
                           time_avg_psip_list_states(:,pos) = qs%psip_list%states(:,i)
                           time_avg_psip_list_pops(pos) = (real(qs%psip_list%pops(1,i))/qs%psip_list%pop_real_factor)
                           time_avg_psip_list_sq(1,pos) = qs%psip_list%states(1,i)
                           time_avg_psip_list_sq(2,pos) = (real(qs%psip_list%pops(1,i))/qs%psip_list%pop_real_factor)**2
                           nstates_sq = nstates_sq + 1
                           ! Extract the real sign from the encoded sign.
                           !real_population = real(time_avg_psip_list%pops(:,pos))/time_avg_psip_list%pop_real_factor
                           !time_avg_psip_list%nparticles = time_avg_psip_list%nparticles + abs(real_population)
                       end if
                end do
                ! Update time average.
                !time_avg_psip_list%nparticles = time_avg_psip_list%nparticles/(iter-avg_start+1)
                time_avg_psip_list_pops(:nstates_sq) =  time_avg_psip_list_pops(:nstates_sq)/(-avg_start+iter+1)

                time_avg_psip_list_sq(2,:nstates_sq) =  time_avg_psip_list_sq(2,:nstates_sq)/(-avg_start+iter+1)
                end if
                call end_mc_cycle(nspawn_events, ndeath_nc, qs%psip_list%pop_real_factor, nattempts_spawn, qs%spawn_store%rspawn)
            end do

            update_tau = bloom_stats%nblooms_curr > 0

            error = qs%spawn_store%spawn%error .or. qs%psip_list%error

            qs%estimators%D0_population = real(qs%estimators%D0_population_comp,p)
            qs%estimators%proj_energy = real(qs%estimators%proj_energy_comp,p)
            if (uccmc_in%variational_energy) then
                call var_energy_uccmc(sys, time_avg_psip_list_ci_states, time_avg_psip_list_ci_pops, nstates_ci, var_energy, real(D0_normalisation,p))
                !qs%estimators%var_energy = var_energy
            end if 
            if (debug) call write_logging_select_ccmc(logging_info, iter, selection_data)
          
            call end_report_loop(io_unit, qmc_in, iter, update_tau, qs, nparticles_old, nspawn_events, &
                                 -1, semi_stoch_it, soft_exit=soft_exit, &
                                 load_bal_in=load_bal_in, bloom_stats=bloom_stats, comp=sys%read_in%comp, &
                                 error=error, vary_shift_reference=uccmc_in%vary_shift_reference)
            if (error) exit

            call cpu_time(t2)
            
            if (parent) then
                if (bloom_stats%nblooms_curr > 0) call bloom_stats_warning(bloom_stats, io_unit=io_unit)
                !if (all(qs%vary_shift)) then
                !write (io_unit, '(1X, "Cluster populations",/)')
                !do i = 1, qs%psip_list%nstates
                !    call write_qmc_var(io_unit, qs%psip_list%states(1,i))
                !    call write_qmc_var(io_unit, qs%psip_list%pops(1, i))
                !    write (io_unit,'()')
                !end do
                !else
                call write_qmc_report(qmc_in, qs, ireport, nparticles_old, t2-t1, .false., .false., &
                                        io_unit=io_unit, cmplx_est=sys%read_in%comp, rdm_energy=uccmc_in%density_matrices, &
                                        nattempts=.true.)
                !end if
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

        if (parent) then
            write (io_unit, '(1X, "Time-averaged cluster populations",/)')
            do i = 1, nstates_sq
                call write_qmc_var(io_unit, time_avg_psip_list_states(1,i))
                call write_qmc_var(io_unit, time_avg_psip_list_states(2,i))
                call write_qmc_var(io_unit, time_avg_psip_list_states(3,i))
                call write_qmc_var(io_unit, time_avg_psip_list_pops(i))
                write (io_unit,'()')
            end do
            write (io_unit, '(1X, "Time-averaged cluster populations squared",/)')
            do i = 1, nstates_sq
                call write_qmc_var(io_unit, int(time_avg_psip_list_sq(1,i)))
                call write_qmc_var(io_unit, time_avg_psip_list_sq(2,i))
                write (io_unit,'()')
            end do
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
            call var_energy_uccmc(sys, time_avg_psip_list_ci_states, time_avg_psip_list_ci_pops, nstates_ci, var_energy, real(D0_normalisation,p))
            print*, 'Variational energy: ', var_energy
        end if 
        
        if (debug) call end_logging(logging_info)
        if (debug) call end_selection_data(selection_data)

        call dealloc_contrib(contrib, uccmc_in%linked)
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
        ! [review] - Verena: consider further deallocations.
    
    end subroutine do_trot_uccmc

    subroutine select_trot_ucc_cluster(rng, sys, psip_list, f0, ex_level, nattempts, normalisation, &
                              initiator_pop, D0_pos, logging_info, cdet, cluster, excit_gen_data)
! [review] - AJWT: Note that this is based on select_cluster
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
        !    logging_info: derived type containing information on currently logging status
        !    excit_gen_data: information about excitation generators

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
        use ccmc_utils, only: convert_excitor_to_determinant, collapse_cluster
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
        type(dSFMT_t), intent(inout) :: rng
        type(det_info_t), intent(inout) :: cdet
        type(cluster_t), intent(inout) :: cluster
        type(logging_t), intent(in) :: logging_info
        type(excit_gen_data_t), intent(in) :: excit_gen_data

        real(dp) :: rand
        real(p) :: pexcit
        complex(p) :: cluster_population, excitor_pop
        integer :: i, j, ierr
        logical :: allowed, conjugate
        integer :: order(12)
        real(p) :: pop_real, ref_real

        !order = (/3, 12, 48, 192, 36, 24, 132, 72, 18, 66, 33, 129/)
        ! We shall accumulate the factors which comprise cluster%pselect as we go.
        !   cluster%pselect = n_sel * p_excit(1) * p_excit (2) * ...
        ! where
        !   n_sel   is the number of cluster selections made;
        !   p_excit is the probability of choosing/not choosing a particular excitor to be in the cluster;
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

        ![todo] psize?

        cluster%pselect = real(nattempts*nprocs, p)
        cluster%nexcitors = 0

        ! Initiator approximation.
        ! This is sufficiently quick that we'll just do it in all cases, even
        ! when not using the initiator approximation.  This matches the approach
        ! used by Alex Thom in 'Initiator Stochastic Coupled Cluster Theory'
        ! (unpublished).
        ! Assume all excitors in the cluster are initiators (initiator_flag=0)
        ! until proven otherwise (initiator_flag=1).
        cdet%initiator_flag = 0
        cdet%f = f0
        cluster_population = 1
        ! Assume cluster is allowed unless collapse_cluster finds out otherwise
        ! when collapsing/combining excitors or if it could never have been
        ! valid

! [review] - AJWT: It may be more sensible eventually to consider a 'base probability' being the probability of not selecting any excips,
! [review] - AJWT: and then for each excip we do choose, multiply it by pexcip/(1-pexcip).  This would allow us to use a similar algorithm to
! [review] - AJWT: select_cluster.
        allowed = .true. 
        ref_real = real(psip_list%pops(1, D0_pos))/real(psip_list%pop_real_factor)
        !do j =1, 12
        do i = 1, psip_list%nstates
         !   if(psip_list%states(1,i) == order(j)) then
            if (i /= D0_pos) then
                rand = get_rand_close_open(rng)
                pop_real = real(psip_list%pops(1, i))/real(psip_list%pop_real_factor)
                !The probability of selecting a given excitor is given by abs(N_i/(N_i+N_0))
                pexcit = abs(sin(pop_real/ref_real))/(abs(cos(pop_real/ref_real))+ abs(sin(pop_real/ref_real)))
                if (cluster%nexcitors == 0) then
                   ! print*, 'still 0', rand, pexcit
                    if (rand < pexcit) then
                       !print*, 'adding first excitor'
                       cluster%nexcitors = cluster%nexcitors + 1
                       cluster%pselect = cluster%pselect*pexcit
                       excitor_pop = sin(pop_real/ref_real)
                       ! First excitor 'seeds' the cluster:
                       cdet%f = psip_list%states(:,i)
                       cdet%data => psip_list%dat(:,i) ! Only use if cluster is non-composite!
                       cluster_population = cluster_population*excitor_pop
                       ! Counter the additional *nprocs above.
                       cluster%pselect = cluster%pselect/nprocs
                       !print*, 'now cdet%f is', cdet%f(1)
                    else
                       excitor_pop = cos(pop_real/ref_real)
                       cluster%pselect = cluster%pselect*(1-pexcit)
                       cluster_population = cluster_population* excitor_pop
                    end if
                else
                   conjugate = .false.
                    if (all(iand(ieor(f0(:sys%basis%bit_string_len),psip_list%states(:sys%basis%bit_string_len,i)), &
                           ieor(f0(:sys%basis%bit_string_len),cdet%f(:sys%basis%bit_string_len))) == &
                           ieor(f0(:sys%basis%bit_string_len),psip_list%states(:sys%basis%bit_string_len,i)))) then
          !             !print*, 'conjugate'
                       conjugate = .true.
                       if (rand < pexcit) then
          !                !print*, 'applying conjugate'
                          excitor_pop = -sin(pop_real/ref_real)
                          cluster%nexcitors = cluster%nexcitors + 1
                          cluster%pselect = cluster%pselect*pexcit
                          call trot_collapse_cluster(sys%basis, f0, psip_list%states(:,i), excitor_pop, cdet%f, &
                                          cluster_population, allowed, conjugate)
                          cluster%excitors(cluster%nexcitors)%f => psip_list%states(:,i)
                          if (abs(excitor_pop) <= initiator_pop) cdet%initiator_flag = 3
          !             !print*, 'now cdet%f is', cdet%f(1)
                       else
                          excitor_pop = cos(pop_real/ref_real)
                          cluster%pselect = cluster%pselect*(1-pexcit)
                          cluster_population = cluster_population*excitor_pop
                       end if
                    else if(all(iand(ieor(f0(:sys%basis%bit_string_len),psip_list%states(:sys%basis%bit_string_len,i)), &
                           ieor(f0(:sys%basis%bit_string_len),cdet%f(:sys%basis%bit_string_len))) == 0)) then
                       
                       if (rand < pexcit) then
                          !print*, 'applying normal'
                          excitor_pop = sin(pop_real/ref_real)
                          cluster%nexcitors = cluster%nexcitors + 1
                          cluster%pselect = cluster%pselect*pexcit
                          call trot_collapse_cluster(sys%basis, f0, psip_list%states(:,i), excitor_pop, cdet%f, &
                                          cluster_population, allowed, conjugate)
                          cluster%excitors(cluster%nexcitors)%f => psip_list%states(:,i)
                          if (abs(excitor_pop) <= initiator_pop) cdet%initiator_flag = 3
                       !print*, 'now cdet%f is', cdet%f(1)
                       else
                          excitor_pop = cos(pop_real/ref_real)
                          cluster%pselect = cluster%pselect*(1-pexcit)
                          cluster_population = cluster_population*excitor_pop
                       end if
                    end if
                    ! Each excitor spends the same amount of time on each processor on
                    ! average.  If this excitor is different from the previous excitor,
                    ! then the probability this excitor is on the same processor as the
                    ! previous excitor is 1/nprocs.  (Note choosing the same excitor
                    ! multiple times is valid in linked CC.)
                    cluster%pselect = cluster%pselect/nprocs
                   ! If the excitor's population is below the initiator threshold, we remove the
                   ! initiator status for the cluster
                end if
            end if
         !   end if
        !end do
        end do
        if (cluster%nexcitors == 0) then
            call create_null_cluster(sys, f0, cluster%pselect, normalisation*cluster_population, initiator_pop, &
                                    cdet, cluster, excit_gen_data)
        else 

            if (allowed) then
                cluster%excitation_level = get_excitation_level(f0(:sys%basis%bit_string_len), cdet%f(:sys%basis%bit_string_len))
                ! To contribute the cluster must be within a double excitation of
                ! the maximum excitation included in the CC wavefunction.
                allowed = cluster%excitation_level <= ex_level+2
            end if

            if (allowed) then

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
        end if
        !print*, 'final cdet', cdet%f(1),allowed

       ![todo] write logging to account for missing parameters
       !if (debug) call write_logging_stoch_selection(logging_info, cluster%nexcitors, cluster%excitation_level, pop, &
        !        max_size, cluster%pselect, cluster%amplitude, allowed)

    end subroutine select_trot_ucc_cluster

    subroutine direct_annihilation_trot(sys, rng, reference, annihilation_flags, psip_list, spawn, &
                                   nspawn_events, determ)

        ! Annihilation algorithm.
! [review] - AJWT: Based on direct_annihilation routine.
        ! Spawned walkers are added to the main list, by which new walkers are
        ! introduced to the main list and existing walkers can have their
        ! populations either enhanced or diminished.

        ! This is a wrapper around various utility functions which perform the
        ! different parts of the annihilation process.

        ! In:
        !    sys: system being studied.
        !    reference: current reference determinant.
        !    annihilation_flags: calculation specific annihilation flags.
        ! In/Out:
        !    rng: random number generator.
        !    psip_list: particle_t object containing psip information after the
        !       death step for the current iteration; combined with the spawned
        !       particles on exit.
        !    spawn: spawn_t object containing the set of spawned particles.
        !    determ (optional): Derived type containing information on the
        !       semi-stochastic part of the simulation.
        ! Out:
        !    nspawn_events (optional): number of successful spawning events on
        !       the processor.

        use parallel, only: iproc
        use spawn_data, only: spawn_t, annihilate_wrapper_spawn_t, calc_events_spawn_t, memcheck_spawn_t
        use system, only: sys_t
        use dSFMT_interface, only: dSFMT_t
        use qmc_data, only: semi_stoch_t, particle_t, annihilation_flags_t, semi_stoch_separate_annihilation
        use reference_determinant, only: reference_t
        

        type(sys_t), intent(in) :: sys
        type(dSFMT_t), intent(inout) :: rng
        type(reference_t), intent(in) :: reference
        type(annihilation_flags_t), intent(in) :: annihilation_flags
        type(particle_t), intent(inout) :: psip_list
        type(spawn_t), intent(inout) :: spawn
        integer, optional, intent(out) :: nspawn_events
        type(semi_stoch_t), intent(inout), optional :: determ
        logical :: doing_semi_stoch

        doing_semi_stoch = .false.
        if (present(determ)) doing_semi_stoch = determ%doing_semi_stoch
        if (present(nspawn_events)) nspawn_events = calc_events_spawn_t(spawn)

        call memcheck_spawn_t(spawn, dont_warn=spawn%warned)

        call annihilate_wrapper_spawn_t_single_trot(spawn, annihilation_flags%initiator_approx)
        call annihilate_main_list_wrapper_trot(sys, rng, reference, annihilation_flags, psip_list, spawn)

    end subroutine direct_annihilation_trot

    subroutine annihilate_wrapper_spawn_t_single_trot(spawn, tinitiator, determ_size)

! [review] - AJWT: based on annihilate_wrapper_spawn_t_single - only difference is a different sort routine_called
        ! Helper procedure for performing annihilation within a spawn_t object.

        ! In:
        !    tinitiator: true if the initiator approximation is being used.
        !    determ_size (optional): The size of the deterministic space in
        !       use, on this process. If input then the deterministic states
        !       received from the various processes will be combined in a
        !       separate call to compress_determ_repeats.
        ! In/Out:
        !    spawn: spawn_t object containing spawned particles.  On output, the
        !        spawned particles are sent to the processor which 'owns' the
        !        determinant they are located on and annihilation is performed
        !        internally, so each determinant appears (at most) once in the
        !        spawn%sdata array.

        use parallel, only: nthreads, nprocs
        use basis_types, only: basis_t
        use spawn_data, only: compress_threaded_spawn_t, comm_spawn_t, compress_determ_repeats, annihilate_spawn_t, &
                        annihilate_spawn_t_initiator, spawn_t
        use sort, only: qsort

        type(spawn_t), intent(inout) :: spawn
        logical, intent(in) :: tinitiator
        integer, intent(in), optional :: determ_size

        integer :: nstates_received(0:nprocs-1)
        integer, parameter :: thread_id = 0
        integer :: i

        ! Compress the successful spawning events from each thread so the
        ! spawned list being sent to each processor contains no gaps.
        if (nthreads > 1) call compress_threaded_spawn_t(spawn)
        if (nprocs > 1) then
            ! Send spawned walkers to the processor which "owns" them and
            ! receive the walkers "owned" by this processor.
            call comm_spawn_t(spawn, nstates_received)

            ! Compress the repeats of the various deterministic states, each of
            ! which is received once from each process.
            if (present(determ_size)) call compress_determ_repeats(spawn, nstates_received, determ_size)
        end if
        if (spawn%head(thread_id,0) > 0) then
            ! Have spawned walkers on this processor.
            call qsort_i0_list_trot(spawn%sdata, spawn%head(thread_id,0), spawn%bit_str_len)
            !call qsort(spawn%sdata, spawn%head(thread_id,0), spawn%bit_str_len)
            ! Annihilate within spawned walkers list.
            ! Compress the remaining spawned walkers list.

            if (tinitiator) then
                call annihilate_spawn_t_initiator(spawn)
            else
                call annihilate_spawn_t(spawn)
            end if
        end if

    end subroutine annihilate_wrapper_spawn_t_single_trot

    subroutine annihilate_main_list_wrapper_trot(sys, rng, reference, annihilation_flags, psip_list, spawn, &
                                            lower_bound, determ_flags)
! [review] - AJWT: based on annihilate_main_list_wrapper
        ! This is a wrapper around various utility functions which perform the
        ! different parts of the annihilation process during non-blocking
        ! communications.

        ! In:
        !    sys: system being studied.
        !    reference: current reference determinant.
        !    annihilation_flags: calculation specific annihilation flags.
        ! In/Out:
        !    rng: random number generator.
        !    psip_list: particle_t object containing psip information after the
        !       death step for the current iteration; combined with the spawned
        !       particles on exit.
        !    spawn: spawn_t object containing spawned particles. For non-blocking
        !       communications a subsection of the spawned walker list will be annihilated
        !       with the main list, otherwise the entire list will be annihilated and merged.
        !    determ_flags (optional): A list of flags specifying whether determinants in
        !        the corresponding particle_t%states are deterministic or not.
        ! In (optional):
        !     lower_bound: starting point we annihiliate from in spawn_t object.

        use system, only: sys_t
        use spawn_data, only: spawn_t
        use dSFMT_interface, only: dSFMT_t
        use qmc_data, only: particle_t, annihilation_flags_t
        use reference_determinant, only: reference_t
        use annihilation, only: remove_unoccupied_dets, round_low_population_spawns

        type(sys_t), intent(in) :: sys
        type(dSFMT_t), intent(inout) :: rng
        type(reference_t), intent(in) :: reference
        type(annihilation_flags_t), intent(in) :: annihilation_flags
        integer, optional, intent(in) :: lower_bound
        type(particle_t), intent(inout) :: psip_list
        type(spawn_t), intent(inout) :: spawn
        integer, intent(inout), optional :: determ_flags(:)

        integer, parameter :: thread_id = 0
        integer :: spawn_start

        if (present(lower_bound)) then
            spawn_start = lower_bound
        else
            spawn_start = 1
        end if

        if (spawn%head(thread_id,0) >= spawn_start) then
            ! Have spawned walkers on this processor.

            call annihilate_main_list_trot(psip_list, spawn, sys, reference, lower_bound)

            ! Remove determinants with zero walkers on them from the main
            ! walker list.
            call remove_unoccupied_dets(rng, psip_list, annihilation_flags%real_amplitudes, determ_flags)

            ! Remove low-population spawned walkers by stochastically
            ! rounding their population up to one or down to zero.
            if (annihilation_flags%real_amplitudes) then
                call round_low_population_spawns(rng, psip_list%pop_real_factor, spawn, lower_bound)
            end if

            ! Insert new walkers into main walker list.
            call insert_new_walkers_trot(sys, psip_list, reference, annihilation_flags, spawn, determ_flags, lower_bound)

            call remove_unoccupied_dets(rng, psip_list, annihilation_flags%real_amplitudes, determ_flags)
        else

            ! No spawned walkers so we only have to check to see if death has
            ! killed the entire population on a determinant.
            call remove_unoccupied_dets(rng, psip_list, annihilation_flags%real_amplitudes, determ_flags)

        end if

    end subroutine annihilate_main_list_wrapper_trot

    subroutine annihilate_main_list_trot(psip_list, spawn, sys, ref, lower_bound)
! [review] - AJWT: based on annihilate_main_list
        ! Annihilate particles in the main walker list with those in the spawned
        ! walker list.

        ! In:
![review] -  AJWT: why not pass in tensor_label_len?
        !    sys: system being studied.
        !    ref: current reference determinant.
        ! In/Out:
        !    psip_list: particle_t object containing psip information.
        !       On exit the particles spawned onto occupied sites have been
        !       annihilated with the existing populations.
        !    spawn: spawn_t obeject containing spawned particles to be annihilated with main
        !       list.  On exit contains particles spawned onto non-occupied
        !       sites.
        ! In (optional):
        !    lower_bound: starting point we annihiliate from in spawn_t object.
        !       Default: 1.

        use search, only: binary_search
        use spawn_data, only: spawn_t
        use qmc_data, only: particle_t
        use system, only: sys_t
        use reference_determinant, only: reference_t

        type(particle_t), intent(inout) :: psip_list
        type(spawn_t), intent(inout) :: spawn
        type(sys_t), intent(in) :: sys
        type(reference_t), intent(in) :: ref
        integer, intent(in), optional :: lower_bound

        integer :: i, pos, k, istart, iend, nannihilate, spawn_start
        integer(int_p) :: old_pop(psip_list%nspaces)
        integer(i0) :: f(sys%basis%tensor_label_len)

        logical :: hit
        integer, parameter :: thread_id = 0

        nannihilate = 0
        if (present(lower_bound)) then
            spawn_start = lower_bound
        else
            spawn_start = 1
        end if
        istart = 1
        iend = psip_list%nstates

        do i = spawn_start, spawn%head(thread_id,0)
            f = int(spawn%sdata(:sys%basis%tensor_label_len,i), i0)
            call binary_search_i0_list_trot(psip_list%states, f, istart, iend, hit, pos)
            !call binary_search(psip_list%states, f, istart, iend, hit, pos)
            if (hit) then
                ! Annihilate!
                old_pop = psip_list%pops(:,pos)
                psip_list%pops(:,pos) = psip_list%pops(:,pos) + &
                    int(spawn%sdata(spawn%bit_str_len+1:spawn%bit_str_len+spawn%ntypes,i), int_p)
                nannihilate = nannihilate + 1
                ! The change in the number of particles is a bit subtle.
                ! We need to take into account:
                !   i) annihilation enhancing the population on a determinant.
                !  ii) annihilation diminishing the population on a determinant.
                ! iii) annihilation changing the sign of the population (i.e.
                !      killing the population and then some).
                associate(nparticles=>psip_list%nparticles)
                    nparticles = nparticles + real(abs(psip_list%pops(:,pos)) - abs(old_pop),p)/psip_list%pop_real_factor
                end associate
                ! Next spawned walker cannot annihilate any determinant prior to
                ! this one as the lists are sorted.
                istart = pos + 1
            else
                ! Compress spawned list.
                k = i - nannihilate
                spawn%sdata(:,k) = spawn%sdata(:,i)
            end if
        end do

        spawn%head(thread_id,0) = spawn%head(thread_id,0) - nannihilate

    end subroutine annihilate_main_list_trot

    subroutine insert_new_walkers_trot(sys, psip_list, ref, annihilation_flags, spawn, determ_flags, lower_bound)
! [review] - AJWT: based on insert_new_walkers - just with a different search.
        ! Insert new walkers into the main walker list from the spawned list.
        ! This is done after all particles have been annihilated, so the spawned
        ! list contains only new walkers.

        ! In:
        !    sys: system being studied.
        !    ref: reference determinant --- the diagonal matrix elements are required.
        !    annihilation_flags: calculation specific annihilation flags.
        ! In/Out:
        !    psip_list: psip information.
        !    spawn: spawn_t object containing list of particles spawned onto
        !        previously occupied sites.
        ! In (optional):
        !    determ_flags: A list of flags specifying whether determinants in
        !        psip_list%states are deterministic or not.
        !    lower_bound: starting point we annihiliate from in spawn_t object.

        use errors, only: stop_all
        use parallel, only: iproc
        use qmc_data, only: particle_t, annihilation_flags_t
        use reference_determinant, only: reference_t
        use search, only: binary_search
        use spawn_data, only: spawn_t
        use system, only: sys_t
        use utils, only: int_fmt
        use annihilation, only: insert_new_walker

        use, intrinsic :: iso_fortran_env, only: error_unit

        type(sys_t), intent(in) :: sys
        type(particle_t), intent(inout) :: psip_list
        type(reference_t), intent(in) :: ref
        type(annihilation_flags_t), intent(in) :: annihilation_flags
        type(spawn_t), intent(inout) :: spawn
        integer, intent(inout), optional :: determ_flags(:)
        integer, intent(in), optional :: lower_bound

        integer :: i, istart, iend, j, k, pos, spawn_start, disp
        real(p) :: real_population(psip_list%nspaces)

        logical :: hit
        integer, parameter :: thread_id = 0
        real :: fill_fraction

        ! Merge new walkers into the main list.

        ! Both the main list and the spawned list are sorted: the spawned list
        ! is sorted explicitly and the main list is sorted by construction via
        ! merging.

        ! 1. Find the position where the spawned walker should go.
        ! 2. Move all walkers above it to create the vacant slot for the new
        ! walker.  As we know how many elements we are inserting, we only need
        ! move a given walker at most once.
        ! 3. Insert the new walker at the bottom of the shifted block so it
        ! doesn't have to be moved again to accommodate other new walkers.

        ! We can make the search faster as we iterate through the spawned
        ! walkers in descending order, so once we know where one walker goes, we
        ! know that the next new walker has to go below it, allowing us to
        ! search through an ever-decreasing number of elements.

        if (present(lower_bound)) then
            spawn_start = lower_bound
            disp = lower_bound - 1
        else
            spawn_start = 1
            disp = 0
        end if

        ! Don't bother to perform these checks and print more error messages
        ! if we've run out of memory already.
        if (.not. psip_list%error) then
            fill_fraction = real(psip_list%nstates+(spawn%head(thread_id,0)-spawn_start+1))/size(psip_list%states,2)
            if (fill_fraction > 1.00) then
                write (error_unit,'(1X,"# Error: No space left in main particle array on processor",'//int_fmt(iproc,1)//',".")') &
                              iproc
                write (error_unit,'(1X,"# Error: HANDE will exit at the end of this report loop.")')
                write (error_unit,'(1X,"# Error: Note that spawning until the end of the report loop will be affected and&
                              & so results from this final loop may be slightly incorrect.")')
                write (error_unit,'(1X,"# Error: Some reconvergence time should be allowed if continuing from a&
                              & subsequent restart file.")')

                psip_list%error = .true.
            else if (fill_fraction > 0.95) then
                if (psip_list%warn) then
                    write (error_unit,'(1X,"# Warning: filled over 95% of main particle array on processor",'//int_fmt(iproc,1)&
                              //',".")') iproc
                    write (error_unit,'(1x,"This warning only prints once")')
                    psip_list%warn = .false.
                end if
                psip_list%warning_count = psip_list%warning_count + 1
            end if
        end if

        if (.not. psip_list%error) then
            istart = 1
            iend = psip_list%nstates
            do i = spawn%head(thread_id,0), spawn_start, -1

                ! spawned det is not in the main walker list.
                call binary_search_i0_list_trot(psip_list%states, int(spawn%sdata(:sys%basis%tensor_label_len,i), i0), &
                !call binary_search(psip_list%states, int(spawn%sdata(:sys%basis%tensor_label_len,i), i0), &
                                   istart, iend, hit, pos)
                ! f should be in slot pos.  Move all determinants above it.
                do j = iend, pos, -1
                    ! i is the number of determinants that will be inserted below j.
                    k = j + i - disp
                    psip_list%states(:,k) = psip_list%states(:,j)
                    psip_list%pops(:,k) = psip_list%pops(:,j)
                    psip_list%dat(:,k) = psip_list%dat(:,j)
                    if (present(determ_flags)) determ_flags(k) = determ_flags(j)
                end do

                ! Insert new walker into pos and shift it to accommodate the number
                ! of elements that are still to be inserted below it.
                k = pos + i - 1 - disp
                ! The encoded spawned walker sign.
                associate(spawned_population => spawn%sdata(spawn%bit_str_len+1:spawn%bit_str_len+spawn%ntypes, i), &
                        tbl=>sys%basis%tensor_label_len)
                    call insert_new_walker(sys, psip_list, annihilation_flags, k, int(spawn%sdata(:tbl,i), i0), &
                                           int(spawned_population, int_p), ref)
                    ! Extract the real sign from the encoded sign.
                    real_population = real(spawned_population,p)/psip_list%pop_real_factor
                    psip_list%nparticles = psip_list%nparticles + abs(real_population)
                end associate

                ! A deterministic state can never leave the main list so cannot be
                ! in the spawned list at this point. So set the flag to specify
                ! that this state is not deterministic.
                if (present(determ_flags)) determ_flags(k) = 1

                ! Next walker will be inserted below this one.
                iend = pos - 1
            end do

            ! Update psip_list%nstates.
            psip_list%nstates = psip_list%nstates + spawn%head(thread_id,0) - disp
        end if

    end subroutine insert_new_walkers_trot

    pure subroutine qsort_i0_list_trot(list, head, nsort)
! [review] - AJWT: based on sort routines in sort.f90
! [review] - AJWT: Should probably be in sort.f90 and called qsort_i0_list_rev

        ! Sort a 2D array of int_64 integers.

        ! list(:,i) is regarded as greater than list(:,j) if the first
        ! non-identical element between list(:,i) and list(:,j) is lower in
        ! list(:,i).

        ! In/Out:
        !    list: 2D array of int_64 integers.  Sorted on output.
        ! In:
        !    head (optional): sort list up to and including list(:,:head) and
        !        leave the rest of the array untouched.  Default: sort the
        !        entire array.
        !    nsort (optional): sort list only using the first nsort elements in
        !        each 1D slice to compare entries (ie compare list(:nsort,i) and
        !        list(:nsort,j), so list is sorted according to list(:nsort,:)).
        !        Default: use entire slice.

        use bit_utils, only: operator(.bitstrge.), operator(.bitstrgt.)

        integer(i0), intent(inout) :: list(:,:)
        integer, intent(in), optional :: head, nsort

        ! Threshold.  When a sublist gets to this length, switch to using
        ! insertion sort to sort the sublist.
        integer, parameter :: switch_threshold = 7

        ! sort needs auxiliary storage of length 2*log_2(n).
        integer, parameter :: stack_max = 50

        integer :: pivot, lo, hi, i, j, ns
        integer(int_64) :: tmp(ubound(list,dim=1))

        ! Stack.  This is the auxilliary memory required by quicksort.
        integer :: stack(2,stack_max), nstack

        if (present(nsort)) then
            ns = nsort
        else
            ns = ubound(list, dim=1)
        end if

        nstack = 0
        lo = 1
        if (present(head)) then
            hi = head
        else
            hi = ubound(list, dim=2)
        end if
        do
            ! If the section/partition we are looking at is smaller than
            ! switch_threshold then perform an insertion sort.
            if (hi - lo < switch_threshold) then
                do j = lo + 1, hi
                    tmp = list(:,j)
                    do i = j - 1, 1, -1
                        if (.not.(tmp(1:ns) .bitstrgt. list(1:ns,i))) exit
                        list(:,i+1) = list(:,i)
                    end do
                    list(:,i+1) = tmp
                end do

                if (nstack == 0) exit
                hi = stack(2,nstack)
                lo = stack(1,nstack)
                nstack = nstack - 1

            else
                ! Otherwise start partitioning with quicksort.

                ! Pick the pivot element to be the median of list(:,lo), list(:,hi)
                ! and list(:,(lo+hi)/2).
                ! This largely overcomes a major problem with quicksort, where it
                ! degrades if the pivot is always the smallest element.
                pivot = (lo + hi)/2
                call swap_sublist(list(:,pivot), list(:,lo + 1))
                if (.not.(list(1:ns,lo) .bitstrge. list(1:ns,hi))) then
                    call swap_sublist(list(:,lo), list(:,hi))
                end if
                if (.not.(list(1:ns,lo+1) .bitstrge. list(1:ns,hi))) then
                    call swap_sublist(list(:,lo+1), list(:,hi))
                end if
                if (.not.(list(1:ns,lo) .bitstrge. list(1:ns,lo+1))) then
                    call swap_sublist(list(:,lo), list(:,lo+1))
                end if

                i = lo + 1
                j = hi
                tmp = list(:,lo + 1) ! a is the pivot value
                do while (.true.)
                    ! Scan down list to find element > a.
                    i = i + 1
                    do while (.not.(tmp(1:ns) .bitstrge. list(1:ns,i)))
                        i = i + 1
                    end do

                    ! Scan down list to find element < a.
                    j = j - 1
                    do while (.not.(list(1:ns,j) .bitstrge. tmp(1:ns)))
                        j = j - 1
                    end do

                    ! When the pointers crossed, partitioning is complete.
                    if (j < i) exit

                    ! Swap the elements, so that all elements < a end up
                    ! in lower indexed variables.
                    call swap_sublist(list(:,i), list(:,j))
                end do

                ! Insert partitioning element
                list(:,lo + 1) = list(:,j)
                list(:,j) = tmp

                ! Push the larger of the partitioned sections onto the stack
                ! of sections to look at later.
                ! --> need fewest stack elements.
                nstack = nstack + 1

                ! With a stack_max of 50, we can sort arrays of length
                ! 1125899906842624.  It is safe to say this will never be
                ! exceeded, and so this test can be skipped.
!                if (nstack > stack_max) call stop_all('qsort_int_64_list', "parameter stack_max too small")

                if (hi - i + 1 >= j - lo) then
                    stack(2,nstack) = hi
                    stack(1,nstack) = i
                    hi = j - 1
                else
                    stack(2,nstack) = j - 1
                    stack(1,nstack) = lo
                    lo = i
                end if

            end if
        end do

    contains

        pure subroutine swap_sublist(s1,s2)

            integer(int_64), intent(inout) :: s1(:), s2(:)
            integer(int_64) :: tmp(ubound(s1,dim=1))

            tmp = s1
            s1 = s2
            s2 = tmp

        end subroutine swap_sublist

    end subroutine qsort_i0_list_trot

    subroutine regenerate_trot_info_psip_list(basis, nel, qs)

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
        !   qmc_state: information on current state of calculation. We update and
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
             call qsort_psip_info_trot(pl%nstates, pl%states, pl%pops, pl%dat)
        end associate

    end subroutine regenerate_trot_info_psip_list

    pure function bit_str_cmp_trot(b1, b2) result(cmp)

        ! In:
        !    b1(:), b2(:): bit string.
        ! Returns:
        !    0 if b1 and b2 are identical;
        !    1 if the most significant non-identical element in b1 is bitwise
        !      greater than the corresponding element in b2;
        !    -1 if the most significant non-identical element in b1 is bitwise
        !      less than the corresponding element in b2;

        integer :: cmp
        integer(i0), intent(in) :: b1(:), b2(:)

        integer :: i

        cmp = 0
        do i = ubound(b1, dim=1), 1, -1
            if (blt(b1(i),b2(i))) then
                cmp = -1
                exit
            else if (bgt(b1(i),b2(i))) then
                cmp = 1
                exit
            end if
        end do

    end function bit_str_cmp_trot

    pure subroutine binary_search_i0_list_trot(list, item, istart, iend, hit, pos)
! [review] - AJWT: based on binary_search_i0_list with a different comparator operator.

        ! Find where an item resides in a list of such items.
        ! Only elements between istart and iend are examined (use the
        ! array boundaries in the worst case).
        !
        ! In:
        !    list: a sorted i0 integer 2D list/array; the first dimension
        !        corresponds to 1D arrays to compare to item.
        !    item: an i0 integer 1D list/array.
        !    istart: first position to examine in the list.
        !    iend: last position to examine in the list.
        ! Out:
        !    hit: true if found item in list.
        !    pos: the position corresponding to item in list.
        !        If hit is true, then the element in this position is the same
        !        as item, else this is where item should go to keep the list
        !        sorted.

        use const, only: i0

        integer(i0), intent(in) :: list(:,:), item(:)
        integer, intent(in) :: istart, iend
        logical, intent(out) :: hit
        integer, intent(out) :: pos

        integer :: hi, lo, compare

        if (istart > iend) then

            ! Already know the element has to be appended to the list.
            ! This should only occur if istart = iend + 1.
            pos = istart
            hit = .false.

        else

            ! Search range.
            lo = istart
            hi = iend

            ! Assume item doesn't exist in the list initially.
            hit = .false.

            do while (hi /= lo)
                ! Narrow the search range down in steps.

                ! Mid-point.
                ! We shift one of the search limits to be the mid-point.
                ! The successive dividing the search range by 2 gives a O[log N]
                ! search algorithm.
                pos = (hi+lo)/2

                compare = bit_str_cmp_trot(list(:,pos), item)
                select case(compare)
                case (0)
                    ! hit!
                    hit = .true.
                    exit
                case(1)
                    ! list(:,pos) is "smaller" than item.
                    ! The lowest position item can take is hence pos + 1 (i.e. if
                    ! item is greater than pos by smaller than pos + 1).
                    lo = pos + 1
                case(-1)
                    ! list(:,pos) is "greater" than item.
                    ! The highest position item can take is hence pos (i.e. if item is
                    ! smaller than pos but greater than pos - 1).  This is why
                    ! we differ slightly from a standard binary search (where lo
                    ! is set to be pos+1 and hi to be pos-1 accordingly), as
                    ! a standard binary search assumes that the element you are
                    ! searching for actually appears in the array being
                    ! searched...
                    hi = pos
                end select

            end do

            ! If hi == lo, then we have narrowed the search down to one position but
            ! not checked if that position is the item we're hunting for.
            ! Because list can expand (i.e. we might be searching for an
            ! element which doesn't exist yet) the binary search can find either
            ! the element before or after where item should be placed.
            if (hi == lo) then
                compare = bit_str_cmp_trot(list(:,hi), item)
                select case(compare)
                case (0)
                    ! hit!
                    hit = .true.
                    pos = hi
                case(1)
                    ! list(:,pos) is "smaller" than item.
                    ! item should be placed in the next slot.
                    pos = hi + 1
                case(-1)
                    ! list(:,pos) is "greater" than item.
                    ! item should ber placed here.
                    pos = hi
                end select
            end if

        end if

    end subroutine binary_search_i0_list_trot

    pure subroutine qsort_psip_info_trot(nstates, states, pops, dat)
! [review] - AJWT: based on qsort_psip_info with alternative ordering.
        ! Sort a set of psip information (states, populations and data) in order according
        ! to the state labels.

        ! states(:,i) is regarded as greater than states(:,j) if the first
        ! non-identical element between states(:,i) and states(:,j) is smaller in
        ! states(:,i).

        ! In:
        !    nstates: number of occupied states.
        ! In/Out:
        !    states: 2D array of i0 integers containing the state label for each occupied state.
        !        Sorted on output.
        !    pops, dat: population and data arrays for each state.  Sorted by states on output.

        ! Note: the size of the first dimension of states is immaterial, as we do a comparison
        ! based on the entire slice.  The second dimensions of pops, dat and states must be >= nstates.

        use const, only: int_p, i0, p
        use bit_utils, only: operator(.bitstrge.), operator(.bitstrgt.)

        integer, intent(in) :: nstates
        integer(i0), intent(inout) :: states(:,:)
        integer(int_p), intent(inout) :: pops(:,:)
        real(p), intent(inout) :: dat(:,:)

        ! Threshold.  When a substates gets to this length, switch to using
        ! insertion sort to sort the substates.
        integer, parameter :: switch_threshold = 7

        ! sort needs auxiliary storage of length 2*log_2(n).
        integer, parameter :: stack_max = 50

        integer :: pivot, lo, hi, i, j
        integer(i0) :: tmp_state(ubound(states,dim=1))
        integer(int_p) :: tmp_pop(ubound(pops,dim=1))
        real(p) :: tmp_dat(ubound(dat,dim=1))

        ! Stack.  This is the auxilliary memory required by quicksort.
        integer :: stack(2,stack_max), nstack

        hi = nstates

        nstack = 0
        lo = 1

        do
            ! If the section/partition we are looking at is smaller than
            ! switch_threshold then perform an insertion sort.
            if (hi - lo < switch_threshold) then
                do j = lo + 1, hi
                    tmp_state = states(:,j)
                    tmp_pop = pops(:,j)
                    tmp_dat = dat(:,j)
                    do i = j - 1, 1, -1
                        if (.not.(tmp_state .bitstrgt. states(:,i))) exit
                        states(:,i+1) = states(:,i)
                        pops(:,i+1) = pops(:,i)
                        dat(:,i+1) = dat(:,i)
                    end do
                    states(:,i+1) = tmp_state
                    pops(:,i+1) = tmp_pop
                    dat(:,i+1) = tmp_dat
                end do

                if (nstack == 0) exit
                hi = stack(2,nstack)
                lo = stack(1,nstack)
                nstack = nstack - 1

            else
                ! Otherwise start partitioning with quicksort.

                ! Pick the pivot element to be the median of states(:,lo), states(:,hi)
                ! and states(:,(lo+hi)/2).
                ! This largely overcomes a major problem with quicksort, where it
                ! degrades if the pivot is always the smallest element.
                pivot = (lo + hi)/2
                call swap_states(states(:,pivot), pops(:,pivot), dat(:,pivot), states(:,lo+1), pops(:,lo+1), dat(:,lo+1))
                if (.not.(states(:,lo) .bitstrge. states(:,hi))) then
                    call swap_states(states(:,lo), pops(:,lo), dat(:,lo), states(:,hi), pops(:,hi), dat(:,hi))
                end if
                if (.not.(states(:,lo+1) .bitstrge. states(:,hi))) then
                    call swap_states(states(:,lo+1), pops(:,lo+1), dat(:,lo+1), states(:,hi), pops(:,hi), dat(:,hi))
                end if
                if (.not.(states(:,lo) .bitstrge. states(:,lo+1))) then
                    call swap_states(states(:,lo), pops(:,lo), dat(:,lo), states(:,lo+1), pops(:,lo+1), dat(:,lo+1))
                end if

                i = lo + 1
                j = hi
                tmp_state = states(:,lo+1) ! a is the pivot value
                tmp_pop = pops(:,lo+1)
                tmp_dat = dat(:,lo+1)
                do while (.true.)
                    ! Scan down states to find element > a.
                    i = i + 1
                    do while (.not.(tmp_state .bitstrge. states(:,i)))
                        i = i + 1
                    end do

                    ! Scan down states to find element < a.
                    j = j - 1
                    do while (.not.(states(:,j) .bitstrge. tmp_state))
                        j = j - 1
                    end do

                    ! When the pointers crossed, partitioning is complete.
                    if (j < i) exit

                    ! Swap the elements, so that all elements < a end up
                    ! in lower indexed variables.
                    call swap_states(states(:,i), pops(:,i), dat(:,i), states(:,j), pops(:,j), dat(:,j))
                end do

                ! Insert partitioning element
                states(:,lo + 1) = states(:,j)
                pops(:,lo + 1) = pops(:,j)
                dat(:,lo + 1) = dat(:,j)
                states(:,j) = tmp_state
                pops(:,j) = tmp_pop
                dat(:,j) = tmp_dat

                ! Push the larger of the partitioned sections onto the stack
                ! of sections to look at later.
                ! --> need fewest stack elements.
                nstack = nstack + 1

                ! With a stack_max of 50, we can sort arrays of length
                ! 1125899906842624.  It is safe to say this will never be
                ! exceeded, and so this test can be skipped.
!                if (nstack > stack_max) call stop_all('qsort_int_64_states', "parameter stack_max too small")

                if (hi - i + 1 >= j - lo) then
                    stack(2,nstack) = hi
                    stack(1,nstack) = i
                    hi = j - 1
                else
                    stack(2,nstack) = j - 1
                    stack(1,nstack) = lo
                    lo = i
                end if

            end if
        end do

    contains

        pure subroutine swap_states(s1,p1,d1,s2,p2,d2)

            integer(i0), intent(inout) :: s1(:), s2(:)
            integer(int_p), intent(inout) :: p1(:), p2(:)
            real(p), intent(inout) :: d1(:), d2(:)
            integer(i0) :: tmp_state(ubound(s1,dim=1))
            integer(int_p) :: tmp_pop(ubound(p1,dim=1))
            real(p) :: tmp_dat(ubound(d1,dim=1))

            tmp_state = s1
            s1 = s2
            s2 = tmp_state

            tmp_pop = p1
            p1 = p2
            p2 = tmp_pop

            tmp_dat = d1
            d1 = d2
            d2 = tmp_dat

        end subroutine swap_states

    end subroutine qsort_psip_info_trot

    subroutine find_D0_trot(psip_list, f0, D0_pos)

        ! Find the reference determinant in the list of walkers

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
            call binary_search_i0_list_trot(psip_list%states, f0, 1, psip_list%nstates, hit, D0_pos)
            !call binary_search(psip_list%states, f0, 1, psip_list%nstates, hit, D0_pos)
        else
            D0_pos_old = D0_pos
            select case(bit_str_cmp_trot(f0, psip_list%states(:,D0_pos)))
            case(0)
                ! D0 hasn't moved.
                hit = .true.
            case(1)
                ! D0 < psip_list%states(:,D0_pos) -- it has moved to earlier in
                ! the list and the old D0_pos is an upper bound.
                call binary_search_i0_list_trot(psip_list%states, f0, 1, D0_pos_old, hit, D0_pos)
                !call binary_search(psip_list%states, f0, 1, D0_pos_old, hit, D0_pos)
            case(-1)
                ! D0 > psip_list%states(:,D0_pos) -- it has moved to later in
                ! the list and the old D0_pos is a lower bound.
                call binary_search_i0_list_trot(psip_list%states, f0, D0_pos_old, psip_list%nstates, hit, D0_pos)
                !call binary_search(psip_list%states, f0, D0_pos_old, psip_list%nstates, hit, D0_pos)
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

    subroutine get_probabilities(pops, nstates, D0_pop, D0_pos, real_factor, probs, pref, pavg, pcompavg)


    integer(i0), intent(in) :: pops(:,:), real_factor
    real(p), intent(in) :: D0_pop
    real(p), intent(out) :: probs(:), pref, pavg, pcompavg
    integer, intent(in) :: nstates, D0_pos

    integer :: i
    real(p) :: prob, pop

    pref = 1
    pavg = 1
    pcompavg = 1

!    print*, pops(:, :nstates)
    if (D0_pos == 1) then
        probs(1) = 0
        do i = 2, nstates
            pop = abs(real(pops(1,i),p)/real(real_factor,p))
!            print*, 'pop', pop
            prob = pop/(pop + D0_pop)
!            print*, 'prob', prob
            probs(i) = prob/(1-prob)
            pref = pref * (1 - prob)
!            print*, 'pref', pref
            pavg = pavg * prob  
!            print*, 'pavg', pavg
        end do
    else
        do i = 1, nstates
            if (i == D0_pos) then
                probs(i) = 0
            else
                pop = abs(real(pops(1,i),p)/real(real_factor,p))
                prob = pop/(pop + D0_pop)
                probs(i) = prob/(1.0-prob)
                pref = pref * (1.0 - prob)
                pavg = pavg * prob  
            end if
        end do
    end if
    
    if (nstates > 1) then
        pavg = pavg ** (1.0/real(nstates-1,p))
 !       print*, 'pavg', pavg
        pcompavg = pref ** (1.0/real(nstates - 1,p))
  !      print*, 'pcavg', pcompavg
    else
       pref = 1
       pavg = 0
       pcompavg = 1
    end if

    end subroutine

    function binomial_coeff(a,b) result (n)

    integer, intent(in) :: a, b
    real(dp) :: n
    integer :: i

    n = 1.0_dp
    do i = 1, b
        n = n / real(i,dp)
        n = n * real(a - i+1, dp)
    end do
    end function
    
    function binomial_sum(prob, pcomp, max_size, nstates) result(binsum)

    integer, intent(in) :: max_size, nstates
    real(p), intent(in) :: prob, pcomp
    real(dp) :: binsum
    integer :: i

    binsum = 0.0_dp
    do i = 0, max_size 
!       print*, i, nstates, prob, pcomp, binomial_coeff(nstates, i)
       binsum = binsum + (prob**i)*(pcomp**(nstates-i))*binomial_coeff(nstates,i)
!       print*, binsum
    end do
    end function

    function binomial_sum_repeat(prob, pcomp, max_size, nstates) result(binsum)

    integer, intent(in) :: max_size, nstates
    real(p), intent(in) :: prob, pcomp
    real(dp) :: binsum
    integer :: i

    binsum = 0.0_dp
    do i = 0, max_size 
!       print*, i, nstates, prob, pcomp, binomial_coeff(nstates, i)
       binsum = binsum + (prob**i)*(pcomp**(nstates-i))*binomial_coeff(nstates,i)/prob_no_repeat(nstates, i)
!       print*, binsum
    end do
    end function

    function prob_no_repeat(n,i) result(prob)

    integer, intent(in) :: n, i
    real(dp) :: prob
    integer :: j 

    prob = 1.0_dp
    do j = 1, i-1
        prob = prob * real(n-j,dp)/real(n, dp)
    end do
    end function
    function exp_sum(prob, max_size) result(expsum)

    use utils, only: factorial
    integer, intent(in) :: max_size 
    real(p), intent(in) :: prob
    real(dp) :: expsum
    integer :: i

    expsum = 0.0_dp
    do i = 0, max_size 
!       print*, i, nstates, prob, pcomp, binomial_coeff(nstates, i)
       expsum = expsum + prob**i/factorial(i)
!       print*, binsum
    end do
    end function
   subroutine trot_collapse_cluster(basis, f0, excitor, excitor_population, cluster_excitor, cluster_population, allowed, conjugate)

        ! Collapse two excitors.  The result is returned in-place.

        ! In:
        !    basis: information about the single-particle basis.
        !    f0: bit string representation of the reference determinant.
        !    excitor: bit string of the Slater determinant formed by applying
        !        the excitor, e1, to the reference determinant.
        !    excitor_population: number of excips on the excitor e1.
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

        type(basis_t), intent(in) :: basis
        integer(i0), intent(in) :: f0(basis%tot_string_len)
        integer(i0), intent(in) :: excitor(basis%tot_string_len)
        complex(p), intent(in) :: excitor_population
        integer(i0), intent(inout) :: cluster_excitor(basis%tot_string_len)
        complex(p), intent(inout) :: cluster_population
        logical,  intent(out) :: allowed
        logical, intent(in) :: conjugate
        integer(i0) :: excitor_loc(basis%tot_string_len)
! [review] - AJWT:  Is this copy needed?
        integer(i0) :: f0_loc(basis%tot_string_len)

        integer :: ibasis, ibit
        integer(i0) :: excitor_excitation(basis%tot_string_len)
        integer(i0) :: excitor_annihilation(basis%tot_string_len)
        integer(i0) :: excitor_creation(basis%tot_string_len)
        integer(i0) :: cluster_excitation(basis%tot_string_len)
        integer(i0) :: cluster_annihilation(basis%tot_string_len)
        integer(i0) :: cluster_creation(basis%tot_string_len)
        integer(i0) :: permute_operators(basis%tot_string_len)

        excitor_loc = excitor
        f0_loc = f0
        call reset_extra_info_bit_string(basis, excitor_loc)
        call reset_extra_info_bit_string(basis, f0_loc)

        ! Apply excitor to the cluster of excitors.

        ! orbitals involved in excitation from reference
        excitor_excitation = ieor(f0_loc, excitor_loc)
        cluster_excitation = ieor(f0_loc, cluster_excitor)
        ! annihilation operators (relative to the reference)
        if (conjugate) then
            excitor_annihilation = excitor_excitation - iand(excitor_excitation, f0_loc)
        else
            excitor_annihilation = iand(excitor_excitation, f0_loc)
        end if
        ! creation operators (relative to the reference)
        if (conjugate) then
            excitor_creation = iand(excitor_excitation, f0_loc)
        else
            excitor_creation = iand(excitor_excitation, excitor_loc)
        end if
        ! annihilation operators (relative to the reference)
        cluster_annihilation = iand(cluster_excitation, f0_loc)
        ! creation operators (relative to the reference)
        cluster_creation = iand(cluster_excitation, cluster_excitor)

        if(conjugate) then
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
                        if (btest(cluster_excitor(ibasis),ibit)) then
                           allowed = .false.
                           exit
                        end if
                    end if
                end do
           end do
        
           if (allowed) then
                cluster_population = cluster_population*excitor_population

                ! Now apply the excitor to the cluster (which is, in its own right,
                ! an excitor).
                ! Consider a cluster, e.g. t_i^a = a^+_a a_i (which corresponds to i->a).
                ! We wish to collapse two excitors, e.g. t_i^a t_j^b, to a single
                ! excitor, e.g. t_{ij}^{ab} = a^+_a a^+_b a_j a_i (where i<j and
                ! a<b).  However t_i^a t_j^b = a^+_a a_i a^+_b a_j.  We thus need to
                ! permute the creation and annihilation operators.  Each permutation
                ! incurs a sign change.
! [review] - AJWT: Even if we can't reuse the whole of collapse_cluster, this part looks teh same, so can that be factored out and called?
                do ibasis = 1, basis%bit_string_len
                    do ibit = 0, i0_end
                        if (btest(excitor_excitation(ibasis),ibit)) then
                            if (.not. btest(f0_loc(ibasis),ibit)) then
                                ! Exciting from this orbital.
                                cluster_excitor(ibasis) = ibclr(cluster_excitor(ibasis),ibit)
                                ! We need to swap it with every annihilation
                                ! operator and every creation operator referring to
                                ! an orbital with a higher index already in the
                                ! cluster.
                                ! Note that an orbital cannot be in the list of
                                ! annihilation operators and the list of creation
                                ! operators.
                                ! First annihilation operators:
                                permute_operators = iand(basis%excit_mask(:,basis%basis_lookup(ibit,ibasis)),cluster_annihilation)
                                ! Now add the creation operators:
                                permute_operators = ior(permute_operators,cluster_creation)
                            else
                                ! Exciting into this orbital.
                                cluster_excitor(ibasis) = ibset(cluster_excitor(ibasis),ibit)
                                ! Need to swap it with every creation operator with
                                ! a lower index already in the cluster.
                                permute_operators = iand(not(basis%excit_mask(:,basis%basis_lookup(ibit,ibasis))),cluster_creation)
                                permute_operators(ibasis) = ibclr(permute_operators(ibasis),ibit)
                            end if
                            if (mod(sum(count_set_bits(permute_operators)),2) == 1) &
                                cluster_population = -cluster_population
                        end if
                    end do
                end do
           end if  
        else
        ! First, let's find out if the excitor is valid...
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

                ! Now apply the excitor to the cluster (which is, in its own right,
                ! an excitor).
                ! Consider a cluster, e.g. t_i^a = a^+_a a_i (which corresponds to i->a).
                ! We wish to collapse two excitors, e.g. t_i^a t_j^b, to a single
                ! excitor, e.g. t_{ij}^{ab} = a^+_a a^+_b a_j a_i (where i<j and
                ! a<b).  However t_i^a t_j^b = a^+_a a_i a^+_b a_j.  We thus need to
                ! permute the creation and annihilation operators.  Each permutation
                ! incurs a sign change.

                do ibasis = 1, basis%bit_string_len
                    do ibit = 0, i0_end
                        if (btest(excitor_excitation(ibasis),ibit)) then
                            if (btest(f0(ibasis),ibit)) then
                                ! Exciting from this orbital.
                                cluster_excitor(ibasis) = ibclr(cluster_excitor(ibasis),ibit)
                                ! We need to swap it with every annihilation
                                ! operator and every creation operator referring to
                                ! an orbital with a higher index already in the
                                ! cluster.
                                ! Note that an orbital cannot be in the list of
                                ! annihilation operators and the list of creation
                                ! operators.
                                ! First annihilation operators:
                                permute_operators = iand(basis%excit_mask(:,basis%basis_lookup(ibit,ibasis)),cluster_annihilation)
                                ! Now add the creation operators:
                                permute_operators = ior(permute_operators,cluster_creation)
                            else
                                ! Exciting into this orbital.
                                cluster_excitor(ibasis) = ibset(cluster_excitor(ibasis),ibit)
                                ! Need to swap it with every creation operator with
                                ! a lower index already in the cluster.
                                permute_operators = iand(not(basis%excit_mask(:,basis%basis_lookup(ibit,ibasis))),cluster_creation)
                                permute_operators(ibasis) = ibclr(permute_operators(ibasis),ibit)
                            end if
                            if (mod(sum(count_set_bits(permute_operators)),2) == 1) &
                                cluster_population = -cluster_population
                        end if
                    end do
                end do

            end if
        end if

    end subroutine trot_collapse_cluster
end module

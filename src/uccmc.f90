module uccmc

use const, only: i0, int_p, int_64, p, dp, debug

implicit none

contains

! [review] - AJWT: Put in a brief summary of the differences between this and conventional CC.
    subroutine do_uccmc(sys, qmc_in, uccmc_in, restart_in, load_bal_in, reference_in, &
                        logging_in, io_unit, qs, qmc_state_restart)

! [review] - AJWT: Perhaps note that this is derived from do_ccmc (and similarly in do_ccmc that this derives from it).  Are there any commonalities which can be sensibly shared as a called function?
        ! Run the UCCMC algorithm starting from the initial walker distribution
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
        use ccmc_data
        use ccmc_selection, only: create_null_cluster
        use ccmc_selection, only: init_selection_data, update_selection_probabilities, &
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
        use spawn_data, only: calc_events_spawn_t, write_memcheck_report
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
 
        real(p), allocatable :: time_avg_psip_list_ci_pops(:), time_avg_psip_list_pops(:),time_avg_psip_list_sq(:,:)
        integer(i0), allocatable :: time_avg_psip_list_ci_states(:,:), time_avg_psip_list_states(:,:)
        integer :: semi_stoch_it, pos, j, k, nstates_ci, nstates_sq
        logical :: hit
        integer(i0), allocatable :: state(:)
        real(p) :: population
        real(p) :: real_population, var_energy

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
            call check_uccmc_opts(sys, uccmc_in, qmc_in)
        end if

        ! Initialise data.
        call init_qmc(sys, qmc_in, restart_in, load_bal_in, reference_in, io_unit, annihilation_flags, qs, &
                      uuid_restart, restart_version_restart, qmc_state_restart=qmc_state_restart, &
                      regenerate_info=regenerate_info, uccmc_in=uccmc_in)

        allocate(state(sys%basis%bit_string_len))

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
            if (allocated(err_msg)) call stop_all('do_uccmc', 'Failed to reset RNG state: '//err_msg)
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
            call init_report_loop(qs, bloom_stats)

            do icycle = 1, qmc_in%ncycles
                
                iter = qs%mc_cycles_done + (ireport-1)*qmc_in%ncycles + icycle

                if(uccmc_in%variational_energy) then
                          time_avg_psip_list_ci_pops(:nstates_ci) =  time_avg_psip_list_ci_pops(:nstates_ci)*(iter-1)
                end if
                time_avg_psip_list_pops(:nstates_sq) =  time_avg_psip_list_pops(:nstates_sq)*(iter - 1)
                time_avg_psip_list_sq(2,:nstates_sq) =  time_avg_psip_list_sq(2,:nstates_sq)*(iter - 1)

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

                call ucc_cumulative_population(qs%psip_list%pops, qs%psip_list%states(sys%basis%tot_string_len,:), &
                                           qs%psip_list%nstates, D0_proc, D0_pos, qs%psip_list%pop_real_factor, &
                                           sys%read_in%comp, cumulative_abs_real_pops, &
                                           tot_abs_real_pop)

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
                !$omp        nattempts_spawn, D0_population_ucc_cycle,&
                !$omp        uccmc_in, nprocs, ms_stats, ps_stats, qmc_in, load_bal_in, &
                !$omp        ndeath_nc,   &
                !$omp        nparticles_change, ndeath, logging_info, time_avg_psip_list_ci)

                it = get_thread_id()
                iexcip_pos = 0
                seen_D0 = .false.
                proj_energy_cycle = cmplx(0.0, 0.0, p)
                D0_population_cycle = cmplx(0.0, 0.0, p)
                D0_population_ucc_cycle = 0.0_p

                !$omp do schedule(dynamic,200) reduction(+:D0_population_cycle,proj_energy_cycle, D0_population_ucc_cycle, nattempts_spawn,ndeath)
                do iattempt = 1, selection_data%nclusters
                    ! For OpenMP scalability, have this test inside a single loop rather
                    ! than attempt to parallelise over three separate loops.
                    call select_ucc_cluster(rng(it), sys, qs%psip_list, qs%ref%f0, qs%ref%max_ex_level, &
                                            selection_data%nstochastic_clusters, D0_normalisation, qmc_in%initiator_pop, D0_pos, &
                                            cumulative_abs_real_pops, tot_abs_real_pop, min_cluster_size, max_cluster_size, &
                                            logging_info, contrib(it)%cdet, contrib(it)%cluster, qs%excit_gen_data)
                    !print*, 'selected',  contrib(it)%cdet%f

                    if (uccmc_in%variational_energy .and. .not. all(contrib(it)%cdet%f==0) .and. contrib(it)%cluster%excitation_level <= qs%ref%ex_level)  then
                       state = contrib(it)%cdet%f 
                       call binary_search(time_avg_psip_list_ci_states, state, 1, nstates_ci, hit, pos)
                       population = contrib(it)%cluster%amplitude*contrib(it)%cluster%cluster_to_det_sign/contrib(it)%cluster%pselect
                       if (hit) then
                          time_avg_psip_list_ci_pops(pos) = time_avg_psip_list_ci_pops(pos) + population 
                       else
                           do j = nstates_ci, pos, -1
                               ! i is the number of determinants that will be inserted below j.
                               k = j + 1 
                               time_avg_psip_list_ci_states(:,k) = time_avg_psip_list_states(:,j)
                               time_avg_psip_list_ci_pops(k) = time_avg_psip_list_pops(j)
                           end do

                           time_avg_psip_list_ci_states(:,pos) = state
                           time_avg_psip_list_pops(pos) = population
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
                                                    D0_population_cycle, proj_energy_cycle,uccmc_in, ref_det, rdm, selection_data, D0_population_ucc_cycle)

                        call do_stochastic_uccmc_propagation(rng(it), sys, qs, &
                                                                uccmc_in, logging_info, ms_stats(it), bloom_stats, &
                                                                contrib(it), nattempts_spawn, ndeath, ps_stats(it))
                    end if
                end do
                !$omp end do
                ndeath_nc=0
                !$omp end parallel

                ! Add the accumulated ps_stats data to qs%excit_gen_data%p_single_double.
                if (qs%excit_gen_data%p_single_double%vary_psingles) then
                    call ps_stats_reduction_update(qs%excit_gen_data%p_single_double%rep_accum, ps_stats)
                end if

                if (uccmc_in%density_matrices .and. qs%vary_shift(1) .and. parent .and. .not. sys%read_in%comp) then
                    ! Add in diagonal contribution to RDM (only once per cycle not each time reference
                    ! is selected as this is O(N^2))
                    call update_rdm(sys, ref_det, ref_det, real(D0_normalisation,p), 1.0_p, 1.0_p, rdm)
                end if


                qs%psip_list%nparticles = qs%psip_list%nparticles + nparticles_change
                qs%estimators%D0_population_comp = qs%estimators%D0_population_comp + D0_population_cycle
                qs%estimators%proj_energy_comp = qs%estimators%proj_energy_comp + proj_energy_cycle
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
                    call direct_annihilation(sys, rng(0), qs%ref, annihilation_flags, pl, spawn)
                end associate
                if (debug) call write_logging_calc_ccmc(logging_info, iter, nspawn_events, ndeath + ndeath_nc, &
                                                        selection_data%nD0_select, &
                                                        selection_data%nclusters, selection_data%nstochastic_clusters, &
                                                        selection_data%nsingle_excitors)

                if(uccmc_in%variational_energy) then
                          time_avg_psip_list_ci_pops(:nstates_ci) =  time_avg_psip_list_ci_pops(:nstates_ci)/(iter)
                end if

! [review] - AJWT: Comment on what this block is doing.
                do i = 1, qs%psip_list%nstates
                    state = qs%psip_list%states(:,i) 
                    call binary_search(time_avg_psip_list_states, state, 1, nstates_sq, hit, pos)
                    if (hit) then
                          time_avg_psip_list_pops(pos) = time_avg_psip_list_pops(pos) + (real(qs%psip_list%pops(1,i))/qs%psip_list%pop_real_factor)
                          time_avg_psip_list_sq(2,pos) = time_avg_psip_list_sq(2,pos) + (real(qs%psip_list%pops(1,i))/qs%psip_list%pop_real_factor)**2 
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
                       end if
                end do
                ! Update time average.
                time_avg_psip_list_pops(:nstates_sq) =  time_avg_psip_list_pops(:nstates_sq)/(iter)

                time_avg_psip_list_sq(2,:nstates_sq) =  time_avg_psip_list_sq(2,:nstates_sq)/(iter)
            
                call end_mc_cycle(nspawn_events, ndeath_nc, qs%psip_list%pop_real_factor, nattempts_spawn, qs%spawn_store%rspawn)
            end do

            update_tau = bloom_stats%nblooms_curr > 0

            if (uccmc_in%density_matrices .and. qs%vary_shift(1)) then
                call calc_rdm_energy(sys, qs%ref, rdm, qs%estimators(1)%rdm_energy, qs%estimators(1)%rdm_trace)
            end if

            error = qs%spawn_store%spawn%error .or. qs%psip_list%error

            qs%estimators%D0_population = real(qs%estimators%D0_population_comp,p)
            qs%estimators%proj_energy = real(qs%estimators%proj_energy_comp,p)
 
            if (debug) call write_logging_select_ccmc(logging_info, iter, selection_data)
            call end_report_loop(io_unit, qmc_in, iter, update_tau, qs, nparticles_old, nspawn_events, &
                                 -1, semi_stoch_it, soft_exit=soft_exit, &
                                 load_bal_in=load_bal_in, bloom_stats=bloom_stats, comp=sys%read_in%comp, &
                                 error=error, vary_shift_reference=uccmc_in%vary_shift_reference)
            if (error) exit

            call cpu_time(t2)
            if (parent) then
                if (bloom_stats%nblooms_curr > 0) call bloom_stats_warning(bloom_stats, io_unit=io_unit)
                call write_qmc_report(qmc_in, qs, ireport, nparticles_old, t2-t1, .false., .false., &
                                        io_unit=io_unit, cmplx_est=sys%read_in%comp, rdm_energy=uccmc_in%density_matrices, &
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

        if (parent) then
            write (io_unit, '(1X, "Time-averaged cluster populations",/)')
            do i = 1, nstates_sq
                call write_qmc_var(io_unit, time_avg_psip_list_states(1,i))
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
            call var_energy_uccmc(sys, time_avg_psip_list_ci_states,time_avg_psip_list_ci_pops,nstates_ci, var_energy)
            print*, 'Variational energy: ', var_energy
        end if 
        
        if (debug) call end_logging(logging_info)
        if (debug) call end_selection_data(selection_data)

        if (uccmc_in%density_matrices) then
            call write_final_rdm(rdm, sys%nel, sys%basis%nbasis, uccmc_in%density_matrix_file, io_unit)
            call calc_rdm_energy(sys, qs%ref, rdm, qs%estimators(1)%rdm_energy, qs%estimators(1)%rdm_trace)
            if (parent) &
                write (io_unit,'(1x,"# Final energy from RDM",2x,es17.10, /)') qs%estimators%rdm_energy/qs%estimators%rdm_trace
            deallocate(rdm, stat=ierr)
            call check_deallocate('rdm',ierr)
            call dealloc_det_info_t(ref_det)
        end if

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

    end subroutine do_uccmc

    subroutine ucc_set_cluster_selections(selection_data, nattempts, min_cluster_size)
! [review] - AJWT: based on set_cluster_selections with even_selection and full_nc removed.
! [revuew] - AJWT:  Could this then actually use the original routine (assuming both of those are set to .false.)?
        ! Function to set total number of selections of different cluster
        ! types within CCMC. This effectively controls the relative sampling
        ! of different clusters within the CC expansion.

        ! In/Out:
        !   selection_data: derived type containing information on various aspects
        !       of cluster selection.
        !   nattempts: total number of selection attempts to make this iteration.
        !       Initially set to default value in original algorithm, but updated
        !       otherwise. Used to calculate spawning rate and included in output
        !       file.
        ! Out:
        !   min_cluster_size: minimum cluster size to select stochastically.

        use ccmc_data, only: selection_data_t

        type(selection_data_t), intent(inout) :: selection_data
        integer(int_64), intent(inout) :: nattempts
        integer, intent(out) :: min_cluster_size
        integer(int_64) :: nselections

        min_cluster_size = 0
        selection_data%nD0_select = 0 ! instead of this number of deterministic selections, these are chosen stochastically
        selection_data%nstochastic_clusters = nattempts*10**3
        selection_data%nsingle_excitors = 0

        selection_data%nclusters = selection_data%nD0_select + selection_data%nsingle_excitors &
                            + selection_data%nstochastic_clusters

    end subroutine ucc_set_cluster_selections

    subroutine select_ucc_cluster(rng, sys, psip_list, f0, ex_level, nattempts, normalisation, &
                              initiator_pop, D0_pos, cumulative_excip_pop, tot_excip_pop, min_size, max_size, &
                              logging_info, cdet, cluster, excit_gen_data)

! [review] - AJWT: based on select_cluster (without the linked_cluster parts) and with information about de-excitors.
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
        !    cumulative_excip_population: running cumulative excip population on
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
! [review] - AJWT: collapse_cluster is no longer used.
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
        real(p), intent(in) :: cumulative_excip_pop(:), tot_excip_pop
        integer :: min_size, max_size
        type(dSFMT_t), intent(inout) :: rng
        type(det_info_t), intent(inout) :: cdet
        type(cluster_t), intent(inout) :: cluster
        type(logging_t), intent(in) :: logging_info
        type(excit_gen_data_t), intent(in) :: excit_gen_data

        real(dp) :: rand
        real(p) :: psize
        complex(p) :: cluster_population, excitor_pop
        integer :: i, pos, prev_pos, deexcit_count, ierr
        real(p) :: pop(max_size)
        logical :: hit, allowed
        logical, allocatable :: deexcitation(:)
       
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
            ! Rather than selecting one excitor at a time and adding it to the
            ! cluster, select all excitors and then find their locations and
            ! apply them.  This allows us to sort by population first (as the
            ! number of excitors is small) and hence allows for a more efficient
            ! searching of the cumulative population list.

            do i = 1, cluster%nexcitors
                ! Select a position in the excitors list.
                pop(i) = get_rand_close_open(rng)*tot_excip_pop
                rand = get_rand_close_open(rng)
                !Decide if this is to be an excitation or deexcitation operator
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
            call insert_sort(pop(:cluster%nexcitors))
            prev_pos = 1
            do i = 1, cluster%nexcitors
                call binary_search(cumulative_excip_pop, pop(i), prev_pos, psip_list%nstates, hit, pos)
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
            else
                ! Simply set excitation level to a too high (fake) level to avoid
                ! this cluster being used.
                cluster%excitation_level = huge(0)
            end if
        end select

        if(cluster%nexcitors>0)  then
            deallocate(deexcitation, stat=ierr)
            call check_deallocate('deexcitation',ierr)
        end if 
        if (debug) call write_logging_stoch_selection(logging_info, cluster%nexcitors, cluster%excitation_level, pop, &
                max_size, cluster%pselect, cluster%amplitude, allowed)

    end subroutine select_ucc_cluster

    subroutine ucc_collapse_cluster(basis, f0, excitor, excitor_population, cluster_excitor, cluster_population, &
                    allowed, conjugate)
! [review] - AJWT: based on collapse_clustera

! [review] - AJWT: Document 'conjugate'

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

        integer(i0) :: excitor_loc(basis%tot_string_len), cluster_loc(basis%tot_string_len), f0_loc(basis%tot_string_len)

        integer :: ibasis, ibit
        integer(i0) :: excitor_excitation(basis%tot_string_len)
        integer(i0) :: excitor_annihilation(basis%tot_string_len)
        integer(i0) :: excitor_creation(basis%tot_string_len)
        integer(i0) :: cluster_excitation(basis%tot_string_len)
        integer(i0) :: cluster_annihilation(basis%tot_string_len)
        integer(i0) :: cluster_creation(basis%tot_string_len)
        integer(i0) :: permute_operators(basis%tot_string_len)
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
    end subroutine ucc_collapse_cluster

    subroutine ucc_cumulative_population(pops, ex_lvls, nactive, D0_proc, D0_pos, real_factor, complx, &
                                    cumulative_pops, tot_pop)
! [review] - AJWT: based on cumulative_population
! [review] - AJWT: is this the same just without calc_dist? and can that be set to .false. and cumulative_population called?

        ! Calculate the cumulative population, i.e. the number of psips/excips
        ! residing on a determinant/an excitor and all determinants/excitors which
        ! occur before it in the determinant/excitor list.

        ! This is primarily so in CCMC we can select clusters of excitors with each
        ! excip being equally likely to be part of a cluster.  (If we just select
        ! each occupied excitor with equal probability, then we get wildy
        ! fluctuating selection probabilities and hence large population blooms.)
        ! As 'excips' on the reference cannot be part of a cluster, then the
        ! population on the reference is treated as 0 if required.

        ! In:
        !    pops: list of populations on each determinant/excitor.  Must have
        !       minimum length of nactive.
        !    nactive: number of occupied determinants/excitors (ie pops(:,1:nactive)
        !       contains the population(s) on each currently "active"
        !       determinant/excitor.
        !    D0_proc: processor on which the reference resides.
        !    D0_pos: position in the pops list of the reference.  Only relevant if
        !       1<=D0_pos<=nactive and the processor holds the reference.
        !    real_factor: the encoding factor by which the stored populations are multiplied
        !       to enable non-integer populations.
        ! Out:
        !    cumulative_pops: running total of excitor population, i.e.
        !        cumulative_pops(i) = sum(abs(pops(1:i))), excluding the
        !        population on the reference if appropriate.
        !    tot_pop: total population (possibly excluding the population on the
        !       reference).

        ! NOTE: currently only the populations in the first psip/excip space are
        ! considered.  This should be changed if we do multiple simulations at
        ! once/Hellmann-Feynman sampling/etc.

        use parallel, only: iproc
        use ccmc_utils, only: get_pop_contrib

        integer(int_p), intent(in) :: pops(:,:), real_factor
        integer, intent(in) :: nactive, D0_proc, D0_pos
        real(p), allocatable, intent(inout) :: cumulative_pops(:)
        real(p), intent(out) :: tot_pop
        logical, intent(in) :: complx
        integer(i0), intent(in) :: ex_lvls(:)

        integer :: i
        integer(i0) :: j, ex_lvl
        ! Need to combine spaces if doing complex; we choose combining in quadrature.
        cumulative_pops(1) = get_pop_contrib(pops(:,1), real_factor, complx)
        if (D0_proc == iproc) then
            ! Let's be a bit faster: unroll loops and skip over the reference
            ! between the loops.
            do i = 2, d0_pos-1
                cumulative_pops(i) = cumulative_pops(i-1) + &
                                        get_pop_contrib(pops(:,i), real_factor, complx)
            end do
            ! Set cumulative on the reference to be the running total merely so we
            ! can continue accessing the running total from the i-1 element in the
            ! loop over excitors in slots above the reference.
            if (d0_pos == 1) then
                cumulative_pops(d0_pos) = 0
            end if
            if (d0_pos > 1) cumulative_pops(d0_pos) = cumulative_pops(d0_pos-1)
            do i = d0_pos+1, nactive
                cumulative_pops(i) = cumulative_pops(i-1) + &
                                        get_pop_contrib(pops(:,i), real_factor, complx)
            end do
        else
            ! V simple on other processors: no reference to get in the way!
            do i = 2, nactive
                cumulative_pops(i) = cumulative_pops(i-1) + &
                                        get_pop_contrib(pops(:,i), real_factor, complx)
            end do
        end if
        if (nactive > 0) then
            tot_pop = cumulative_pops(nactive)
        else
            tot_pop = 0.0_p
        end if

    end subroutine ucc_cumulative_population

    subroutine do_uccmc_accumulation(sys, qs, cdet, cluster, logging_info, D0_population_cycle, proj_energy_cycle, &
                                     uccmc_in, ref_det, rdm, selection_data, D0_population_ucc_cycle)
! [review] - AJWT: based on do_ccmc_accumulation
! [review] - AJWT: If uccmc_in is unfolded from ccmc_in, then can the above be reused by changing what the
! [review] - AJWT: update_proj_energy_ptr points to?



        ! Performs all accumulation of values required for given ccmc clusters.
        ! Updates projected energy and any RDMs with the contribution from the
        ! current cluster.

        ! In:
        !   sys: information on system under consideration.
        !   qs: information on current state of qmc calculation.
        !   ccmc_in: options relating to ccmc passed in to calculation.
        !   cdet: information on determinant resulting from collapsing current
        !       cluster.
        !   ref_det: reference determinant information.
        !   cluster: information on cluster currently under consideration.
        !   logging_info: current logging settings in use.

        ! In/Out:
        !   D0_population_cycle: running total of reference population.
        !   proj_energy_cycle: running total of projected energy contributions.
        !   rdm: array containing reduced density matrix.
        !   selection_data: info on clsuter selection.

        use system, only: sys_t
        use qmc_data, only: qmc_state_t, uccmc_in_t, estimators_t
        use determinants, only: det_info_t
        use ccmc_data, only: cluster_t, selection_data_t
        use ccmc_selection, only: update_selection_data
        use qmc_data, only: zero_estimators_t

        use excitations, only: excit_t, get_excitation
        use hamiltonian_data, only: hmatel_t
        use replica_rdm, only: update_rdm
        use logging, only: logging_t
        use proc_pointers, only: update_proj_energy_ptr

        type(sys_t), intent(in) :: sys
        type(qmc_state_t), intent(in) :: qs
        type(uccmc_in_t), intent(in) :: uccmc_in
        type(det_info_t), intent(in) :: cdet, ref_det
        type(cluster_t), intent(in) :: cluster
        type(logging_t), intent(in) :: logging_info
        complex(p), intent(inout) :: D0_population_cycle, proj_energy_cycle
        real(p), intent(inout), optional :: D0_population_ucc_cycle
        real(p), allocatable, intent(inout) :: rdm(:,:)
        type(selection_data_t), intent(inout) :: selection_data
        type(excit_t) :: connection
        type(hmatel_t) :: hmatel
        type(estimators_t) :: estimators_cycle

        if (debug) call update_selection_data(selection_data, cluster, logging_info)

        if (cluster%excitation_level /= huge(0)) then
            ! FCIQMC calculates the projected energy exactly.  To do
            ! so in CCMC would involve enumerating over all pairs of
            ! single excitors, which is rather painful and slow.
            ! Instead, as we are randomly sampling clusters in order
            ! to evolve the excip population anyway, we can just use
            ! the random clusters to *sample* the projected
            ! estimator.  See comments in spawning.F90 for why we
            ! must divide through by the probability of selecting
            ! the cluster.
            call zero_estimators_t(estimators_cycle)
            connection = get_excitation(sys%nel, sys%basis, cdet%f, qs%ref%f0)
            !if (uccmc_in%trot) then
            !    call update_proj_energy_ptr(sys, qs%ref%f0, qs%trial%wfn_dat, cdet, &
            !         [real(cluster%amplitude,p),aimag(cluster%amplitude)]*&
            !         cluster%cluster_to_det_sign/cluster%pselect, &
            !         estimators_cycle, connection, hmatel)
            !else
                call update_proj_energy_mol_ucc(sys, qs%ref%f0, qs%trial%wfn_dat, cdet, &
                     [real(cluster%amplitude,p),aimag(cluster%amplitude)]*&
                     cluster%cluster_to_det_sign/cluster%pselect, &
                     estimators_cycle, connection, hmatel, cluster%nexcitors)
            !end if

            D0_population_cycle = D0_population_cycle + estimators_cycle%D0_population
            proj_energy_cycle = proj_energy_cycle + estimators_cycle%proj_energy
            !if (.not. uccmc_in%trot)  
            D0_population_ucc_cycle = D0_population_ucc_cycle + estimators_cycle%D0_population_ucc
        end if

        if (uccmc_in%density_matrices .and. cluster%excitation_level <= 2 .and. qs%vary_shift(1) &
            .and. cluster%excitation_level /= 0 .and. .not. sys%read_in%comp) then
            ! Add contribution to density matrix
            ! d_pqrs = <HF|a_p^+a_q^+a_sa_r|CC>
            !$omp critical
            call update_rdm(sys, cdet, ref_det, &
                            real(cluster%amplitude, p)*cluster%cluster_to_det_sign, &
                            1.0_p, cluster%pselect, rdm)
            !$omp end critical
        end if

    end subroutine do_uccmc_accumulation

    subroutine do_stochastic_uccmc_propagation(rng, sys, qs, uccmc_in, &
                                            logging_info, ms_stats, bloom_stats, &
                                            contrib, nattempts_spawn_tot, ndeath, ps_stat)
! [review] - AJWT: based on do_stochastic_ccmc_propagation.
! [review] - AJWT: seems to be much the same with multiref .false. and calling  do_ucc_spawning_death.  Could this be function-pointered?

        ! Perform stochastic propogation of a cluster in an appropriate manner
        ! for the given inputs. For multireference systems, it allows death for
        ! clusters within the desired truncation level of either reference.

        ! In:
        !   sys: information on system under consideration.
        !   uccmc_in: options relating to uccmc passed in to calculation.
        !   logging_info: logging_t object with info about current logging
        !        when debug is true. 
        ! In/Out:
        !   rng: random number generator.
        !   qs: qmc_state_t type, contains information about calculation.
        !   ms_stats: statistics on multispawn performance.
        !   bloom_stats: statistics on blooms during calculation.
        !   contrib: derived type containing information on the current
        !       wavefunction contribution being considered.
        !   nattempts_spawn_tot: running total of number of spawning attempts
        !       made during this mc cycle.
        !   ndeath: total number of particles created via death.
        !   ps_stat: Accumulating the following (and more) on this OpenMP thread:
        !       h_pgen_singles_sum: total of |Hij|/pgen for single excitations attempted.
        !       h_pgen_doubles_sum: total on |Hij|/pgen for double excitations attempted.
        !       excit_gen_singles: counter on number of single excitations attempted.
        !       excit_gen_doubles: counter on number of double excitations attempted.


        use dSFMT_interface, only: dSFMT_t
        use system, only: sys_t
        use qmc_data, only: qmc_state_t, uccmc_in_t
        use ccmc_data, only: multispawn_stats_t, ms_stats_update, wfn_contrib_t
        use bloom_handler, only: bloom_stats_t, accumulate_bloom_stats
        use logging, only: logging_t
        use excit_gens, only: p_single_double_coll_t

        use excitations, only: get_excitation_level
        
        type(sys_t), intent(in) :: sys
        type(dSFMT_T), intent(inout) :: rng
        type(qmc_state_t), intent(inout) :: qs
        type(uccmc_in_t), intent(in) :: uccmc_in
        type(wfn_contrib_t), intent(inout) :: contrib
        type(bloom_stats_t), intent(inout) :: bloom_stats
        type(logging_t), intent(in) :: logging_info

        integer(int_64), intent(inout) :: nattempts_spawn_tot
        integer(int_p), intent(inout) :: ndeath
        type(multispawn_stats_t), intent(inout) :: ms_stats
        type(p_single_double_coll_t), intent(inout) :: ps_stat

        integer :: nspawnings_cluster
        logical :: attempt_death

        ! Spawning
        ! This has the potential to create blooms, so we allow for multiple
        ! spawning events per cluster.
        ! The number of spawning events is decided by the value
        ! of cluster%amplitude/cluster%pselect.  If this is
        ! greater than cluster_multispawn_threshold, then nspawnings is
        ! increased to the ratio of these.
        nspawnings_cluster=max(1,ceiling((abs(contrib%cluster%amplitude)&
                                     /contrib%cluster%pselect)/ uccmc_in%cluster_multispawn_threshold))
        call ms_stats_update(nspawnings_cluster, ms_stats)
        nattempts_spawn_tot = nattempts_spawn_tot + nspawnings_cluster
        attempt_death = (contrib%cluster%excitation_level <= qs%ref%ex_level) 
        call do_ucc_spawning_death(rng, sys, qs, uccmc_in, &
                             logging_info, ms_stats, bloom_stats, contrib, &
                             nattempts_spawn_tot, ndeath, ps_stat, nspawnings_cluster, attempt_death)

    end subroutine do_stochastic_uccmc_propagation  

    subroutine do_ucc_spawning_death(rng, sys, qs, uccmc_in, &
                                 logging_info, ms_stats, bloom_stats, contrib, &
                                 nattempts_spawn_tot, ndeath, ps_stat, nspawnings_cluster, attempt_death)

! [review] - AJWT: based on do_spawning_death, but calling perform_uccmc_spawning_attempt instead.  Function pointer?
        ! For stochastically selected clusters this
        ! attempts spawning and death, adding any created particles to the
        ! spawned list. For deterministically selected clusters spawning is
        ! performed, while death is performed in-place separately.

        ! This is currently non-specific to a any given type of selected
        ! cluster, though this may change in future.

        ! In:
        !   sys: information on system under consideration.
        !   ccmc_in: options relating to ccmc passed in to calculation.
        !   logging_info: logging_t object with info about current logging
        !        when debug is true. 
        !   attempt_death = logical variable that encodes whether the selected 
        !        cluster is within the desired CC truncation. 
        ! In/Out:
        !   rng: random number generator.
        !   qs: qmc_state_t type, contains information about calculation.
        !   ms_stats: statistics on multispawn performance.
        !   bloom_stats: statistics on blooms during calculation.
        !   contrib: derived type containing information on the current
        !       wavefunction contribution being considered.
        !   nattempts_spawn_tot: running total of number of spawning attempts
        !       made during this mc cycle.
        !   ndeath: total number of particles created via death.
        !   ps_stat: Accumulating the following (and more) on this OpenMP thread:
        !       h_pgen_singles_sum: total of |Hij|/pgen for single excitations attempted.
        !       h_pgen_doubles_sum: total on |Hij|/pgen for double excitations attempted.
        !       excit_gen_singles: counter on number of single excitations attempted.
        !       excit_gen_doubles: counter on number of double excitations attempted.

        use dSFMT_interface, only: dSFMT_t
        use system, only: sys_t
        use qmc_data, only: qmc_state_t, uccmc_in_t
        use ccmc_data, only: multispawn_stats_t, ms_stats_update, wfn_contrib_t

        use ccmc_death_spawning, only: stochastic_ccmc_death
        use bloom_handler, only: bloom_stats_t, accumulate_bloom_stats
        use logging, only: logging_t
        use excit_gens, only: p_single_double_coll_t

        use excitations, only: get_excitation_level
        
        type(sys_t), intent(in) :: sys
        type(dSFMT_T), intent(inout) :: rng
        type(qmc_state_t), intent(inout) :: qs
        type(uccmc_in_t), intent(in) :: uccmc_in
        type(wfn_contrib_t), intent(inout) :: contrib
        type(bloom_stats_t), intent(inout) :: bloom_stats
        type(logging_t), intent(in) :: logging_info

        integer(int_64), intent(inout) :: nattempts_spawn_tot
        integer(int_p), intent(inout) :: ndeath
        type(multispawn_stats_t), intent(inout) :: ms_stats
        type(p_single_double_coll_t), intent(inout) :: ps_stat
        integer, intent(in) :: nspawnings_cluster
        logical, intent(in) :: attempt_death

        integer :: i

        do i = 1, nspawnings_cluster
            call perform_uccmc_spawning_attempt(rng, sys, qs, uccmc_in, logging_info, bloom_stats, contrib, &
                                               nspawnings_cluster, ps_stat)
        end do

        ! Does the cluster collapsed onto D0 produce
        ! a determinant is in the truncation space?  If so, also
        ! need to attempt a death/cloning step.
        ! optimisation: call only once per iteration for clusters of size 0 or 1 for ccmc_in%full_nc.
        if (attempt_death) then
           ! Clusters above size 2 can't die in linked ccmc.
            if ((.not. uccmc_in%linked) .or. contrib%cluster%nexcitors <= 2) then
                call stochastic_ccmc_death(rng, qs%spawn_store%spawn, uccmc_in%linked, .false., sys, &
                                           qs, contrib%cdet, contrib%cluster, logging_info, ndeath)
            end if
        end if
    end subroutine do_ucc_spawning_death

    subroutine perform_uccmc_spawning_attempt(rng, sys, qs, uccmc_in, logging_info, bloom_stats, contrib, nspawnings_total, &
                                        ps_stat)

! [review] - AJWT: based on perform_ccmc_spawning_death.
! [review] - AJWT: but copes with Trotterized spawning. Could that be function-pointered?

        ! Performs a single ccmc spawning attempt, as appropriate for a
        ! given setting combination of linked, complex or none of the above.

        ! This could be replaced by a procedure pointer with a little
        ! refactoring if desired.

        ! In:
        !   sys: the system being studied.
        !   ccmc_in: input settings related to ccmc.
        !   logging_info: input settings related to logging.
        !   nspawnings_total: number of spawning attempts made
        !       from the same cluster if using multispawn.
        ! In/Out:
        !   rng: random number generator.
        !   qs: information on current state of calculation.
        !   contrib: information on contribution to wavefunction
        !       currently under consideration. Contains both the
        !       cluster selected and the determinant formed on
        !       collapsing, as well as scratch spaces for partitoning
        !       within linked.
        !   ps_stat: Accumulating the following (and more) on this OpenMP thread:
        !       h_pgen_singles_sum: total of |Hij|/pgen for single excitations attempted.
        !       h_pgen_doubles_sum: total on |Hij|/pgen for double excitations attempted.
        !       excit_gen_singles: counter on number of single excitations attempted.
        !       excit_gen_doubles: counter on number of double excitations attempted.
        use dSFMT_interface, only: dSFMT_t
        use system, only: sys_t
        use qmc_data, only: qmc_state_t, uccmc_in_t
        use ccmc_data, only: wfn_contrib_t

        use excitations, only: excit_t
        use proc_pointers, only: gen_excit_ptr

        use ccmc_death_spawning, only: spawner_ccmc, linked_spawner_ccmc
        use ccmc_death_spawning, only: spawner_complex_ccmc, create_spawned_particle_ccmc
        use bloom_handler, only: bloom_stats_t, accumulate_bloom_stats
        use logging, only: logging_t
        use excit_gens, only: p_single_double_coll_t

        type(sys_t), intent(in) :: sys
        type(dSFMT_T), intent(inout) :: rng
        type(qmc_state_t), intent(inout) :: qs
        type(uccmc_in_t), intent(in) :: uccmc_in
        type(wfn_contrib_t), intent(inout) :: contrib
        type(excit_t) :: connection
        type(bloom_stats_t), intent(inout) :: bloom_stats
        type(logging_t), intent(in) :: logging_info

        integer, intent(in) :: nspawnings_total
        type(p_single_double_coll_t), intent(inout) :: ps_stat
        integer(int_p) :: nspawned
        integer(i0) :: fexcit(sys%basis%tot_string_len)
 
        call spawner_ccmc(rng, sys, qs, qs%spawn_store%spawn%cutoff, &
                      uccmc_in%linked, contrib%cdet, contrib%cluster, gen_excit_ptr, logging_info, &
                      nspawned, connection, nspawnings_total, ps_stat)
        if (nspawned /= 0_int_p) then
            if(uccmc_in%trot) then
                ! Must decide if particle is labelled in any way.
                call create_spawned_particle_uccmc_trot(sys, qs%ref, contrib%cdet, connection, &
                                            nspawned, 1, contrib%cluster%excitation_level, &
                                            .true., fexcit, qs%spawn_store%spawn, bloom_stats)
            else
                call create_spawned_particle_ccmc(sys%basis, qs%ref, contrib%cdet, connection, &
                                            nspawned, 1, contrib%cluster%excitation_level, &
                                            .false., fexcit, qs%spawn_store%spawn, bloom_stats)
            end if
        end if
    end subroutine perform_uccmc_spawning_attempt

    subroutine create_spawned_particle_uccmc_trot(sys, ref, cdet, connection, nspawned, ispace, &
                                            parent_cluster_ex_level, trot, fexcit, spawn, bloom_stats)

! [review] - AJWT: How does this differ (from what?) for the uccmc_trotterization?

        ! Function to create spawned particle in spawned list for ccmc
        ! calculations. Performs required manipulations of bit string
        ! beforehand and accumulation on blooming.

        ! In:
        !   basis: info on current basis functions.
        !   reference: info on current reference state.
        !   cdet: determinant representing state currently spawning
        !       spawning from.
        !   connection: connection from state cdet particle has been spawned
        !       from.
        !   nspawned: number of (encoded) particles to be created via this spawning.
        !   ispace: index of space particles are to be added to.
        !   parent_cluster_ex_level: excitation level of parent cluster.
        !   fexcit: bit string for state spawned to. Only available for linked
        !       ccmc, otherwise generated using cdet+connection.
        !   ex_lvl_sort: true if require states to be sorted by excitation
        !       level within walker list, false otherwise.
        ! In/Out:
        !   spawn: spawn_t type containing information on particles created
        !       via spawning this iteration. Spawned particles will be added
        !       to this on exit.
        !   bloom_stats: information on blooms within a calculation. Will be
        !       updated if a bloom has occurred.

        use system, only: sys_t
        use reference_determinant, only: reference_t
        use spawn_data, only: spawn_t
        use determinant_data, only: det_info_t
        use bloom_handler, only: bloom_stats_t, accumulate_bloom_stats
        use excitations, only: excit_t, create_excited_det
        use proc_pointers, only: create_spawned_particle_ptr
        use ccmc_utils, only: add_ex_level_bit_string_calc

        type(sys_t), intent(in) :: sys
        type(reference_t), intent(in) :: ref
        type(spawn_t), intent(inout) :: spawn
        type(det_info_t), intent(in) :: cdet
        type(bloom_stats_t), intent(inout) :: bloom_stats
        type(excit_t), intent(in) :: connection

        integer(int_p), intent(in) :: nspawned
        integer, intent(in) :: ispace, parent_cluster_ex_level
        integer(i0), intent(in) :: fexcit(:)
        logical, intent(in) :: trot 
        integer(i0) :: fexcit_loc(lbound(fexcit,dim=1):ubound(fexcit,dim=1))

        
        if (parent_cluster_ex_level /= huge(0)) then
            call create_excited_det(sys%basis, cdet%f, connection, fexcit_loc)
        else
            fexcit_loc = fexcit
        end if

        if (trot) call add_info_str_trot(sys%basis, ref%f0, sys%nel, fexcit_loc)

        
        call create_spawned_particle_ptr(sys%basis, ref, cdet, connection, nspawned, &
                                        ispace, spawn, fexcit_loc)
        call accumulate_bloom_stats(bloom_stats, nspawned)

    end subroutine create_spawned_particle_uccmc_trot
    
    pure function earliest_unset(f, f0, nel, basis) result (early)
! [review] - AJWT: probably best in ccmc_utils or the like
        
         ! Function to find earliest unset bit in a determinant bit string.

         ! In:
         !    f: bit_string encoding determinant
         !    nel: number of electrons in the system
         !    basis: basis_t object with information on one-electron basis in use.
     
         use basis_types, only: basis_t
         use bit_utils, only: count_set_bits

         type(basis_t), intent(in) :: basis
         integer(i0), intent(in) :: f(basis%tot_string_len), f0(basis%tot_string_len)
         integer, intent(in) :: nel
         integer :: i, early
         integer(i0) :: diff
         integer(i0) :: f0_loc
         
         early = 0
         i = 0
          
         diff = ieor(f(1),f0(1))
         if (diff /= 0) then
             do 
                 if (btest(diff,i)) then
                     early = i
                     exit
                 end if
                 i = i+1
             end do
         else
             f0_loc = f0(1)
             do while (f0_loc/=0)
                 f0_loc = ibclr(f0_loc, early)
                 early = early + 1
             end do
         end if
    end function

    pure function latest_unset(f, f0, nel, basis) result (late)
! [review] - AJWT: probably best in ccmc_utils or the like
        
         ! Function to find earliest unset bit in a determinant bit string.

         ! In:
         !    f: bit_string encoding determinant
         !    nel: number of electrons in the system
         !    basis: basis_t object with information on one-electron basis in use.
     
         use basis_types, only: basis_t
         use bit_utils, only: count_set_bits

         type(basis_t), intent(in) :: basis
         integer(i0), intent(in) :: f(basis%tot_string_len), f0(basis%tot_string_len)
         integer, intent(in) :: nel
         integer :: late
         integer(i0) :: diff
         integer(i0) :: f0_loc
         
         late = 0
          
         diff = ieor(f(1),f0(1))
         diff = iand(diff, f0(1))
         if (diff /= 0) then
             do 
                 diff = ibclr(diff, late)
                 if (diff /= 0) then
                     late = late + 1 
                 else
                     exit
                 end if
             end do
         else
             f0_loc = f0(1)
             do while (f0_loc/=0)
                 f0_loc = ibclr(f0_loc, late)
                 late = late + 1
             end do
         end if
    end function
    subroutine add_info_str_trot(basis, f0, nel, f)

        ! Sets bits within bit string to give excitation level at end of bit strings.
        ! This routine sets ex level from provided reference.

        use basis_types, only: basis_t
        use excitations, only: get_excitation_level
        type(basis_t), intent(in) :: basis
        integer(i0), intent(inout) :: f(:)
        integer(i0), intent(in) :: f0(:)
        integer, intent(in) :: nel

        integer(i0) :: counter(basis%tot_string_len)
       

        if (basis%info_string_len/=0) then

            f(basis%bit_string_len+2) = latest_unset(f, f0, nel, basis)     
            f(basis%bit_string_len+1) = nel - get_excitation_level(f(:basis%bit_string_len), f0(:basis%bit_string_len)) 

        end if
 
    end subroutine add_info_str_trot

    subroutine add_info_str_trot_alt(basis, f0, nel, f)

        ! Sets bits within bit string to give excitation level at end of bit strings.
        ! This routine sets ex level from provided reference.

        use basis_types, only: basis_t
        use excitations, only: get_excitation_level
        type(basis_t), intent(in) :: basis
        integer(i0), intent(inout) :: f(:)
        integer(i0), intent(in) :: f0(:)
        integer, intent(in) :: nel

        integer(i0) :: counter(basis%tot_string_len)
       

        if (basis%info_string_len/=0) then

            select case(f(1))

            case(3)
            f(basis%bit_string_len+2) = 1
            case(6)
            f(basis%bit_string_len+2) = 2
            case(9)
            f(basis%bit_string_len+2) = 3
            case(258) 
            f(basis%bit_string_len+2) = 4
            case(513)
            f(basis%bit_string_len+2) = 5
            case(12)
            f(basis%bit_string_len+2) = 6
            case(48)
            f(basis%bit_string_len+2) = 7
            case(192)
            f(basis%bit_string_len+2) = 8
            case(768)
            f(basis%bit_string_len+2) = 9
            case(516)
            f(basis%bit_string_len+2) = 10
            case(264)
            f(basis%bit_string_len+2) = 11
            
            end select
        end if
 
    end subroutine add_info_str_trot_alt

    subroutine update_proj_energy_mol_ucc(sys, f0, wfn_dat, cdet, pop, estimators, excitation, hmatel, cluster_size)

! [review] - AJWT: based on update_proj_energy_mol - this might be able to be folded back into this eventually.
        ! Add the contribution of the current determinant to the projected
        ! energy.
        ! The correlation energy given by the projected energy is:
        !   \sum_{i \neq 0} <D_i|H|D_0> N_i/N_0
        ! where N_i is the population on the i-th determinant, D_i,
        ! and 0 refers to the reference determinant.
        ! During a MC cycle we store N_0 and \sum_{i \neq 0} <D_i|H|D_0> N_i.
        ! For UCCMC, we store separately the population on the empty cluster D0 and
        ! the sum of all contributions which land on the reference determinant, as the
        ! latter is needed to correctly compute E_proj.
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
        use energy_evaluation, only: estimators_t, hmatel_t

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f0(:)
        real(p), intent(in) :: wfn_dat(:)
        type(det_info_t), intent(in) :: cdet
        real(p), intent(in) :: pop(:)
        type(estimators_t), intent(inout) :: estimators
        type(excit_t), intent(inout) :: excitation
        type(hmatel_t), intent(out) :: hmatel
        integer, intent(in) :: cluster_size

        integer :: ij_sym, ab_sym

        hmatel%r = 0.0_p

        select case(excitation%nexcit)
        case (0)
            ! Have reference determinant.
            estimators%D0_population = estimators%D0_population + pop(1)
            if (cluster_size == 0) then
                estimators%D0_population_ucc = estimators%D0_population_ucc + pop(1)
            end if
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

    end subroutine update_proj_energy_mol_ucc

    subroutine var_energy_uccmc(sys, time_avg_psip_list_states, time_avg_psip_list_pops, nstates,var_energy)

! [review] - AJWT: Documentation...         
       use excitations, only: excit_t, get_excitation
       use hamiltonian, only: get_hmatel
       use energy_evaluation, only: hmatel_t
       use system, only: sys_t
       use qmc_data, only: particle_t
       use read_in_symmetry, only: cross_product_basis_read_in
       use determinants, only: decode_det

       type(sys_t), intent(in) :: sys
       integer(i0), intent(in) :: time_avg_psip_list_states(:,:)
       integer, intent(in) :: nstates
       real(p), intent(in) :: time_avg_psip_list_pops(:)
       real(p), intent(out) :: var_energy
       real(p) :: normalisation

       type(excit_t) :: excitation
       type(hmatel_t) :: hmatel
       integer :: occ_list(sys%nel)

       integer :: i, j
       integer :: ij_sym, ab_sym

       normalisation = 0.0_p
       var_energy = 0.0_p
! [review] - AJWT:  Consider real(x, p) which is usually used elsewhere.
       do i = 1, nstates
           normalisation = normalisation + (time_avg_psip_list_pops(i)/time_avg_psip_list_pops(1))**2
           do j = 1, nstates
               hmatel = get_hmatel(sys, time_avg_psip_list_states(:,i), time_avg_psip_list_states(:,j))
               var_energy = var_energy + hmatel%r*time_avg_psip_list_pops(i)/time_avg_psip_list_pops(1)*time_avg_psip_list_pops(j)/time_avg_psip_list_pops(1)
           end do
       end do
       var_energy = var_energy/normalisation
   end subroutine var_energy_uccmc

end module

module trotterized_uccmc

use const, only: i0, int_p, int_64, p, dp, debug

implicit none

contains

    subroutine do_trot_uccmc(sys, qmc_in, uccmc_in, restart_in, load_bal_in, reference_in, &
                        logging_in, io_unit, qs, qmc_state_restart)

        ! Run the UCCMC algorithm starting from the initial walker distribution
        ! using the timestep algorithm.

        ! See notes about the implementation of this using function pointers
        ! in fciqmc_main.

        ! In:
        !    sys: system being studied.
        !    uccmc_in: input options relating to UCCMC.
        !    restart_in: input options for HDF5 restart files.
        !    reference_in: current reference determinant.  If not set (ie
        !       components allocated) then a best guess is made based upon the
        !       desired spin/symmetry.
        !    load_bal_in: input options for load balancing.
        !    qmc_in: input options relating to QMC methods.
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
        use ccmc_selection, only: init_selection_data, update_selection_probabilities, set_cluster_selections, &
                                  init_amp_psel_accumulation
        use ccmc_utils, only: get_D0_info, init_contrib, dealloc_contrib, cumulative_population, & 
                              regenerate_ex_levels_psip_list
        use determinants, only: alloc_det_info_t, dealloc_det_info_t, sum_sp_eigenvalues_occ_list, &
                                sum_sp_eigenvalues_bit_string, decode_det
        use determinant_data, only: det_info_t
        use excitations, only: excit_t, get_excitation_level, get_excitation
        use qmc_io, only: write_qmc_report, write_qmc_report_header
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

        real(p), allocatable :: rdm(:,:)

        integer :: iunit, restart_version_restart
        integer :: date_values(8)
        character(:), allocatable :: err_msg
 
        type(particle_t) :: time_avg_psip_list_ci, time_avg_psip_list
        integer :: semi_stoch_it, pos, j, k
        logical :: hit
        integer(i0), allocatable :: state(:)
        integer(int_p), allocatable :: population(:)
        real(p), allocatable :: real_population(:), var_energy(:)

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


        allocate(state(sys%basis%bit_string_len))
        allocate(real_population(qs%psip_list%nspaces))

        if(uccmc_in%variational_energy) then
             allocate(state(sys%basis%bit_string_len))
             allocate(real_population(qs%psip_list%nspaces))
             allocate(population(qs%psip_list%nspaces))
             allocate(var_energy(qs%psip_list%nspaces))
             population(:) = 0
             call init_particle_t(qmc_in%walker_length, 1, sys%basis%tensor_label_len, qmc_in%real_amplitudes, &
                             qmc_in%real_amplitude_force_32, time_avg_psip_list_ci, io_unit=io_unit)
             time_avg_psip_list_ci%pops(:,1) = qs%psip_list%pops(:,1)
             time_avg_psip_list_ci%states(:,1) = qs%psip_list%states(:,1)
             time_avg_psip_list_ci%nparticles = qs%psip_list%pops(:,1)/qs%psip_list%pop_real_factor
             time_avg_psip_list_ci%nstates = 1
             ![todo] deal with restarting
        end if

        call init_particle_t(qmc_in%walker_length, 1, sys%basis%tensor_label_len, qmc_in%real_amplitudes, &
                             qmc_in%real_amplitude_force_32, time_avg_psip_list, io_unit=io_unit)
             time_avg_psip_list%pops(:,1) = qs%psip_list%pops(:,1)
             time_avg_psip_list%states(:,1) = qs%psip_list%states(:,1)
             time_avg_psip_list%nparticles = qs%psip_list%pops(:,1)/qs%psip_list%pop_real_factor
             time_avg_psip_list%nstates = 1
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
            call init_report_loop(qs, bloom_stats)

            do icycle = 1, qmc_in%ncycles
                
                iter = qs%mc_cycles_done + (ireport-1)*qmc_in%ncycles + icycle

                if(uccmc_in%variational_energy) then
                          time_avg_psip_list_ci%nparticles = time_avg_psip_list_ci%nparticles*(iter-1)
                          time_avg_psip_list_ci%pops(:,:time_avg_psip_list_ci%nstates) =  time_avg_psip_list_ci%pops(:,:time_avg_psip_list_ci%nstates)*(iter-1)
                end if
                time_avg_psip_list%nparticles = time_avg_psip_list%nparticles*(iter-1)
                time_avg_psip_list%pops(:,:time_avg_psip_list%nstates) =  time_avg_psip_list%pops(:,:time_avg_psip_list%nstates)*(iter-1)

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

                call trot_ucc_set_cluster_selections(selection_data, qs%estimators(1)%nattempts, min_cluster_size)
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
                !$omp        nattempts_spawn, &
                !$omp        uccmc_in, nprocs, ms_stats, ps_stats, qmc_in, load_bal_in, &
                !$omp        ndeath_nc,   &
                !$omp        nparticles_change, ndeath, logging_info, time_avg_psip_list_ci)

                it = get_thread_id()
                iexcip_pos = 0
                seen_D0 = .false.
                proj_energy_cycle = cmplx(0.0, 0.0, p)
                D0_population_cycle = cmplx(0.0, 0.0, p)

                !$omp do schedule(dynamic,200) reduction(+:D0_population_cycle,proj_energy_cycle, nattempts_spawn,ndeath)
                do iattempt = 1, selection_data%nclusters
                    ! For OpenMP scalability, have this test inside a single loop rather
                    ! than attempt to parallelise over three separate loops.
                    call select_trot_ucc_cluster(rng(it), sys, qs%psip_list, qs%ref%f0, qs%ref%max_ex_level, &
                                            selection_data%nstochastic_clusters, D0_normalisation, qmc_in%initiator_pop, D0_pos, &
                                            logging_info, contrib(it)%cdet, contrib(it)%cluster, qs%excit_gen_data)

                    if (uccmc_in%variational_energy .and. .not. all(contrib(it)%cdet%f==0) .and. contrib(it)%cluster%excitation_level <= qs%ref%ex_level)  then
                       state = contrib(it)%cdet%f 
                       call binary_search(time_avg_psip_list_ci%states, state, 1, time_avg_psip_list_ci%nstates, hit, pos)
                       population(1) = int(contrib(it)%cluster%amplitude*contrib(it)%cluster%cluster_to_det_sign*time_avg_psip_list_ci%pop_real_factor, int_p)
                       if(contrib(it)%cluster%nexcitors>0) then
                       end if
                       if (hit) then
                          time_avg_psip_list_ci%nparticles = time_avg_psip_list_ci%nparticles - abs(real(time_avg_psip_list_ci%pops(:,pos))/time_avg_psip_list_ci%pop_real_factor)
                          time_avg_psip_list_ci%pops(:,pos) = time_avg_psip_list_ci%pops(:,pos) + population 
                          time_avg_psip_list_ci%nparticles = time_avg_psip_list_ci%nparticles + abs(real(time_avg_psip_list_ci%pops(:,pos))/time_avg_psip_list_ci%pop_real_factor)
                       else
                           do j = time_avg_psip_list_ci%nstates, pos, -1
                               ! i is the number of determinants that will be inserted below j.
                               k = j + 1 
                               time_avg_psip_list_ci%states(:,k) = time_avg_psip_list_ci%states(:,j)
                               time_avg_psip_list_ci%pops(:,k) = time_avg_psip_list_ci%pops(:,j)
                               time_avg_psip_list_ci%dat(:,k) = time_avg_psip_list_ci%dat(:,j)
                           end do

                           call insert_new_walker(sys, time_avg_psip_list_ci, annihilation_flags, pos, &
                               contrib(it)%cdet%f, population, qs%ref)
                           ! Extract the real sign from the encoded sign.
                           time_avg_psip_list_ci%nstates = time_avg_psip_list_ci%nstates+1
                           real_population = real(time_avg_psip_list_ci%pops(:,pos))/time_avg_psip_list_ci%pop_real_factor
                           time_avg_psip_list_ci%nparticles = time_avg_psip_list_ci%nparticles + abs(real_population)

                       end if
                    end if

                    if (contrib(it)%cluster%excitation_level <= qs%ref%max_ex_level+2) then
                            ! cluster%excitation_level == huge(0) indicates a cluster
                            ! where two excitors share an elementary operator
                        if (qs%propagator%quasi_newton) contrib(it)%cdet%fock_sum = &
                                            sum_sp_eigenvalues_occ_list(sys, contrib(it)%cdet%occ_list) - qs%ref%fock_sum

                        call do_trot_uccmc_accumulation(sys, qs, contrib(it)%cdet, contrib(it)%cluster, logging_info, &
                                                    D0_population_cycle, proj_energy_cycle, uccmc_in, ref_det, rdm, selection_data)

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

                if(uccmc_in%variational_energy) then
                          time_avg_psip_list_ci%nparticles = time_avg_psip_list_ci%nparticles/(iter)
                          time_avg_psip_list_ci%pops(:,:time_avg_psip_list_ci%nstates) =  time_avg_psip_list_ci%pops(:,:time_avg_psip_list_ci%nstates)/(iter)
                end if
                do i = 1, qs%psip_list%nstates
                    state = qs%psip_list%states(:,i) 
                    call binary_search_trot(time_avg_psip_list%states, state, 1, time_avg_psip_list%nstates, hit, pos, qs%ref%f0, sys%nel, sys%basis)
                    if (hit) then
                          time_avg_psip_list%nparticles = time_avg_psip_list%nparticles - abs(real(time_avg_psip_list%pops(:,pos))/time_avg_psip_list%pop_real_factor)
                          time_avg_psip_list%pops(:,pos) = time_avg_psip_list%pops(:,pos) + qs%psip_list%pops(:,i) 
                          time_avg_psip_list%nparticles = time_avg_psip_list%nparticles + abs(real(time_avg_psip_list%pops(:,pos))/time_avg_psip_list%pop_real_factor)
                       else
                           do j = time_avg_psip_list%nstates, pos, -1
                               ! i is the number of determinants that will be inserted below j.
                               k = j + 1 
                               time_avg_psip_list%states(:,k) = time_avg_psip_list%states(:,j)
                               time_avg_psip_list%pops(:,k) = time_avg_psip_list%pops(:,j)
                               time_avg_psip_list%dat(:,k) = time_avg_psip_list%dat(:,j)
                           end do
                           call insert_new_walker(sys, time_avg_psip_list, annihilation_flags, pos, &
                               state, qs%psip_list%pops(:,i), qs%ref)
                           time_avg_psip_list%nstates = time_avg_psip_list%nstates+1
                           ! Extract the real sign from the encoded sign.
                           real_population = real(time_avg_psip_list%pops(:,pos))/time_avg_psip_list%pop_real_factor
                           time_avg_psip_list%nparticles = time_avg_psip_list%nparticles + abs(real_population)
                       end if
                end do

                time_avg_psip_list%nparticles = time_avg_psip_list%nparticles/(iter)
                time_avg_psip_list%pops(:,:time_avg_psip_list%nstates) =  time_avg_psip_list%pops(:,:time_avg_psip_list%nstates)/(iter)

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

        print*, 'end coeffs'
        print*, time_avg_psip_list%states(1,:qs%psip_list%nstates)
        print*, time_avg_psip_list%pops(1,:qs%psip_list%nstates)
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
            call var_energy_uccmc(sys, time_avg_psip_list_ci, var_energy)
            !print*, 'Variational energy: ', var_energy
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

        if(uccmc_in%variational_energy) then 
             deallocate(state, stat=ierr)
             call check_deallocate('state', ierr)
             deallocate(real_population, stat=ierr)
             call check_deallocate('real_population', ierr)
             deallocate(population, stat=ierr)
             call check_deallocate('population', ierr)
        end if 
    
    end subroutine do_trot_uccmc

    subroutine trot_ucc_set_cluster_selections(selection_data, nattempts, min_cluster_size)

        ! Function to set total number of selections of different cluster
        ! types within CCMC. This effectively controls the relative sampling
        ! of different clusters within the CC expansion.

        ! In:
        !   max_size: maximum cluster size we need to select.
        !   nstates: total number of occupied states within calculation.
        !   D0_normalisation: total population on the reference this iteration.
        !   tot_abs_pop: sum of absolute excip populations on all excitors.
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
        selection_data%nstochastic_clusters = nattempts
        selection_data%nsingle_excitors = 0

        selection_data%nclusters = selection_data%nD0_select + selection_data%nsingle_excitors &
                            + selection_data%nstochastic_clusters

    end subroutine trot_ucc_set_cluster_selections

    subroutine select_trot_ucc_cluster(rng, sys, psip_list, f0, ex_level, nattempts, normalisation, &
                              initiator_pop, D0_pos, logging_info, cdet, cluster, excit_gen_data)

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
        real(p) :: psize
        complex(p) :: cluster_population, excitor_pop
        integer :: i, ierr
        logical :: allowed
       
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

        ! Assume cluster is allowed unless collapse_cluster finds out otherwise
        ! when collapsing/combining excitors or if it could never have been
        ! valid
        allowed = .true. 
        do i = 1, psip_list%nstates
            if (.not.(i == D0_pos)) then
                rand = get_rand_close_open(rng)
                psize = real(abs(psip_list%pops(1,i)),p)/real((abs(psip_list%pops(1,D0_pos))+abs(psip_list%pops(1,i))),p)
                if (rand < psize) then
                   cluster%nexcitors = cluster%nexcitors + 1
                   cluster%pselect = cluster%pselect*psize
                   excitor_pop = real(psip_list%pops(1,i),p)/psip_list%pop_real_factor
                   if (cluster%nexcitors == 1) then
                       ! First excitor 'seeds' the cluster:
                       cdet%f = psip_list%states(:,i)
                       cdet%data => psip_list%dat(:,i) ! Only use if cluster is non-composite!
                       cluster_population = excitor_pop
                       ! Counter the additional *nprocs above.
                       cluster%pselect = cluster%pselect/nprocs
                   else
                       call trot_ucc_collapse_cluster(sys%basis, f0, psip_list%states(:,i), excitor_pop, cdet%f, &
                                          cluster_population, allowed)
                       if (.not. allowed) then
                           cluster%excitation_level = huge(0)
                           exit
                       end if
                    ! Each excitor spends the same amount of time on each processor on
                    ! average.  If this excitor is different from the previous excitor,
                    ! then the probability this excitor is on the same processor as the
                    ! previous excitor is 1/nprocs.  (Note choosing the same excitor
                    ! multiple times is valid in linked CC.)
                       cluster%pselect = cluster%pselect/nprocs
                   end if
                   ! If the excitor's population is below the initiator threshold, we remove the
                   ! initiator status for the cluster
                   if (abs(excitor_pop) <= initiator_pop) cdet%initiator_flag = 3
                   cluster%excitors(cluster%nexcitors)%f => psip_list%states(:,i)
                else
                   cluster%pselect = cluster%pselect*(1-psize)
                end if
            end if
        end do

        if (cluster%nexcitors == 0) then
            call create_null_cluster(sys, f0, cluster%pselect, normalisation, initiator_pop, &
                                    cdet, cluster, excit_gen_data)
        else 

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
               !cluster%pselect = cluster%pselect*factorial(cluster%nexcitors)

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
        end if

       !if (debug) call write_logging_stoch_selection(logging_info, cluster%nexcitors, cluster%excitation_level, pop, &
        !        max_size, cluster%pselect, cluster%amplitude, allowed)

    end subroutine select_trot_ucc_cluster

    pure subroutine trot_ucc_collapse_cluster(basis, f0, excitor, excitor_population, cluster_excitor, cluster_population, &
                    allowed)

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

        integer(i0) :: excitor_loc(basis%tot_string_len)

        integer :: ibasis, ibit
        integer(i0) :: excitor_excitation(basis%tot_string_len)
        integer(i0) :: excitor_annihilation(basis%tot_string_len)
        integer(i0) :: excitor_creation(basis%tot_string_len)
        integer(i0) :: cluster_excitation(basis%tot_string_len)
        integer(i0) :: cluster_annihilation(basis%tot_string_len)
        integer(i0) :: cluster_creation(basis%tot_string_len)
        integer(i0) :: permute_operators(basis%tot_string_len)
        excitor_loc = excitor
        call reset_extra_info_bit_string(basis, excitor_loc)

        ! Apply excitor to the cluster of excitors.

        ! orbitals involved in excitation from reference
        excitor_excitation = ieor(f0, excitor_loc)
        cluster_excitation = ieor(f0, cluster_excitor)
        ! annihilation operators (relative to the reference)
        excitor_annihilation = iand(excitor_excitation, f0)
        cluster_annihilation = iand(cluster_excitation, f0)
        ! creation operators (relative to the reference)
        excitor_creation = iand(excitor_excitation, excitor_loc)
        cluster_creation = iand(cluster_excitation, cluster_excitor)

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

    end subroutine trot_ucc_collapse_cluster

    subroutine do_trot_uccmc_accumulation(sys, qs, cdet, cluster, logging_info, D0_population_cycle, proj_energy_cycle, &
                                    uccmc_in, ref_det, rdm, selection_data)



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
        use proc_pointers, only: update_proj_energy_ptr
        use replica_rdm, only: update_rdm
        use logging, only: logging_t

        type(sys_t), intent(in) :: sys
        type(qmc_state_t), intent(in) :: qs
        type(uccmc_in_t), intent(in) :: uccmc_in
        type(det_info_t), intent(in) :: cdet, ref_det
        type(cluster_t), intent(in) :: cluster
        type(logging_t), intent(in) :: logging_info
        complex(p), intent(inout) :: D0_population_cycle, proj_energy_cycle
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
            call update_proj_energy_trot(sys, qs%ref%f0, qs%trial%wfn_dat, cdet, &
                     [real(cluster%amplitude,p),aimag(cluster%amplitude)]*&
                     cluster%cluster_to_det_sign/cluster%pselect, &
                     estimators_cycle, connection, hmatel)
            if (sys%read_in%comp) then
                D0_population_cycle = D0_population_cycle + estimators_cycle%D0_population_comp
                proj_energy_cycle = proj_energy_cycle + estimators_cycle%proj_energy_comp
            else
                D0_population_cycle = D0_population_cycle + estimators_cycle%D0_population
                proj_energy_cycle = proj_energy_cycle + estimators_cycle%proj_energy
            end if
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

    end subroutine do_trot_uccmc_accumulation

    subroutine do_stochastic_uccmc_propagation(rng, sys, qs, uccmc_in, &
                                            logging_info, ms_stats, bloom_stats, &
                                            contrib, nattempts_spawn_tot, ndeath, ps_stat)

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
        integer(int_p) :: nspawned, nspawned_im
        integer(i0) :: fexcit(sys%basis%tot_string_len)
 
        if (contrib%cluster%excitation_level == huge(0)) then
            ! When sampling e^-T H e^T, the cluster operators in e^-T
            ! and e^T can excite to/from the same orbital, requiring
            ! a different spawning routine
            call linked_spawner_ccmc(rng, sys, qs, qs%spawn_store%spawn%cutoff, &
                      contrib%cluster, gen_excit_ptr, nspawned, connection, nspawnings_total, &
                      fexcit, contrib%cdet, contrib%ldet, contrib%rdet, contrib%left_cluster, contrib%right_cluster, ps_stat)
            nspawned_im = 0_int_p
        !else if (sys%read_in%comp) then
        !    call spawner_complex_ccmc(rng, sys, qs, qs%spawn_store%spawn%cutoff, &
        !              uccmc_in%linked, contrib%cdet, contrib%cluster, gen_excit_ptr, logging_info,  nspawned, nspawned_im, &
        !              connection, nspawnings_total, ps_stat)
        else
            call spawner_ccmc(rng, sys, qs, qs%spawn_store%spawn%cutoff, &
                      uccmc_in%linked, contrib%cdet, contrib%cluster, gen_excit_ptr, logging_info, &
                      nspawned, connection, nspawnings_total, ps_stat)
            nspawned_im = 0_int_p
        end if
        if (nspawned /= 0_int_p) call create_spawned_particle_ccmc(sys%basis, qs%ref, contrib%cdet, connection, &
                                            nspawned, 1, contrib%cluster%excitation_level, &
                                            .false., fexcit, qs%spawn_store%spawn, bloom_stats)
        !if (nspawned_im /= 0_int_p) call create_spawned_particle_ccmc(sys%basis, qs%ref, contrib%cdet, connection,&
        !                                    nspawned_im, 2, contrib%cluster%excitation_level, &
        !                                    .false., fexcit, qs%spawn_store%spawn, bloom_stats)

    end subroutine perform_uccmc_spawning_attempt

    subroutine update_proj_energy_mol_ucc(sys, f0, wfn_dat, cdet, pop, estimators, excitation, hmatel, cluster_size)

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

    subroutine var_energy_uccmc(sys, time_avg_psip_list, var_energy)
         
       use excitations, only: excit_t, get_excitation
       use hamiltonian, only: get_hmatel
       use energy_evaluation, only: hmatel_t
       use system, only: sys_t
       use qmc_data, only: particle_t
       use read_in_symmetry, only: cross_product_basis_read_in
       use determinants, only: decode_det

       type(sys_t), intent(in) :: sys
       type(particle_t), intent(in) :: time_avg_psip_list
       real(p), intent(out) :: var_energy(:)
       real(p) :: normalisation(time_avg_psip_list%nspaces)

       type(excit_t) :: excitation
       type(hmatel_t) :: hmatel
       integer :: occ_list(sys%nel)

       integer :: i, j
       integer :: ij_sym, ab_sym

       normalisation = 0.0_p
       var_energy(:) = 0.0_p
       do i = 1, time_avg_psip_list%nstates
           normalisation = normalisation +(real(time_avg_psip_list%pops(:,i))/real(time_avg_psip_list%pops(:,1)))**2
           
           do j = 1, time_avg_psip_list%nstates
               hmatel = get_hmatel(sys, time_avg_psip_list%states(:,i), time_avg_psip_list%states(:,j))
               var_energy = var_energy + hmatel%r*real(time_avg_psip_list%pops(:,i))/real(time_avg_psip_list%pops(:,1))*real(time_avg_psip_list%pops(:,j))/real(time_avg_psip_list%pops(:,1))
           end do
           
       end do
       var_energy = var_energy/normalisation
   end subroutine var_energy_uccmc

   pure function earliest_unset(f, nel, basis) result (early)
    
       use basis_types, only: basis_t

       type(basis_t), intent(in) :: basis
       integer(i0), intent(in) :: f(basis%tot_string_len)
       integer, intent(in) :: nel
       integer :: i, early
     
       early = nel
       do i = 0, nel-1
           if (.not. btest(f(1),i)) then
               early = i
               exit
           end if
       end do
   end function

   pure function count_unset(f, f0, basis) result (nunset)

   use basis_types, only: basis_t
   use bit_utils, only: count_set_bits

   type(basis_t), intent(in) :: basis
   integer(i0), intent(in) :: f(basis%tot_string_len), f0(basis%tot_string_len)
   integer(i0) :: annihilation(basis%tot_string_len)
   integer :: nunset(basis%bit_string_len)

   annihilation = iand(ieor(f(:basis%bit_string_len),f0(:basis%bit_string_len)),f0(:basis%bit_string_len))
   nunset = count_set_bits(annihilation)
    
   end function

   subroutine binary_search_trot(list, item, istart, iend, hit, pos, f0, nel, basis)

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
        use bit_utils, only: bit_str_cmp
        use basis_types, only: basis_t

        integer(i0), intent(in) :: list(:,:), item(:)
        integer, intent(in) :: istart, iend
        integer(i0), intent(in) :: f0(:)
        integer, intent(in) :: nel
        type(basis_t), intent(in) :: basis
        logical, intent(out) :: hit
        integer, intent(out) :: pos

        integer :: hi, lo, compare
        integer :: pos_nunset(basis%bit_string_len), item_nunset(basis%bit_string_len), item_earliest_unset, pos_earliest_unset

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

            item_earliest_unset = earliest_unset(item, nel, basis)
            item_nunset = count_unset(item, f0, basis)
            
            if (item_earliest_unset == nel) then
                hi = 1
                lo = 1
            end if

            do while (hi /= lo)
                ! Narrow the search range down in steps.

                ! Mid-point.
                ! We shift one of the search limits to be the mid-point.
                ! The successive dividing the search range by 2 gives a O[log N]
                ! search algorithm.
                pos = (hi+lo)/2

                pos_earliest_unset = earliest_unset(list(:,pos),nel,basis)
                pos_nunset = count_unset(list(:,pos), f0, basis)
               
                if (pos_earliest_unset > item_earliest_unset) then
                    lo = pos + 1
                else if (pos_earliest_unset < item_earliest_unset) then
                    hi = pos
                else
                    if (pos_nunset(1) > item_nunset(1)) then
                        lo = pos + 1
                    else if (pos_nunset(1) < item_nunset(1)) then
                        hi = pos
                    else
                        compare = bit_str_cmp(list(:,pos), item)
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
                    end if
                end if
            end do

            ! If hi == lo, then we have narrowed the search down to one position but
            ! not checked if that position is the item we're hunting for.
            ! Because list can expand (i.e. we might be searching for an
            ! element which doesn't exist yet) the binary search can find either
            ! the element before or after where item should be placed.
            if (hi == lo) then
                pos_earliest_unset = earliest_unset(list(:,hi),nel,basis)
                pos_nunset = count_unset(list(:,hi), f0, basis)
               
                if (pos_earliest_unset > item_earliest_unset) then
                    pos = hi + 1
                else if (pos_earliest_unset < item_earliest_unset) then
                    pos = hi
                else
                    if (pos_nunset(1) > item_nunset(1)) then
                        pos = hi + 1
                    else if (pos_nunset(1) < item_nunset(1)) then
                        pos = hi
                    else
                        compare = bit_str_cmp(list(:,hi), item)
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
            end if

        end if

    end subroutine binary_search_trot

    subroutine direct_annihilation_trot(sys, rng, reference, annihilation_flags, psip_list, spawn, &
                                   nspawn_events, determ)

        ! Annihilation algorithm.
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

        ! If performing a semi-stochastic calculation then the annihilation
        ! process is slightly different, so call the correct routines depending
        ! on the situation.
        call annihilate_wrapper_spawn_t(spawn, annihilation_flags%initiator_approx)
        call annihilate_main_list_wrapper_trot(sys, rng, reference, annihilation_flags, psip_list, spawn)

    end subroutine direct_annihilation_trot

    subroutine annihilate_main_list_wrapper_trot(sys, rng, reference, annihilation_flags, psip_list, spawn, &
                                            lower_bound, determ_flags)

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

        else

            ! No spawned walkers so we only have to check to see if death has
            ! killed the entire population on a determinant.
            call remove_unoccupied_dets(rng, psip_list, annihilation_flags%real_amplitudes, determ_flags)

        end if

    end subroutine annihilate_main_list_wrapper_trot

    subroutine annihilate_main_list_trot(psip_list, spawn, sys, ref, lower_bound)

        ! Annihilate particles in the main walker list with those in the spawned
        ! walker list.

        ! In:
        !    tensor_label_len: number of elements in the bit array describing the position
        !       of the particle in the space (i.e.  determinant label in vector/pair of
        !       determinants label in array).
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
            call binary_search_trot(psip_list%states, f, istart, iend, hit, pos, ref%f0, sys%nel, sys%basis)
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
                !istart = pos + 1
            else
                ! Compress spawned list.
                k = i - nannihilate
                spawn%sdata(:,k) = spawn%sdata(:,i)
            end if
        end do

        spawn%head(thread_id,0) = spawn%head(thread_id,0) - nannihilate

    end subroutine annihilate_main_list_trot

    subroutine insert_new_walkers_trot(sys, psip_list, ref, annihilation_flags, spawn, determ_flags, lower_bound)

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
                call binary_search_trot(psip_list%states, int(spawn%sdata(:sys%basis%tensor_label_len,i), i0), &
                                   istart, iend, hit, pos, ref%f0, sys%nel, sys%basis)
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

    subroutine update_proj_energy_trot(sys, f0, wfn_dat, cdet, pop, estimators, excitation, hmatel)

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
        use energy_evaluation, only: hmatel_t, estimators_t

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

    end subroutine update_proj_energy_trot
end module

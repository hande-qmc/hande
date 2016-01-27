module ccmc

! Module  for performing coupled cluster Monte Carlo (CCMC) calculations).

! Due to the similarities with FCIQMC, we can re-use lots of the same routines
! (especially the spawning, death and annihilation).  As a result, the structure
! of do_ccmc is remarkably similar to the other do_*mc routines.

! NOTE:
!     FCIQMC routines work entirely in a determinant framework.  This is
!     somewhat inconsistent with coupled cluster which is most naturally defined
!     in terms of excitors.  Thus when considering the interface between CC
!     routines and FCI routines, one should mentally translate 'Slater
!     determinant' into 'the Slater determinant formed by applying the excitor
!     to the reference determinant'.

! Numerous discussions with Alex Thom (AJWT) and access to drafts/preprints of
! AJWT's CCMC papers gratefully acknowledged.  The following comments are some
! useful implementation notes that fill in some technical details not (yet)
! explicitly described in the CCMC papers.

! CCMC
! ====
!
! In the following, a capital index refers to combined index (ie a single label
! for a determinant or excitor) and a lower case index referis to a spin
! orbital.  Hence |D_I> can be equivalently represented as |D_{ij...}^{ab..}.
!
! Analogue with FCIQMC
! --------------------
!
! In FCIQMC we essentially use first-order finite-difference approximation to the
! imaginary-time Schroedinger equation::
!
!     c_I(t+dt) = c_I(t) - < D_I | H - S | \Psi(t) > dt  (1)
!
! where::
!
!     \Psi(t) = \sum_J c_J(t) | D_J >.  (2)
!
! CCMC involves propogating a very similar equation::
!
!     < D_I | a_I D_0 > t_I(t+dt) = < D_I | a_I D_0 > t_I(t) - < D_I | H - S | \Psi_{CC}(t) > dt  (3)
!
! where::
!
!     \Psi_{CC}(t) = e^{T(t)} | D_0 >,  (4)
!
! and the cluster operator, T, is::
!
!     T(t) = \sum_{ia} t_i^a(t) a_i^a + 1/2!^2 \sum_{ijab} t_{ij}^{ab}(t) a_{ij}^{ab} + ...  (5)
!
! (I apologise for using a in two different senses...)
!
! The operators, {a_I=a_{ij...}^{ab...}}, known as excitors, cause electrons to
! be excited from occupied orbitals in the reference determinant to virtual
! orbitals and have corresponding amplitudes, {t_I}, which evolve in time.  In
! the infinite-time limit, the amplitudes converge onto the true CC
! wavefunction.
!
! Note that, depending upon the definition of a_I, < D_I | a_I D_0 > can be
! either +1 or -1.  Alex Thom defines a_I such that this term is always +1 in his
! CCMC papers but this definition is not convenient for computation.  See
! comments in collapse_cluster and convert_excitor_to_determinant.
!
! We can only handle uniquely defined excitors, i.e. we consider t_{ij}^{ab},
! t_{ji}^{ba}, t_{ji}^{ab} and t_{ij}^{ba} to be one entity.  This is fine as
! t_{ij}^{ab} a_{ij}^{ab} = t_{ji}^{ba} a_{ji}^{ba} = -t_{ji}^{ab} = -t_{ij}^{ba}.
! As a consequence, we must be very careful with the above expansion in (5).  We
! can simplify (5) such that each unique cluster only occurs once::
!
!     T(t) = \sum_{ia} t_i^a(t) a_i^a + \sum_{i<j,a<b} t_{ij}^{ab}(t) a_{ij}^{ab} + ...  (6)
!
! Evaluating e^{T(t)} is, as in traditional CC, performed using a series
! expansion::
!
!     e^{T(t)} = e^{\sum_I t_I a_I}                                           (7)
!              = 1 + sum_I t_I a_I + 1/2! sum_{IJ} t_I t_J a_I a_J + ...      (8)
!
! Alternatively, as pointed out by Trygve Helgaker, one can consider::
!
!     e^{T(t)} = e^{t_I(t) a_I} e^{t_J(t) a_J} e^{t_K(t) a_K} ...             (9)
!              = (1+t_I(t) a_I) (1+t_J(t) a_J) (1+t_K(t) a_K) ...             (10)
!
! which expands out to (8), except for the absence of the factorial prefactors,
! and uses the property that a_I a_I = 0.  However, there is no computational
! advantage to using this and, due to the need to consider an order for {I}, it
! is actually trickier to code than using (8).
!
! Normalisation
! -------------
!
! The convention that the overlap with the reference is set to unity is rather
! hard to maintain in CCMC.  Instead we follow the prescription of Alex Thom and
! introduce an additional variable, N_0, to act as a normalisation constant::
!
!     |\Psi_{CC}> = N_0 e^{T/N_0} | D_0 >    (11)
!
! The number of particles on the reference (not strictly excips though we shall
! refer to them by this name for the sake of conciseness and generality) is
! hence N_0 and can be varied stochastically to maintain the required
! (relative) normalisation.
!
! Dynamics
! --------
!
! Dynamics in FCIQMC involve stochastically sampling the action of the
! Hamiltonian.  As the projection onto the CC wavefunction in (3) involves many
! terms (higher-order 'clusters' of excitors), CCMC also requires one to
! stochastically sample the CC wavefunction itself.  This results in a process
! similar (but more complicated) than the FCIQMC dynamics.  The similarities are
! such that we can reuse much of the FCIQMC routines, including the excitation
! generators and annihilation procedures.
!
! The CC wavefunction can be sampled by selected a random cluster of excitors.
! Applying this cluster to the reference determinant generates a determinant,
! D_I.  The action of the Hamiltonian is sampled in an identical fashion to
! FCIQMC:
!
! #. Generate D_J, a random (single or double) excitation of D_I.  Create a new
!    particle on D_J with a probability proportional to |H_{IJ}| dt.
! #. Create a particle on D_I with a probability proportional to |H_{II}-S| dt.
!
! As with FCIQMC, because we do not allow every possible event to take place each
! iteration, the probabilities must be weighted accordingly by the probability
! that the event was selected.
!
! After we have sampled the CC wavefunction and sampled its evolution the desired
! number of times, we then annihilate all particles (excips) on the same excitor
! with opposite signs.
!
! There are two main differences between FCIQMC and CCMC evolution.  As the FCI
! wavefunction is a linear combination of determinants, we can simply loop over
! all particles (psips) on the determinants and allow them to spawn and die.  As
! the CC wavefunction ansatz involves an exponentiation, we must consider
! combinations of excitors and hence combinations of excips.  The stochastic
! sampling of the wavefunction is achieved by selecting random clusters;
! hopefully the comments in select_cluster provide suitable illumination.
!
! The other main difference (and certainly the hardest to get right) is that the
! signs of the excips and their progeny are *far* more subtle and involved than
! in FCIQMC.  In addition to the explicit minus sign in the propogation equation
! (3), the other sources of sign changes are:
!
! #. the sign of the excip from which a new excip is spawned/cloned/killed (also
!    in FCIQMC).  Note that the amplitude of the excips on a cluster of excitors
!    is the product of the amplitudes of the excips on the individual excitors.
! #. the sign of the Hamiltonian matrix element (also in FCIQMC).
! #. combining excitors to form a 'cluster' excitor requires reordering of the
!    set of annihilation and creation operators.  Each permutation causes a sign
!    change and thus an overall odd number of permutations causes a sign change
!    in the excip population.  See collapse_cluster.
! #. Applying the I-th excitor to D_0 can result in a sign change due to the
!    different definitions of an excitor and determinant which also requires
!    permutations of creation and annihilation operators.  See
!    convert_excitor_to_determinant.
! #. As we are creating an excip on the J-th excitor from the I-th excitor, we
!    must be consistent with our definition of the excitor; in particular they
!    might produce determinants of different signs when applied to the reference
!    determinant.  Thus if we consider creating an excip on t_J from t_I (where
!    J can be I or any connected excitor), we *must* take into account the signs
!    of the excitors relative to the reference determinant.
!
! In order for the correct dynamics to be performed, we must carefully
! accumulate all negative signs and take them into account when
! determining the sign of the child excips.
!   
! Shift
! -----
!
! The shift is a free parameter introduced to provide control of the total population.
! The value of the shift should not affect the convergence of the wavefunction, only
! the overall growth.  Once the distribution of excips is (on average) correctly
! representing the wavefunction, this should remain true, regardless of the shift.
! But consider applying equation (3) above to the exact wavefunction.  The change in
! population from one iteration to the next is given by:
!
!   \Delta t_I = - dt < D_I | H - S | \Psi_{CC}(t) >                        (12)
!              = - dt (E_{CC} - S) < D_I | \Psi_{CC}(t) >
!
! This is not, in general, proportional to t_I, unless S = E_{CC}, so the amplitudes
! are not uniformly scaled and the new wavefunction is incorrect.  This can be
! remedied by instead using:
!
!   \Delta t_I = - dt < D_I | H - E | \Psi_{CC}(t) > - dt (E - S) t_I.      (13)
!
! As the exact energy E_{CC} is not known during the calculation, we use the projected
! energy estimator.  Once the shift has converged to E_{CC} (within stochastic 
! fluctuations), the equations (3) and (13) are equivalent and either may be used.
!
! This new expression can be thought of as separating the two roles of the dynamics: 
! the first term optimises the amplitudes to solve the coupled cluster equations, 
! and the second provides control of the total population.
!
! Two sorts of clusters may be distinguished: composite clusters, that contain a product
! of two or more excitors, and non-composite that have only one or no excitors.  There is
! no difference in the dynamics of non-composite clusters between the use of equations (3)
! and (13), but on non-composite clusters the projected energy rather than the shift is
! used for the diagonal death step.  Death on the composite clusters does not necessarily
! remove particles, as all the population resides on non-composite clusters, but
! potentially increases the population on a (possibly previously unoccupied) excitor, so
! is not effective for population control. From this perspective it makes sense to not
! have the shift involved in such steps.
!
! Linked CCMC
! ===========
!
! Instead of sampling the amplitude equations (13), the equivalent equations::
!
!   t_I(t+dt) = t_I(t) - dt (< D_I | e^{-T(t)} (H - E) e^{T(t)} | D_0 > - (E - S) t_I(t))  (14)
!
! may be used. Using the identity::
!
!   e^T H e^{-T} = H + [H,T]_c + 1/2 [[H,T],T]_c + 1/3! [[[H,T],T],T]_c + 1/4! [[[[H,T],T],T],T]_c,  (15)
!
! where the subscript c indicates that only terms in the commutators coming from
! linked diagrams need to be included, gives a very similar form to the original
! equations. These equations, however, only include terms of at most fourth order
! in T regardless of the truncation level. The equations (14) can be sampled in
! very much the same way as (13), but require some modifications due to the
! presence of the commutators instead of a simple product of operators:
!
! #. Clusters that include two excitors that excite from (to) the same orbital
!    give a contribution to the equations, in contrast to the original form
!    where they do not as the product of the excitors is 0. This corresponds to
!    one excitor coming from the e^T and the other from the e^{-T} in (14). See
!    comments in select_cluster and linked_spawner_ccmc for details.
! #. Excitations not linked to the cluster being spawned from can be rejected.
!    See linked_excitation.
! #. Matrix elements used for spawning and death probabilities can have more than
!    one term from the commutator contributing. See stochastic_ccmc_death and
!    unlinked_commutator.
!
! Clusters
! ========
!
! We use clusters in two slightly different but related senses.  In both cases the cluster
! contains a set of excitors in a specific order.  The difference really only arises when we
! come to evaluate them.
!
! unlinked CC
!   A cluster is the product of the excitors and arises from sampling (8).
! linked CC
!   A cluster is the set of excitors and arises from sampling (15).  As the selected set
!   of excitors arise from the sampling of a sum of commutators, we must evaluate the
!   commutator by considering all possible partitions, where a partition is a defined by
!   a subset of excitors before the Hamiltonian and the remaining subset after the
!   Hamiltonian.  As a concrete example, the cluster {t_i, t_j} has partitions H t_i t_j,
!   t_i H t_j, t_j H t_i and t_j t_i H (with appropriate signs).
!
!   Helpfully, we also refer to the product of excitors either side of the Hamiltonian as
!   clusters in procedures which act upon products (ie in the same sense as unlinked CC)
!   but do not require (subsequent) application of the Hamiltonian.
!
! Note that we cannot separate the sampling the action of the Hamiltonian and the sampling
! of the wavefunction in linked CC in the same way that we can in unlinked CC.

use const, only: i0, int_p, int_64, p

implicit none

contains

    subroutine do_ccmc(sys, qmc_in, ccmc_in, semi_stoch_in, restart_in, load_bal_in, reference_in, qs, qmc_state_restart)

        ! Run the CCMC algorithm starting from the initial walker distribution
        ! using the timestep algorithm.

        ! See notes about the implementation of this using function pointers
        ! in fciqmc_main.

        ! In:
        !    sys: system being studied.
        !    ccmc_in: input options relating to CCMC.
        !    semi_stoch_in: Input options for the semi-stochastic adaptation.
        !    restart_in: input options for HDF5 restart files.
        !    reference_in: current reference determinant.  If not set (ie
        !       components allocated) then a best guess is made based upon the
        !       desired spin/symmetry.
        !    load_bal_in: input options for load balancing.
        !    qmc_in: input options relating to QMC methods.
        ! In/Out:
        !    qmc_state_restart (optional): if present, restart from a previous fciqmc calculation.
        !       Deallocated on exit.
        ! Out:
        !    qs: qmc_state for use if restarting the calculation

        use checking, only: check_allocate, check_deallocate
        use dSFMT_interface, only: dSFMT_t, dSFMT_init
        use errors, only: stop_all
        use parallel
        use restart_hdf5, only: dump_restart_hdf5, restart_info_t, init_restart_info_t, dump_restart_file_wrapper

        use annihilation, only: direct_annihilation
        use bloom_handler, only: init_bloom_stats_t, bloom_stats_t, bloom_mode_fractionn, &
                                 accumulate_bloom_stats, write_bloom_report, bloom_stats_warning
        use ccmc_data
        use determinants, only: det_info_t, dealloc_det_info_t
        use excitations, only: excit_t, get_excitation_level, get_excitation
        use fciqmc_data, only: write_fciqmc_report, &
                               write_fciqmc_report_header
        use qmc, only: init_qmc
        use qmc_common, only: initial_fciqmc_status, cumulative_population, load_balancing_report, &
                              init_report_loop, init_mc_cycle, end_report_loop, end_mc_cycle,      &
                              redistribute_particles, rescale_tau
        use proc_pointers
        use spawning, only: assign_particle_processor
        use system, only: sys_t, sys_t_json
        use spawn_data, only: calc_events_spawn_t, write_memcheck_report

        use qmc_data, only: qmc_in_t, ccmc_in_t, semi_stoch_in_t, restart_in_t, reference_t
        use qmc_data, only: load_bal_in_t, qmc_state_t, annihilation_flags_t
        use qmc_data, only: qmc_in_t_json, ccmc_in_t_json, semi_stoch_in_t_json, restart_in_t_json, reference_t_json
        use check_input, only: check_qmc_opts, check_ccmc_opts
        use json_out, only: json_out_t, json_object_init, json_object_end

        type(sys_t), intent(in) :: sys
        type(qmc_in_t), intent(in) :: qmc_in
        type(ccmc_in_t), intent(in) :: ccmc_in
        type(semi_stoch_in_t), intent(in) :: semi_stoch_in
        type(restart_in_t), intent(in) :: restart_in
        type(load_bal_in_t), intent(in) :: load_bal_in
        type(reference_t), intent(in) :: reference_in
        type(qmc_state_t), target, intent(out) :: qs
        type(qmc_state_t), intent(inout), optional :: qmc_state_restart

        integer :: i, ireport, icycle, iter, semi_stoch_iter, it
        integer(int_64) :: iattempt, nattempts, nclusters, nstochastic_clusters, nsingle_excitors, nD0_select
        integer(int_64) :: nattempts_spawn
        real(p), allocatable :: nparticles_old(:), nparticles_change(:)
        type(det_info_t), allocatable :: cdet(:)
        type(det_info_t), allocatable :: ldet(:), rdet(:)

        integer(int_p) :: nspawned, ndeath
        integer :: nspawn_events, ierr
        type(excit_t) :: connection
        type(cluster_t), allocatable, target :: cluster(:)
        type(cluster_t), allocatable :: left_cluster(:), right_cluster(:)
        type(multispawn_stats_t), allocatable :: ms_stats(:)
        type(dSFMT_t), allocatable :: rng(:)
        real(p) :: junk, bloom_threshold
        type(json_out_t) :: js
        type(qmc_in_t) :: qmc_in_loc

        logical :: soft_exit, dump_restart_shift, restarting

        integer(int_p), allocatable :: cumulative_abs_nint_pops(:)
        integer :: D0_proc, D0_pos, nD0_proc, min_cluster_size, max_cluster_size, iexcip_pos, slot
        integer(int_p) :: tot_abs_nint_pop
        real(p) :: D0_normalisation
        type(bloom_stats_t) :: bloom_stats
        type(annihilation_flags_t) :: annihilation_flags
        type(restart_info_t) :: ri, ri_shift

        real :: t1, t2

        logical :: update_tau, error

        integer :: nspawnings_total

        integer(i0) :: fexcit(sys%basis%string_len)
        logical :: seen_D0
        real(p) :: D0_population_cycle, proj_energy_cycle, proj_energy_old

        if (parent) then
            write (6,'(1X,"CCMC")')
            write (6,'(1X,"----",/)')
        end if

        ! Check input options.
        if (parent) then
            restarting = present(qmc_state_restart) .or. restart_in%read_restart
            call check_qmc_opts(qmc_in, .not.present(qmc_state_restart), restarting)
            call check_ccmc_opts(sys, ccmc_in)
        end if

        ! Initialise data.
        call init_qmc(sys, qmc_in, restart_in, load_bal_in, reference_in, annihilation_flags, qs, &
                      qmc_state_restart=qmc_state_restart)

        if (parent) then
            call json_object_init(js, tag=.true.)
            call sys_t_json(js, sys)
            ! The default values of pattempt_* are not in qmc_in
            qmc_in_loc = qmc_in
            qmc_in_loc%pattempt_single = qs%pattempt_single
            qmc_in_loc%pattempt_double = qs%pattempt_double
            call qmc_in_t_json(js, qmc_in_loc)
            call ccmc_in_t_json(js, ccmc_in)
            call semi_stoch_in_t_json(js, semi_stoch_in)
            call restart_in_t_json(js, restart_in)
            call reference_t_json(js, qs%ref, sys, .true.)
            call json_object_end(js, terminal=.true., tag=.true.)
            write (js%io, '()')
        end if

        allocate(nparticles_old(qs%psip_list%nspaces), stat=ierr)
        call check_allocate('nparticles_old', qs%psip_list%nspaces, ierr)
        allocate(nparticles_change(qs%psip_list%nspaces), stat=ierr)
        call check_allocate('nparticles_change', qs%psip_list%nspaces, ierr)

        ! Initialise bloom_stats components to the following parameters.
        call init_bloom_stats_t(bloom_stats, mode=bloom_mode_fractionn, encoding_factor=qs%psip_list%pop_real_factor)

        if (qs%ref%ex_level+2 > 12 .and. .not. ccmc_in%linked) then
            call stop_all('do_ccmc', 'CCMC can currently only handle clusters up to size 12 due to&
                                     &integer overflow in factorial routines for larger clusters.  &
                                     &Please implement better factorial routines or use linked-CCMC.')
        end if

        ! Allocate and initialise per thread...
        allocate(rng(0:nthreads-1), stat=ierr)
        call check_allocate('rng', nthreads, ierr)
        allocate(ms_stats(0:nthreads-1), stat=ierr)
        call check_allocate('ms_stats', nthreads, ierr)

        if (ccmc_in%linked) then
            call init_cluster(sys, 4, cdet, cluster)
            call init_cluster(sys, 4, ldet, left_cluster)
            call init_cluster(sys, 4, rdet, right_cluster)
        else
            call init_cluster(sys, qs%ref%ex_level+2, cdet, cluster)
        end if

        do i = 0, nthreads-1
            ! Initialise and allocate RNG store.
            call dSFMT_init(qmc_in%seed+iproc+i*nprocs, 50000, rng(i))
        end do

        ! Whilst cluster data can be accessed from cdet, I recommend explicitly
        ! passing it as an argument rather than accessing cdet%cluster both for
        ! the sake of brevity and clarity.  In particular, I wish to encourage
        ! not using cdet%cluster in order to maintain (where possible and
        ! relevant) generality in routines applicable to FCIQMC and CCMC.
        do i = 0, nthreads-1
            cdet(i)%cluster => cluster(i)
        end do

        ! ...and scratch space for calculative cumulative probabilities.
        allocate(cumulative_abs_nint_pops(size(qs%psip_list%states,dim=2)), stat=ierr)
        call check_allocate('cumulative_abs_nint_pops', size(qs%psip_list%states, dim=2), ierr)

        nparticles_old = qs%psip_list%tot_nparticles

        ! Initialise D0_pos to be somewhere (anywhere) in the list.
        D0_pos = 1

        ! Main fciqmc loop.
        if (parent) call write_fciqmc_report_header(qs%psip_list%nspaces)
        call initial_fciqmc_status(sys, qmc_in, qs)
        ! Initialise timer.
        call cpu_time(t1)

        associate(spawn=>qs%spawn_store%spawn)
            ! Initialise hash shift if restarting...
            spawn%hash_shift = qs%mc_cycles_done
            ! Hard code how frequently (ie 2^10) a determinant can move.
            spawn%move_freq = ccmc_in%move_freq
        end associate

        ! The iteration on which to start performing semi-stochastic.
        semi_stoch_iter = qs%mc_cycles_done + semi_stoch_in%start_iter

        ! Should we dump a restart file just before the shift is turned on?
        dump_restart_shift = restart_in%write_restart_shift
        call init_restart_info_t(ri, write_id=restart_in%write_id)
        call init_restart_info_t(ri_shift, write_id=restart_in%write_shift_id)

        do ireport = 1, qmc_in%nreport

            ! Projected energy from last report loop to correct death
            proj_energy_old = qs%estimators%proj_energy/qs%estimators%D0_population

            call init_report_loop(qs, bloom_stats)

            do icycle = 1, qmc_in%ncycles

                iter = qs%mc_cycles_done + (ireport-1)*qmc_in%ncycles + icycle

                associate(spawn=>qs%spawn_store%spawn, pm=>qs%spawn_store%spawn%proc_map)
                    call assign_particle_processor(qs%ref%f0, spawn%bit_str_nbits, spawn%hash_seed, spawn%hash_shift, &
                                                   spawn%move_freq, nprocs, D0_proc, slot, pm%map, pm%nslots)

                    ! Update the shift of the excitor locations to be the end of this
                    ! current iteration.
                    spawn%hash_shift = spawn%hash_shift + 1
                end associate

                if (iproc == D0_proc) then

                    ! Population on reference determinant.
                    ! As we might select the reference determinant multiple times in
                    ! a cycle, the running total of D0_population is incorrect (by
                    ! a factor of the number of times it was selected).
                    call find_D0(qs%psip_list, qs%ref%f0, D0_pos)
                    D0_normalisation = real(qs%psip_list%pops(1,D0_pos),p)/qs%psip_list%pop_real_factor

                    nD0_proc = 1

                else

                    ! Can't find D0 on this processor.  (See how D0_pos is used
                    ! in select_cluster.)
                    D0_pos = -1
                    nD0_proc = 0 ! No reference excitor on the processor.

                end if

                if (ccmc_in%linked) then
                    ! The BCH expansion of the Hamiltonian terminates at fourth
                    ! order in T so at most four excitors needed in the cluster
                    if (qs%psip_list%nstates == nD0_proc) then
                        ! All excips are on the reference, so no possible clusters
                        ! In linked CCMC we can select the same excip multiple times.
                        max_cluster_size = 0
                    else
                        max_cluster_size = 4
                    end if
                else
                    ! Maximum possible cluster size that we can generate.
                    ! Usually this is either the number of electrons or the
                    ! truncation level + 2 but we must handle the case where we are
                    ! growing the initial population from a single/small number of
                    ! excitors.
                    ! Can't include the reference in the cluster, so -1 from the
                    ! total number of excitors.
                    max_cluster_size = min(sys%nel, qs%ref%ex_level+2, qs%psip_list%nstates-nD0_proc)
                end if

#ifdef PARALLEL
                call mpi_bcast(D0_normalisation, 1, mpi_preal, D0_proc, MPI_COMM_WORLD, ierr)
#endif

                ! Note that 'death' in CCMC creates particles in the spawned
                ! list, so the number of deaths not in the spawned list is
                ! always 0.
                call init_mc_cycle(qs%psip_list, qs%spawn_store%spawn, nattempts, ndeath, &
                                   min_attempts=nint(D0_normalisation,int_64))
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
                call cumulative_population(qs%psip_list%pops, qs%psip_list%nstates, D0_proc, D0_pos, qs%psip_list%pop_real_factor, &
                                           cumulative_abs_nint_pops, tot_abs_nint_pop)

                associate(bs=>bloom_stats, nstates_active=>qs%psip_list%nstates)
                    bloom_threshold = ceiling(max(nattempts, int(nstates_active,int_64))*bs%prop)*real(bs%encoding_factor,p)
                end associate

                ! Two options for evolution:

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
                if (ccmc_in%full_nc) then
                    ! Note that nattempts /= tot_abs_nint_pop+D0_normalisation if the
                    ! reference is not on the current processor.  Instead work
                    ! out how many clusters of each type we will sample
                    ! explicitly.
                    min_cluster_size = 2
                    nD0_select = nint(D0_normalisation)
                    nclusters = 2*tot_abs_nint_pop + nD0_select
                    nstochastic_clusters = tot_abs_nint_pop
                    nsingle_excitors = tot_abs_nint_pop
                else
                    min_cluster_size = 0
                    nclusters = nattempts
                    nD0_select = 0 ! instead of this number of deterministic selections, these are chosen stochastically
                    nstochastic_clusters = nattempts
                end if

                ! OpenMP chunk size determined completely empirically from a single
                ! test.  Please feel free to improve...
                ! NOTE: we can't refer to procedure pointers in shared blocks so
                ! can't use default(none).  I *strongly* recommend turning
                ! default(none) on when making changes and ensure that the only
                ! errors relate to the procedure pointers...
                !$omp parallel &
                ! --DEFAULT(NONE) DISABLED-- !$omp default(none) &
                !$omp private(it, iexcip_pos, nspawned, connection, junk,       &
                !$omp         nspawnings_total, fexcit, i,     &
                !$omp         seen_D0) &
                !$omp shared(nattempts, rng, cumulative_abs_nint_pops, tot_abs_nint_pop,  &
                !$omp        max_cluster_size, cdet, cluster, &
                !$omp        D0_normalisation, D0_pos, nD0_select, qs,               &
                !$omp        sys, bloom_threshold, bloom_stats,                      &
                !$omp        proj_energy_cycle, min_cluster_size,       &
                !$omp        nclusters, nstochastic_clusters, nattempts_spawn,       &
                !$omp        nsingle_excitors, ccmc_in, ldet, rdet, left_cluster,    &
                !$omp        right_cluster, nprocs, ms_stats, qmc_in, load_bal_in,   &
                !$omp        nparticles_change, ndeath, D0_population_cycle)
                it = get_thread_id()
                iexcip_pos = 0
                seen_D0 = .false.
                D0_population_cycle = 0.0_p
                proj_energy_cycle = 0.0_p
                !$omp do schedule(dynamic,200) reduction(+:D0_population_cycle,proj_energy_cycle,nattempts_spawn)
                do iattempt = 1, nclusters

                    ! For OpenMP scalability, have this test inside a single loop rather
                    ! than attempt to parallelise over three separate loops.
                    if (iattempt <= nstochastic_clusters) then
                        call select_cluster(rng(it), sys, qs%psip_list, qs%ref%f0, qs%ref%ex_level, ccmc_in%linked, &
                                            nstochastic_clusters, D0_normalisation, qmc_in%initiator_pop, D0_pos, &
                                            cumulative_abs_nint_pops, tot_abs_nint_pop, min_cluster_size, max_cluster_size, &
                                            cdet(it), cluster(it))
                    else if (iattempt <= nstochastic_clusters+nD0_select) then
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
                            call create_null_cluster(sys, qs%ref%f0, nprocs*real(nD0_select,p), D0_normalisation, &
                                                     qmc_in%initiator_pop, cdet(it), cluster(it))
                        end if
                    else
                        ! Deterministically select each excip as a non-composite cluster.
                        call select_cluster_non_composite(sys, qs%psip_list, qs%ref%f0, iattempt-nstochastic_clusters-nD0_select, &
                                                          iexcip_pos, nsingle_excitors, qmc_in%initiator_pop, D0_pos, &
                                                          cumulative_abs_nint_pops, tot_abs_nint_pop, cdet(it), cluster(it))
                    end if

                    if (cluster(it)%excitation_level <= qs%ref%ex_level+2 .or. &
                            (ccmc_in%linked .and. cluster(it)%excitation_level == huge(0))) then
                        ! cluster%excitation_level == huge(0) indicates a cluster
                        ! where two excitors share an elementary operator

                        if (cluster(it)%excitation_level /= huge(0)) then
                            ! FCIQMC calculates the projected energy exactly.  To do
                            ! so in CCMC would involve enumerating over all pairs of
                            ! single excitors, which is rather painful and slow.
                            ! Instead, as we are randomly sampling clusters in order
                            ! to evolve the excip population anyway, we can just use
                            ! the random clusters to *sample* the projected
                            ! estimator.  See comments in spawning.F90 for why we
                            ! must divide through by the probability of selecting
                            ! the cluster.
                            connection = get_excitation(sys%nel, sys%basis, cdet(it)%f, qs%ref%f0)
                            call update_proj_energy_ptr(sys, qs%ref%f0, qs%trial%wfn_dat, cdet(it), &
                                     cluster(it)%cluster_to_det_sign*cluster(it)%amplitude/cluster(it)%pselect, &
                                     D0_population_cycle, proj_energy_cycle, connection, junk)
                        end if

                        ! Spawning
                        ! This has the potential to create blooms, so we allow for multiple
                        ! spawning events per cluster.
                        ! The number of spawning events is decided by the value
                        ! of cluster%amplitude/cluster%pselect.  If this is
                        ! greater than cluster_multispawn_threshold, then nspawnings is
                        ! increased to the ratio of these.
                        nspawnings_total=max(1,ceiling( abs(cluster(it)%amplitude/cluster(it)%pselect)/ &
                                                         ccmc_in%cluster_multispawn_threshold))
                        call ms_stats_update(nspawnings_total, ms_stats(it))
                        nattempts_spawn = nattempts_spawn + nspawnings_total

                        do i = 1, nspawnings_total
                            if (cluster(it)%excitation_level == huge(0)) then
                                ! When sampling e^-T H e^T, the cluster operators in e^-T
                                ! and e^T can excite to/from the same orbital, requiring
                                ! a different spawning routine
                                call linked_spawner_ccmc(rng(it), sys, qmc_in, qs, qs%spawn_store%spawn%cutoff, &
                                          cluster(it), gen_excit_ptr, nspawned, connection, nspawnings_total, &
                                          fexcit, ldet(it), rdet(it), left_cluster(it), right_cluster(it))
                            else
                                call spawner_ccmc(rng(it), sys, qs, qs%spawn_store%spawn%cutoff, &
                                          ccmc_in%linked, cdet(it), cluster(it), gen_excit_ptr, nspawned, connection, &
                                          nspawnings_total)
                            end if

                           if (nspawned /= 0_int_p) then
                               if (cluster(it)%excitation_level == huge(0)) then
                                   call create_spawned_particle_ptr(sys%basis, qs%ref, cdet(it), connection, nspawned, &
                                                                    1, qs%spawn_store%spawn, fexcit)
                               else
                                   call create_spawned_particle_ptr(sys%basis, qs%ref, cdet(it), connection, nspawned, 1, &
                                                                    qs%spawn_store%spawn)
                               end if
                               if (abs(nspawned) > bloom_threshold) call accumulate_bloom_stats(bloom_stats, nspawned)
                           end if
                        end do

                        ! Does the cluster collapsed onto D0 produce
                        ! a determinant is in the truncation space?  If so, also
                        ! need to attempt a death/cloning step.
                        ! optimisation: call only once per iteration for clusters of size 0 or 1 for ccmc_in%full_nc.
                        if (cluster(it)%excitation_level <= qs%ref%ex_level) then
                            ! Clusters above size 2 can't die in linked ccmc.
                            if ((.not. ccmc_in%linked) .or. cluster(it)%nexcitors <= 2) then
                                ! Do death for non-composite clusters directly and in a separate loop
                                if (cluster(it)%nexcitors >= 2 .or. .not. ccmc_in%full_nc) then
                                    call stochastic_ccmc_death(rng(it), qs%spawn_store%spawn, ccmc_in%linked, sys, &
                                                               qs, cdet(it), cluster(it), proj_energy_old)
                                end if
                            end if
                        end if

                    end if

                end do
                !$omp end do

                if (ccmc_in%full_nc .and. qs%psip_list%nstates > 0) then
                    ! Do death exactly and directly for non-composite clusters
                    !$omp do schedule(dynamic,200) reduction(+:ndeath,nparticles_change)
                    do iattempt = 1, qs%psip_list%nstates
                        ! Note we use the (encoded) population directly in stochastic_ccmc_death_nc
                        ! (unlike the stochastic_ccmc_death) to avoid unnecessary decoding/encoding
                        ! steps (cf comments in stochastic_death for FCIQMC).
                        call stochastic_ccmc_death_nc(rng(it), ccmc_in%linked, qs, iattempt==D0_pos, &
                                              qs%psip_list%dat(1,iattempt), proj_energy_old, qs%psip_list%pops(1, iattempt), &
                                              nparticles_change(1), ndeath)
                    end do
                    !$omp end do
                end if
                !$omp end parallel

                qs%psip_list%nparticles = qs%psip_list%nparticles + nparticles_change
                qs%estimators%D0_population = qs%estimators%D0_population + D0_population_cycle
                qs%estimators%proj_energy = qs%estimators%proj_energy + proj_energy_cycle

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

                call end_mc_cycle(nspawn_events, ndeath, qs%psip_list%pop_real_factor, nattempts_spawn, qs%spawn_store%rspawn)

            end do

            update_tau = bloom_stats%nblooms_curr > 0

            error = qs%spawn_store%spawn%error .or. qs%psip_list%error

            call end_report_loop(qmc_in, iter, update_tau, qs, nparticles_old, nspawn_events, &
                                 semi_stoch_in%shift_iter, semi_stoch_iter, soft_exit, &
                                 load_bal_in, bloom_stats=bloom_stats, error=error)
            if (error) exit

            call cpu_time(t2)
            if (parent) then
                if (bloom_stats%nblooms_curr > 0) call bloom_stats_warning(bloom_stats)
                call write_fciqmc_report(qmc_in, qs, ireport, nparticles_old, t2-t1, .false., .false.)
            end if

            ! Update the time for the start of the next iteration.
            t1 = t2

            call dump_restart_file_wrapper(qs, dump_restart_shift, restart_in%write_freq, nparticles_old, ireport, &
                                           qmc_in%ncycles, sys%basis%nbasis, ri, ri_shift, .false.)

            qs%psip_list%tot_nparticles = nparticles_old

            if (soft_exit) exit

            if (update_tau) call rescale_tau(qs%tau)

        end do

        if (parent) write (6,'()')
        call write_bloom_report(bloom_stats)
        call multispawn_stats_report(ms_stats)
        call load_balancing_report(qs%psip_list%nparticles, qs%psip_list%nstates, qmc_in%use_mpi_barriers,&
                                   qs%spawn_store%spawn%mpi_time)
        call write_memcheck_report(qs%spawn_store%spawn)

        if (soft_exit .or. error) then
            qs%mc_cycles_done = qs%mc_cycles_done + qmc_in%ncycles*ireport
        else
            qs%mc_cycles_done = qs%mc_cycles_done + qmc_in%ncycles*qmc_in%nreport
        end if

        if (restart_in%write_restart) then
            call dump_restart_hdf5(ri, qs, qs%mc_cycles_done, nparticles_old, sys%basis%nbasis, .false.)
            if (parent) write (6,'()')
        end if

        ! TODO: deallocation...
!        call dealloc_det_info_t(cdet)
!        cdet%cluster => NULL()
!        deallocate(cluster%excitors, stat=ierr)
!        call check_deallocate('cluster%excitors', ierr)

    end subroutine do_ccmc

    subroutine init_cluster(sys, cluster_size, cdet, cluster)

        ! Allocates cdet and cluster, and their components.

        ! In:
        !    sys: system being studied
        !    cluster_size: the maximum number of excitors allowed in a cluster
        ! Out:
        !    cdet: Array of det_info_t variables, one for each thread, with
        !       components allocated
        !    cluster: Array of cluster_t variables, one for each thread, with
        !       components allocated

        use parallel, only: nthreads
        use determinants, only: det_info_t, alloc_det_info_t
        use ccmc_data, only: cluster_t
        use system, only: sys_t
        use checking, only: check_allocate

        type(sys_t), intent(in) :: sys
        integer, intent(in) :: cluster_size
        type(det_info_t), allocatable, intent(out) :: cdet(:)
        type(cluster_t), allocatable, intent(out) :: cluster(:)

        integer :: i, ierr

        ! Allocate arrays
        allocate(cdet(0:nthreads-1), stat=ierr)
        call check_allocate('cdet', nthreads, ierr)
        allocate(cluster(0:nthreads-1), stat=ierr)
        call check_allocate('cluster', nthreads, ierr)

        do i = 0, nthreads-1
            ! Allocate det_info_t and cluster_t components
            call alloc_det_info_t(sys, cdet(i))
            allocate(cluster(i)%excitors(cluster_size), stat=ierr)
            call check_allocate('cluster%excitors', cluster_size, ierr)
        end do

    end subroutine init_cluster

    subroutine find_D0(psip_list, f0, D0_pos)

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
            call binary_search(psip_list%states, f0, 1, psip_list%nstates, hit, D0_pos)
        else
            D0_pos_old = D0_pos
            select case(bit_str_cmp(f0, psip_list%states(:,D0_pos)))
            case(0)
                ! D0 hasn't moved.
                hit = .true.
            case(1)
                ! D0 < psip_list%states(:,D0_pos) -- it has moved to earlier in
                ! the list and the old D0_pos is an upper bound.
                call binary_search(psip_list%states, f0, 1, D0_pos_old, hit, D0_pos)
            case(-1)
                ! D0 > psip_list%states(:,D0_pos) -- it has moved to later in
                ! the list and the old D0_pos is a lower bound.
                call binary_search(psip_list%states, f0, D0_pos_old, psip_list%nstates, hit, D0_pos)
            end select
        end if
        if (.not.hit) call stop_all('find_D0', 'Cannot find reference!')

    end subroutine find_D0

    subroutine select_cluster(rng, sys, psip_list, f0, ex_level, linked_ccmc, nattempts, normalisation, &
                              initiator_pop, D0_pos, cumulative_excip_pop, tot_excip_pop, min_size, max_size, cdet, cluster)

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
        !    linked_ccmc: if true then only sample linked clusters.
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

        use determinants, only: det_info_t
        use ccmc_data, only: cluster_t
        use excitations, only: get_excitation_level
        use dSFMT_interface, only: dSFMT_t, get_rand_close_open
        use qmc_data, only: particle_t
        use proc_pointers, only: decoder_ptr
        use utils, only: factorial
        use search, only: binary_search
        use sort, only: insert_sort
        use parallel, only: nprocs
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        type(particle_t), intent(in), target :: psip_list
        integer(i0), intent(in) :: f0(sys%basis%string_len)
        integer, intent(in) :: ex_level
        integer(int_64), intent(in) :: nattempts
        logical, intent(in) :: linked_ccmc
        integer, intent(in) :: D0_pos
        real(p), intent(in) :: normalisation, initiator_pop
        integer(int_p), intent(in) :: cumulative_excip_pop(:), tot_excip_pop
        integer :: min_size, max_size
        type(dSFMT_t), intent(inout) :: rng
        type(det_info_t), intent(inout) :: cdet
        type(cluster_t), intent(inout) :: cluster

        real(p) :: rand, psize, cluster_population, excitor_pop
        integer :: i, pos, prev_pos
        integer(int_p) :: pop(max_size)
        logical :: hit, allowed, all_allowed

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
        cluster%pselect = nattempts*nprocs

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
        do i = 0, max_size-min_size-1
            psize = psize + 1.0_p/2**(i+1)
            if (rand < psize) then
                ! Found size!
                cluster%nexcitors = i+min_size
                cluster%pselect = cluster%pselect/2**(i+1)
                exit
            end if
        end do
        ! If not set, then must be the largest possible cluster
        if (cluster%nexcitors == -1) then
            cluster%nexcitors = max_size
            cluster%pselect = cluster%pselect*(1.0_p - psize)
        end if

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
        all_allowed = allowed

        select case(cluster%nexcitors)
        case(0)
            call create_null_cluster(sys, f0, cluster%pselect, normalisation, initiator_pop, cdet, cluster)
        case default
            ! Select cluster from the excitors on the current processor with
            ! probability for choosing an excitor proportional to the excip
            ! population on that excitor.  (For convenience, we use a probability
            ! proportional to the nint(pop), as it makes finding the right excitor
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
                pop(i) = int(get_rand_close_open(rng)*tot_excip_pop, int_p) + 1
            end do
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
                if (pos == D0_pos) pos = pos - 1
                excitor_pop = real(psip_list%pops(1,pos),p)/psip_list%pop_real_factor
                if (i == 1) then
                    ! First excitor 'seeds' the cluster:
                    cdet%f = psip_list%states(:,pos)
                    cdet%data => psip_list%dat(:,pos) ! Only use if cluster is non-composite!
                    cluster_population = excitor_pop
                    ! Counter the additional *nprocs above.
                    cluster%pselect = cluster%pselect/nprocs
                else
                    call collapse_cluster(sys%basis, f0, psip_list%states(:,pos), excitor_pop, cdet%f, &
                                          cluster_population, allowed)
                    if (.not.allowed) then
                        if (.not. linked_ccmc) exit
                        all_allowed = .false.
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
                if (abs(excitor_pop) <= initiator_pop) cdet%initiator_flag = 1
                ! Probability of choosing this excitor = nint(pop)/tot_pop.
                cluster%pselect = (cluster%pselect*nint(abs(excitor_pop)))/tot_excip_pop
                cluster%excitors(i)%f => psip_list%states(:,pos)
                prev_pos = pos
            end do

            if (allowed) cluster%excitation_level = get_excitation_level(f0, cdet%f)
            ! To contribute the cluster must be within a double excitation of
            ! the maximum excitation included in the CC wavefunction.
            if (cluster%excitation_level > ex_level+2) allowed = .false.

            if (allowed.or.linked_ccmc) then
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
                call decoder_ptr(sys, cdet%f, cdet)

                ! Normalisation factor for cluster%amplitudes...
                cluster%amplitude = cluster_population/(real(normalisation,p)**(cluster%nexcitors-1))
            else
                ! Simply set excitation level to a too high (fake) level to avoid
                ! this cluster being used.
                cluster%excitation_level = huge(0)
            end if

            if (.not.all_allowed) cluster%excitation_level = huge(0)

        end select

    end subroutine select_cluster

    subroutine create_null_cluster(sys, f0, prob, D0_normalisation, initiator_pop, cdet, cluster)

        ! Create a cluster with no excitors in it, and set it to have
        ! probability of generation prob.

        ! In:
        !    sys: system being studied
        !    f0: bit string of the reference
        !    prob: The probability we set in it of having been generated
        !    D0_normalisation:  The number of excips at the reference (which
        !        will become the amplitude of this cluster)
        !    initiator_pop: the population above which a determinant is an initiator.

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
        use determinants, only: det_info_t
        use ccmc_data, only: cluster_t
        use proc_pointers, only: decoder_ptr

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f0(sys%basis%string_len)
        real(p), intent(in) :: prob, D0_normalisation, initiator_pop
        type(det_info_t), intent(inout) :: cdet
        type(cluster_t), intent(inout) :: cluster

        ! Note only one null cluster to choose => p_clust = 1.
        cluster%pselect = prob

        cluster%nexcitors = 0

        ! Initiator approximation.
        ! This is sufficiently quick that we'll just do it in all cases, even
        ! when not using the initiator approximation.  This matches the approach
        ! used by Alex Thom in 'Initiator Stochastic Coupled Cluster Theory'
        ! (unpublished).
        ! Surely the reference has an initiator population?
        cdet%initiator_flag = 0

        ! Must be the reference.
        cdet%f = f0
        cluster%excitation_level = 0
        cluster%amplitude = D0_normalisation
        cluster%cluster_to_det_sign = 1
        if (cluster%amplitude <= initiator_pop) then
             ! Something has gone seriously wrong and the CC
             ! approximation is (most likely) not suitably for this system.
             ! Let the user be an idiot if they want to be...
             cdet%initiator_flag = 1
        end if

        call decoder_ptr(sys, cdet%f, cdet)

    end subroutine create_null_cluster

    subroutine select_cluster_non_composite(sys, psip_list, f0, iexcip, iexcip_pos, nattempts, initiator_pop,  D0_pos, &
                                            cumulative_excip_pop, tot_excip_pop, cdet, cluster)

        ! Select (deterministically) the non-composite cluster containing only
        ! the single excitor iexcitor and set the same information as select_cluster.

        ! In:
        !    sys: system being studied
        !    psip_list: particle_t object containing current excip distribution on
        !       this processor.
        !    f0: bit string of the reference
        !    iexcip: the index (in range [1,tot_excip_pop]) of the excip to select.
        !    nattempts: the number of times (on this processor) a random cluster
        !        of excitors is generated in the current timestep.
        !    initiator_pop: the population above which a determinant is an initiator.
        !    D0_pos: position in the excip list of the reference.
        !    cumulative_excip_population: running cumulative excip population on
        !        all excitors; i.e. cumulative_excip_population(i) = sum(particle_t%pops(1:i)).
        !    tot_excip_pop: total excip population.

        ! NOTE: cumulative_excip_pop and tot_excip_pop ignore the population on the
        ! reference as excips on the reference cannot form a cluster and the rounds the
        ! population on all other excitors to the nearest integer (for convenience--see
        ! comments in do_ccmc).  Both these quantities should be generated by
        ! cumulative_population (or be in the same format).

        ! In/Out:
        !    iexcip_pos: on output position of iexcip in the
        !        cumulative_excip_pop list.  Set to 0 on the initial call and use
        !        the previous return value (or a smaller number) on subsequent
        !        calls.  WARNING: we assume that this is a minimum value for the
        !        position of iexcip (hence loop over excips in order or reset
        !        iexcip pos each time).
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
        use determinants, only: det_info_t
        use ccmc_data, only: cluster_t
        use excitations, only: get_excitation_level
        use qmc_data, only: particle_t
        use search, only: binary_search
        use proc_pointers, only: decoder_ptr

        type(sys_t), intent(in) :: sys
        type(particle_t), intent(in), target :: psip_list
        integer(i0), intent(in) :: f0(sys%basis%string_len)
        integer(int_64), intent(in) :: iexcip, nattempts
        integer, intent(inout) :: iexcip_pos
        integer, intent(in) :: D0_pos
        real(p), intent(in) :: initiator_pop
        integer(int_p), intent(in) :: cumulative_excip_pop(:), tot_excip_pop
        type(det_info_t), intent(inout) :: cdet
        type(cluster_t), intent(inout) :: cluster
        real(p) :: excitor_pop

        integer :: iexcip_last
        ! It is more convenient to find the excitor on which the iexcip-th excip resides
        ! rather than looping over all excips explicitly (as in fciqmc) as it enables us
        ! to use the same control loop for all CCMC spawning which is then simpler and
        ! more performant for OpenMP parallelisation.

        ! Note that whilst we select each excip in turn, we actually set the amplitude to
        ! the population of excips on that excitor and the selection probability to the
        ! ratio of the population and the total population.  These factors cancel out in
        ! the spawning attempt (see spawner_ccmc) but doing so means that cluster is set
        ! here is an identical fashion to select_cluster.

        iexcip_last = iexcip_pos
        if (iexcip_pos == 0) iexcip_pos = 1

        ! Most of the time the excip is either on the current position or the
        ! next one, so special case to avoid the loop overhead.
        if (cumulative_excip_pop(iexcip_pos) >= iexcip) then
            ! Do nothing---already on the right position
        else if (cumulative_excip_pop(iexcip_pos+1) >= iexcip) then
            ! In the next slot...
            iexcip_pos = iexcip_pos + 1
        else
            ! Need to hunt for it (ie the reference position is in the way).
            iexcip_pos = iexcip_pos + 1
            do
                if (cumulative_excip_pop(iexcip_pos) >= iexcip .and. cumulative_excip_pop(iexcip_pos-1) < iexcip) exit
                iexcip_pos = iexcip_pos + 1
            end do
        end if
        ! Adjust for reference---cumulative_excip_pop(D0_pos) = cumulative_excip_pop(D0_pos-1).
        if (iexcip_pos == D0_pos) iexcip_pos = iexcip_pos - 1

        ! cdet and cluster only need to be set the first time the cluster is selected. On subsequent
        ! spawning attempts the same values can be reused, saving calls to decode_det_* and
        ! convert_excitor_to_determinant
        if (iexcip_last /= iexcip_pos) then
            ! We shall accumulate the factors which comprise cluster%pselect as we go.
            !   cluster%pselect = n_sel p_size p_clust
            ! where
            !   n_sel   is the number of cluster selections made (by this
            !           processor);
            !   p_size  is the probability of choosing a cluster of that size (1 in this case);
            !   p_clust is the probability of choosing a specific cluster given
            !           the choice of size.

            ! This excitor can only be selected on this processor and only one excitor is
            ! selected in the cluster, so unlike selecting the reference or composite
            ! clusters, there are no additional factors of nprocs or 1/nprocs to include.
            cluster%pselect = nattempts

            cluster%nexcitors = 1

            ! Initiator approximation.
            ! This is sufficiently quick that we'll just do it in all cases, even
            ! when not using the initiator approximation.  This matches the approach
            ! used by Alex Thom in 'Initiator Stochastic Coupled Cluster Theory'
            ! (unpublished).
            ! Assume all excitors in the cluster are initiators (initiator_flag=0)
            ! until proven otherwise (initiator_flag=1).
            cdet%initiator_flag = 0

            cdet%f = psip_list%states(:,iexcip_pos)
            cdet%data => psip_list%dat(:,iexcip_pos)
            cluster%excitors(1)%f => psip_list%states(:,iexcip_pos)
            excitor_pop = real(psip_list%pops(1,iexcip_pos),p)/psip_list%pop_real_factor
            if (abs(excitor_pop) <= initiator_pop) cdet%initiator_flag = 1
            ! pclust = |population|/total_population, as just a single excitor in the cluster..
            cluster%pselect = (cluster%pselect*nint(abs(excitor_pop)))/tot_excip_pop
            cluster%excitation_level = get_excitation_level(f0, cdet%f)
            cluster%amplitude = excitor_pop

            ! Sign change due to difference between determinant
            ! representation and excitors and excitation level.
            call convert_excitor_to_determinant(cdet%f, cluster%excitation_level, cluster%cluster_to_det_sign, f0)
            call decoder_ptr(sys, cdet%f, cdet)
        end if

    end subroutine select_cluster_non_composite

    subroutine spawner_ccmc(rng, sys, qs, spawn_cutoff, linked_ccmc, cdet, cluster, &
                            gen_excit_ptr, nspawn, connection, nspawnings_total)

        ! Attempt to spawn a new particle on a connected excitor with
        ! probability
        !     \tau |<D'|H|D_s> A_s|
        !   -------------------------
        !   n_sel p_s p_clust p_excit
        ! where |D_s> is the determinant formed by applying the excitor to the
        ! reference determinant, A_s is the amplitude and D' is the determinant
        ! formed from applying a connected excitor to the reference determinant.
        ! See comments in select_cluster about n_sel, p_s and p_clust.  p_excit
        ! is the probability of choosing D' given D_s.

        ! This is just a thin wrapper around a system-specific excitation
        ! generator and a utility function.  We need to modify the spawning
        ! probability compared to the FCIQMC algorithm as we spawn from multiple
        ! excips at once (in FCIQMC we allow each psip to spawn individually)
        ! and have additional probabilities to take into account.

        ! This routine will only attempt one spawning event, but needs to know
        ! the total number attempted for this cluster which is passed into
        ! nspawnings_total.

        ! In linked CCMC, the probability of spawning is modified by changing
        ! the matrix elements <D'|H|D_s> so that only linked diagrams are used
        ! for spawning.

        ! In:
        !    sys: system being studied.
        !    qs: qmc_state_t object. The timestep and reference determinant are used.
        !    spawn_cutoff: The size of the minimum spawning event allowed, in
        !        the encoded representation. Events smaller than this will be
        !        stochastically rounded up to this value or down to zero.
        !    linked_ccmc: if true then only sample linked clusters.
        !    cdet: info on the current excitor (cdet) that we will spawn
        !        from.
        !    cluster: information about the cluster which forms the excitor.  In
        !        particular, we use the amplitude, cluster_to_det_sign and pselect
        !        (i.e. n_sel.p_s.p_clust) attributes in addition to any used in
        !        the excitation generator.
        !    gen_excit_ptr: procedure pointer to excitation generators.
        !        gen_excit_ptr%full *must* be set to a procedure which generates
        !        a complete excitation.
        ! In/Out:
        !    rng: random number generator.
        !    nspawnings_total: The total number of spawnings attemped by the current cluster
        !        in the current timestep.
        ! Out:
        !    nspawn: number of particles spawned, in the encoded representation.
        !        0 indicates the spawning attempt was unsuccessful.
        !    connection: excitation connection between the current excitor
        !        and the child excitor, on which progeny are spawned.

        use ccmc_data, only: cluster_t
        use determinants, only: det_info_t
        use dSFMT_interface, only: dSFMT_t
        use excitations, only: excit_t, create_excited_det, get_excitation_level
        use proc_pointers, only: gen_excit_ptr_t
        use spawning, only: attempt_to_spawn
        use system, only: sys_t
        use const, only: depsilon
        use qmc_data, only: qmc_in_t, qmc_state_t

        type(sys_t), intent(in) :: sys
        type(qmc_state_t), intent(in) :: qs
        integer(int_p), intent(in) :: spawn_cutoff
        logical, intent(in) :: linked_ccmc
        type(det_info_t), intent(in) :: cdet
        type(cluster_t), intent(in) :: cluster
        type(dSFMT_t), intent(inout) :: rng
        integer, intent(in) :: nspawnings_total
        type(gen_excit_ptr_t), intent(in) :: gen_excit_ptr
        integer(int_p), intent(out) :: nspawn
        type(excit_t), intent(out) :: connection

        ! We incorporate the sign of the amplitude into the Hamiltonian matrix
        ! element, so we 'pretend' to attempt_to_spawn that all excips are
        ! actually spawned by positive excips.
        integer(int_p), parameter :: parent_sign = 1_int_p
        real(p) :: hmatel, pgen
        integer(i0) :: fexcit(sys%basis%string_len), funlinked(sys%basis%string_len)
        integer :: excitor_sign, excitor_level
        logical :: linked, single_unlinked, allowed_excitation

        ! 1. Generate random excitation.
        ! Note CCMC is not (yet, if ever) compatible with the 'split' excitation
        ! generators of the sys%lattice%lattice models.  It is trivial to implement and (at
        ! least for now) is left as an exercise to the interested reader.
        call gen_excit_ptr%full(rng, sys, qs%pattempt_single, cdet, pgen, connection, hmatel, allowed_excitation)

        if (linked_ccmc .and. allowed_excitation) then
            ! For Linked Coupled Cluster we reject any spawning where the
            ! Hamiltonian is not linked to every cluster operator
            ! The matrix element to be evaluated is not <D_j|H a_i|D0> but <D_j|[H,a_i]|D0>
            ! (and similarly for composite clusters)
            if (cluster%nexcitors > 0) then
                ! Check whether this is an unlinked diagram - if so, the matrix element is 0 and
                ! no spawning is attempted
                call linked_excitation(sys%basis, qs%ref%f0, connection, cluster, linked, single_unlinked, funlinked)
                if (.not. linked) then
                    hmatel = 0.0_p
                else if (single_unlinked) then
                    ! Single excitation: need to modify the matrix element
                    ! Subtract off the matrix element from the cluster without
                    ! the unlinked a_i operator
                    hmatel = hmatel - unlinked_commutator(sys, qs%ref%f0, connection, cluster, cdet%f, funlinked)
                end if
            end if
        end if

        ! 2, Apply additional factors.
        hmatel = hmatel*cluster%amplitude*cluster%cluster_to_det_sign
        pgen = pgen*cluster%pselect*nspawnings_total

        ! 3. Attempt spawning.
        nspawn = attempt_to_spawn(rng, qs%tau, spawn_cutoff, qs%psip_list%pop_real_factor, hmatel, pgen, parent_sign)

        if (nspawn /= 0_int_p) then
            ! 4. Convert the random excitation from a determinant into an
            ! excitor.  This might incur a sign change and hence result in
            ! a change in sign to the sign of the progeny.
            ! This is the same process as excitor to determinant and hence we
            ! can reuse code...
            call create_excited_det(sys%basis, cdet%f, connection, fexcit)
            excitor_level = get_excitation_level(qs%ref%f0, fexcit)
            call convert_excitor_to_determinant(fexcit, excitor_level, excitor_sign, qs%ref%f0)
            if (excitor_sign < 0) nspawn = -nspawn
        end if

    end subroutine spawner_ccmc

    subroutine stochastic_ccmc_death(rng, spawn, linked_ccmc, sys, qs, cdet, cluster, proj_energy)

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

        ! In:
        !    sys: system being studied.
        !    qs: qmc_state_t containing information about the reference and estimators.
        !    linked_ccmc: if true then only sample linked clusters.
        !    cdet: info on the current excitor (cdet) that we will spawn
        !        from.
        !    cluster: information about the cluster which forms the excitor.
        !    proj_energy: projected energy.  This should be the average value from the last
        !        report loop, not the running total in qs%estimators.
        ! In/Out:
        !    rng: random number generator.
        !    spawn: spawn_t object to which the spanwed particle will be added.

        use ccmc_data, only: cluster_t
        use determinants, only: det_info_t
        use excitations, only: excit_t
        use proc_pointers, only: sc0_ptr, create_spawned_particle_ptr
        use dSFMT_interface, only: dSFMT_t, get_rand_close_open
        use spawn_data, only: spawn_t
        use system, only: sys_t
        use qmc_data, only: qmc_state_t

        type(sys_t), intent(in) :: sys
        type(qmc_state_t), intent(in) :: qs
        logical, intent(in) :: linked_ccmc
        type(det_info_t), intent(in) :: cdet
        type(cluster_t), intent(in) :: cluster
        real(p), intent(in) :: proj_energy
        type(dSFMT_t), intent(inout) :: rng
        type(spawn_t), intent(inout) :: spawn

        real(p) :: pdeath, KiiAi
        integer(int_p) :: nkill
        type(excit_t), parameter :: null_excit = excit_t( 0, [0,0], [0,0], .false.)

        ! Spawning onto the same excitor so no change in sign due to
        ! a difference in the sign of the determinant formed from applying the
        ! parent excitor to the qs%ref and that formed from applying the
        ! child excitor.
        if (linked_ccmc) then
            select case (cluster%nexcitors)
            case(0)
                ! Death on the reference is unchanged
                KiiAi = (-qs%shift(1))*cluster%amplitude
            case(1)
                ! Evaluating the commutator gives
                ! <D1|[H,a1]|D0> = <D1|H|D1> - <D0|H|D0>
                KiiAi = (cdet%data(1) + proj_energy - qs%shift(1))*cluster%amplitude
            case(2)
                ! Evaluate the commutator
                ! The cluster operators are a1 and a2 (with a1 D0 = D1, a2 D0 = D2,
                ! a1 a2 D0 = D3) so the commutator gives:
                ! <D3|[[H,a1],a2]|D0> = <D3|H|D3> - <D2|H|D2> - <D1|H|D1> + <D0|H|D0>
                KiiAi = (sc0_ptr(sys, cdet%f) - sc0_ptr(sys, cluster%excitors(1)%f) &
                    - sc0_ptr(sys, cluster%excitors(2)%f) + qs%ref%H00)*cluster%amplitude
            case default
                ! At most two cluster operators can be linked to the diagonal
                ! part of H so this must be an unlinked cluster
                KiiAi = 0.0_p
            end select
        else
            select case (cluster%nexcitors)
            case(0)
                KiiAi = (-qs%shift(1))*cluster%amplitude
            case(1)
                KiiAi = (cdet%data(1) - qs%shift(1))*cluster%amplitude
            case default
                KiiAi = (sc0_ptr(sys, cdet%f) - qs%ref%H00 - proj_energy)*cluster%amplitude
            end select
        end if

        ! Amplitude is the decoded value.  Scale here so death is performed exactly (bar precision).
        ! See comments in stochastic_death.
        KiiAi = qs%psip_list%pop_real_factor*KiiAi
        pdeath = qs%tau*abs(KiiAi)/cluster%pselect

        if (pdeath < spawn%cutoff) then
            ! Calling death once per excip (and hence with a low pselect) without any
            ! stochastic rounding leads to a large number of excips being spawned with low
            ! weight (with the death then performed during annihilation).  This is not
            ! good for performance of the communication and an annihilation algorithms, so
            ! stochastically round here to overcome this.  Unlike in FCIQMC, death in CCMC
            ! can spawn new particles on basis functions that are not yet occupied so
            ! treating it on the same footing as death is not the worst idea...
            ! Note death is performed exactly for non-composite clusters with the non-composite
            ! algorithm in stochastic_ccmc_death_nc, as in FCIQMC.
            if (pdeath > get_rand_close_open(rng)*spawn%cutoff) then
                nkill = spawn%cutoff
            else
                nkill = 0_int_p
            end if
        else
            ! Number that will definitely die
            nkill = int(pdeath,int_p)
            ! Stochastic death...
            pdeath = pdeath - nkill
            if (pdeath > get_rand_close_open(rng)) nkill = nkill + 1
        end if

        if (nkill /= 0) then
            ! Create nkill excips with sign of -K_ii A_i
            if (KiiAi > 0) nkill = -nkill
!            cdet%initiator_flag=0  !All death is allowed
            ! The excitor might be a composite cluster so we'll just create
            ! excips in the spawned list and allow the annihilation process to take
            ! care of the rest.
            ! Pass through a null excitation so that we create a spawned particle on
            ! the current excitor.
            call create_spawned_particle_ptr(sys%basis, qs%ref, cdet, null_excit, nkill, 1, spawn)
        end if

    end subroutine stochastic_ccmc_death

    subroutine stochastic_ccmc_death_nc(rng, linked_ccmc, qs, D0, Hii, proj_energy, population, tot_population, ndeath)

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
        !    D0: true if the current excip is the null (reference) excitor
        !    Hii: the diagonal matrix element of the determinant formed by applying the excip to the
        !       reference
        !    proj_energy: projected energy.  This should be the average value from the last
        !        report loop, not the running total in qs%estimators.
        ! In/Out:
        !    rng: random number generator.
        !    ndeath: running (encoded) total of number of particles killed/cloned.
        !    population: the (encoded) population on the current excip
        !    tot_population: total number of particles.

        use dSFMT_interface, only: dSFMT_t, get_rand_close_open
        use qmc_data, only: qmc_state_t

        logical, intent(in) :: linked_ccmc
        type(qmc_state_t), intent(in) :: qs
        logical, intent(in) :: D0
        real(p), intent(in) :: Hii, proj_energy
        type(dSFMT_t), intent(inout) :: rng
        integer(int_p), intent(inout) :: population, ndeath
        real(p), intent(inout) :: tot_population

        real(p) :: pdeath, KiiAi
        integer(int_p) :: nkill, old_pop

        ! Spawning onto the same excitor so no change in sign due to
        ! a difference in the sign of the determinant formed from applying the
        ! parent excitor to the reference and that formed from applying the
        ! child excitor.

        if (D0) then
            KiiAi = (-qs%shift(1))*population
        else
            if (linked_ccmc) then
                KiiAi = (Hii + proj_energy - qs%shift(1))*population
            else
                KiiAi = (Hii - qs%shift(1))*population
            end if
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

    end subroutine stochastic_ccmc_death_nc

    pure subroutine collapse_cluster(basis, f0, excitor, excitor_population, cluster_excitor, cluster_population, allowed)

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

        use basis_types, only: basis_t

        use bit_utils, only: count_set_bits
        use const, only: i0_end

        type(basis_t), intent(in) :: basis
        integer(i0), intent(in) :: f0(basis%string_len)
        integer(i0), intent(in) :: excitor(basis%string_len)
        real(p), intent(in) :: excitor_population
        integer(i0), intent(inout) :: cluster_excitor(basis%string_len)
        real(p), intent(inout) :: cluster_population
        logical,  intent(out) :: allowed

        integer :: ibasis, ibit
        integer(i0) :: excitor_excitation(basis%string_len)
        integer(i0) :: excitor_annihilation(basis%string_len)
        integer(i0) :: excitor_creation(basis%string_len)
        integer(i0) :: cluster_excitation(basis%string_len)
        integer(i0) :: cluster_annihilation(basis%string_len)
        integer(i0) :: cluster_creation(basis%string_len)
        integer(i0) :: permute_operators(basis%string_len)

        ! Apply excitor to the cluster of excitors.

        ! orbitals involved in excitation from reference
        excitor_excitation = ieor(f0, excitor)
        cluster_excitation = ieor(f0, cluster_excitor)
        ! annihilation operators (relative to the reference)
        excitor_annihilation = iand(excitor_excitation, f0)
        cluster_annihilation = iand(cluster_excitation, f0)
        ! creation operators (relative to the reference)
        excitor_creation = iand(excitor_excitation, excitor)
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

            do ibasis = 1, basis%string_len
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

    end subroutine collapse_cluster

    pure subroutine convert_excitor_to_determinant(excitor, excitor_level, excitor_sign, f)

        ! We usually consider an excitor as a bit string representation of the
        ! determinant formed by applying the excitor (a group of annihilation
        ! and creation operators) to the reference determinant; indeed the
        ! determinant form is required when constructing Hamiltonian matrix
        ! elements.  However, the resulting determinant might actually contain
        ! a sign change, which needs to be absorbed into the (signed) population
        ! of excips on the excitor.

        ! This results from the fact that a determinant, |D>, is defined as:
        !   |D> = a^+_i a^+_j ... a^+_k |0>,
        ! where |0> is the vacuum, a^+_i creates an electron in the i-th
        ! orbital, i<j<...<k and |0> is the vacuum.  An excitor is defined as
        !   t_{ij...k}^{ab...c} = a^+_a a^+_b ... a^+_c a_k ... a_j a_i
        ! where i<j<...<k and a<b<...<c.  (This definition is somewhat
        ! arbitrary; the key thing is to be consistent.)  Hence applying an
        ! excitor to the reference might result in a change of sign, i.e.
        ! t_{ij}^{ab} |D_0> = -|D_{ij}^{ab}>.  As a more concrete example,
        ! consider a set of spin states (as the extension to fermions is
        ! irrelevant to the argument), with the reference:
        !   |D_0> = | 1 2 3 >
        ! and the excitor
        !   t_{13}^{58}
        ! Thus, using |0> to denote the vacuum:
        !   t_{13}^{58} | 1 2 3 > = + a^+_5 a^+_8 a_3 a_1 a^+_1 a^+_2 a^+_3 |0>
        !                         = + a^+_5 a^+_8 a_3 a^+_2 a^+_3 |0>
        !                         = - a^+_5 a^+_8 a_3 a^+_3 a^+_2 |0>
        !                         = - a^+_5 a^+_8 a^+_2 |0>
        !                         = + a^+_5 a^+_2 a^+_8 |0>
        !                         = - a^+_2 a^+_5 a^+_8 |0>
        !                         = - | 2 5 8 >
        ! Similarly
        !   t_{12}^{58} | 1 2 3 > = + a^+_5 a^+_8 a_2 a_1 a^+_1 a^+_2 a^+_3 |0>
        !                         = + a^+_5 a^+_8 a_2 a^+_2 a^+_3 |0>
        !                         = + a^+_5 a^+_8 a^+_3 |0>
        !                         = - a^+_5 a^+_3 a^+_8 |0>
        !                         = + a^+_3 a^+_5 a^+_8 |0>
        !                         = + | 3 5 8 >

        ! This potential sign change must be taken into account; we do so by
        ! absorbing the sign into the signed population of excips on the
        ! excitor.

        ! Essentially taken from Alex Thom's original implementation.

        ! In:
        !    excitor: bit string of the Slater determinant formed by applying
        !        the excitor to the reference determinant.
        !    excitor_level: excitation level, relative to the determinant f,
        !        of the excitor.  Equal to the number of
        !        annihilation (or indeed creation) operators in the excitor.
        !    f: bit string of determinant to which excitor is
        !       applied to generate a new determinant.
        ! Out:
        !    excitor_sign: sign due to applying the excitor to the
        !       determinant f to form a Slater determinant, i.e. < D_i | a_i D_f >,
        !       which is +1 or -1, where D_i is the determinant formed from
        !       applying the cluster of excitors, a_i, to D_f

        use const, only: i0_end

        integer(i0), intent(in) :: excitor(:)
        integer, intent(in) :: excitor_level
        integer, intent(inout) :: excitor_sign
        integer(i0), intent(in) :: f(:)

        integer(i0) :: excitation(size(excitor))
        integer :: ibasis, ibit, ncreation, nannihilation

        ! Bits involved in the excitation from the reference determinant.
        excitation = ieor(f, excitor)

        nannihilation = excitor_level
        ncreation = excitor_level

        excitor_sign = 1

        ! Obtain sign change by (implicitly) constructing the determinant formed
        ! by applying the excitor to the reference determinant.
        do ibasis = 1, size(excitor)
            do ibit = 0, i0_end
                if (btest(f(ibasis),ibit)) then
                    ! Occupied orbital in reference.
                    if (btest(excitation(ibasis),ibit)) then
                        ! Orbital excited from...annihilate electron.
                        ! This amounts to one fewer operator in the cluster through
                        ! which the other creation operators in the determinant must
                        ! permute.
                        nannihilation = nannihilation - 1
                    else
                        ! Orbital is occupied in the reference and once the
                        ! excitor has been applied.
                        ! Permute the corresponding creation operator through
                        ! the remaining creation and annihilation operators of
                        ! the excitor (which operate on orbitals with a higher
                        ! index than the current orbital).
                        ! If the permutation is odd, then we incur a sign
                        ! change.
                        if (mod(nannihilation+ncreation,2) == 1) &
                            excitor_sign = -excitor_sign
                    end if
                else if (btest(excitation(ibasis),ibit)) then
                    ! Orbital excited to...create electron.
                    ! This amounts to one fewer operator in the cluster through
                    ! which the creation operators in the determinant must
                    ! permute.
                    ! Note due to definition of the excitor, it is guaranteed
                    ! that this is created in the correct place, ie there are no
                    ! other operators in the excitor it needs to be interchanged
                    ! with.
                    ncreation = ncreation - 1
                end if
            end do
        end do

    end subroutine convert_excitor_to_determinant

    pure subroutine linked_excitation(basis, f0, connection, cluster, linked, single_unlinked, excitor)

        ! For Linked Coupled Cluster, the only terms of H that need to
        ! be sampled are those which are connected to (ie have a
        ! creation/annihilation operator in common with) each excitor in the
        ! cluster (by Wick's Theorem). This routine tests which excitors are
        ! connected to the Hamiltonian.

        ! In:
        !    basis: information about the single-particle basis
        !    f0: bit string of the reference
        !    connection: the excitation connecting the current excitor and the
        !       child excitor
        !    cluster: the cluster of excitation operators
        ! Out:
        !    linked: true if the Hamiltonian is connected to the cluster
        !    single_unlinked: true if connection is a single excitation and
        !       cluster contains one excitor that does not share an orbital with
        !       connection. (It is instead connected to the Hamiltonian by the
        !       dummy index in the two-body term \sum_j <aj||ij> a^+_a a^+_j a_j a_i)
        !       Any other excitors in the cluster must share an orbital with connection.
        !    excitor: if single_unlinked is true, this is the bit string of the
        !       excitor not linked to connection (otherwise 0)

        use excitations, only: excit_t, create_excited_det
        use basis_types, only: basis_t
        use ccmc_data, only: cluster_t

        type(basis_t), intent(in) :: basis
        integer(i0), intent(in) :: f0(basis%string_len)
        type(excit_t), intent(in) :: connection
        type(cluster_t), intent(in) :: cluster
        logical, intent(out) :: linked, single_unlinked
        integer(i0), intent(out) :: excitor(basis%string_len)

        integer :: i, orb, bit_pos, bit_element, unconnected
        integer(i0) :: excitor_excitation(basis%string_len)
        integer(i0) :: h_excitation(basis%string_len)

        single_unlinked = .false.
        unconnected = 0
        excitor = 0

        ! Get bit string of orbitals in H excitation
        ! (modified from create_excited_det)
        h_excitation = 0
        do i = 1, connection%nexcit
            ! i/j orbital
            orb = connection%from_orb(i)
            bit_pos = basis%bit_lookup(1,orb)
            bit_element = basis%bit_lookup(2,orb)
            h_excitation(bit_element) = ibset(h_excitation(bit_element), bit_pos)
            ! a/b orbital
            orb = connection%to_orb(i)
            bit_pos = basis%bit_lookup(1,orb)
            bit_element = basis%bit_lookup(2,orb)
            h_excitation(bit_element) = ibset(h_excitation(bit_element), bit_pos)
        end do

        do i = 1, cluster%nexcitors
            ! check each cluster operator shares an index with connection
            ! orbitals involved in cluster operator excitation (from reference)
            excitor_excitation = ieor(cluster%excitors(i)%f, f0)
            if (all(iand(h_excitation, excitor_excitation) == 0)) then
                ! no orbitals in common between H and cluster
                unconnected = unconnected + 1
                excitor = cluster%excitors(i)%f
            end if
        end do

        select case(connection%nexcit)
        case(1)
            ! For a single excitation, H contains the sum
            ! \sum_j <ij||aj>a^+j^+ji so can connect to one excitor that has no
            ! orbitals in common with the excitation
            linked = (unconnected <= 1)
            single_unlinked = (unconnected == 1)
        case(2)
            ! Double excitation, H only has the term <ab||ij>a^+b^+ji so must
            ! have an orbital in common with all cluster operators
            linked = (unconnected == 0)
        end select

        if (.not. linked) excitor = 0

    end subroutine linked_excitation

    pure function unlinked_commutator(sys, f0, connection, cluster, cdet, funlinked) result(hmatel)

        ! When H is a single excitation, it can be connected to one of the
        ! excitors in a cluster (a1) by the creation / annihilation operators
        ! corresponding to the dummy index in the two-body term
        !   \sum_j <aj||ij> a^+ j^+ i j.
        ! In that case, both terms of the commutator <D_j|[H,a1]|D_i> (where D_i
        ! is the determinant obtained by applying the other excitors in the
        ! cluster to the reference) are non-zero. <D_j|H a1|D_i> has already been
        ! calculated in the excitation generator; this function evaluates <D_j|a1 H|D_i>.

        ! In:
        !    sys: the system being studied
        !    f0: bit string of the reference
        !    connection: excitation connection between the current excitor
        !        and the child excitor, on which progeny are spawned.
        !    cluster: the cluster of excitors
        !    cdet: the determinant formed by the cluster
        !    funlinked: the excitor in the cluster that is not linked to H
        ! Returns:
        !    <Dj|a H|Di> (to be subtracted from <Dj|H a|Di> calculated when the
        !           excitation was chosen to give the commutator)

        use system, only: sys_t
        use excitations, only: excit_t, create_excited_det, get_excitation_level
        use ccmc_data, only: cluster_t
        use hamiltonian, only: get_hmatel

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f0(sys%basis%string_len)
        type(excit_t), intent(in) :: connection
        type(cluster_t), intent(in) :: cluster
        integer(i0), intent(in) :: cdet(sys%basis%string_len)
        integer(i0), intent(in) :: funlinked(sys%basis%string_len)
        real(p) :: hmatel

        integer(i0) :: deti(sys%basis%string_len), detj(sys%basis%string_len)
        integer(i0) :: temp(sys%basis%string_len)
        integer :: i, found, excitor_level, excitor_sign
        logical :: allowed
        real(p) :: population

        ! Find the determinant obtained by applying all of the cluster operators
        ! linked to H to D0
        population = 0.0_p ! The population doesn't matter as the commutator does not change the amplitude
        found = 0
        if (cluster%nexcitors > 1) then
            do i = 1, cluster%nexcitors
                if (any(cluster%excitors(i)%f /= funlinked)) then
                    ! Linked excitor, needed in cluster
                    found = found + 1
                    if (found == 1) then
                        deti = cluster%excitors(i)%f
                    else
                        call collapse_cluster(sys%basis, f0, cluster%excitors(i)%f, 1.0_p, &
                            deti, population, allowed)
                    end if
                end if
            end do
        else
            ! Only the unlinked excitor
            deti = f0
        end if

        ! Now we want to evaluate <D_i^a|H_i^a|D> ...
        call create_excited_det(sys%basis, deti, connection, detj)
        ! [todo] - general case call is slow.  Improvements: Slater--Condon procedure for
        ! [todo] - the relevant excitation level and system-specific procedures.
        hmatel = get_hmatel(sys, deti, detj)

        ! hmatel will be multiplied by cluster%amplitude and cluster%cluster_to_det_sign which
        ! potentially introduce unwanted sign changes, so we deal with them here
        hmatel = hmatel*cluster%cluster_to_det_sign
        ! Multiplying excitors can give a sign change, which is absorbed into cluster%amplitude
        population = 1.0_p
        temp = deti
        call collapse_cluster(sys%basis, f0, funlinked, 1.0_p, temp, population, allowed)
        hmatel = population*hmatel

        ! Possible sign changes from <D|a|D_0> ...
        if (cluster%nexcitors > 1) then
            excitor_level = get_excitation_level(f0, deti)
            call convert_excitor_to_determinant(deti, excitor_level, excitor_sign, f0)
            if (excitor_sign < 0) hmatel = -hmatel
        end if

        ! ... and <D_k|a_unlinked|D_i^a>
        call create_excited_det(sys%basis, cdet, connection, deti)
        excitor_level = get_excitation_level(deti, detj)
        call convert_excitor_to_determinant(deti, excitor_level, excitor_sign, detj)
        if (excitor_sign < 0) hmatel = -hmatel

    end function unlinked_commutator

    subroutine linked_spawner_ccmc(rng, sys, qmc_in, qs, spawn_cutoff, cluster, gen_excit_ptr, nspawn, &
                            connection, nspawnings_total, fexcit, ldet, rdet, left_cluster, right_cluster)

        ! When sampling e^-T H e^T, clusters need to be considered where two
        ! operators excite from/to the same orbital (one in the "left cluster"
        ! for e^-T and one in the "right cluster" for e^T). This makes the
        ! selection of an excitation more complicated as all possible
        ! partitionings of the cluster need to be accounted for: just because
        ! one partition has a zero contribution doesn't mean all will.

        ! See comments in spawner_ccmc for more details about spawning

        ! In:
        !    sys: system being studied.
        !    qmc_in: input options relating to QMC methods.
        !    qs: qmc_state_t object. ref and tau are used.
        !    spawn_cutoff: The size of the minimum spawning event allowed, in
        !        the encoded representation. Events smaller than this will be
        !        stochastically rounded up to this value or down to zero.
        !    cluster: information about the cluster which forms the excitor.  In
        !        particular, we use the amplitude, cluster_to_det_sign and pselect
        !        (i.e. n_sel.p_s.p_clust) attributes in addition to any used in
        !        the excitation generator.
        !    gen_excit_ptr: procedure pointer to excitation generators.
        !        gen_excit_ptr%full *must* be set to a procedure which generates
        !        a complete excitation.
        !    nspawnings_total: The total number of spawnings attemped by the current cluster
        !        in the current timestep.
        ! In/Out:
        !    rng: random number generator.
        !    ldet, rdet, left_cluster, right_cluster: used to store temporary information for
        !        selecting an excitor
        ! Out:
        !    nspawn: number of particles spawned, in the encoded representation.
        !        0 indicates the spawning attempt was unsuccessful.
        !    connection: excitation connection between the current excitor
        !        and the child excitor, on which progeny are spawned.
        !    fexcit: the bitstring of the determinant to spawn on to (as it is
        !        not necessarily easy to get from the connection)

        use ccmc_data, only: cluster_t
        use determinants, only: det_info_t
        use dSFMT_interface, only: dSFMT_t
        use excitations, only: excit_t, create_excited_det, get_excitation_level
        use proc_pointers, only: gen_excit_ptr_t, decoder_ptr
        use spawning, only: attempt_to_spawn
        use system, only: sys_t
        use const, only: depsilon
        use hamiltonian, only: get_hmatel
        use bit_utils, only: count_set_bits
        use qmc_data, only: qmc_in_t, qmc_state_t

        type(sys_t), intent(in) :: sys
        type(qmc_in_t), intent(in) :: qmc_in
        type(qmc_state_t), intent(in) :: qs
        integer(int_p), intent(in) :: spawn_cutoff
        type(cluster_t), intent(in) :: cluster
        type(dSFMT_t), intent(inout) :: rng
        integer, intent(in) :: nspawnings_total
        type(gen_excit_ptr_t), intent(in) :: gen_excit_ptr
        integer(int_p), intent(out) :: nspawn
        type(excit_t), intent(out) :: connection
        integer(i0), intent(out) :: fexcit(sys%basis%string_len)
        type(det_info_t), intent(inout) :: ldet, rdet
        type(cluster_t), intent(inout) :: left_cluster, right_cluster

        ! We incorporate the sign of the amplitude into the Hamiltonian matrix
        ! element, so we 'pretend' to attempt_to_spawn that all excips are
        ! actually spawned by positive excips.
        integer(int_p), parameter :: parent_sign = 1_int_p
        integer :: excitor_sign, excitor_level

        integer :: i, j, npartitions, orb, bit_pos, bit_element
        real(p) :: ppart, pgen, hmatel, pop, delta_h
        logical :: allowed, sign_change, linked, single_unlinked
        integer(i0) :: new_det(sys%basis%string_len)
        integer(i0) :: excitor(sys%basis%string_len)


        ! 1) Choose an order for the excitors
        ! The number of clusters with disallowed partitions is relatively small (20% in
        ! Ne pVDZ CCSDTQ) and should decrease with increasing number of electrons and/or
        ! basis functions. Rather than attempt to find an allowed partition (expensive)
        ! we select one at random  to excite from and only evaluate the full contribution
        ! (from all partitions) if the randomly-selected partition has a non-zero contribution,
        ! with an appropriately rescaled pgen.
        call partition_cluster(rng, sys, qs%ref%f0, cluster, left_cluster, right_cluster, ppart, ldet%f, &
                               rdet%f, allowed, sign_change)
        pop = 1

        ! 2) Choose excitation from right_cluster|D_0>
        if (allowed) then
            call decoder_ptr(sys, rdet%f, rdet)
            call gen_excit_ptr%full(rng, sys, qs%pattempt_single, rdet, pgen, connection, hmatel, allowed)
        end if

        if (allowed) then
            ! check that left_cluster can be applied to the resulting excitor to
            ! give a cluster to spawn on to
            call create_excited_det(sys%basis, rdet%f, connection, fexcit)
            call collapse_cluster(sys%basis, qs%ref%f0, ldet%f, 1.0_p, fexcit, pop, allowed)
        end if

        if (allowed) then
            ! If the excitation is not linked to the cluster, the matrix element
            ! is 0 and we can reject the spawning attempt.
            call linked_excitation(sys%basis, qs%ref%f0, connection, cluster, linked, single_unlinked, excitor)
        else
            linked = .false.
        end if

        if (linked) then
            ! 3) Evaluate commutator and pgen
            ! pgen and hmatel need recalculating to account for other permutations
            npartitions = nint(1.0/ppart)
            hmatel = 0.0_p
            pgen = 0.0_p
            do i = 1, npartitions
                ! Iterate over all allowed partitions and get contribution to
                ! hmatel and pgen
                call partition_cluster(rng, sys, qs%ref%f0, cluster, left_cluster, right_cluster, ppart, &
                                       ldet%f, rdet%f, allowed, sign_change, i)
                if (allowed) then
                    ! need to check that the excitation is valid!
                    do j = 1, connection%nexcit
                        ! i/j orbital should be occupied
                        orb = connection%from_orb(j)
                        bit_pos = sys%basis%bit_lookup(1,orb)
                        bit_element = sys%basis%bit_lookup(2,orb)
                        if (.not. btest(rdet%f(bit_element), bit_pos)) then
                            allowed = .false.
                            exit
                        end if
                        ! a/b orbital should be unoccupied
                        orb = connection%to_orb(j)
                        bit_pos = sys%basis%bit_lookup(1,orb)
                        bit_element = sys%basis%bit_lookup(2,orb)
                        if (btest(rdet%f(bit_element), bit_pos)) then
                            allowed = .false.
                            exit
                        end if
                    end do
                end if

                if (allowed) then

                    ! It's possible to get the same excitation from different partitionings
                    ! of the cluster so they all need to be accounted for in pgen
                    call create_excited_det(sys%basis, rdet%f, connection, new_det)
                    pgen = pgen + calc_pgen(sys, qmc_in%excit_gen, qs, rdet%f, connection, rdet)

                    ! Sign of the term in the commutator depends on the number of Ts in left_cluster
                    ! also need to account for possible sign change on going from excitor to determinant
                    ! for each of the determinants and when collapsing the clusters
                    ! Need <D_right|right_cluster|D0> and <D_spawn|left_cluster|D>
                    delta_h = get_hmatel(sys, rdet%f, new_det)

                    if (mod(left_cluster%nexcitors,2) /= 0) delta_h = -delta_h

                    excitor_level = get_excitation_level(qs%ref%f0, rdet%f)
                    call convert_excitor_to_determinant(rdet%f, excitor_level, excitor_sign, qs%ref%f0)
                    if (excitor_sign < 0) delta_h = -delta_h

                    excitor_level = get_excitation_level(fexcit, new_det)
                    call convert_excitor_to_determinant(fexcit, excitor_level, excitor_sign, new_det)
                    if (excitor_sign < 0) delta_h = -delta_h

                    if (sign_change) delta_h = -delta_h

                    hmatel = hmatel + delta_h
                end if
            end do

            ! apply additional factors to pgen
            pgen = pgen*cluster%pselect*nspawnings_total/npartitions

            ! correct hmatel for cluster amplitude
            hmatel = hmatel*cluster%amplitude
            excitor_level = get_excitation_level(fexcit, qs%ref%f0)
            call convert_excitor_to_determinant(fexcit, excitor_level, excitor_sign, qs%ref%f0)
            if (excitor_sign < 0) hmatel = -hmatel

            ! 4) Attempt to spawn
            nspawn = attempt_to_spawn(rng, qs%tau, spawn_cutoff, qs%psip_list%pop_real_factor, hmatel, pgen, parent_sign)

        else
            nspawn = 0
        end if

    end subroutine linked_spawner_ccmc

    subroutine partition_cluster(rng, sys, f0, cluster, left_cluster, right_cluster, ppart, &
                                 ldet, rdet, allowed, sign_change, part_number)

        ! Divides a cluster into two halves such that any excitors that share an
        ! orbital are in different halves

        ! In:
        !    cluster: the cluster of excitors
        !    sys: the system being studied
        !    f0: bit string of the reference
        !    part_number: (optional) use to enumerate partitions rather than
        !            choosing a random one
        ! In/Out:
        !    rng: random number generator
        ! Out:
        !    left_cluster: the excitors to be applied after the Hamiltonian
        !    right_cluster: the excitors to be applied before the Hamiltonian
        !    ppart: the probability of choosing this partition
        !    ldet: the determinant formed by applying left_cluster to the reference
        !    rdet: the determinant formed by applying right_cluster to the reference
        !    allowed: are left_cluster and right_cluster both valid clusters
        !    sign_change: is there a net sign change on collapsing the two clusters

        use ccmc_data, only: cluster_t
        use system, only: sys_t
        use dSFMT_interface, only: dSFMT_t, get_rand_close_open

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f0(sys%basis%string_len)
        type(cluster_t), intent(in) :: cluster
        type(cluster_t), intent(inout) :: left_cluster, right_cluster
        type(dSFMT_t), intent(inout) :: rng
        real(p), intent(out) :: ppart
        integer(i0), intent(inout) :: ldet(sys%basis%string_len), rdet(sys%basis%string_len)
        logical, intent(out) :: allowed, sign_change
        integer, intent(in), optional :: part_number

        integer :: i, in_left, in_right, side
        real(p) :: rand, population

        in_left = 0
        in_right = 0
        ppart = 1.0_p
        allowed = .true.
        population = 1.0_p

        do i = 1, cluster%nexcitors
            if (present(part_number)) then
                ! Use integers 1 to 2^n as bitstrings to consider all
                ! partititions in turn
                if (btest(part_number, i-1)) then
                    side = 1
                else
                    side = -1
                end if
            else
                ! for simplicity of generation probabilities just assign everything
                ! randomly to left or right
                ppart = ppart * 0.5_p
                rand = get_rand_close_open(rng)
                if (rand < 0.5_p) then
                    side = -1
                else
                    side = 1
                end if
            end if
            if (side == -1) then
                ! add to left_cluster
                in_left = in_left + 1
                left_cluster%excitors(in_left)%f => cluster%excitors(i)%f
                if (in_left == 1) then
                    ldet = cluster%excitors(i)%f
                else
                    call collapse_cluster(sys%basis, f0, cluster%excitors(i)%f, 1.0_p, ldet, population, allowed)
                    if (.not.allowed) exit
                end if
            else
                ! add to right
                in_right = in_right + 1
                right_cluster%excitors(in_right)%f => cluster%excitors(i)%f
                if (in_right == 1) then
                    rdet = cluster%excitors(i)%f
                else
                    call collapse_cluster(sys%basis, f0, cluster%excitors(i)%f, 1.0_p, rdet, population, allowed)
                    if (.not.allowed) exit
                end if
            end if
        end do

        left_cluster%nexcitors = in_left
        right_cluster%nexcitors = in_right

        sign_change = (population < 0)

    end subroutine partition_cluster

    function calc_pgen(sys, excit_gen, qmc_state, f, connection, parent_det) result(pgen)

        ! Calculate the probability of an excitation being selected.
        ! Wrapper round system specific functions.

        ! In:
        !    sys: the system being studied
        !    excit_gen: which excitation generator is being used.  Should correspond to a value in the excit_gen_* enum in qmc_data.
        !    qmc_state: input options relating to QMC methods.
        !    f: bit string representation of parent excitor
        !    connection: excitation connection between the current excitor
        !        and the child excitor, on which progeny are spawned.
        !    parent_det: information on the parent excitor
        ! Returns:
        !    The probability of generating the excitation specified by
        !    connection, assuming the excitation is valid

        use bit_utils, only: count_set_bits
        use errors, only: stop_all
        use system
        use excitations, only: excit_t
        use excit_gen_mol, only: calc_pgen_single_mol_no_renorm, calc_pgen_double_mol_no_renorm, &
                                 calc_pgen_single_mol, calc_pgen_double_mol
        use excit_gen_ueg, only: calc_pgen_ueg_no_renorm
        use excit_gen_ringium, only: calc_pgen_ringium
        use point_group_symmetry, only: cross_product_pg_basis, pg_sym_conj
        use determinants, only: det_info_t
        use qmc_data, only: qmc_state_t, excit_gen_no_renorm

        type(sys_t), intent(in) :: sys
        integer, intent(in) :: excit_gen
        type(qmc_state_t), intent(in) :: qmc_state
        integer(i0), intent(in) :: f(sys%basis%string_len)
        type(excit_t), intent(in) :: connection
        type(det_info_t), intent(in) :: parent_det
        real(p) :: pgen

        integer :: spin, ij_sym, max_na, ij_k(3)
        integer(i0) :: poss_a(sys%basis%string_len)

        associate(a=>connection%to_orb(1), b=>connection%to_orb(2), i=>connection%from_orb(1), j=>connection%from_orb(2))
            select case(sys%system)
            case(read_in)
                if (excit_gen == excit_gen_no_renorm) then
                    if (connection%nexcit == 1) then
                        pgen = qmc_state%pattempt_single * calc_pgen_single_mol_no_renorm(sys, a)
                    else
                        spin = sys%basis%basis_fns(a)%ms + sys%basis%basis_fns(b)%ms
                        pgen = qmc_state%pattempt_double * calc_pgen_double_mol_no_renorm(sys, a, b, spin)
                    end if
                else
                    if (connection%nexcit == 1) then
                        pgen = qmc_state%pattempt_single * calc_pgen_single_mol(sys, sys%read_in%pg_sym%gamma_sym, &
                                                                                parent_det%occ_list, parent_det%symunocc, a)
                    else
                        spin = sys%basis%basis_fns(a)%ms + sys%basis%basis_fns(b)%ms
                        ij_sym = pg_sym_conj(sys%read_in%pg_sym, &
                                             cross_product_pg_basis(sys%read_in%pg_sym, a, b, sys%basis%basis_fns))
                        pgen = qmc_state%pattempt_double * calc_pgen_double_mol(sys, ij_sym, a, b, spin, parent_det%symunocc)
                    end if
                end if
            case(ueg)
                spin = sys%basis%basis_fns(i)%ms + sys%basis%basis_fns(j)%ms
                ij_k = 0
                ij_k(1:sys%lattice%ndim) = sys%basis%basis_fns(i)%l + sys%basis%basis_fns(j)%l
                if (spin == -2) then
                    poss_a = iand(not(f), ishft(sys%ueg%ternary_conserve(1:,ij_k(1),ij_k(2),ij_k(3)),1))
                else
                    poss_a = iand(not(f), sys%ueg%ternary_conserve(1:,ij_k(1),ij_k(2),ij_k(3)))
                end if
                max_na = sum(count_set_bits(poss_a))
                pgen = calc_pgen_ueg_no_renorm(sys, max_na, spin)
            case(ringium)
                pgen = calc_pgen_ringium(sys)
            case default
                call stop_all('calc_pgen', 'Linked CCMC is not implemented for this system.')
            end select
        end associate
    end function calc_pgen

end module ccmc

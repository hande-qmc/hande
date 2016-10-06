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

use const, only: i0, int_p, int_64, p, dp, debug

implicit none

contains

    subroutine do_ccmc(sys, qmc_in, ccmc_in, semi_stoch_in, restart_in, load_bal_in, reference_in, &
                        logging_in, blocking_in, io_unit, qs, qmc_state_restart)

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
        !    io_unit: input unit to write all output to.
        ! In/Out:
        !    qmc_state_restart (optional): if present, restart from a previous fciqmc calculation.
        !       Deallocated on exit.
        ! Out:
        !    qs: qmc_state for use if restarting the calculation

        use checking, only: check_allocate, check_deallocate
        use dSFMT_interface, only: dSFMT_t, dSFMT_init, dSFMT_end
        use errors, only: stop_all
        use parallel
        use restart_hdf5, only: dump_restart_hdf5, restart_info_t, init_restart_info_t, dump_restart_file_wrapper

        use annihilation, only: direct_annihilation
        use bloom_handler, only: init_bloom_stats_t, bloom_stats_t, bloom_mode_fractionn, bloom_mode_fixedn, &
                                 write_bloom_report, bloom_stats_warning, update_bloom_threshold_prop
        use ccmc_data
        use ccmc_selection, only: select_cluster, create_null_cluster, select_nc_cluster, select_cluster_truncated
        use ccmc_selection, only: init_selection_data, update_selection_probabilities, set_cluster_selections, &
                                  init_amp_psel_accumulation
        use ccmc_death_spawning, only: stochastic_ccmc_death_nc
        use ccmc_utils, only: init_contrib, dealloc_contrib, find_D0, cumulative_population, init_ex_lvl_dist_t, &
                              update_ex_lvl_dist, regenerate_ex_levels_psip_list
        use determinants, only: det_info_t, alloc_det_info_t, dealloc_det_info_t, sum_sp_eigenvalues_occ_list, &
                                sum_sp_eigenvalues_bit_string, decode_det
        use excitations, only: excit_t, get_excitation_level, get_excitation
        use qmc_io, only: write_qmc_report, write_qmc_report_header
        use qmc, only: init_qmc
        use qmc_common, only: initial_fciqmc_status, load_balancing_report, &
                              init_report_loop, init_mc_cycle, end_report_loop, end_mc_cycle,      &
                              redistribute_particles, rescale_tau
        use proc_pointers
        use spawning, only: assign_particle_processor
        use system, only: sys_t, sys_t_json
        use spawn_data, only: calc_events_spawn_t, write_memcheck_report
        use replica_rdm, only: update_rdm, calc_rdm_energy, write_final_rdm

        use qmc_data, only: qmc_in_t, ccmc_in_t, semi_stoch_in_t, restart_in_t
        use qmc_data, only: blocking_in_t
        use qmc_data, only: load_bal_in_t, qmc_state_t, annihilation_flags_t, estimators_t, blocking_t
        use qmc_data, only: qmc_in_t_json, ccmc_in_t_json, semi_stoch_in_t_json, restart_in_t_json
        use qmc_data, only: blocking_in_t_json
        use reference_determinant, only: reference_t, reference_t_json
        use check_input, only: check_qmc_opts, check_ccmc_opts
        use json_out, only: json_out_t, json_object_init, json_object_end
        use hamiltonian_data
        use energy_evaluation, only: get_sanitized_projected_energy, get_sanitized_projected_energy_cmplx
        use blocking
        use utils, only: get_free_unit

        use logging, only: init_logging, end_logging, prep_logging_mc_cycle, write_logging_calc_ccmc
        use logging, only: logging_in_t, logging_t, logging_in_t_json, logging_t_json, write_logging_select_ccmc

        type(sys_t), intent(in) :: sys
        type(qmc_in_t), intent(in) :: qmc_in
        type(ccmc_in_t), intent(in) :: ccmc_in
        type(blocking_in_t), intent(in) :: blocking_in
        type(semi_stoch_in_t), intent(in) :: semi_stoch_in
        type(restart_in_t), intent(in) :: restart_in
        type(load_bal_in_t), intent(in) :: load_bal_in
        type(reference_t), intent(in) :: reference_in
        type(logging_in_t), intent(in) :: logging_in
        type(qmc_state_t), target, intent(out) :: qs
        type(qmc_state_t), intent(inout), optional :: qmc_state_restart
        integer, intent(in) :: io_unit

        integer :: i, ireport, icycle, iter, semi_stoch_iter, it
        integer(int_64) :: iattempt
        integer(int_64) :: nattempts_spawn
        real(dp), allocatable :: nparticles_old(:), nparticles_change(:)
        type(det_info_t) :: ref_det

        integer(int_p) :: ndeath, ndeath_nc
        integer(int_p) :: ndeath_im, ndeath_nc_im
        integer :: nspawn_events, ierr
        type(wfn_contrib_t), allocatable :: contrib(:)
        type(multispawn_stats_t), allocatable :: ms_stats(:)
        type(dSFMT_t), allocatable :: rng(:)
        type(json_out_t) :: js
        type(qmc_in_t) :: qmc_in_loc
        type(logging_t) :: logging_info
        type(selection_data_t) :: selection_data

        logical :: soft_exit, dump_restart_shift, restarting

        real(p), allocatable :: cumulative_abs_real_pops(:)
        integer :: D0_proc, D0_pos, nD0_proc, min_cluster_size, max_cluster_size, iexcip_pos, slot
        real(p) :: tot_abs_real_pop
        complex(p) :: D0_normalisation
        type(bloom_stats_t) :: bloom_stats
        type(annihilation_flags_t) :: annihilation_flags
        type(restart_info_t) :: ri, ri_shift
        character(36) :: uuid_restart
        type(ex_lvl_dist_t) :: ex_lvl_dist

        real :: t1, t2

        logical :: update_tau, error

        logical :: seen_D0, regenerate_info
        real(p) :: dfock
        complex(p) :: D0_population_cycle, proj_energy_cycle

        real(p), allocatable :: rdm(:,:)

        type(blocking_t) :: bl
        integer :: iunit

        if (parent) then
            write (io_unit,'(1X,"CCMC")')
            write (io_unit,'(1X,"----",/)')
        end if

        ! Check input options.
        if (parent) then
            restarting = present(qmc_state_restart) .or. restart_in%read_restart
            call check_qmc_opts(qmc_in, sys, .not.present(qmc_state_restart), restarting)
            call check_ccmc_opts(sys, ccmc_in)
        end if

        ! Initialise data.
        call init_qmc(sys, qmc_in, restart_in, load_bal_in, reference_in, io_unit, annihilation_flags, qs, &
                      uuid_restart, qmc_state_restart=qmc_state_restart, regenerate_info=regenerate_info)

        if (ccmc_in%even_selection .and. regenerate_info) then
            call regenerate_ex_levels_psip_list(sys%basis, qs)
        else if (regenerate_info) then
            call stop_all('do_ccmc', &
            'Asked to regenerate extra information after restart but no extra information expected.')
        end if

        if (debug) call init_logging(logging_in, logging_info, qs%ref%ex_level)

        if (parent) then
            call json_object_init(js, tag=.true., io=io_unit)
            call sys_t_json(js, sys)
            ! The default values of pattempt_* are not in qmc_in
            qmc_in_loc = qmc_in
            qmc_in_loc%pattempt_single = qs%excit_gen_data%pattempt_single
            qmc_in_loc%pattempt_double = qs%excit_gen_data%pattempt_double
            call qmc_in_t_json(js, qmc_in_loc)
            call ccmc_in_t_json(js, ccmc_in)
            call semi_stoch_in_t_json(js, semi_stoch_in)
            call restart_in_t_json(js, restart_in, uuid_restart)
            call reference_t_json(js, qs%ref, sys)
            call blocking_in_t_json(js, blocking_in)
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
        if (ccmc_in%even_selection) then
            call init_bloom_stats_t(bloom_stats, mode=bloom_mode_fixedn, encoding_factor=qs%psip_list%pop_real_factor)
        else
            call init_bloom_stats_t(bloom_stats, mode=bloom_mode_fractionn, encoding_factor=qs%psip_list%pop_real_factor)
        end if

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

        call init_contrib(sys, qs%ref%ex_level+2, ccmc_in%linked, contrib)

        do i = 0, nthreads-1
            ! Initialise and allocate RNG store.
            call dSFMT_init(qmc_in%seed+iproc+i*nprocs, 50000, rng(i))
        end do

        ! ...and scratch space for calculative cumulative probabilities.
        allocate(cumulative_abs_real_pops(size(qs%psip_list%states,dim=2)), stat=ierr)
        call check_allocate('cumulative_abs_real_pops', size(qs%psip_list%states, dim=2), ierr)

        if (ccmc_in%even_selection) then
            if (ccmc_in%linked) then
                call init_selection_data(qs%ref%ex_level, 4, selection_data)
            else
                call init_selection_data(qs%ref%ex_level, qs%ref%ex_level+2, selection_data)
            end if
            call init_ex_lvl_dist_t(qs%ref%ex_level, ex_lvl_dist)
        end if
        if (debug) call init_amp_psel_accumulation(qs%ref%ex_level+2, logging_info, ccmc_in%linked, selection_data)

        nparticles_old = qs%psip_list%tot_nparticles

        ! Initialise D0_pos to be somewhere (anywhere) in the list.
        D0_pos = 1

        ! Main fciqmc loop.
        if (parent) call write_qmc_report_header(qs%psip_list%nspaces, cmplx_est=sys%read_in%comp, &
                                            rdm_energy=ccmc_in%density_matrices, nattempts=.true., io_unit=io_unit)
        call initial_fciqmc_status(sys, qmc_in, qs, doing_ccmc=.true., io_unit=io_unit)
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

        if (ccmc_in%density_matrices) then
            associate(nbasis=>sys%basis%nbasis)
                allocate(rdm(nbasis*(nbasis-1)/2,nbasis*(nbasis-1)/2), stat=ierr)
                call check_allocate('rdm', nbasis**2*(nbasis-1)**2/4, ierr)
                rdm = 0.0_p
            end associate
            call alloc_det_info_t(sys, ref_det)
            ref_det%f = qs%ref%f0
            call decode_det(sys%basis, ref_det%f, ref_det%occ_list)
        end if

        if (parent .and. blocking_in%blocking_on_the_fly) then
            iunit = get_free_unit()
            open(iunit, file=blocking_in%filename, status='unknown')
            call write_blocking_report_header(iunit)
        end if

        if (blocking_in%blocking_on_the_fly) call allocate_blocking(qmc_in, blocking_in, bl)
        do ireport = 1, qmc_in%nreport

            ! Projected energy from last report loop to correct death
            if (sys%read_in%comp) then
                qs%estimators%proj_energy_old = get_sanitized_projected_energy_cmplx(qs)
            else
                qs%estimators%proj_energy_old = get_sanitized_projected_energy(qs)
            end if
            call init_report_loop(qs, bloom_stats)

            do icycle = 1, qmc_in%ncycles

                iter = qs%mc_cycles_done + (ireport-1)*qmc_in%ncycles + icycle

                if (debug) call prep_logging_mc_cycle(iter, logging_in, logging_info, sys%read_in%comp, &
                                                        min(sys%nel, qs%ref%ex_level+2))

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
                    if (sys%read_in%comp) then
                        D0_normalisation = cmplx(qs%psip_list%pops(1,D0_pos), qs%psip_list%pops(2,D0_pos), p)&
                                                    /qs%psip_list%pop_real_factor
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
                call mpi_bcast(D0_normalisation, 1, mpi_pcomplex, D0_proc, MPI_COMM_WORLD, ierr)
#endif

                ! Note that 'death' in CCMC creates particles in the spawned
                ! list, so the number of deaths not in the spawned list is
                ! always 0.
                call init_mc_cycle(qs%psip_list, qs%spawn_store%spawn, qs%estimators(1)%nattempts, ndeath, &
                                   min_attempts=nint(abs(D0_normalisation), kind=int_64), &
                                   ndeath_im = ndeath_im, &
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
                                           tot_abs_real_pop, ex_lvl_dist)

                if (.not.ccmc_in%even_selection) call update_bloom_threshold_prop(bloom_stats, nparticles_old(1))

                if (ccmc_in%even_selection) then
                    call update_ex_lvl_dist(ex_lvl_dist)
                    call update_selection_probabilities(cumulative_abs_real_pops, ex_lvl_dist, &
                                                    abs(D0_normalisation), tot_abs_real_pop, selection_data)
                end if


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
                call set_cluster_selections(selection_data, qs%estimators(1)%nattempts, min_cluster_size, max_cluster_size, &
                                            D0_normalisation, tot_abs_real_pop, qs%psip_list%nstates, ccmc_in%full_nc, &
                                            ccmc_in%even_selection)
                ! OpenMP chunk size determined completely empirically from a single
                ! test.  Please feel free to improve...
                ! NOTE: we can't refer to procedure pointers in shared blocks so
                ! can't use default(none).  I *strongly* recommend turning
                ! default(none) on when making changes and ensure that the only
                ! errors relate to the procedure pointers...
                !$omp parallel default(none) &
                !$omp private(it, iexcip_pos, i, seen_D0) &
                !$omp shared(rng, cumulative_abs_real_pops, tot_abs_real_pop,  &
                !$omp        max_cluster_size, contrib, D0_normalisation, D0_pos, rdm,    &
                !$omp        qs, sys, bloom_stats, min_cluster_size, ref_det,             &
                !$omp        proj_energy_cycle, D0_population_cycle, selection_data,      &
                !$omp        nattempts_spawn, ex_lvl_dist, &
                !$omp        ccmc_in, nprocs, ms_stats, qmc_in, load_bal_in, ndeath_nc,   &
                !$omp        nparticles_change, ndeath, ndeath_nc_im, logging_info)
                it = get_thread_id()
                iexcip_pos = 0
                seen_D0 = .false.
                proj_energy_cycle = cmplx(0.0, 0.0, p)
                D0_population_cycle = cmplx(0.0, 0.0, p)
                !$omp do schedule(dynamic,200) reduction(+:D0_population_cycle,proj_energy_cycle,nattempts_spawn,ndeath)
                do iattempt = 1, selection_data%nclusters
                    if (iattempt <= selection_data%nsingle_excitors) then
                        ! As noncomposite clusters can't be above truncation level or linked-only all can accumulate +
                        ! propagate. Only need to check not selecting the reference as we treat it separately.
                        if (iattempt /= D0_pos) then
                            ! Deterministically select each excip as a non-composite cluster.
                            call select_nc_cluster(sys, qs%psip_list, qs%ref%f0, &
                                        iattempt, qmc_in%initiator_pop, ccmc_in%even_selection, &
                                        contrib(it)%cdet, contrib(it)%cluster)

                            if (qs%propagator%quasi_newton) contrib(it)%cdet%fock_sum = &
                                            sum_sp_eigenvalues_occ_list(sys, contrib(it)%cdet%occ_list) - qs%ref%fock_sum

                            call do_ccmc_accumulation(sys, qs, contrib(it)%cdet, contrib(it)%cluster, logging_info, &
                                                    D0_population_cycle, proj_energy_cycle, ccmc_in, ref_det, rdm, selection_data)
                            call do_nc_ccmc_propagation(rng(it), sys, qs, ccmc_in, logging_info, bloom_stats, &
                                                                contrib(it), nattempts_spawn)
                        end if

                    ! For OpenMP scalability, have this test inside a single loop rather
                    ! than attempt to parallelise over three separate loops.
                    else if (iattempt <= selection_data%nsingle_excitors + selection_data%nstochastic_clusters) then
                        if (ccmc_in%even_selection) then
                            call select_cluster_truncated(rng(it), sys, qs%psip_list, qs%ref%f0, &
                                                        ccmc_in%linked, selection_data%nstochastic_clusters, D0_normalisation, &
                                                        qmc_in%initiator_pop, selection_data, cumulative_abs_real_pops, &
                                                        qs%ref%ex_level, min_cluster_size, max_cluster_size, &
                                                        ex_lvl_dist, contrib(it)%cluster, contrib(it)%cdet)

                        else
                            call select_cluster(rng(it), sys, qs%psip_list, qs%ref%f0, qs%ref%ex_level, ccmc_in%linked, &
                                            selection_data%nstochastic_clusters, D0_normalisation, qmc_in%initiator_pop, D0_pos, &
                                            cumulative_abs_real_pops, tot_abs_real_pop, min_cluster_size, max_cluster_size, &
                                            logging_info, contrib(it)%cdet, contrib(it)%cluster)
                        end if

                        if (contrib(it)%cluster%excitation_level <= qs%ref%ex_level+2 .or. &
                                (ccmc_in%linked .and. contrib(it)%cluster%excitation_level == huge(0))) then
                            ! cluster%excitation_level == huge(0) indicates a cluster
                            ! where two excitors share an elementary operator

                            if (qs%propagator%quasi_newton) contrib(it)%cdet%fock_sum = &
                                            sum_sp_eigenvalues_occ_list(sys, contrib(it)%cdet%occ_list) - qs%ref%fock_sum

                            call do_ccmc_accumulation(sys, qs, contrib(it)%cdet, contrib(it)%cluster, logging_info, &
                                                    D0_population_cycle, proj_energy_cycle, ccmc_in, ref_det, rdm, selection_data)
                            call do_stochastic_ccmc_propagation(rng(it), sys, qs, &
                                                                ccmc_in, logging_info, ms_stats(it), bloom_stats, &
                                                                contrib(it), nattempts_spawn, ndeath)
                        end if
                    else
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
                                                     qmc_in%initiator_pop, contrib(it)%cdet, contrib(it)%cluster)
                        end if

                        if (qs%propagator%quasi_newton) contrib(it)%cdet%fock_sum = &
                                        sum_sp_eigenvalues_occ_list(sys, contrib(it)%cdet%occ_list) - qs%ref%fock_sum

                        call do_ccmc_accumulation(sys, qs, contrib(it)%cdet, contrib(it)%cluster, logging_info, &
                                                D0_population_cycle, proj_energy_cycle, ccmc_in, ref_det, rdm, selection_data)
                        nattempts_spawn = nattempts_spawn + 1
                        call perform_ccmc_spawning_attempt(rng(it), sys, qs, ccmc_in, logging_info, bloom_stats, contrib(it), 1)
                    end if
                end do
                !$omp end do

                ndeath_nc = 0
                ndeath_nc_im = 0
                if (ccmc_in%full_nc .and. qs%psip_list%nstates > 0) then
                    ! Do death exactly and directly for non-composite clusters
                    !$omp do schedule(dynamic,200) private(dfock) reduction(+:ndeath_nc,nparticles_change,ndeath_nc_im)
                    do iattempt = 1, qs%psip_list%nstates
                        ! Note we use the (encoded) population directly in stochastic_ccmc_death_nc
                        ! (unlike the stochastic_ccmc_death) to avoid unnecessary decoding/encoding
                        ! steps (cf comments in stochastic_death for FCIQMC).
                        if (qs%propagator%quasi_newton) then
                            dfock = sum_sp_eigenvalues_bit_string(sys, qs%psip_list%states(:,iattempt)) - qs%ref%fock_sum
                        end if
                        call stochastic_ccmc_death_nc(rng(it), ccmc_in%linked, sys, qs, iattempt==D0_pos, dfock, &
                                          qs%psip_list%dat(1,iattempt), qs%estimators(1)%proj_energy_old, &
                                          qs%psip_list%pops(1, iattempt), nparticles_change(1), ndeath_nc, &
                                          logging_info)
                        if (sys%read_in%comp) then
                            call stochastic_ccmc_death_nc(rng(it), ccmc_in%linked, sys, qs, iattempt==D0_pos, dfock, &
                                              qs%psip_list%dat(1,iattempt), qs%estimators(2)%proj_energy_old, &
                                              qs%psip_list%pops(2, iattempt), nparticles_change(2), ndeath_nc, &
                                              logging_info)
                        end if
                    end do
                    !$omp end do
                end if
                !$omp end parallel

                if (ccmc_in%density_matrices .and. qs%vary_shift(1) .and. parent .and. .not. sys%read_in%comp) then
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

                    call direct_annihilation(sys, rng(0), qs%ref, annihilation_flags, pl, spawn)
                end associate
                if (debug) call write_logging_calc_ccmc(logging_info, iter, nspawn_events, ndeath + ndeath_nc, &
                                                        selection_data%nD0_select, &
                                                        selection_data%nclusters, selection_data%nstochastic_clusters, &
                                                        selection_data%nsingle_excitors)
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

            call end_report_loop(qmc_in, iter, update_tau, qs, nparticles_old, nspawn_events, &
                                 semi_stoch_in%shift_iter, semi_stoch_iter, soft_exit, &
                                 load_bal_in, bloom_stats=bloom_stats, comp=sys%read_in%comp, &
                                 error=error, vary_shift_reference=ccmc_in%vary_shift_reference)
            if (error) exit

            call cpu_time(t2)
            if (parent) then
                if (bloom_stats%nblooms_curr > 0) call bloom_stats_warning(bloom_stats, io_unit=io_unit)
                call write_qmc_report(qmc_in, qs, ireport, nparticles_old, t2-t1, .false., .false., &
                                        io_unit=io_unit, cmplx_est=sys%read_in%comp, rdm_energy=ccmc_in%density_matrices, &
                                        nattempts=.true.)
                if (blocking_in%blocking_on_the_fly) call do_blocking(bl, qs, qmc_in, ireport, iter, iunit, blocking_in)
            end if

            ! Update the time for the start of the next iteration.
            t1 = t2

            call dump_restart_file_wrapper(qs, dump_restart_shift, restart_in%write_freq, nparticles_old, ireport, &
                                           qmc_in%ncycles, sys%basis%nbasis, ri, ri_shift, .false., sys%basis%info_string_len)

            qs%psip_list%tot_nparticles = nparticles_old

            if (soft_exit) exit

            if (update_tau) call rescale_tau(qs%tau)

        end do

        if (blocking_in%blocking_on_the_fly) call deallocate_blocking(bl)
        if (blocking_in%blocking_on_the_fly) close(iunit, status='keep')

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
            call dump_restart_hdf5(ri, qs, qs%mc_cycles_done, nparticles_old, sys%basis%nbasis, .false., sys%basis%info_string_len)
            if (parent) write (io_unit,'()')
        end if

        if (debug) call end_logging(logging_info)
        if (debug .or. ccmc_in%even_selection) call end_selection_data(selection_data)

        if (ccmc_in%density_matrices) then
            call write_final_rdm(rdm, sys%nel, sys%basis%nbasis, ccmc_in%density_matrix_file, io_unit)
            call calc_rdm_energy(sys, qs%ref, rdm, qs%estimators(1)%rdm_energy, qs%estimators(1)%rdm_trace)
            if (parent) write (io_unit,'(1x,"# Final energy from RDM",2x,es17.10)') qs%estimators%rdm_energy/qs%estimators%rdm_trace
            deallocate(rdm, stat=ierr)
            call check_deallocate('rdm',ierr)
            call dealloc_det_info_t(ref_det)
        end if

        call dealloc_contrib(contrib, ccmc_in%linked)
        do i = 0, nthreads-1
            call dSFMT_end(rng(i))
        end do

    end subroutine do_ccmc

    subroutine do_ccmc_accumulation(sys, qs, cdet, cluster, logging_info, D0_population_cycle, proj_energy_cycle, &
                                    ccmc_in, ref_det, rdm, selection_data)



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
        use qmc_data, only: qmc_state_t, ccmc_in_t, estimators_t
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
        type(ccmc_in_t), intent(in) :: ccmc_in
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
            call update_proj_energy_ptr(sys, qs%ref%f0, qs%trial%wfn_dat, cdet, &
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

        if (ccmc_in%density_matrices .and. cluster%excitation_level <= 2 .and. qs%vary_shift(1) &
            .and. cluster%excitation_level /= 0 .and. .not. sys%read_in%comp) then
            ! Add contribution to density matrix
            ! d_pqrs = <HF|a_p^+a_q^+a_sa_r|CC>
            !$omp critical
            call update_rdm(sys, cdet, ref_det, &
                            real(cluster%amplitude, p)*cluster%cluster_to_det_sign, &
                            1.0_p, cluster%pselect, rdm)
            !$omp end critical
        end if

    end subroutine do_ccmc_accumulation

    subroutine do_stochastic_ccmc_propagation(rng, sys, qs, &
                                            ccmc_in, logging_info, ms_stats, bloom_stats, &
                                            contrib, nattempts_spawn_tot, ndeath)

        ! Perform stochastic propogation of a cluster in an appropriate manner
        ! for the given inputs. For stochastically selected clusters this
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


        use dSFMT_interface, only: dSFMT_t
        use system, only: sys_t
        use qmc_data, only: qmc_state_t, ccmc_in_t
        use ccmc_data, only: multispawn_stats_t, ms_stats_update, wfn_contrib_t

        use ccmc_death_spawning, only: spawner_ccmc, linked_spawner_ccmc, stochastic_ccmc_death
        use ccmc_death_spawning, only: stochastic_ccmc_death_nc, spawner_complex_ccmc
        use bloom_handler, only: bloom_stats_t, accumulate_bloom_stats
        use logging, only: logging_t

        type(sys_t), intent(in) :: sys
        type(dSFMT_T), intent(inout) :: rng
        type(qmc_state_t), intent(inout) :: qs
        type(ccmc_in_t), intent(in) :: ccmc_in
        type(wfn_contrib_t), intent(inout) :: contrib
        type(bloom_stats_t), intent(inout) :: bloom_stats
        type(logging_t), intent(in) :: logging_info

        integer(int_64), intent(inout) :: nattempts_spawn_tot
        integer(int_p), intent(inout) :: ndeath
        type(multispawn_stats_t), intent(inout) :: ms_stats

        integer :: i, nspawnings_cluster

        ! Spawning
        ! This has the potential to create blooms, so we allow for multiple
        ! spawning events per cluster.
        ! The number of spawning events is decided by the value
        ! of cluster%amplitude/cluster%pselect.  If this is
        ! greater than cluster_multispawn_threshold, then nspawnings is
        ! increased to the ratio of these.
        nspawnings_cluster=max(1,ceiling((abs(contrib%cluster%amplitude)&
                                     /contrib%cluster%pselect)/ ccmc_in%cluster_multispawn_threshold))

        call ms_stats_update(nspawnings_cluster, ms_stats)
        nattempts_spawn_tot = nattempts_spawn_tot + nspawnings_cluster

        do i = 1, nspawnings_cluster
            call perform_ccmc_spawning_attempt(rng, sys, qs, ccmc_in, logging_info, bloom_stats, contrib, nspawnings_cluster)
        end do

        ! Does the cluster collapsed onto D0 produce
        ! a determinant is in the truncation space?  If so, also
        ! need to attempt a death/cloning step.
        ! optimisation: call only once per iteration for clusters of size 0 or 1 for ccmc_in%full_nc.
        if (contrib%cluster%excitation_level <= qs%ref%ex_level) then
            ! Clusters above size 2 can't die in linked ccmc.
            if ((.not. ccmc_in%linked) .or. contrib%cluster%nexcitors <= 2) then
                ! Do death for non-composite clusters directly and in a separate loop
                if (contrib%cluster%nexcitors >= 2 .or. .not. ccmc_in%full_nc) then
                    call stochastic_ccmc_death(rng, qs%spawn_store%spawn, ccmc_in%linked, ccmc_in%even_selection, sys, &
                                               qs, contrib%cdet, contrib%cluster, logging_info, ndeath)
                end if
            end if
        end if

    end subroutine do_stochastic_ccmc_propagation

    subroutine do_nc_ccmc_propagation(rng, sys, qs, ccmc_in, logging_info, bloom_stats, &
                                            contrib, nattempts_spawn_tot)

        ! Perform stochastic propogation of a cluster selected deterministically
        ! in full non-composite selection. This performs all spawning attempts
        ! for a given noncomposite cluster, adding any created particles to the
        ! spawned list, while death is performed in-place later.

        ! This is in many ways similar to the approach used in fciqmc except in
        ! the case of solid-state calculations. In this case, noncomposite CCMC
        ! spawns from multiple chunks of magnitude one and phase equal to the
        ! excip population. FCIQMC instead spawns from purely real or imaginary
        ! parents.
        ! The CCMC approach is consistent with that used throughout complex CCMC,
        ! and requires fewer attempts per iteration to sample the same population,
        ! though the efficacy of the two approaches has not yet been investigated.

        ! In:
        !   sys: information on system under consideration.
        !   ccmc_in: options relating to ccmc passed in to calculation.
        !   logging_info: logging_t object with info about current logging
        !        when debug is true.
        ! In/Out:
        !   rng: random number generator.
        !   qs: qmc_state_t type, contains information about calculation.
        !   bloom_stats: statistics on blooms during calculation.
        !   contrib: derived type containing information on the current
        !       wavefunction contribution being considered.
        !   nattempts_spawn: running total of number of spawning attempts
        !       made during this mc cycle.

        use dSFMT_interface, only: dSFMT_t
        use system, only: sys_t
        use qmc_data, only: qmc_state_t, ccmc_in_t
        use ccmc_data, only: multispawn_stats_t, ms_stats_update, wfn_contrib_t
        use qmc_common, only: decide_nattempts

        use ccmc_death_spawning, only: spawner_ccmc, linked_spawner_ccmc, stochastic_ccmc_death
        use ccmc_death_spawning, only: stochastic_ccmc_death_nc, spawner_complex_ccmc
        use bloom_handler, only: bloom_stats_t, accumulate_bloom_stats
        use logging, only: logging_t

        type(sys_t), intent(in) :: sys
        type(dSFMT_T), intent(inout) :: rng
        type(qmc_state_t), intent(inout) :: qs
        type(ccmc_in_t), intent(in) :: ccmc_in
        type(wfn_contrib_t), intent(inout) :: contrib
        type(bloom_stats_t), intent(inout) :: bloom_stats
        type(logging_t), intent(in) :: logging_info

        integer(int_64), intent(inout) :: nattempts_spawn_tot

        integer :: i, nspawnings_cluster

        ! Spawning
        ! We select each non-composite cluster only once, then perform
        ! a number of spawning attempts proportional to the excip population
        ! on the corresponding excitor. This is similar to the previous
        ! approach but removes the requirement to have a deterministic
        ! method of deciding upon the number of attempts to make on each
        ! non-composite cluster previously imposed.

        ! We now use this to decide the number of attempts on each cluster
        ! stochastically, as in fciqmc.

        nspawnings_cluster = decide_nattempts(rng, abs(contrib%cluster%amplitude)/contrib%cluster%pselect)

        nattempts_spawn_tot = nattempts_spawn_tot + nspawnings_cluster

        contrib%cluster%amplitude = contrib%cluster%amplitude / abs(contrib%cluster%amplitude)

        do i = 1, nspawnings_cluster
            call perform_ccmc_spawning_attempt(rng, sys, qs, ccmc_in, logging_info, bloom_stats, contrib, 1)
        end do

    end subroutine do_nc_ccmc_propagation

    subroutine perform_ccmc_spawning_attempt(rng, sys, qs, ccmc_in, logging_info, bloom_stats, contrib, nspawnings_total)

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

        use dSFMT_interface, only: dSFMT_t
        use system, only: sys_t
        use qmc_data, only: qmc_state_t, ccmc_in_t
        use ccmc_data, only: wfn_contrib_t

        use excitations, only: excit_t
        use proc_pointers, only: gen_excit_ptr

        use ccmc_death_spawning, only: spawner_ccmc, linked_spawner_ccmc
        use ccmc_death_spawning, only: spawner_complex_ccmc, create_spawned_particle_ccmc
        use bloom_handler, only: bloom_stats_t, accumulate_bloom_stats
        use logging, only: logging_t

        type(sys_t), intent(in) :: sys
        type(dSFMT_T), intent(inout) :: rng
        type(qmc_state_t), intent(inout) :: qs
        type(ccmc_in_t), intent(in) :: ccmc_in
        type(wfn_contrib_t), intent(inout) :: contrib
        type(excit_t) :: connection
        type(bloom_stats_t), intent(inout) :: bloom_stats
        type(logging_t), intent(in) :: logging_info

        integer, intent(in) :: nspawnings_total

        integer(int_p) :: nspawned, nspawned_im
        integer(i0) :: fexcit(sys%basis%tot_string_len)

        if (contrib%cluster%excitation_level == huge(0)) then
            ! When sampling e^-T H e^T, the cluster operators in e^-T
            ! and e^T can excite to/from the same orbital, requiring
            ! a different spawning routine
            call linked_spawner_ccmc(rng, sys, qs, qs%spawn_store%spawn%cutoff, &
                      contrib%cluster, gen_excit_ptr, nspawned, connection, nspawnings_total, &
                      fexcit, contrib%cdet, contrib%ldet, contrib%rdet, contrib%left_cluster, contrib%right_cluster)
            nspawned_im = 0_int_p
        else if (sys%read_in%comp) then
            call spawner_complex_ccmc(rng, sys, qs, qs%spawn_store%spawn%cutoff, &
                      contrib%cdet, contrib%cluster, gen_excit_ptr, nspawned, nspawned_im, &
                      connection, nspawnings_total)
        else
            call spawner_ccmc(rng, sys, qs, qs%spawn_store%spawn%cutoff, &
                      ccmc_in%linked, contrib%cdet, contrib%cluster, gen_excit_ptr, logging_info, nspawned, &
                      connection, nspawnings_total)
            nspawned_im = 0_int_p
        end if

        if (nspawned /= 0_int_p) call create_spawned_particle_ccmc(sys%basis, qs%ref, contrib%cdet, connection, &
                                            nspawned, 1, contrib%cluster%excitation_level, &
                                            ccmc_in%even_selection, fexcit, qs%spawn_store%spawn, bloom_stats)
        if (nspawned_im /= 0_int_p) call create_spawned_particle_ccmc(sys%basis, qs%ref, contrib%cdet, connection,&
                                            nspawned_im, 2, contrib%cluster%excitation_level, &
                                            ccmc_in%even_selection, fexcit, qs%spawn_store%spawn, bloom_stats)

    end subroutine perform_ccmc_spawning_attempt

end module ccmc

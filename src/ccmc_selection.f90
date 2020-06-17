module ccmc_selection

! Module containing all ccmc selection subroutines.
! Three options for selection:

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

! A general overview of implementation details for each approach is given below.
! Guide to general terms: N_0 = reference population
!                         N_ex = absolute population on non-reference excips
!                         N_attempts = number of attempts made to generate a
!                                   cluster on a given iteration.

! Original CCMC algorithm
! -----------------------
! N_attempts = N_0 + N_ex
! p_size(N) = 1/2^(N+1)
! All excitors for clusters selected at random from whole list with probability
! proportional to absolute walker population on excitor.

! Full non-composite algorithm
! ----------------------------
! Overall N_attempts = N_0 + 2 N_ex
! N_0 attempts on the reference
! N_ex selections of non-composite clusters (made deterministically).
! N_ex selections of composite clusters. p_size(N) = 1/2^(N-1) for N>=2.
! Excitors for composite clusters selected at random from whole list with
! probability proportional to absolute walker population on excitor.

! Even Selection
! --------------
! Based upon the premise that all components of the wavefunction should be
! sampled evenly, regardless of provenance. In practice, this involves
! selecting individual clusters with probability proportional to the
! product of all constituent absolute excip populations.

! ie. for a cluster of size N
!       p_clust \propto \product_i=1^N |N_i|

! And weighting p_size according to the expected amplitude contribution of all
! possible clusters at a given cluster size. This is proportional to the product
! of the fraction of excips at each level used within a cluster. It also includes
! additional factors of N_0^(-N+1) and 1/(N!) due to our wavefunction anzatz and
! expansion.

! Once we obtain this p_size distribution, N_attempts can be set by requiring
!   N_attempts * p_size(1) = N_ex
! As we are utilising full non-composite selection for non-composite clusters.

! Thus described we have a functional method, but the combinatorially growing
! number of possible excitor combinations with increasing cluster sizes leads
! to a huge number of selections being required at these higher levels. This
! is not a problem in and of itself, but noting that the Hamiltonian is a
! two-electron operator we observe that the vast majority of these higher-size
! selections made in anything but CCSD will have excitation level greater than our
! truncation level + 2. In practice, this means these selections cannot affect
! our stored coefficients via death or spawning, and so are simply discarded.

! Remedying this situation is possible, by instead considering at each size
! all possible combinations of excitors at different excitation levels that
! results in a cluster of excitation level less than our truncation level + 2.
! Selecting a known combination of excitation levels requires easy selection of
! an excitor with known excitation level when acting upon the reference. This
! is achieved by sorting excitors according to excitation level, requiring
! extension of excitor bit string with the excitation level. This is achieved by
! modification of sys%basis%info_string_len in lua_hande_calc.f90 where needed.

! With these improvements, we observe a reduction in attempts required, especially
! in combination with the initiator approximation. This is an area of current active
! research.

use const, only: i0, int_p, int_32, int_64, debug, p, dp, depsilon

implicit none

contains

    subroutine select_cluster(rng, sys, psip_list, f0, ex_level, linked_ccmc, nattempts, normalisation, &
                              initiator_pop, cumulative_excip_pop, tot_excip_pop, min_size, max_size, &
                              logging_info, cdet, cluster, excit_gen_data)

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

        type(sys_t), intent(in) :: sys
        type(particle_t), intent(in), target :: psip_list
        integer(i0), intent(in) :: f0(sys%basis%tot_string_len)
        integer, intent(in) :: ex_level
        integer(int_64), intent(in) :: nattempts
        logical, intent(in) :: linked_ccmc
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
        integer :: i, pos, prev_pos
        real(p) :: pop(max_size)
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
        all_allowed = allowed

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
                ! If contain multiple spaces we can have this in a more general
                ! case, where an excitor has population in another space but not
                ! that which we're currently concerned with. More general test
                ! should account for this.
                do
                    if (pos == 1) exit
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
                call decoder_ptr(sys, cdet%f, cdet, excit_gen_data)

                ! Normalisation factor for cluster%amplitudes...
                cluster%amplitude = cluster_population/(normalisation**(cluster%nexcitors-1))
            else
                ! Simply set excitation level to a too high (fake) level to avoid
                ! this cluster being used.
                cluster%excitation_level = huge(0)
            end if

            if (.not.all_allowed) cluster%excitation_level = huge(0)

        end select

        if (debug) call write_logging_stoch_selection(logging_info, cluster%nexcitors, cluster%excitation_level, pop, &
                min(sys%nel, ex_level+2), cluster%pselect, cluster%amplitude, allowed)

    end subroutine select_cluster

    subroutine create_null_cluster(sys, f0, prob, D0_normalisation, initiator_pop, cdet, cluster, excit_gen_data)

        ! Create a cluster with no excitors in it, and set it to have
        ! probability of generation prob.

        ! In:
        !    sys: system being studied
        !    f0: bit string of the reference
        !    prob: The probability we set in it of having been generated
        !    D0_normalisation:  The number of excips at the reference (which
        !        will become the amplitude of this cluster)
        !    initiator_pop: the population above which a determinant is an initiator.
        !    excit_gen_data: info for excitation generators.

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
        use ccmc_data, only: cluster_t
        use proc_pointers, only: decoder_ptr
        use excit_gens, only: excit_gen_data_t

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f0(sys%basis%tot_string_len)
        real(p), intent(in) :: prob, initiator_pop
        complex(p), intent(in) :: D0_normalisation
        type(excit_gen_data_t), intent(in) :: excit_gen_data
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
        ! If not initiator something has gone seriously wrong and the CC
        ! approximation is (most likely) not suitably for this system.
        ! Let the user be an idiot if they want to be...
        if (abs(D0_normalisation) <= initiator_pop) cdet%initiator_flag = 3

        call decoder_ptr(sys, cdet%f, cdet, excit_gen_data)

    end subroutine create_null_cluster

    subroutine select_nc_cluster(sys, psip_list, f0, iexcitor, initiator_pop, ex_lvl_sort, &
                                            cdet, cluster, excit_gen_data)

        ! Select (deterministically) the non-composite cluster containing only
        ! the single excitor iexcitor and set the same information as select_cluster.

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
        use search, only: binary_search
        use proc_pointers, only: decoder_ptr
        use excit_gens, only: excit_gen_data_t

        type(sys_t), intent(in) :: sys
        type(particle_t), intent(in), target :: psip_list
        integer(i0), intent(in) :: f0(sys%basis%tot_string_len)
        integer(int_64), intent(in) :: iexcitor
        real(p), intent(in) :: initiator_pop
        logical, intent(in) :: ex_lvl_sort
        type(excit_gen_data_t), intent(in) :: excit_gen_data
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

        cluster%excitation_level = get_excitation_level(f0, cdet%f)
        cluster%amplitude = excitor_pop

        ! Sign change due to difference between determinant
        ! representation and excitors and excitation level.
        call convert_excitor_to_determinant(cdet%f, cluster%excitation_level, cluster%cluster_to_det_sign, f0)
        call decoder_ptr(sys, cdet%f, cdet, excit_gen_data)

    end subroutine select_nc_cluster

    subroutine select_cluster_truncated(rng, sys, psip_list, f0, linked_ccmc, nattempts, normalisation, initiator_pop, &
                        selection_data, cumulative_excip_pop, ex_level, min_size, max_size, ex_lvl_dist, cluster, cdet, &
                        excit_gen_data)

        ! Select (stochastically) cluster of size 2+ with restriction to ensure only
        ! select clusters with overall excitation level within truncation level+2.

        ! In:
        !    sys: system being studied
        !    psip_list: particle_t object containing current excip distribution on
        !       this processor.
        !    f0: bit string of the reference
        !    linked_ccmc: are we doing linked CCMC?
        !    nattempts: the number of times (on this processor) a random cluster
        !        of excitors is generated in the current timestep.
        !    normalisation: intermediate normalisation factor, N_0, where we use the
        !       wavefunction ansatz |\Psi_{CC}> = N_0 e^{T/N_0} | D_0 >.
        !    initiator_pop: the population above which a determinant is an initiator.
        !    selection_data: contains weightings used for selection.
        !    cumulative_excip_pop: running cumulative excip population on
        !        all excitors; i.e. cumulative_excip_population(i) = sum(particle_t%pops(1:i)).
        !    ex_level: maximum excitor excitation level to store coefficients for.
        !    min_size: minimum size of cluster to select stocastically
        !       (should always be 2).
        !    max_size: maximum size of cluster to select stochastically.
        !    ex_lvl_dist: derived type containing information on distribution of
        !       wavefunction between different excitation levels.
        !    excit_gen_data: data for excitation generators.

        ! NOTE: cumulative_excip_pop and tot_excip_pop ignore the population on the
        ! reference as excips on the reference cannot form a cluster.  Both these
        ! quantities should be generated by cumulative_population (or be in the same
        ! format).

        ! In/Out:
        !    rng: random number generator.
        !    cluster:
        !        Additional information about the cluster of excitors.  On
        !        input this is a bare cluster_t variable with the excitors array
        !        allocated to the maximum number of excitors in a cluster.  On
        !        output all fields in cluster have been set.
        !    cdet: information about the cluster of excitors applied to the
        !        reference determinant.  This is a bare det_info variable on input
        !        with only the relevant fields allocated.  On output the
        !        appropriate (system-specific) fields have been filled by
        !        decoding the bit string of the determinant formed from applying
        !        the cluster to the reference determinant.

        use ccmc_data, only: cluster_t, selection_data_t, ex_lvl_dist_t
        use ccmc_utils, only: convert_excitor_to_determinant
        use determinant_data, only: det_info_t
        use excitations, only: get_excitation_level
        use qmc_data, only: particle_t
        use dSFMT_interface, only: dSFMT_t, get_rand_close_open
        use proc_pointers, only: decoder_ptr
        use system, only: sys_t
        use utils, only: factorial
        use excit_gens, only: excit_gen_data_t
        use bit_utils, only: count_set_bits
        use errors, only: stop_all

        type(particle_t), intent(in), target :: psip_list
        type(sys_t), intent(in) :: sys
        type(cluster_t), intent(inout) :: cluster
        type(det_info_t), intent(inout) :: cdet
        real(p), intent(in) :: cumulative_excip_pop(:)
        type(dSFMT_t), intent(inout) :: rng
        integer(i0), intent(in) :: f0(sys%basis%tot_string_len)
        real(p), intent(in) :: initiator_pop
        complex(p), intent(in) :: normalisation
        integer(int_64), intent(in) :: nattempts
        integer, intent(in) :: min_size, max_size, ex_level
        logical, intent(in) :: linked_ccmc
        type(excit_gen_data_t), intent(in) :: excit_gen_data

        type(selection_data_t), intent(in) :: selection_data
        type(ex_lvl_dist_t), intent(in) :: ex_lvl_dist

        real(dp) :: rand, cumulative
        complex(p) :: cluster_population
        integer :: i, num, choice, iex
        logical :: first, allowed, all_allowed

        cluster%pselect = real(nattempts, p)

        ! Select cluster size according to probability distribution set by
        ! update_selection_probabilites.
        rand = get_rand_close_open(rng)
        cluster%nexcitors = -1
        do i = min_size, max_size - 1
            ! Add equality for edge case.
            if (rand < selection_data%cumulative_stoch_size_weighting(i)) then
                cluster%nexcitors = i
                cluster%pselect = cluster%pselect * selection_data%stoch_size_weighting(i)
                exit
            end if
        end do

        if (cluster%nexcitors == -1) then
                cluster%nexcitors = max_size
                cluster%pselect = cluster%pselect * selection_data%stoch_size_weighting(max_size)
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
        ! valid.
        allowed = (min_size <= max_size)

        ! For linked coupled cluster we keep building the cluster after a
        ! disallowed excitation so need to know if there has been a disallowed
        ! excitation at all.
        all_allowed = allowed

        ! Now choose the combination of excitor-levels, given the number of excitors in the cluster.
        associate(select_info => selection_data%cluster_sizes_info(cluster%nexcitors)%v, &
                select_proportion => selection_data%cluster_sizes_proportion(cluster%nexcitors)%v)
            rand = get_rand_close_open(rng)

            choice = -1
            cumulative = 0.0_dp
            do num = 1, size(select_proportion, dim=1) - 1
                cumulative = cumulative + select_proportion(num)
                if (rand < cumulative) then
                    choice = num
                    exit
                end if
            end do

            if (choice == -1) then
                choice = size(select_proportion, dim=1)
            end if

            cluster%pselect = cluster%pselect * select_proportion(choice)

            ! choice indexes the combinations, and select_info(choice, i) is the number of
            ! excitors with level i needed in this particular combination.
            first = .true.
            iex = 1
            accumulate_excitors : do i = 1, size(select_info, dim=2)
                if (select_info(choice, i) /= 0) then

                    call add_excitors_given_level(rng, sys, psip_list, f0, linked_ccmc, initiator_pop, first, i, &
                            select_info(choice, i), iex, cumulative_excip_pop, ex_lvl_dist, cluster_population, &
                            cluster, cdet, allowed)
                    if (.not.allowed) then
                        if (.not.linked_ccmc) exit accumulate_excitors
                        all_allowed = .false.
                    end if
                    first = .false.
                    ! We chose excitors with a probability proportional to their
                    ! occupation.  However, because (for example) the cluster t_X t_Y
                    ! and t_Y t_X collapse onto the same excitor (where X and Y each
                    ! label an excitor), the probability of selecting a given cluster is
                    ! proportional to the number of ways the cluster could have been
                    ! formed. Obviously this now only applies to clusters of the same
                    ! excitation level.
                    ! If two excitors in the cluster are the same, the factorial
                    ! overcounts the number of ways the cluster could have been formed
                    ! but the extra factor(s) of 2 are cancelled by a similar
                    ! overcounting in the calculation of hmatel.
                    cluster%pselect = cluster%pselect*factorial(select_info(choice, i))
                end if
            end do accumulate_excitors

            if (allowed) then
                cluster%excitation_level = get_excitation_level(f0, cdet%f)
                ! To contribute the cluster must be within a double excitation of
                ! the maximum excitation included in the CC wavefunction.
                allowed = cluster%excitation_level <= ex_level+2
            end if

            if (allowed .or. linked_ccmc) then
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

            if (.not.all_allowed) cluster%excitation_level = huge(0)

        end associate

    end subroutine select_cluster_truncated

    subroutine add_excitors_given_level(rng, sys, psip_list, f0, linked_ccmc, initiator_pop, first, ex_level, nexcitors, iexcitor, &
                                        cumulative_excip_pop, ex_lvl_dist, cluster_population, cluster, cdet, allowed)
        ! Adds excitors of a given level to cluster being created.
        ! In:
        !   sys: system being studied.
        !   psip_list: particle_t object containing current excip distribution on this
        !       processor.
        !   f0: bit string of the reference.
        !   linked_ccmc: are we doing linked CCMC?
        !   initiator_pop: the population above which a determinant is an initiator.
        !   first: if true first cluster selected will seed the cluster
        !   ex_level: excitation level from the reference of the excitors to add
        !   nexcitors: number of excitors of this excitation level to add to cluster.
        !   cumulative_excip_pop: running excip population on all excitors.
        !   ex_lvl_dist: derived type containing information on excip distribution between
        !       excitation levels.
        ! In/Out:
        !   rng: random number generator.
        !   iexcitor: number of excitor being added - this is incremented with each added excitor
        !   cluster_population: current population of cluster.
        !   cluster: cluster of excitors currently being accumulated.
        !   cdet: information anout the cluster of excitors applied to the reference determinant.
        ! Out:
        !   allowed: false if collapsed cluster is invalid (ie. tried to excite to/from
        !       same spinorbital twice).
        use ccmc_data, only: cluster_t, ex_lvl_dist_t
        use ccmc_utils, only: collapse_cluster
        use determinant_data, only: det_info_t
        use qmc_data, only: particle_t
        use dSFMT_interface, only: dSFMT_t, get_rand_close_open
        use sort, only: insert_sort
        use search, only: binary_search
        use system, only: sys_t
        use basis_types, only: reset_extra_info_bit_string
        use parallel, only: nprocs

        logical, intent(in) :: first, linked_ccmc
        type(particle_t), intent(in), target :: psip_list
        integer, intent(in) :: ex_level, nexcitors
        integer, intent(inout) :: iexcitor
        type(sys_t), intent(in) :: sys
        type(cluster_t), intent(inout) :: cluster
        type(det_info_t), intent(inout) :: cdet
        real(p), intent(in) :: cumulative_excip_pop(:)
        type(dSFMT_t), intent(inout) :: rng
        complex(p), intent(inout) :: cluster_population
        integer(i0), intent(in) :: f0(sys%basis%tot_string_len)
        real(p), intent(in) :: initiator_pop
        type(ex_lvl_dist_t), intent(in) :: ex_lvl_dist

        integer :: i, pos, prev_pos
        real(p) :: pops(1:nexcitors), tot_level_pop
        logical :: hit
        logical, intent(inout) :: allowed
        complex(p) :: excitor_pop

        ! Calculate total population of excips on states with this excit level.
        tot_level_pop = ex_lvl_dist%pop_ex_lvl(ex_level)

        if (nexcitors /= 0) then
            ! Want to generate populations in range:
            !       cumulative_pop_ex_lvl(ex_lvl-1) < pop <= cumulative_pop_ex_lvl(ex_lvl)
            do i = 1, nexcitors
                pops(i) = get_rand_close_open(rng)*(tot_level_pop-depsilon) + depsilon + &
                                ex_lvl_dist%cumulative_pop_ex_lvl(ex_level - 1)
            end do
            call insert_sort(pops)
            ! Lowest possible position is first state of chosen ex lvl ie. cumulative
            ! states to one lower ex lvl plus 1.
            prev_pos = ex_lvl_dist%cumulative_nstates_ex_lvl(ex_level - 1) + 1
            do i = 1, nexcitors
                call binary_search(cumulative_excip_pop, pops(i), prev_pos, &
                            ex_lvl_dist%cumulative_nstates_ex_lvl(ex_level), hit, pos)
                if (sys%read_in%comp) then
                    excitor_pop = cmplx(psip_list%pops(1,pos),psip_list%pops(2,pos), p)/psip_list%pop_real_factor
                else
                    excitor_pop = cmplx(psip_list%pops(1,pos), 0.0_p, p)/psip_list%pop_real_factor
                end if
                if (i == 1 .and. first) then
                    ! First excitor 'seeds' the cluster:
                    cdet%f = psip_list%states(:,pos)
                    call reset_extra_info_bit_string(sys%basis, cdet%f)
                    cdet%data => psip_list%dat(:,pos) ! Only use if cluster is non-composite!
                    cluster_population = excitor_pop
                else
                    call collapse_cluster(sys%basis, f0, psip_list%states(:,pos), excitor_pop, cdet%f, &
                                          cluster_population, allowed)

                    if (.not.allowed .and. .not. linked_ccmc) exit
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
                ! Probability of choosing this excitor = abs(pop)/tot_pop.
                cluster%pselect = (cluster%pselect*abs(excitor_pop))/tot_level_pop
                cluster%excitors(iexcitor)%f => psip_list%states(:,pos)
                prev_pos = pos
                iexcitor = iexcitor + 1
            end do
        end if

    end subroutine add_excitors_given_level

!---- nselections update function ----

    subroutine set_cluster_selections(selection_data, nattempts, min_cluster_size, max_size, D0_normalisation, tot_abs_pop, &
                                    nstates, full_nc, even_selection)

        ! Function to set total number of selections of different cluster
        ! types within CCMC. This effectively controls the relative sampling
        ! of different clusters within the CC expansion.

        ! In:
        !   max_size: maximum cluster size we need to select.
        !   nstates: total number of occupied states within calculation.
        !   D0_normalisation: total population on the reference this iteration.
        !   tot_abs_pop: sum of absolute excip populations on all excitors.
        !   full_nc: whether using full non-composite algorithm.
        !   even_selection: whether using even selection algorithm for
        !       stochastic clusters.
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
        integer, intent(in) :: max_size, nstates
        complex(p), intent(in) :: D0_normalisation
        real(p), intent(in) :: tot_abs_pop
        logical, intent(in) :: full_nc, even_selection

        if (full_nc .or. even_selection) then
            ! Note that nattempts /= tot_abs_nint_pop+D0_normalisation if the
            ! reference is not on the current processor.  Instead work
            ! out how many clusters of each type we will sample
            ! explicitly.
            min_cluster_size = 2_int_32
            selection_data%nD0_select = nint(abs(D0_normalisation),kind=int_64)
            selection_data%nstochastic_clusters = ceiling(tot_abs_pop,kind=int_64)
            selection_data%nsingle_excitors = int(nstates,kind=int_64)
            if (even_selection) then
                if (max_size > 1) then
                    ! Set total selections so that expected proportion of selections of noncomposite gives
                    ! correct number of selections.
                    nselections = ceiling(tot_abs_pop / selection_data%size_weighting(1), kind=int_64)
                    selection_data%nD0_select = ceiling(nselections * selection_data%size_weighting(0), kind=int_64)
                    selection_data%nstochastic_clusters = ceiling(nselections * sum(selection_data%size_weighting(2:)), kind=int_64)
                else
                    selection_data%nstochastic_clusters = 0_int_64
                end if
            end if
            nattempts = nint(tot_abs_pop, kind=int_64) + selection_data%nD0_select + selection_data%nstochastic_clusters
        else
            min_cluster_size = 0
            selection_data%nD0_select = 0 ! instead of this number of deterministic selections, these are chosen stochastically
            selection_data%nstochastic_clusters = nattempts
            selection_data%nsingle_excitors = 0
        end if

        selection_data%nclusters = selection_data%nD0_select + selection_data%nsingle_excitors &
                            + selection_data%nstochastic_clusters

    end subroutine set_cluster_selections

!---- Cluster information accumulation ---

    subroutine update_selection_data(selection_data, cluster, logging_info)

        ! Updates selection_data derived type with information relating to amp/pselect. Only required
        ! for use with selection logging.
        ! Information accumulated is the mean and mean square values, used to obtain the average and
        ! variance/standard deviation of the amplitude distribution.

        ! In:
        !   cluster: info on currently selected cluster to add to accumulated information.
        !   logging_info: info on logging settings currently in use.
        ! In/Out:
        !   selection_data: information on cluster selection within calculation. On output will have
        !       contribution from cluster added into mean and mean square value accumulation.

        use ccmc_data, only: cluster_t, selection_data_t
        use logging, only: logging_t

        type(selection_data_t), intent(inout) :: selection_data
        type(logging_t), intent(in) :: logging_info
        type(cluster_t), intent(in) :: cluster

        if (logging_info%write_amp_psel) then
            associate( nex=>cluster%nexcitors, amp=>real(cluster%amplitude,kind=dp)/real(cluster%pselect,kind=dp), &
                aa=>selection_data%average_amplitude, va=>selection_data%variance_amplitude, &
                nsuccess=>selection_data%nsuccessful)

                if (cluster%excitation_level /= huge(0)) then
                    aa(nex) = aa(nex) * (real(nsuccess(nex),kind=dp)/real(nsuccess(nex)+1_int_64,kind=dp)) &
                        + abs(amp)/ real(nsuccess(nex)+1_int_64, kind=dp)
                    va(nex) = va(nex) * (real(nsuccess(nex),kind=dp)/real(nsuccess(nex)+1_int_64,kind=dp)) &
                        + amp**2/ real(nsuccess(nex)+1_int_64, kind=dp)
                    nsuccess(nex) = nsuccess(nex) + 1_int_64
                end if
            end associate
        end if

    end subroutine update_selection_data

! --- From on here all functions are only used for truncated selection ---

!---- p_comb and p_size Probability update functions ----

    subroutine update_selection_probabilities(ex_lvl_dist, abs_D0_normalisation, tot_abs_pop, cluster_selection)

        ! Updates all probabilities for selecting different excitation level combinations within
        ! cluster_selection object in accordance with cumulative population distribution and
        ! number of states per excitation level given.
        ! In:
        !    ex_lvl_dist: derived type containing information on distribution of excip population
        !       between excitation levels.
        !    abs_D0_normalisation: absolute magnitude of D0 normalisation.
        !    tot_abs_pop: absolute sum of cumulative population.
        ! In/Out:
        !    cluster_selection: selection_data_t object containing all information required for
        !       truncated selection (all valid excitation level combinations). On output
        !       cluster_sizes_proportion will be updated.

        use ccmc_data, only: selection_data_t, ex_lvl_dist_t
        use ccmc_utils, only: update_cumulative_dist_real
        use parallel, only: nprocs

        type(selection_data_t), intent(inout) :: cluster_selection
        type(ex_lvl_dist_t), intent(in) :: ex_lvl_dist
        real(p), intent(in) :: tot_abs_pop
        real(p), intent(in) :: abs_D0_normalisation
        integer :: i

        cluster_selection%size_weighting(0) = abs_D0_normalisation
        cluster_selection%size_weighting(1) = real(tot_abs_pop, kind=dp)

        ! To sample between all combinations effectively requires selecting:
        !   - individual excitors with probability proportional to absolute excip population.
        !   - different excitation level combinations with probability proportional to the proportion
        !       of absolute excip population on constituent excitation levels, accounting for the
        !       number of possible ways of selecting the same cluster.
        !   - different sizes in proportion to the total chance of selecting all permitted excitation
        !       level combinations of a given size.

        ! This is achieved below via update_selection_block_probability calculating the weighting for
        ! all combinations of a given size given the current excitation distribution, then applying
        ! reweighting for the reference normalisation and the number of processors running to the
        ! probability of selecting each size. After normalisation this gives a normalised probability of
        ! selecting each cluster size and each combination of excitation levels.

        do i = lbound(cluster_selection%cluster_sizes_info, dim=1), ubound(cluster_selection%cluster_sizes_info, dim=1)
            associate(select_info => cluster_selection%cluster_sizes_info(i)%v, &
                      select_proportion => cluster_selection%cluster_sizes_proportion(i)%v)

                call update_selection_block_probability(select_info, ex_lvl_dist%pop_ex_lvl, &
                                                    select_proportion, cluster_selection%size_weighting(i))
                ! Must include factor of D0_normalisation from wavefunction anzatx expansion and factor of
                ! 1/(nprocs^(i-1)) to account for the probability of two excitors being on the same
                ! proccessor on any given iteration.
                cluster_selection%size_weighting(i) = cluster_selection%size_weighting(i) * &
                                                        ((nprocs/abs_D0_normalisation) ** (i-1))
            end associate
        end do

        cluster_selection%size_weighting = cluster_selection%size_weighting / sum(cluster_selection%size_weighting, dim=1)

        call update_cumulative_dist_real(cluster_selection%size_weighting, cluster_selection%cumulative_size_weighting, &
                                        normalised=.true.)
        ! Now need to update lists storing only sizes selected stochastically.
        cluster_selection%stoch_size_weighting(2:) = cluster_selection%size_weighting(2:)/sum(cluster_selection%size_weighting(2:))
        call update_cumulative_dist_real(cluster_selection%stoch_size_weighting, &
                                        cluster_selection%cumulative_stoch_size_weighting, normalised=.true.)

    end subroutine update_selection_probabilities

!---- Initialisation routines for improved selection information ----

    subroutine init_selection_data(ex_level, max_cluster_size, selection_data)

        ! Take cluster selection object and initialise all data required for calculation.

        ! In:
        !   ex_level: maximum excitation level allowed for stored coefficients in calculation.
        !   max_cluster_size: maximum allowed cluster size.
        ! In/Out:
        !   cluster_selection: selection_data_t object. On output cluster_sizes_info components
        !       will be allocated and set as appropriate, and cluster_sizes_proportion allocated.

        use ccmc_data, only: selection_data_t
        use checking, only: check_allocate

        integer, intent(in) :: ex_level, max_cluster_size
        type(selection_data_t), intent(inout) :: selection_data

        call init_possible_clusters(ex_level, max_cluster_size, selection_data)

        call init_psize_data(max_cluster_size, selection_data)

    end subroutine init_selection_data

    subroutine init_possible_clusters(ex_level, max_cluster_size, selection_data)

        ! Take cluster selection object and initialise all possible size combinations of clusters.

        ! In:
        !   ex_level: maximum excitation level allowed in calculation.
        !   max_cluster_size: maximum allowed cluster size.
        ! In/Out:
        !   cluster_selection: selection_data_t object. On output cluster_sizes_info components
        !       will be allocated and set as appropriate, and cluster_sizes_proportion allocated.

        use ccmc_data, only: selection_data_t
        use checking, only: check_allocate
        use parallel, only: parent
        integer, intent(in) :: ex_level, max_cluster_size
        type(selection_data_t), intent(inout) :: selection_data

        integer :: ncombinations, cluster_size, ierr, dummy
        integer :: temporary(1:ex_level), i, j

        allocate(selection_data%cluster_sizes_info(2:max_cluster_size), stat=ierr)
        call check_allocate('cluster_selection%cluster_sizes_info',max_cluster_size-2, ierr)
        allocate(selection_data%cluster_sizes_proportion(2:max_cluster_size), stat=ierr)
        call check_allocate('cluster_selection%cluster_sizes_proportion',max_cluster_size-2, ierr)

        if (parent) then
            write (6, '(1X,"Truncated Selection Initialisation")')
            write (6, '(1X,"----------------------------------",/)')
            write (6, '(1X,"Setting up required data storage to sample all composite clusters of size <= ",i0,", ")') &
                                        max_cluster_size
            write (6, '(1X,"cluster excitation level <= ",i0," using excitors of excitation level <= ",i0,".",/)') &
                                        ex_level+2, ex_level
        end if

        ! Go through possible sizes...
        do cluster_size = 2, max_cluster_size
            ! Calc number possible clusters for each
            ncombinations = calc_available_perms(cluster_size, 1, 0, ex_level)
            ! Allocate storage of appropriate size to store all cluster info
            allocate(selection_data%cluster_sizes_info(cluster_size)%v(1:ncombinations,1:ex_level), stat=ierr)
            call check_allocate('cluster_selection_component', ncombinations*ex_level,ierr)

            allocate(selection_data%cluster_sizes_proportion(cluster_size)%v(1:ncombinations), stat=ierr)
            call check_allocate('cluster_selection_component', ncombinations,ierr)

            selection_data%cluster_sizes_info(cluster_size)%v(:,:) = 0
            selection_data%cluster_sizes_proportion(cluster_size)%v(:) = 0
            dummy = 0
            temporary(:) = 0
            call find_available_perms(dummy, selection_data%cluster_sizes_info(cluster_size)%v, temporary, &
                                    cluster_size, 1, 0, ex_level)
            if (parent) then
                write(6, '(1X,"Found ",i0," possible excitation level combinations for a cluster of size ",i0,".")') &
                                                            ncombinations, cluster_size
                write(6, '(1X,"Combinations are:",/)')
                write(6, '(12X,"|",5X,a30)') "N_excitors @ excitation level:"
                write(6, '(a11,1X,"|",41("-"))') "Combo"
                write(6, '(a11,1X,"|",1X)', advance='no') "Number"
                do j = 1, ex_level
                    write(6, '(1X,a9,i2,1X)', advance='no') "ex level=", j
                end do
                write(6, '(1X)')
                write(6, '(4X,50("-"))')
                do i = 1, ncombinations
                    write (6, '(9X,i2,1X,"|",1X)',advance='no') i
                    do j = 1, ex_level
                        write(6, '(1X,5X,i2,5X)', advance='no') selection_data%cluster_sizes_info(cluster_size)%v(i, j)
                    end do
                    write(6, '(1x)')
                end do
                write(6, '(1x)')
            end if
        end do

    end subroutine init_possible_clusters

    subroutine init_psize_data(max_cluster_size, selection_data)

        ! Take cluster selection object and initialise all data required for psize variation.

        ! In:
        !   max_cluster_size: maximum allowed cluster size.
        ! In/Out:
        !   cluster_selection: selection_data_t object. On output cluster_sizes_info components
        !       will be allocated and set as appropriate, and cluster_sizes_proportion allocated.

        use ccmc_data, only: selection_data_t
        use checking, only: check_allocate

        integer, intent(in) :: max_cluster_size
        type(selection_data_t), intent(inout) :: selection_data
        integer :: ierr, i

        allocate(selection_data%size_weighting(0:max_cluster_size), stat=ierr)
        call check_allocate('size_weighting', max_cluster_size + 1, ierr)
        allocate(selection_data%cumulative_size_weighting(0:max_cluster_size), stat=ierr)
        call check_allocate('cumulative_size_weighting', max_cluster_size + 1, ierr)

        allocate(selection_data%stoch_size_weighting(2:max_cluster_size), stat=ierr)
        call check_allocate('stoch_size_weighting', max_cluster_size - 1, ierr)
        allocate(selection_data%cumulative_stoch_size_weighting(2:max_cluster_size), stat=ierr)
        call check_allocate('cumulative_stoch_size_weighting', max_cluster_size - 1, ierr)

        associate(size_weighting => selection_data%size_weighting, &
                cumulative_size_weighting => selection_data%cumulative_size_weighting)

            size_weighting(0) = 0.5_dp
            cumulative_size_weighting(0) = 0.5_dp

            do i = 1, max_cluster_size -1
                size_weighting(i) = 1.0_dp/(2.0_dp**(i+1))
                cumulative_size_weighting(i) = cumulative_size_weighting(i-1) + size_weighting(i)
            end do

            size_weighting(max_cluster_size) = 1.0_dp - cumulative_size_weighting(max_cluster_size-1)
            cumulative_size_weighting(max_cluster_size) = 1.0_dp

        end associate

    end subroutine init_psize_data

    subroutine init_amp_psel_accumulation(max_cluster_size, logging_info, linked_ccmc, selection_data)

        ! Take selection_data derived type and initialise data strunctures required to accumulate
        ! amp/pselect for different cluster sizes.

        ! In:
        !   max_cluster_size: maximum allowed cluster size.
        ! In/Out:
        !   cluster_selection: selection_data_t object. On output cluster_sizes_info components
        !       will be allocated and set as appropriate.


        use ccmc_data, only: selection_data_t
        use checking, only: check_allocate
        use logging, only: logging_t

        integer, intent(in) :: max_cluster_size
        type(logging_t), intent(in) :: logging_info
        logical, intent(in) :: linked_ccmc
        type(selection_data_t), intent(inout) :: selection_data
        integer :: ierr, max_cluster_size_loc

        max_cluster_size_loc = max_cluster_size
        if (linked_ccmc) max_cluster_size_loc = 4

        if (logging_info%write_amp_psel) then
            allocate(selection_data%average_amplitude(0:max_cluster_size_loc), stat=ierr)
            call check_allocate('average_amplitude', max_cluster_size_loc + 1, ierr)
            allocate(selection_data%variance_amplitude(0:max_cluster_size_loc), stat=ierr)
            call check_allocate('variance_amplitude', max_cluster_size_loc + 1, ierr)
            allocate(selection_data%nsuccessful(0:max_cluster_size_loc), stat=ierr)
            call check_allocate('nsuccessful', max_cluster_size_loc + 1, ierr)

            selection_data%average_amplitude = 0.0_dp
            selection_data%variance_amplitude = 0.0_dp
            selection_data%nsuccessful = 0_int_64
        end if

    end subroutine init_amp_psel_accumulation

!---- Cluster generation functions ----
! Below functions generate all allowed excitation level combinations for clusters satifying given
! initial conditions.
! This is achieved by starting from lowest excitation level of excitor considered, iterating through all
! permissible numbers of said excitors, and depending upon the properties of any suggested combination
! rejecting it, accepting it or calling the same function for the next excitation level to search for
! acceptable combinations.

     pure recursive function calc_available_perms(selections_remain, excitor_excit_level, current_excit_level, &
                                 max_excit_level) result(nfound)
         ! For given number of selections and excitations remaining,
         ! calculate number of possible clusters. Performs this via a recursive
         ! search of all combinations.

         ! In:
         !   selections_remain: total number more excitors to be added to the
         !       cluster to reach chosen total size.
         !   excitor_excit_level: excitation level of excitors currently being
         !       considered for addition to the cluster.
         !   current_excit_level: excitation level of partial cluster as it stands.
         !   max_excit_level: maximum excitation level permitted in calculation.
         ! Returns:
         !   Number of possible size combinations fitting given initial conditions.

         integer, intent(in) :: selections_remain, excitor_excit_level, current_excit_level
         integer, intent(in) :: max_excit_level
         integer :: nfound
         integer :: new_selections_remain, new_current_excit_level
         integer :: i

         nfound = 0

         if (max_excit_level + 2 - current_excit_level >= excitor_excit_level * selections_remain) then
             ! Is still possible to form a valid cluster.
             do i = 0, selections_remain
                 new_selections_remain = selections_remain - i
                 new_current_excit_level = current_excit_level + i * (excitor_excit_level)
                 if (new_selections_remain == 0) then
                     ! Have selected appropriate number of excitors to give cluster
                     ! of required size without exceeding excitation level; increase.
                     nfound = nfound + 1
                 else if (excitor_excit_level < max_excit_level) then
                     ! Haven't exceeded allowed excitation level, so have to
                     ! continue searching at higher levels.
                     nfound = nfound + calc_available_perms(new_selections_remain, excitor_excit_level + 1, &
                                 new_current_excit_level, max_excit_level)
                 end if
             end do
         end if

     end function calc_available_perms

     pure recursive subroutine find_available_perms(nfound, main_storage, temporary, selections_remain, &
                                    excitor_excit_level, current_excit_level, max_excit_level)
        ! For given number of selections and max excitation level calculate
        ! all possible excitation level combinations forming a valid
        ! cluster recursively. Once sufficient excitors are selected the overall
        ! combination is stored in 'main_storage' in the next free slot.

        ! main_storage should be allocated to the correct size before
        ! initialisation; calc_available_perms enables calculations of the
        ! size required.
        ! In:
        !   temporary: 1d integer array containing information on excitation
        !       level combination accrued so far in recursive search.
        !   selections_remain: total number more excitors to be added to the
        !       cluster to reach chosen total size.
        !   excitor_excit_level: excitation level of excitors currently being
        !       considered for addition to the cluster.
        !   current_excit_level: excitation level of cluster as currently
        !       stands.
        !   max_excit_level: maximum excitation level permitted in calculation.
        ! In/Out:
        !   nfound: number of allowed combinations already found.
        !   main_storage: 2d integer array to store all found combinations.
        !       main_storage(j,k) stores number of individual excitors of
        !       excitation level k in combination j.
        integer, intent(inout) :: main_storage(:,:)
        integer, intent(in) :: temporary(:), max_excit_level, excitor_excit_level
        integer, intent(inout) :: nfound
        integer, intent(in) :: selections_remain, current_excit_level
        integer :: new_selections_remain, new_current_excit_level
        integer :: i
        integer :: temporary_loc(1:max_excit_level)

        temporary_loc = temporary

        if (max_excit_level + 2 - current_excit_level >= excitor_excit_level * selections_remain) then
            ! Is still possible to form a valid cluster.
            do i = 0, selections_remain
                new_selections_remain = selections_remain - i
                new_current_excit_level = current_excit_level + i * excitor_excit_level
                if (new_selections_remain == 0) then
                    temporary_loc(excitor_excit_level) = i
                    nfound = nfound + 1
                    main_storage(nfound, :) = temporary_loc(:)
                else if (excitor_excit_level < max_excit_level) then
                    temporary_loc(excitor_excit_level) = i
                    call find_available_perms(nfound, main_storage, temporary_loc, new_selections_remain, &
                            excitor_excit_level + 1, new_current_excit_level, max_excit_level)
                end if
            end do
        end if

    end subroutine find_available_perms

!---- Helper Functions -----

    pure subroutine update_selection_block_probability(cluster_info, &
                        excitation_level_populations, cluster_proportion, size_weighting)
        ! Updates single block of selection_data to have selection probability
        ! proportional to occupancy of excitation levels in given cluster. Each
        ! block corresponds to a single size of cluster.

        ! In:
        !   cluster_info: array containing in element (i,j) the number of
        !       excitors of excitation level j within cluster combination
        !       i, where the array contains all combinations of a given
        !       size.
        !   excitation_level_populations: population at each excitation level.
        ! Out:
        !   cluster_proportion: array containing normalised proportion of selections
        !       of clusters of given size that should be a given combination.
        !   size_weighting: overall value of weighting for this cluster size in our
        !       calculation.

        real(dp), intent(out) :: cluster_proportion(:)
        integer, intent(in) :: cluster_info(:,:)
        real(p), intent(in), allocatable :: excitation_level_populations(:)
        integer :: i
        real(dp), intent(out) :: size_weighting

        do i = lbound(cluster_info, dim=1), ubound(cluster_info, dim=1)
            cluster_proportion(i) = calc_combination_weighting(cluster_info(i, :), excitation_level_populations(1:))
        end do

        size_weighting = sum(cluster_proportion)

        if (sum(cluster_proportion) > depsilon) then
            cluster_proportion = cluster_proportion / sum(cluster_proportion)
        else
            cluster_proportion = 1.0_dp / real(size(cluster_proportion, dim=1))
        end if

    end subroutine update_selection_block_probability

    pure function calc_combination_weighting(ar1, v1) result(res)

        ! Calculates overall weighting of combination of excitation levels
        ! ar1 given population distribution by population level of v1.

        ! Combines two vectors together such that
        !   res = product_{i, v1(i)=/0} v1_{i} ** ar1{i} / ar1{i}!

        use utils, only: factorial

        integer, intent(in) :: ar1(:)
        real(p), intent(in) :: v1(:)
        real(dp) :: res
        real(dp) :: v2(size(v1,dim=1))
        integer :: i
        logical :: mask(size(ar1,dim=1))

        mask = .true.

        do i = lbound(ar1, dim=1), ubound(ar1, dim=1)
           v2(i) = real(v1(i),kind=dp) ** ar1(i) / factorial(ar1(i))
           if (ar1(i) == 0_int_32) mask(i) = .false.
        end do

        res = real(product(v2, dim=1, mask=mask), dp)

    end function calc_combination_weighting

end module ccmc_selection

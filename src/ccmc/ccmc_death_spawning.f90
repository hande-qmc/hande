module ccmc_death_spawning

! Module containing all ccmc death and spawning routines, for all possible propogation methods.
! For full explanation see top of ccmc.F90.

use const, only: i0, int_p, int_64, p, dp, debug

implicit none

contains

    subroutine spawner_ccmc(rng, sys, qs, spawn_cutoff, linked_ccmc, cdet, cluster, &
                            gen_excit_ptr, logging_info, nspawn, connection, nspawnings_total)

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
        !    logging_info: logging_t derived type containing information on logging behaviour.
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
        use ccmc_utils, only: convert_excitor_to_determinant
        use ccmc_linked, only: unlinked_commutator, linked_excitation
        use determinants, only: det_info_t
        use dSFMT_interface, only: dSFMT_t
        use excitations, only: excit_t, create_excited_det, get_excitation_level
        use proc_pointers, only: gen_excit_ptr_t
        use spawning, only: attempt_to_spawn, calc_qn_spawned_weighting
        use system, only: sys_t
        use const, only: depsilon, debug
        use qmc_data, only: qmc_in_t, qmc_state_t
        use hamiltonian_data
        use logging, only: logging_t, write_logging_spawn

        type(sys_t), intent(in) :: sys
        type(qmc_state_t), intent(in) :: qs
        integer(int_p), intent(in) :: spawn_cutoff
        logical, intent(in) :: linked_ccmc
        type(det_info_t), intent(in) :: cdet
        type(cluster_t), intent(in) :: cluster
        type(dSFMT_t), intent(inout) :: rng
        integer, intent(in) :: nspawnings_total
        type(gen_excit_ptr_t), intent(in) :: gen_excit_ptr
        type(logging_t), intent(in) :: logging_info
        integer(int_p), intent(out) :: nspawn
        type(excit_t), intent(out) :: connection

        ! We incorporate the sign of the amplitude into the Hamiltonian matrix
        ! element, so we 'pretend' to attempt_to_spawn that all excips are
        ! actually spawned by positive excips.
        integer(int_p), parameter :: parent_sign = 1_int_p
        type(hmatel_t) :: hmatel, hmatel_save
        real(p) :: pgen
        integer(i0) :: fexcit(sys%basis%string_len), funlinked(sys%basis%string_len)
        integer :: excitor_sign, excitor_level
        logical :: linked, single_unlinked, allowed_excitation
        real(p) :: invdiagel

        ! 1. Generate random excitation.
        ! Note CCMC is not (yet, if ever) compatible with the 'split' excitation
        ! generators of the sys%lattice%lattice models.  It is trivial to implement and (at
        ! least for now) is left as an exercise to the interested reader.
        call gen_excit_ptr%full(rng, sys, qs%excit_gen_data, cdet, pgen, connection, hmatel, allowed_excitation)

        if (allowed_excitation) then
            if (linked_ccmc) then
                ! For Linked Coupled Cluster we reject any spawning where the
                ! Hamiltonian is not linked to every cluster operator
                ! The matrix element to be evaluated is not <D_j|H a_i|D0> but <D_j|[H,a_i]|D0>
                ! (and similarly for composite clusters)
                if (cluster%nexcitors > 0) then
                    ! Check whether this is an unlinked diagram - if so, the matrix element is 0 and
                    ! no spawning is attempted
                    call linked_excitation(sys%basis, qs%ref%f0, connection, cluster, linked, single_unlinked, funlinked)
                    if (.not. linked) then
                        hmatel%r = 0.0_p
                    else if (single_unlinked) then
                        ! Single excitation: need to modify the matrix element
                        ! Subtract off the matrix element from the cluster without
                        ! the unlinked a_i operator
                        hmatel%r = hmatel%r - unlinked_commutator(sys, qs%ref%f0, connection, cluster, cdet%f, funlinked)
                    end if
                end if
            end if
            invdiagel = calc_qn_spawned_weighting(sys, qs, cdet%fock_sum, connection)
        else
            invdiagel = 1
        end if
        ! 2, Apply additional factors.
        hmatel_save = hmatel
        hmatel%r = hmatel%r*cluster%amplitude*invdiagel*cluster%cluster_to_det_sign
        pgen = pgen*cluster%pselect*nspawnings_total

        ! 3. Attempt spawning.
        nspawn = attempt_to_spawn(rng, qs%tau, spawn_cutoff, qs%psip_list%pop_real_factor, hmatel%r, pgen, parent_sign)

        if (debug) call write_logging_spawn(logging_info, hmatel_save, pgen, invdiagel, [nspawn], &
                        cluster%amplitude*qs%psip_list%pop_real_factor, sys%read_in%comp)

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

    subroutine stochastic_ccmc_death(rng, spawn, linked_ccmc, sys, qs, cdet, cluster, proj_energy, logging_info, &
                                    ndeath_tot)

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

        ! Quasinewtwon approaches scale this death step, but doing this naively
        ! would break population control.  Instead, we split H-S into 
        ! (H - E_proj) + (E_proj - S).
        ! The former is scaled and produces the step.  The latter effects the population control.

        ! In:
        !    sys: system being studied.
        !    qs: qmc_state_t containing information about the reference and estimators.
        !    linked_ccmc: if true then only sample linked clusters.
        !    cdet: info on the current excitor (cdet) that we will spawn
        !        from.
        !    cluster: information about the cluster which forms the excitor.
        !    proj_energy: projected energy.  This should be the average value from the last
        !        report loop, not the running total in qs%estimators.
        !    logging_info: logging_t derived type containing information on logging behaviour.
        ! In/Out:
        !    rng: random number generator.
        !    spawn: spawn_t object to which the spanwed particle will be added.

        use const, only: debug
        use ccmc_data, only: cluster_t
        use determinants, only: det_info_t
        use excitations, only: excit_t
        use proc_pointers, only: sc0_ptr, create_spawned_particle_ptr
        use dSFMT_interface, only: dSFMT_t, get_rand_close_open
        use spawn_data, only: spawn_t
        use spawning, only: calc_qn_weighting
        use system, only: sys_t
        use qmc_data, only: qmc_state_t
        use logging, only: logging_t, write_logging_death
        type(sys_t), intent(in) :: sys
        type(qmc_state_t), intent(in) :: qs
        logical, intent(in) :: linked_ccmc
        type(det_info_t), intent(in) :: cdet
        type(cluster_t), intent(in) :: cluster
        real(p), intent(in) :: proj_energy
        type(logging_t), intent(in) :: logging_info
        type(dSFMT_t), intent(inout) :: rng
        type(spawn_t), intent(inout) :: spawn
        integer(int_p), intent(inout) :: ndeath_tot

        real(p) :: pdeath, KiiAi
        integer(int_p) :: nkill
        type(excit_t), parameter :: null_excit = excit_t( 0, [0,0], [0,0], .false.)

        real(p) :: invdiagel

        ! Spawning onto the same excitor so no change in sign due to
        ! a difference in the sign of the determinant formed from applying the
        ! parent excitor to the qs%ref and that formed from applying the
        ! child excitor.
        invdiagel = calc_qn_weighting(qs, cdet%fock_sum)
        if (linked_ccmc) then
            select case (cluster%nexcitors)
            case(0)
                ! Death on the reference has H_ii - E_HF = 0.
                KiiAi = (( - proj_energy)*invdiagel + (proj_energy - qs%shift(1)))*cluster%amplitude
            case(1)
                ! Evaluating the commutator gives
                ! <D1|[H,a1]|D0> = <D1|H|D1> - <D0|H|D0>
                ! (this is scaled for quasinewton approaches)
                KiiAi = (cdet%data(1) * invdiagel + proj_energy - qs%shift(1))*cluster%amplitude
            case(2)
                ! Evaluate the commutator
                ! The cluster operators are a1 and a2 (with a1 D0 = D1, a2 D0 = D2,
                ! a1 a2 D0 = D3) so the commutator gives:
                ! <D3|[[H,a1],a2]|D0> = <D3|H|D3> - <D2|H|D2> - <D1|H|D1> + <D0|H|D0>
                KiiAi = (sc0_ptr(sys, cdet%f) - sc0_ptr(sys, cluster%excitors(1)%f) &
                    - sc0_ptr(sys, cluster%excitors(2)%f) + qs%ref%H00)*cluster%amplitude
                ! (this is scaled for quasinewton approaches)
                KiiAi = KiiAi*invdiagel                                             
            case default
                ! At most two cluster operators can be linked to the diagonal
                ! part of H so this must be an unlinked cluster
                KiiAi = 0.0_p
            end select
        else
            select case (cluster%nexcitors)
            case(0)
                ! Death on the reference has H_ii - E_HF = 0.
                KiiAi = (( - proj_energy)*invdiagel + (proj_energy - qs%shift(1)))*cluster%amplitude
            case(1)
                KiiAi = ((cdet%data(1) - proj_energy)*invdiagel + (proj_energy - qs%shift(1)))*cluster%amplitude
            case default
                KiiAi = ((sc0_ptr(sys, cdet%f) - qs%ref%H00) - proj_energy)*invdiagel *cluster%amplitude
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
            ndeath_tot = ndeath_tot + abs(nkill)
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

        if (debug) call write_logging_death(logging_info, KiiAi, proj_energy, qs%shift(1), invdiagel, &
                                            nkill, pdeath, cluster%amplitude, 0.0_p)

    end subroutine stochastic_ccmc_death

    subroutine stochastic_ccmc_death_nc(rng, linked_ccmc,  sys, qs, isD0, dfock, Hii, proj_energy, population, &
                                        tot_population, ndeath, logging_info)

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
        real(dp), intent(inout) :: tot_population

        real(p) :: pdeath, KiiAi
        integer(int_p) :: nkill, old_pop
        real(p) :: invdiagel

        ! Spawning onto the same excitor so no change in sign due to
        ! a difference in the sign of the determinant formed from applying the
        ! parent excitor to the reference and that formed from applying the
        ! child excitor.

        invdiagel = calc_qn_weighting(qs, dfock)
        if (isD0) then
            KiiAi = ((- proj_energy)*invdiagel + (proj_energy - qs%shift(1)))*population
        else
            if (linked_ccmc) then
                KiiAi = (Hii*invdiagel + proj_energy - qs%shift(1))*population
            else
                KiiAi = ((Hii - proj_energy)*invdiagel + (proj_energy - qs%shift(1)))*population
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

        if (debug) call write_logging_death(logging_info, KiiAi, proj_energy, qs%shift(1), invdiagel, &
                                            nkill, pdeath, real(old_pop,p)/qs%psip_list%pop_real_factor, &
                                            real(population,p)/qs%psip_list%pop_real_factor)

    end subroutine stochastic_ccmc_death_nc

    subroutine linked_spawner_ccmc(rng, sys, qmc_in, qs, spawn_cutoff, cluster, gen_excit_ptr, nspawn, &
                            connection, nspawnings_total, fexcit, cdet, ldet, rdet, left_cluster, right_cluster)

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
        use ccmc_linked, only: calc_pgen, partition_cluster, linked_excitation
        use ccmc_utils, only: collapse_cluster, convert_excitor_to_determinant
        use determinants, only: det_info_t, sum_sp_eigenvalues_bit_string
        use dSFMT_interface, only: dSFMT_t
        use excitations, only: excit_t, create_excited_det, get_excitation_level
        use proc_pointers, only: gen_excit_ptr_t, decoder_ptr
        use spawning, only: attempt_to_spawn, calc_qn_weighting, calc_qn_spawned_weighting
        use system, only: sys_t
        use const, only: depsilon
        use hamiltonian, only: get_hmatel
        use bit_utils, only: count_set_bits
        use qmc_data, only: qmc_in_t, qmc_state_t
        use hamiltonian_data

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
        type(det_info_t), intent(in) :: cdet
        type(cluster_t), intent(inout) :: left_cluster, right_cluster

        ! We incorporate the sign of the amplitude into the Hamiltonian matrix
        ! element, so we 'pretend' to attempt_to_spawn that all excips are
        ! actually spawned by positive excips.
        integer(int_p), parameter :: parent_sign = 1_int_p
        integer :: excitor_sign, excitor_level

        integer :: i, j, npartitions, orb, bit_pos, bit_element
        real(p) :: ppart, pgen, pop
        type(hmatel_t) :: hmatel, delta_h
        logical :: allowed, sign_change, linked, single_unlinked
        integer(i0) :: new_det(sys%basis%string_len)
        integer(i0) :: excitor(sys%basis%string_len)
        real(p) :: invdiagel, fock_sum

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
            call gen_excit_ptr%full(rng, sys, qs%excit_gen_data, rdet, pgen, connection, hmatel, allowed)
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
            hmatel%r = 0.0_p
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
                    pgen = pgen + calc_pgen(sys, qmc_in%excit_gen, qs%excit_gen_data, rdet%f, connection, rdet)

                    ! Sign of the term in the commutator depends on the number of Ts in left_cluster
                    ! also need to account for possible sign change on going from excitor to determinant
                    ! for each of the determinants and when collapsing the clusters
                    ! Need <D_right|right_cluster|D0> and <D_spawn|left_cluster|D>
                    delta_h = get_hmatel(sys, rdet%f, new_det)

                    if (mod(left_cluster%nexcitors,2) /= 0) delta_h%r = -delta_h%r

                    excitor_level = get_excitation_level(qs%ref%f0, rdet%f)
                    call convert_excitor_to_determinant(rdet%f, excitor_level, excitor_sign, qs%ref%f0)
                    if (excitor_sign < 0) delta_h%r = -delta_h%r

                    excitor_level = get_excitation_level(fexcit, new_det)
                    call convert_excitor_to_determinant(fexcit, excitor_level, excitor_sign, new_det)
                    if (excitor_sign < 0) delta_h%r = -delta_h%r

                    if (sign_change) delta_h%r = -delta_h%r

                    hmatel%r = hmatel%r + delta_h%r
                end if
            end do

            ! apply additional factors to pgen
            pgen = pgen*cluster%pselect*nspawnings_total/npartitions

            fock_sum = sum_sp_eigenvalues_bit_string(sys, fexcit)
            invdiagel = calc_qn_weighting(qs, fock_sum - qs%ref%fock_sum)
            ! correct hmatel for cluster amplitude
            hmatel%r = hmatel%r * invdiagel * cluster%amplitude
            excitor_level = get_excitation_level(fexcit, qs%ref%f0)
            call convert_excitor_to_determinant(fexcit, excitor_level, excitor_sign, qs%ref%f0)
            if (excitor_sign < 0) hmatel%r = -hmatel%r

            ! 4) Attempt to spawn
            nspawn = attempt_to_spawn(rng, qs%tau, spawn_cutoff, qs%psip_list%pop_real_factor, hmatel%r, pgen, parent_sign)
        else
            nspawn = 0
        end if

    end subroutine linked_spawner_ccmc

end module ccmc_death_spawning

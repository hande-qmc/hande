module excit_gen_hub_k

! Module for random excitation generators and related routines for the Hubbard
! model in momentum space (i.e. Bloch orbitals).

! See top-level comments in spawning about the overall aim of the spawning step.

use const

implicit none

contains

!--- Split excitation gnerators ---

! Because the off-diagonal matrix elements are so simple in the Hubbard model,
! we don't actually need to generate the entire excitation before we have
! sufficient information to test whether the spawning event is accepted.  This
! allows us to save some time when the spawning attempt is rejected (which is
! the majority of the time).

! These split excitation generators *MUST BE USED IN PAIRS*.

! gen_excit_init_hub_k*:     provides sufficient information to attempt spawning
!                            event.
! gen_excit_finalise_hub_k*: provides rest of excitation required to determine
!                            the sign of the child walkers.

    subroutine gen_excit_init_hub_k(rng, sys, qmc_state, cdet, pgen, connection, abs_hmatel)

        ! Select orbitals from which to excite electrons and hence the
        ! generation probability (which is independent of the virtual orbitals
        ! into which the electrons are excited).

        ! In:
        !    sys: system object being studied.
        !    qmc_state: input options relating to QMC methods.
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.
        ! In/Out:
        !    rng: random number generator.
        ! Out:
        !    pgen: probability of generating the excited determinant from cdet.
        !    connection: from_orb field is set to be the orbitals from which
        !        electrons are excited from the current determinant in the
        !        excitation to the child determinant, on which progeny are spawned.
        !    abs_hmatel: |< D | H | D' >|, the absolute value of the Hamiltonian
        !        matrix element between a determinant and a connected determinant in
        !        the Hubbard model in a Bloch basis.

        use determinants, only: det_info_t
        use excitations, only: excit_t
        use qmc_data, only: qmc_state_t
        use system, only: sys_t
        use dSFMT_interface, only: dSFMT_t

        type(sys_t), intent(in) :: sys
        type(qmc_state_t), intent(in) :: qmc_state
        type(det_info_t), intent(in) :: cdet
        type(dSFMT_t), intent(inout) :: rng
        type(excit_t), intent(out) :: connection
        real(p), intent(out) :: pgen, abs_hmatel

        integer :: ij_sym

        ! Single excitations are not connected determinants within the
        ! momentum space formulation of the Hubbard model.

        ! 1. Select a random pair of spin orbitals to excite from.
        call choose_ij_hub_k(rng, sys, cdet%occ_list_alpha, cdet%occ_list_beta, &
                                 connection%from_orb(1), connection%from_orb(2), ij_sym)

        ! 2. Calculate the generation probability of the excitation.
        ! For one-band systems this depends only upon the orbitals excited from.
        pgen = calc_pgen_hub_k(sys, ij_sym, cdet%f, cdet%unocc_list_alpha)

        ! The hubbard model in momentum space is a special case. Connected
        ! non-identical determinants have the following properties:
        !     a) They differ by two spin-orbitals.
        !     b) In the (i,j)->(a,b) connecting excitation, the spins of i and
        !     j have to be opposite.  This is because
        !     < ij | ab > = U/N_k \delta_{s_i,s_a} \delta_{s_j,s_b}
        !     and so < ij || ab > = 0 if s_i = s_a = s_j = s_b.
        !     In fact:
        !     < ij || ab > = 0          if s_i = s_a = s_j = s_b
        !                    U/N        if s_i = s_a & s_j = s_b & s_i /= s_b
        !                   -U/N        if s_i = s_b & s_j = s_a & s_i /= s_a
        ! The FCIQMC method allows us to only generate connected excitations, so
        ! we can actually test whether we accept the excitation before we finish
        ! completing the excitation.

        ! 3. Matrix element: as we can only excite from (alpha,beta) or (beta,alpha),
        !   H_ij =  < ij | ab >  *or*
        !        = -< ij | ba >
        ! so
        !   |H_ij| = U/\Omega
        ! for allowed excitations.
        abs_hmatel = abs(sys%hubbard%coulomb_k)

    end subroutine gen_excit_init_hub_k

    subroutine gen_excit_finalise_hub_k(rng, sys, qmc_state, cdet, connection, hmatel)

        ! Complete the excitation started in gen_excit_hub_k:
        !    * select the virtual orbitals electrons are excited into.
        !    * find the connecting Hamiltonian matrix element.

        ! In:
        !    sys: system object being studied.
        !    qmc_state: input options relating to QMC methods.
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.
        ! In/Out:
        !    rng: random number generator.
        !    connection: excitation connection between the current determinant
        !        and the child determinant, on which progeny are spawned.
        !        On input the from_orb field must be given.  On output the
        !        remaining fields are also set.
        ! Out:
        !    hmatel: < D | H | D' >, the value of the Hamiltonian matrix element
        !        between a determinant and the connected determinant in the Hubbard
        !        model in a Bloch basis.

        use determinants, only: det_info_t
        use excitations, only: excit_t
        use hamiltonian_hub_k, only: slater_condon2_hub_k_excit
        use momentum_symmetry, only: sym_table
        use qmc_data, only: qmc_state_t
        use system, only: sys_t
        use dSFMT_interface, only: dSFMT_t

        type(sys_t), intent(in) :: sys
        type(qmc_state_t), intent(in) :: qmc_state
        type(det_info_t), intent(in) :: cdet
        type(dSFMT_t), intent(inout) :: rng
        type(excit_t), intent(inout) :: connection
        real(p), intent(out) :: hmatel

        integer :: ij_sym

        ! Continuing from gen_excit_init_hub_k...

        ! Recalculate symmetry of (i,j) pair quickly to save passing it round...
        ij_sym = sym_table((connection%from_orb(1)+1)/2,(connection%from_orb(2)+1)/2)

        ! 4. Well, I suppose we should find out which determinant we're spawning
        ! on...
        connection%nexcit = 2
        call choose_ab_hub_k(rng, sys, cdet%f, cdet%unocc_list_alpha, ij_sym, &
                                 connection%to_orb(1), connection%to_orb(2))

        ! 5. Is connecting matrix element positive (in which case we spawn with
        ! negative walkers) or negative (in which case we spawn with positive
        ! walkers)?
        call slater_condon2_hub_k_excit(sys, cdet%f, connection, hmatel)

    end subroutine gen_excit_finalise_hub_k

    subroutine gen_excit_init_hub_k_no_renorm(rng, sys, qmc_state, cdet, pgen, connection, abs_hmatel)

        ! Create a random excitation from cdet and calculate the probability of
        ! selecting that excitation.

        ! This doesn't exclude the case where, having selected 2 occupied
        ! orbitals (i and j) and the first virtual orbital (a), the final
        ! orbital is already occupied and so that excitation is impossible.
        ! Whilst it is somewhat wasteful (generating excitations which can't be
        ! performed), there is a balance between the cost of generating
        ! forbidden excitations and the O(N) cost of renormalising the
        ! generation probabilities.

        ! In:
        !    sys: system object being studied.
        !    qmc_state: input options relating to QMC methods.
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.
        ! In/Out:
        !    rng: random number generator.
        ! Out:
        !    pgen: probability of generating the excited determinant from cdet.
        !    connection: excitation connection between the current determinant
        !        and the child determinant, on which progeny are spawned.
        !    abs_hmatel: |< D | H | D' >|, the absolute value of the Hamiltonian
        !        matrix element between a determinant and a connected determinant in
        !        the Hubbard model in a Bloch basis.

        use determinants, only: det_info_t
        use dSFMT_interface, only: dSFMT_t, get_rand_close_open
        use excitations, only: excit_t
        use hamiltonian_hub_k, only: slater_condon2_hub_k_excit
        use qmc_data, only: qmc_state_t
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        type(qmc_state_t), intent(in) :: qmc_state
        type(det_info_t), intent(in) :: cdet
        type(dSFMT_t), intent(inout) :: rng
        type(excit_t), intent(out) :: connection
        real(p), intent(out) :: pgen, abs_hmatel

        integer :: ij_sym
        logical :: allowed_excitation

        ! Single excitations are not connected determinants within the
        ! momentum space formulation of the Hubbard model.

        ! 1. Select a random pair of spin orbitals to excite from.
        call choose_ij_hub_k(rng, sys, cdet%occ_list_alpha, cdet%occ_list_beta, connection%from_orb(1), &
                             connection%from_orb(2), ij_sym)

        ! 2. Chose a random pair of spin orbitals to excite to.
        ! WARNING: if the implementation of generating the excitation is
        ! changed, the generation probability will also (probably) need to
        ! be altered.
        call find_ab_hub_k(rng, sys, cdet%f, cdet%unocc_list_alpha, ij_sym, connection%to_orb(1), &
                           connection%to_orb(2), allowed_excitation)
        connection%nexcit = 2

        if (allowed_excitation) then

            ! 3. probability of this excitation.
            !
            ! Usually:
            !   pgen = p(i,j) [ p(a|i,j) p(b|i,j,a) + p(b|i,j) p(a|i,j,b) ]
            ! We must be very, very careful here.
            ! Consider the excitation (i,j) -> (a,b).
            ! As the excitation must be alpha,beta->alpha,beta, we chose to pick
            ! a first *and* require a to be alpha spin.
            !   ***WARNING: this is purely an implementation detail***
            ! As such, we cannot generate an excitation where the alpha orbital
            ! in the virtal pair (a,b) is already occupied (even though the beta
            ! orbital in the (a,b) pair can be).  This asymmetry results in the
            ! probability being simply:
            !   pgen = p(i,j) p(a|i,j) p(b|i,j,a)
            ! Here's an example to show this asymmetry.  Consider a system with
            ! N_a alpha electrons, N_b beta electrons and M sites (so M-N_a
            ! (M-N_b) virtual alpha (beta) orbitals). There are N_a*N_b ways of
            ! choosing (i,j).  However, because we choose alpha first in the
            ! (a,b) pair and because the Hubbard model is a one-band system,
            ! there are N_a ways of choosing (a,b) (as we don't require b to be
            ! unoccupied).  Another way of saying this is that our
            ! implementation choices make it to generate an excitation which
            ! attempts to excite into an occupied alpha orbital.
            !
            !   p(i,j) = 1/(nalpha*nbeta)
            !   p(b|i,j,a) = 1
            !   p(a|i,j) = n_virt_alpha
            pgen = 1.0_p/(sys%nalpha*sys%nbeta*sys%nvirt_alpha)

            ! 4. |H_ij| is constant for this system.
            abs_hmatel = abs(sys%hubbard%coulomb_k)

        else

            ! Forbidden excitation.
            pgen = 1.0_p ! just to avoid division by zero issues
            abs_hmatel = 0.0_p

        end if

    end subroutine gen_excit_init_hub_k_no_renorm

    subroutine gen_excit_finalise_hub_k_no_renorm(rng, sys, qmc_state, cdet, connection, hmatel)

        ! Complete the excitation started in gen_excit_hub_k_no_renorm:
        !    * find the connecting Hamiltonian matrix element.

        ! In:
        !    sys: system object being studied.
        !    qmc_state: input options relating to QMC methods.
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.
        !    connection: excitation connection between the current determinant
        !        and the child determinant, on which progeny are spawned.
        ! In/Out:
        !    rng: random number generator (interface compatibility only).
        ! Out:
        !    hmatel: < D | H | D' >, the value of the Hamiltonian matrix element
        !        between a determinant and the connected determinant in the Hubbard
        !        model in a Bloch basis.

        use determinants, only: det_info_t
        use excitations, only: excit_t
        use hamiltonian_hub_k, only: slater_condon2_hub_k_excit
        use qmc_data, only: qmc_state_t
        use system, only: sys_t
        use dSFMT_interface, only: dSFMT_t

        type(sys_t), intent(in) :: sys
        type(qmc_state_t), intent(in) :: qmc_state
        type(det_info_t), intent(in) :: cdet
        type(dSFMT_t), intent(inout) :: rng
        type(excit_t), intent(inout) :: connection ! inout for interface compatibility
        real(p), intent(out) :: hmatel

        ! Continuing on from gen_excit_init_hub_k_no_renorm now the spawning
        ! event has been accepted.

        ! For no_renorm the excitation has been completely specified, so we just
        ! need to find the exact connecting matrix element.
        call slater_condon2_hub_k_excit(sys, cdet%f, connection, hmatel)

    end subroutine gen_excit_finalise_hub_k_no_renorm

!--- Excitation generation (see also split excitation generators) ---

    subroutine gen_excit_hub_k(rng, sys, qmc_state, cdet, pgen, connection, hmatel)

        ! Create a random excitation from cdet and calculate the probability of
        ! selecting that excitation.

        ! In:
        !    sys: system object being studied.
        !    qmc_state: input options relating to QMC methods.
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.
        ! Out:
        ! In/Out:
        !    rng: random number generator.
        !    pgen: probability of generating the excited determinant from cdet.
        !    connection: excitation connection between the current determinant
        !        and the child determinant, on which progeny are spawned.
        !    hmatel: < D | H | D' >, the Hamiltonian matrix element between a
        !        determinant and a connected determinant in the Hubbard model in
        !        a Bloch basis.

        use determinants, only: det_info_t
        use excitations, only: excit_t
        use qmc_data, only: qmc_state_t
        use system, only: sys_t
        use hamiltonian_hub_k, only: slater_condon2_hub_k_excit
        use dSFMT_interface, only: dSFMT_t

        type(sys_t), intent(in) :: sys
        type(qmc_state_t), intent(in) :: qmc_state
        type(det_info_t), intent(in) :: cdet
        type(dSFMT_t), intent(inout) :: rng
        real(p), intent(out) :: pgen, hmatel
        type(excit_t), intent(out) :: connection

        integer :: ij_sym

        ! See notes in gen_excit_init_hub_k and gen_excit_finalise_hub_k for
        ! more details regarding this algorithm.
        ! The split excitation generators are a somewhat more optimised version
        ! for fciqmc and i-fciqmc, as you can decide whether to spawn (or not)
        ! without actually generating the excitation when working in momentum
        ! space.

        connection%nexcit = 2

        ! 1. Select a random pair of spin orbitals to excite from.
        call choose_ij_hub_k(rng, sys, cdet%occ_list_alpha, cdet%occ_list_beta, &
                             connection%from_orb(1), connection%from_orb(2), ij_sym)

        ! 2. Select a random pair of spin orbitals to excite to.
        call choose_ab_hub_k(rng, sys, cdet%f, cdet%unocc_list_alpha, ij_sym, &
                             connection%to_orb(1), connection%to_orb(2))

        ! 3. Calculate the generation probability of the excitation.
        ! For one-band systems this depends only upon the orbitals excited from.
        pgen = calc_pgen_hub_k(sys, ij_sym, cdet%f, cdet%unocc_list_alpha)

        ! 4. find the connecting matrix element.
        call slater_condon2_hub_k_excit(sys, cdet%f, connection, hmatel)

    end subroutine gen_excit_hub_k

    subroutine gen_excit_hub_k_no_renorm(rng, sys, qmc_state, cdet, pgen, connection, hmatel)

        ! Create a random excitation from cdet and calculate the probability of
        ! selecting that excitation.

        ! This doesn't exclude the case where, having selected 2 occupied
        ! orbitals (i and j) and the first virtual orbital (a), the final
        ! orbital is already occupied and so that excitation is impossible.
        ! Whilst it is somewhat wasteful (generating excitations which can't be
        ! performed), there is a balance between the cost of generating
        ! forbidden excitations and the O(N) cost of renormalising the
        ! generation probabilities.

        ! In:
        !    sys: system object being studied.
        !    qmc_state: input options relating to QMC methods.
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.
        ! In/Out:
        !    rng: random number generator.
        ! Out:
        !    pgen: probability of generating the excited determinant from cdet.
        !    connection: excitation connection between the current determinant
        !        and the child determinant, on which progeny are spawned.
        !    hmatel: < D | H | D' >, the Hamiltonian matrix element between a
        !        determinant and a connected determinant in the Hubbard model in
        !        a Bloch basis.

        use determinants, only: det_info_t
        use excitations, only: excit_t
        use qmc_data, only: qmc_state_t
        use system, only: sys_t
        use hamiltonian_hub_k, only: slater_condon2_hub_k_excit
        use dSFMT_interface, only: dSFMT_t

        type(sys_t), intent(in) :: sys
        type(qmc_state_t), intent(in) :: qmc_state
        type(det_info_t), intent(in) :: cdet
        type(dSFMT_t), intent(inout) :: rng
        real(p), intent(out) :: pgen, hmatel
        type(excit_t), intent(out) :: connection

        integer :: ij_sym
        logical :: allowed_excitation

        ! See notes in gen_excit_init_hub_k and gen_excit_hub_k for more details regarding this algorithm.

        ! 1. Select a random pair of spin orbitals to excite from.
        call choose_ij_hub_k(rng, sys, cdet%occ_list_alpha, cdet%occ_list_beta, connection%from_orb(1), &
                             connection%from_orb(2), ij_sym)

        ! 2. Chose a random pair of spin orbitals to excite to.
        call find_ab_hub_k(rng, sys, cdet%f, cdet%unocc_list_alpha, ij_sym, connection%to_orb(1), &
                           connection%to_orb(2), allowed_excitation)

        if (allowed_excitation) then

            ! 3. probability of this excitation.
            !    ***WARNING: dependent upon (arbitrary) implementation choices***
            !    See comments in gen_excit_init_hub_k_no_renorm.
            ! pgen = p(i,j) p(a|i,j) p(b|i,j,a)
            !      = 1/(nalpha*nbeta) 1/nvirt_alpha
            pgen = 1.0_p/(sys%nalpha*sys%nbeta*sys%nvirt_alpha)

            ! 4. find the connecting matrix element.
            connection%nexcit = 2
            call slater_condon2_hub_k_excit(sys, cdet%f, connection, hmatel)

        else

            ! Generated a forbidden excitation (b is already occupied)
            hmatel = 0.0_p
            pgen = 1.0_p

        end if

    end subroutine gen_excit_hub_k_no_renorm

!--- Select random orbitals involved in a valid double excitation ---

    subroutine choose_ij_hub_k(rng, sys, occ_list_alpha, occ_list_beta, i ,j, ij_sym)

        ! Randomly choose a pair of spin-orbitals.
        !
        ! This is specific to the Hubbard model in momentum space.
        ! Only double excitations which excite from an alpha and a beta
        ! spin orbital are connected, so we return only i,j which meet this
        ! criterion.
        !
        ! In:
        !    sys: system object being studied.
        !    occ_list_alpha: Integer list of occupied alpha spin-orbitals
        !        (min length: sys%nalpha.)
        !    occ_list_beta: Integer list of occupied beta spin-orbitals.
        !        (min length: sys%nbeta.)
        ! In/Out:
        !    rng: random number generator.
        ! Out:
        !    i, j: randomly selected spin-orbitals.
        !    ij_sym: symmetry label of the (i,j) combination.

        use momentum_symmetry, only: sym_table
        use system, only: sys_t
        use dSFMT_interface, only: dSFMT_t, get_rand_close_open

        type(sys_t), intent(in) :: sys
        integer, intent(in) :: occ_list_alpha(:), occ_list_beta(:)
        type(dSFMT_t), intent(inout) :: rng
        integer, intent(out) :: i,j, ij_sym
        integer :: ind
        real(dp) :: r

        ! We use a similar indexing scheme to choose_ij, except our initial
        ! indices refer to an index in the occupied alpha array and in the
        ! occupied beta array.  This means that rather than be a triangular
        ! indexing scheme, we have a rectangular one:
        !   1  2  3    to    1,1  2,1  1,1
        !   3  5  6          1,2  2,2  3,2
        !   7  8  9          1,3  2,3  3,3
        !  10 11 12          1,4  2,4  3,4
        ! total number of possible combinations is nalpha*nbeta.
        ! The indexing scheme is:
        !  p = (i-1)*n_j + j
        ! Hence to invert this, following a similar method to Rifkin:
        !  i = floor( (p-1)/n_j ) + 1
        !  j = p - (i-1)*n_j

        r = get_rand_close_open(rng)

        ! i,j initially refer to the indices in the lists of occupied spin-orbitals
        ! rather than the spin-orbitals.

        ind = int(r*sys%nalpha*sys%nbeta) + 1
        i = int( (ind-1.0_p)/sys%nbeta ) + 1
        j = ind - (i-1)*sys%nbeta

        ! i,j are the electrons we're exciting.  Find the occupied corresponding
        ! spin-orbitals.
        i = occ_list_alpha(i)
        j = occ_list_beta(j)

        ! Symmetry info is a simple lookup...
        ij_sym = sym_table((i+1)/2,(j+1)/2)

    end subroutine choose_ij_hub_k

    subroutine choose_ab_hub_k(rng, sys, f, unocc_list_alpha, ij_sym, a, b)

        ! Choose a random pair of (a,b) unoccupied virtual spin-orbitals into
        ! which electrons are excited.
        ! (a,b) are chosen such that the (i,j)->(a,b) excitation is symmetry-
        ! allowed.
        ! In:
        !    sys: system object being studied.
        !    f(string_len): bit string representation of the Slater
        !        determinant.
        !    unocc_alpha: integer list of the unoccupied alpha spin-orbitals.
        !        (min length: sys%nvirt_alpha.)
        !    ij_sym: symmetry spanned by the (i,j) combination of unoccupied
        !        spin-orbitals into which electrons are excited.
        ! In/Out:
        !    rng: random number generator.
        ! Out:
        !    a,b: virtual spin orbitals involved in the excitation.

        use system, only: sys_t
        use dSFMT_interface, only:  dSFMT_t

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f(sys%basis%string_len)
        integer, intent(in) :: unocc_list_alpha(:)
        integer, intent(in) :: ij_sym
        type(dSFMT_t), intent(inout) :: rng
        integer, intent(out) :: a, b

        logical :: allowed_excitation

        allowed_excitation = .false.

        do while (.not.allowed_excitation)

            ! Until we find an allowed excitation.

            call find_ab_hub_k(rng, sys, f, unocc_list_alpha, ij_sym, a, b, allowed_excitation)

        end do

    end subroutine choose_ab_hub_k

!--- Select random orbitals in double excitations ---

    subroutine find_ab_hub_k(rng, sys, f, unocc_list_alpha, ij_sym, a, b, allowed_excitation)

        ! Choose a random pair of (a,b) spin-orbitals.
        ! (a,b) are chosen such that the (i,j)->(a,b) excitation is symmetry-
        ! allowed and a is a virtual spin-orbital.  As (a,b) must be one alpha
        ! orbital and one beta orbital, we can choose a to be an alpha orbital
        ! without loss of generality.
        ! In:
        !    sys: system object being studied.
        !    f(string_len): bit string representation of the Slater
        !        determinant.
        !    unocc_alpha: integer list of the unoccupied alpha spin-orbitals.
        !        (min length: sys%nvirt_alpha.)
        !    ij_sym: symmetry spanned by the (i,j) combination of unoccupied
        !        spin-orbitals into which electrons are excited.
        ! In/Out:
        !    rng: random number generator.
        ! Out:
        !    a,b: spin orbitals involved in the excitation.
        !    allowed_excitation: is true if the excitation (i,j)->(a,b) is
        !        allowed (i.e. both a and b are indeed spin orbitals).

        use dSFMT_interface, only:  dSFMT_t, get_rand_close_open
        use system, only: sys_t
        use momentum_symmetry, only: sym_table, inv_sym

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f(sys%basis%string_len)
        integer, intent(in) :: unocc_list_alpha(:)
        integer, intent(in) :: ij_sym
        type(dSFMT_t), intent(inout) :: rng
        integer, intent(out) :: a, b
        logical, intent(out) :: allowed_excitation

        integer :: r, b_pos, b_el, ka

        ! The excitation i,j -> a,b is only allowed if k_i + k_j - k_a - k_b = 0
        ! (up to a reciprocal lattice vector).  We store k_i + k_j as ij_sym, so
        ! k_a + k_b must be identical to this.
        ! If we view this in terms of the representations spanned by i,j,a,b
        ! under translational symmetry (which forms an Abelian group) then
        !   \Gamma_i* x \Gamma_j* x \Gamma_a x \Gamma_b = \Gamma_1
        ! is equivalent to the conversation of crystal momentum (where \Gamma_1
        ! is the totally symmetric representation).
        ! Further, as
        !   \Gamma_i* x \Gamma_i = 1
        ! and direct products in Abelian groups commute, it follows that:
        !   \Gamma_b = \Gamma_i x \Gamma_j x \Gamma_a*
        ! Thus k_b is defined by i,j and a.  As we only consider one-band
        ! systems, b is thus defined by spin conservation.

        ! One electron must be in unocc_list_alpha, so we can use the
        ! random number just to find which unoccupied alpha orbital is in
        ! the excitation.

        r = int(get_rand_close_open(rng)*(sys%nvirt_alpha)) + 1

        a = unocc_list_alpha(r)
        ! Find corresponding beta orbital which satisfies conservation
        ! of crystal momentum.
        ka = (a+1)/2
        b = 2*sym_table(ij_sym, inv_sym(ka))

        b_pos = sys%basis%bit_lookup(1,b)
        b_el = sys%basis%bit_lookup(2,b)

        ! If b is unoccupied then have found an allowed excitation.
        allowed_excitation = .not.btest(f(b_el), b_pos)

    end subroutine find_ab_hub_k

!--- Excitation generation probabilities ---

    pure function calc_pgen_hub_k(sys, ab_sym, f, unocc_alpha) result(pgen)

        ! Calculate the generation probability of a given excitation for the
        ! Hubbard model in momentum space.  The Hubbard model is a special case
        ! as it is a one-band system and so the generation probability is
        ! independent of the virtual spin-orbitals into which electrons are
        ! excited and depends only upon the spin-orbitals from which we excite.
        !
        ! Note that all the information required for input should be available
        ! during the FCIQMC algorithm and should not be needed to be calculated.
        !
        ! Further, we assume only allowed excitations are generated.
        !
        ! In:
        !    sys: system object being studied.
        !    ab_sym: symmetry spanned by the (a,b) combination of unoccupied
        !        spin-orbitals into which electrons are excited.
        !    f: bit string representation of the determinant we're exciting
        !        from.
        !    unocc_alpha: integer list of the unoccupied alpha spin-orbitals.
        !        (min length: sys%nvirt_alpha.)
        ! Returns:
        !    pgen: the generation probability of the excitation.  See notes in
        !        spawning.

        use system, only: sys_t
        use momentum_symmetry, only: sym_table, inv_sym

        real(p) :: pgen
        type(sys_t), intent(in) :: sys
        integer, intent(in) :: ab_sym
        integer(i0), intent(in) :: f(sys%basis%string_len)
        integer, intent(in) :: unocc_alpha(:)

        integer :: forbidden_excitations, a, b, a_pos, a_el, b_pos, b_el, ka, kb

        forbidden_excitations = 0

        ! pgen = p(i,j) [ p(a|i,j) p(b|i,j,a) + p(b|i,j) p(a|i,j,b) ]
        !
        ! The number of ways of choosing i,j is
        !
        !  nalpha*nbeta
        !
        ! Due to the requirement that crystal momentum is conserved and that the Hubbard
        ! model is a 2-band system:
        !
        !  p(a|i,j,b) = 1
        !  p(b|i,j,a) = 1
        !
        ! i.e. once three spin-orbitals are selected, the fourth is fixed.
        !
        ! We now consider p(a|i,j).  Not all a are possible, as an a virtual spin-orbital
        ! can have an occupied b spin-orbital, as b is fixed by the choice of i,j and a.
        !
        ! The number of spin-orbitals from which a can be chosen is
        !
        !  nbasis - nel - delta
        !
        ! where delta is the number of a orbitals which are forbidden due to b being occupied.
        ! -> p(a|i,j) = 1/(nbasis - nel - delta).
        ! p(b|i,j) is clearly identical to p(a|i,j).
        ! -> p(a|i,j) p(b|i,j,a) + p(b|i,j) p(a|i,j,b) = 2/(nbasis - nel - delta).
        !
        ! However, we can be a bit faster.
        !  nel = nalpha + nbeta
        !  nvirt_alpha = nbasis/2 - nalpha
        !  delta = delta_alpha + delta_beta
        ! where delta_alpha is the number of alpha orbitals which cannot be
        ! excited into because the correspoinding beta orbital is occupied and
        ! similarly for delta_beta.  Note that delta_alpha /= delta_beta in
        ! general.
        !
        ! However, as the Hubbard model is one-band model, for every possible
        ! alpha virtual orbital, there is exactly one possible beta virtual
        ! orbital.  Therefore:
        !  nvirt_alpha - delta_alpha = nvirt_beta - delta_beta
        ! even when delta_alpha /= delta_beta.
        ! This means we need only calculate delta_alpha or delta_beta.
        !
        ! Hence
        !  p(a|i,j) p(b|i,j,a) + p(b|i,j) p(a|i,j,b) = 1/(nvirt_alpha - delta_alpha)
        !
        !                           1
        ! pgen =  --------------------------------------
        !         nalpha*nbeta*(nvirt_alpha-delta_alpha)

        ! Note that we could also use the implementation choice of always
        ! choosing alpha first in the (a,b) pair to obtain the same result (as
        ! must be done for the no_renorm pgen).

        ! We count the number of a orbitals which cannot be excited into due to
        ! the corresponding b orbital not being available.
        ! The Hubbard model is a 2-band system, which makes this pleasingly
        ! easy. :-)

        ! exciting from alpha, beta orbitals.
        ! alpha orbitals are odd (1,3,5,...)
        ! beta orbitals are odd (2,4,6,...)
        ! [Note that this is not the indexing used in bit strings: see basis
        ! module.]
        ! To convert from the wavevector label, 1,2,3,..., where wavevector
        ! 1 corresponds to orbitals 1 and 2, we do:
        !   2*k-1    for alpha
        !   2*k      for beta
        ! and similarly for the reverse transformation.

        ! a is an alpha orbital
        ! b is a beta orbital.
        do a = 1, sys%nvirt_alpha
            ka = (unocc_alpha(a)+1)/2
            b = 2*sym_table(ab_sym, inv_sym(ka))
            b_pos = sys%basis%bit_lookup(1,b)
            b_el = sys%basis%bit_lookup(2,b)
            ! Are (a,b) both unoccupied?
            if (btest(f(b_el), b_pos)) forbidden_excitations = forbidden_excitations + 1
        end do

        pgen = 1.0_p/(sys%nalpha*sys%nbeta*(sys%nvirt_alpha - forbidden_excitations))

    end function calc_pgen_hub_k

end module excit_gen_hub_k

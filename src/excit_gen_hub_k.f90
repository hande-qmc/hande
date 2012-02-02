module excit_gen_hub_k

! Module for random excitation generators and related routines for the Hubbard
! model in momentum space (i.e. Bloch orbitals).

! See top-level comments in spawning about the overall aim of the spawning step.

use const

implicit none

contains

!--- Top level spawning routines ---

    subroutine spawn_hub_k(cdet, parent_sign, nspawn, connection)

        ! Attempt to spawn a new particle on a connected determinant for the
        ! momentum space formulation of the Hubbard model.
        !
        ! In:
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.
        !    parent_sign: sign of the population on the parent determinant (i.e.
        !        either a positive or negative integer).
        ! Out:
        !    nspawn: number of particles spawned.  0 indicates the spawning
        !        attempt was unsuccessful.
        !    connection: excitation connection between the current determinant
        !        and the child determinant, on which progeny are spawned.

        use determinants, only: det_info
        use dSFMT_interface, only:  genrand_real2
        use excitations, only: excit
        use hamiltonian, only: slater_condon2_hub_k_excit
        use fciqmc_data, only: tau
        use system, only: hub_k_coulomb

        type(det_info), intent(in) :: cdet
        integer, intent(in) :: parent_sign
        integer, intent(out) :: nspawn
        type(excit), intent(out) :: connection

        real(p) :: pgen, psuccess, pspawn, hmatel
        integer :: i, j, a, b, ij_sym

        ! Single excitations are not connected determinants within the
        ! momentum space formulation of the Hubbard model.

        ! 1. Select a random pair of spin orbitals to excite from.
        call choose_ij_hub_k(cdet%occ_list_alpha, cdet%occ_list_beta, i ,j, ij_sym)

        ! 2. Calculate the generation probability of the excitation.
        ! For two-band systems this depends only upon the orbitals excited from.
        pgen = calc_pgen_hub_k(ij_sym, cdet%f, cdet%unocc_list_alpha, cdet%unocc_list_beta)

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

        ! 3. Test that whether the attempted spawning is successful.
        ! As we can only excite from (alpha,beta) or (beta,alpha),
        !   H_ij =  < ij | ab >  *or*
        !        = -< ij | ba >
        ! so
        !   |H_ij| = U/\Omega
        ! for allowed excitations.
        pspawn = tau*abs(hub_k_coulomb)/pgen
        psuccess = genrand_real2()

        ! Need to take into account the possibilty of a spawning attempt
        ! producing multiple offspring...
        ! If pspawn is > 1, then we spawn floor(pspawn) as a minimum and
        ! then spawn a particle with probability pspawn-floor(pspawn).
        nspawn = int(pspawn)
        pspawn = pspawn - nspawn

        if (pspawn > psuccess) nspawn = nspawn + 1

        if (nspawn > 0) then
            ! 4. Well, I suppose we should find out which determinant we're spawning
            ! on...
            call choose_ab_hub_k(cdet%f, cdet%unocc_list_alpha, ij_sym, a, b)

            connection%nexcit = 2
            connection%from_orb(1:2) = (/ i,j /)
            connection%to_orb(1:2) = (/ a,b /)

            ! 5. Is connecting matrix element positive (in which case we spawn with
            ! negative walkers) or negative (in which case we spawn with positive
            ! walkers)?
            call slater_condon2_hub_k_excit(cdet%f, connection, hmatel)

            ! If H_ij is positive, then the spawned walker is of opposite sign
            ! to the parent.
            ! If H_ij is negative, then the spawned walker is of the same sign
            ! as the parent.
            if (hmatel > 0) then
                nspawn = -sign(nspawn, parent_sign)
            else
                nspawn = sign(nspawn, parent_sign)
            end if

        end if

    end subroutine spawn_hub_k

    subroutine spawn_hub_k_no_renorm(cdet, parent_sign, nspawn, connection)

        ! Attempt to spawn a new particle on a connected determinant for the
        ! momentum space formulation of the Hubbard model.
        !
        ! This doesn't use excitation generators which exclude the case where,
        ! having selected 2 occupied orbitals (i and j) and the first virtual
        ! orbital (a), the final orbital is already occupied and so that
        ! excitation is impossible.  Whilst it is somewhat wasteful (generating
        ! excitations which can't be performed), there is a balance between the
        ! cost of generating forbidden excitations and the O(N) cost of
        ! renormalising the generation probabilities.
        !
        ! In:
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.
        !    parent_sign: sign of the population on the parent determinant (i.e.
        !        either a positive or negative integer).
        ! Out:
        !    nspawn: number of particles spawned.  0 indicates the spawning
        !        attempt was unsuccessful.
        !    connection: excitation connection between the current determinant
        !        and the child determinant, on which progeny are spawned.

        use determinants, only: det_info
        use dSFMT_interface, only:  genrand_real2
        use excitations, only: excit
        use hamiltonian, only: slater_condon2_hub_k_excit
        use fciqmc_data, only: tau
        use system, only: hub_k_coulomb, nalpha, nbeta, nvirt

        type(det_info), intent(in) :: cdet
        integer, intent(in) :: parent_sign

        integer, intent(out) :: nspawn
        type(excit), intent(out) :: connection

        real(p) :: pgen, psuccess, pspawn, hmatel
        integer :: ij_sym
        logical :: allowed_excitation

        ! Single excitations are not connected determinants within the
        ! momentum space formulation of the Hubbard model.

        ! 1. Select a random pair of spin orbitals to excite from.
        call choose_ij_hub_k(cdet%occ_list_alpha, cdet%occ_list_beta, connection%to_orb(1), connection%to_orb(2), ij_sym)

        ! 2. Chose a random pair of spin orbitals to excite to.
        call find_ab_hub_k(cdet%f, cdet%unocc_list_alpha, ij_sym, connection%to_orb(1), connection%to_orb(2), allowed_excitation)

        if (allowed_excitation) then

            ! 3. probability of this excitation.
            ! pgen = p(i,j) [ p(a|i,j) p(b|i,j,a) + p(b|i,j) p(a|i,j,b) ]
            ! pgen = 1/(nalpha*nbeta) [ 1/(nbasis-nel) + 1/(basis-nel) ]
            !                    2
            !      =  -------------------------
            !         nalpha*nbeta*(nbasis-nel)
            pgen = 2.0_dp/(nalpha*nbeta*nvirt)

            ! 4. Test that whether the attempted spawning is successful.
            pspawn = tau*abs(hub_k_coulomb)/pgen
            psuccess = genrand_real2()
            nspawn = int(pspawn)
            pspawn = pspawn - nspawn
            if (pspawn > psuccess) nspawn = nspawn + 1

            if (nspawn > 0) then

                connection%nexcit = 2

                ! 5. Is connecting matrix element positive (in which case we spawn with
                ! negative walkers) or negative (in which case we spawn with positive
                ! walkers)?
                call slater_condon2_hub_k_excit(cdet%f, connection, hmatel)

                ! If H_ij is positive, then the spawned walker is of opposite sign
                ! to the parent.
                ! If H_ij is negative, then the spawned walker is of the same sign
                ! as the parent.
                if (hmatel > 0) then
                    nspawn = -sign(nspawn, parent_sign)
                else
                    nspawn = sign(nspawn, parent_sign)
                end if

            end if

        else

            ! Generated a forbidden excitation (b is already occupied)
            nspawn = 0

        end if

    end subroutine spawn_hub_k_no_renorm

!--- Excitation generation ---

    subroutine gen_excit_hub_k(cdet, pgen, connection, hmatel)

        ! Create a random excitation from cdet and calculate the probability of
        ! selecting that excitation.

        ! In:
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.
        ! Out:
        !    pgen: probability of generating the excited determinant from cdet.
        !    connection: excitation connection between the current determinant
        !        and the child determinant, on which progeny are spawned.
        !    hmatel: < D | H | D' >, the Hamiltonian matrix element between a
        !        determinant and a connected determinant in the Hubbard model in
        !        a Bloch basis.

        use determinants, only: det_info
        use excitations, only: excit
        use fciqmc_data, only: tau
        use system, only: hub_k_coulomb
        use hamiltonian, only: slater_condon2_hub_k_excit

        type(det_info), intent(in) :: cdet
        real(p), intent(out) :: pgen, hmatel
        type(excit), intent(out) :: connection

        integer :: ij_sym

        ! See notes in spawn_hub_k for more details regarding this algorithm.
        ! spawn_hub_k contains a somewhat more optimised version for fciqmc and
        ! i-fciqmc, as you can decide whether to spawn (or not) without actually
        ! generating the excitation when working in momentum space.

        ! However, we need pgen and an excitation for use with, e.g. the folded
        ! spectrum algorithm.

        connection%nexcit = 2

        ! 1. Select a random pair of spin orbitals to excite from.
        call choose_ij_hub_k(cdet%occ_list_alpha, cdet%occ_list_beta, connection%from_orb(1), connection%from_orb(2), ij_sym)

        ! 2. Select a random pair of spin orbitals to excite to.
        call choose_ab_hub_k(cdet%f, cdet%unocc_list_alpha, ij_sym, connection%to_orb(1), connection%to_orb(2))

        ! 3. Calculate the generation probability of the excitation.
        ! For two-band systems this depends only upon the orbitals excited from.
        pgen = calc_pgen_hub_k(ij_sym, cdet%f, cdet%unocc_list_alpha, cdet%unocc_list_beta)


        ! 4. find the connecting matrix element.
        call slater_condon2_hub_k_excit(cdet%f, connection, hmatel)

    end subroutine gen_excit_hub_k

!--- Select random orbitals involved in a valid double excitation ---

    subroutine choose_ij_k(occ_list, i ,j, ij_sym, ij_spin)

        ! Randomly choose a pair of spin-orbitals for 1-band systems with
        ! Bloch orbitals.
        !
        ! See choose_ij_hub_k for a specific procedure for the momentum space
        ! formulation of the hubbard model.
        !
        ! In:
        !    occ_list: integer list of occupied spin-orbitals in a determinant.
        ! Out:
        !    i, j: randomly selected spin-orbitals.
        !    ij_sym: symmetry label of the (i,j) combination.
        !    ij_spin: -1 if (i,j) are both beta, +1 if (i,j) are both alpha and
        !        0 if it's a "mixed" excitation, i.e. alpha, beta or beta,
        !        alpha.

        use basis, only: basis_fns
        use momentum_symmetry, only: sym_table
        use system, only: nel
        use dSFMT_interface, only: genrand_real2

        integer, intent(in) :: occ_list(nel)
        integer, intent(out) :: i,j, ij_sym, ij_spin
        integer :: ind, spin_sum
        real(dp) :: r

        ! We use a triangular indexing scheme to compress 2 electron indices
        ! into 1.
        ! For i/=j and (for an arbitrary choice of i>j), a 1D index of
        ! a strictly lower triangular array is:
        !   p = (i-1)(i-2)/2 + j,
        ! where 1<=j<i and 1<=p<=n(n-1)/2
        ! This maps the indexing scheme as:
        !    .                  .
        !   2,1  .              1  .
        !   3,1 3,2  .      to  2  3  .
        !   4,1 4,2 4,3  .      3  4  5  .
        ! We want to do the reverse process in order to pick 2 electron labels
        ! from one random number.
        ! Consider the case where j=1.  i can trivially be determined from the
        ! quadratic equation:
        !   i = 3/2 + \sqrt(2p-1.75)
        ! As j<i and, for a fixed i, p increases monotonically with j, the
        ! integer part of i given by the above equation can never exceed the
        ! correct value.  Hence i can be found for arbitrary j by taking the
        ! floor of the preceeding equation.  j follows trivially.
        !
        ! See (for lower triangular arrays rather than strictly lower):
        ! Decoding the sequential indices of a (lower) triangular array
        ! SIS-75-1783,  S Rifkin, CERN report (CERN-DD-75-7).

        ! This might seem odd, but it enables us to pick the (i,j) pair to
        ! excite with half the calls to the random number generator, which
        ! represents a substantial saving. :-)

        r = genrand_real2()
        ind = int(r*nel*(nel-1)/2) + 1

        ! i,j initially refer to the indices in the lists of occupied spin-orbitals
        ! rather than the spin-orbitals.
        i = int(1.50_p + sqrt(2*ind-1.750_p))
        j = ind - ((i-1)*(i-2))/2

        ! i,j are the electrons we're exciting.  Find the occupied corresponding
        ! spin-orbitals.
        i = occ_list(i)
        j = occ_list(j)

        ! Symmetry info is a simple lookup...
        ij_sym = sym_table((i+1)/2,(j+1)/2)

        ! Is mod faster than lookup?  Not sure...
        spin_sum = basis_fns(i)%Ms + basis_fns(j)%Ms
        select case(spin_sum)
        case(2)
            ! alpha, alpha
            ij_spin = 1
        case(0)
            ! alpha, beta
            ij_spin = 0
        case(-2)
            ! beta, beta
            ij_spin = -1
        end select

    end subroutine choose_ij_k

    subroutine choose_ij_hub_k(occ_list_alpha, occ_list_beta, i ,j, ij_sym)

        ! Randomly choose a pair of spin-orbitals.
        !
        ! This is specific to the Hubbard model in momentum space.
        ! Only double excitations which excite from an alpha and a beta
        ! spin orbital are connected, so we return only i,j which meet this
        ! criterion.
        !
        ! In:
        !    occ_list_alpha: Integer list of occupied alpha spin-orbitals.
        !    occ_list_beta: Integer list of occupied beta spin-orbitals.
        ! Out:
        !    i, j: randomly selected spin-orbitals.
        !    ij_sym: symmetry label of the (i,j) combination.

        use momentum_symmetry, only: sym_table
        use system, only: nalpha, nbeta
        use dSFMT_interface, only: genrand_real2

        integer, intent(in) :: occ_list_alpha(nalpha), occ_list_beta(nbeta)
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

        r = genrand_real2()

        ! i,j initially refer to the indices in the lists of occupied spin-orbitals
        ! rather than the spin-orbitals.

        ind = int(r*nalpha*nbeta) + 1
        i = int( (ind-1.0_p)/nbeta ) + 1
        j = ind - (i-1)*nbeta

        ! i,j are the electrons we're exciting.  Find the occupied corresponding
        ! spin-orbitals.
        i = occ_list_alpha(i)
        j = occ_list_beta(j)

        ! Symmetry info is a simple lookup...
        ij_sym = sym_table((i+1)/2,(j+1)/2)

    end subroutine choose_ij_hub_k

    subroutine choose_ab_hub_k(f, unocc_list_alpha, ij_sym, a, b)

        ! Choose a random pair of (a,b) unoccupied virtual spin-orbitals into
        ! which electrons are excited.
        ! (a,b) are chosen such that the (i,j)->(a,b) excitation is symmetry-
        ! allowed.
        ! In:
        !    f(basis_length): bit string representation of the Slater
        !        determinant.
        !    unocc_alpha: integer list of the unoccupied alpha spin-orbitals.
        !    ij_sym: symmetry spanned by the (i,j) combination of unoccupied
        !        spin-orbitals into which electrons are excited.
        ! Returns:
        !    a,b: virtual spin orbitals involved in the excitation.

        use basis, only: basis_length
        use system, only: nvirt_alpha

        integer(i0), intent(in) :: f(basis_length)
        integer, intent(in) :: unocc_list_alpha(nvirt_alpha)
        integer, intent(in) :: ij_sym
        integer, intent(out) :: a, b

        logical :: allowed_excitation

        allowed_excitation = .false.

        do while (.not.allowed_excitation)

            ! Until we find an allowed excitation.

            call find_ab_hub_k(f, unocc_list_alpha, ij_sym, a, b, allowed_excitation)

        end do

    end subroutine choose_ab_hub_k

!--- Select random orbitals in double excitations ---

    subroutine find_ab_hub_k(f, unocc_list_alpha, ij_sym, a, b, allowed_excitation)

        ! Choose a random pair of (a,b) spin-orbitals.
        ! (a,b) are chosen such that the (i,j)->(a,b) excitation is symmetry-
        ! allowed and a is a virtual spin-orbital.  As (a,b) must be one alpha
        ! orbital and one beta orbital, we can choose a to be an alpha orbital
        ! without loss of generality.
        ! In:
        !    f(basis_length): bit string representation of the Slater
        !        determinant.
        !    unocc_alpha: integer list of the unoccupied alpha spin-orbitals.
        !    ij_sym: symmetry spanned by the (i,j) combination of unoccupied
        !        spin-orbitals into which electrons are excited.
        ! Returns:
        !    a,b: spin orbitals involved in the excitation.
        !    allowed_excitation: is true if the excitation (i,j)->(a,b) is
        !        allowed (i.e. both a and b are indeed spin orbitals).

        use basis, only: basis_length, bit_lookup, nbasis
        use dSFMT_interface, only:  genrand_real2
        use system, only: nvirt_alpha
        use momentum_symmetry, only: sym_table, inv_sym

        integer(i0), intent(in) :: f(basis_length)
        integer, intent(in) :: unocc_list_alpha(nvirt_alpha)
        integer, intent(in) :: ij_sym
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
        ! Thus k_b is defined by i,j and a.  As we only consider two-band
        ! systems, b is thus defined by spin conservation.

        ! One electron must be in unocc_list_alpha, so we can use the
        ! random number just to find which unoccupied alpha orbital is in
        ! the excitation.

        r = int(genrand_real2()*(nvirt_alpha)) + 1

        a = unocc_list_alpha(r)
        ! Find corresponding beta orbital which satisfies conservation
        ! of crystal momentum.
        ka = (a+1)/2
        b = 2*sym_table(ij_sym, inv_sym(ka))

        b_pos = bit_lookup(1,b)
        b_el = bit_lookup(2,b)

        ! If b is unoccupied then have found an allowed excitation.
        allowed_excitation = .not.btest(f(b_el), b_pos)

    end subroutine find_ab_hub_k

!--- Excitation generation probabilities ---

    pure function calc_pgen_hub_k(ab_sym, f, unocc_alpha, unocc_beta) result(pgen)

        ! Calculate the generation probability of a given excitation for the
        ! Hubbard model in momentum space.  The Hubbard model is a special case
        ! as it is a two-band system and so the generation probability is
        ! independent of the virtual spin-orbitals into which electrons are
        ! excited and depends only upon the spin-orbitals from which we excite.
        !
        ! Note that all the information required for input should be available
        ! during the FCIQMC algorithm and should not be needed to be calculated.
        !
        ! Further, we assume only allowed excitations are generated.
        !
        ! In:
        !    ab_sym: symmetry spanned by the (a,b) combination of unoccupied
        !        spin-orbitals into which electrons are excited.
        !    f: bit string representation of the determinant we're exciting
        !        from.
        !    unocc_alpha, unocc_beta: integer list of the unoccupied alpha and
        !        beta (respectively) spin-orbitals.
        ! Returns:
        !    pgen: the generation probability of the excitation.  See notes in
        !        spawning.

        use basis, only: basis_length, bit_lookup, nbasis
        use system, only: nvirt, nvirt_alpha, nvirt_beta, nalpha, nbeta, nel
        use momentum_symmetry, only: sym_table, inv_sym

        real(p) :: pgen
        integer, intent(in) :: ab_sym
        integer(i0), intent(in) :: f(basis_length)
        integer, intent(in) :: unocc_alpha(nvirt_alpha), unocc_beta(nvirt_beta)

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
        !  nbasis - nel - delta_d
        !
        ! where delta_d is the number of a orbitals which are forbidden due to b being occupied.
        ! p(b|i,j) is identical.  Hence:
        !
        ! pgen = 1/(nalpha*nbeta) [ 1/(nbasis-nel-delta_d) + 1/(basis-nel-delta_d) ]
        !                       2
        !      =  ---------------------------------
        !         nalpha*nbeta*(nbasis-nel-delta_d)

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
        do a = 1, nvirt_alpha
            ka = (unocc_alpha(a)+1)/2
            b = 2*sym_table(ab_sym, inv_sym(ka))
            b_pos = bit_lookup(1,b)
            b_el = bit_lookup(2,b)
            ! Are (a,b) both unoccupied?
            if (btest(f(b_el), b_pos)) forbidden_excitations = forbidden_excitations + 1
        end do
        do b = 1, nvirt_beta
            kb = unocc_beta(b)/2
            a = 2*sym_table(ab_sym, inv_sym(kb)) - 1
            a_pos = bit_lookup(1,a)
            a_el = bit_lookup(2,a)
            ! Are (a,b) both unoccupied?
            if (btest(f(a_el), a_pos)) forbidden_excitations = forbidden_excitations + 1
        end do

        pgen = 2.0_p/(nalpha*nbeta*(nvirt - forbidden_excitations))

    end function calc_pgen_hub_k

end module excit_gen_hub_k

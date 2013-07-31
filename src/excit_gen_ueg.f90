module excit_gen_ueg

! Module for random excitation generators and related routines for the uniform
! electron gas.

! NOTE: the electron gas is wonderfully simple and only contains double
! excitations (much like, in fact, the Hubbard model).

! TODO: interface comments.

use const, only: i0, p, dp

implicit none

contains

!--- Excitation generators ---

    subroutine gen_excit_ueg(rng, cdet, pgen, connection, hmatel)

        ! WARNING: this is not complete.  If you wish to use it, please
        ! implementation the calculation of the generation probabilities.

        ! Create a random excitation from cdet and calculate both the probability
        ! of selecting that excitation and the Hamiltonian matrix element.

        ! In:
        !    cdet: info on the current determinant (cdet) that we will gen
        !        from.
        !    parent_sign: sign of the population on the parent determinant (i.e.
        !        either a positive or negative integer).
        ! In/Out:
        !    rng: random number generator.
        ! Out:
        !    pgen: probability of generating the excited determinant from cdet.
        !    connection: excitation connection between the current determinant
        !        and the child determinant, on which progeny are gened.
        !    hmatel: < D | H | D' >, the Hamiltonian matrix element between a
        !    determinant and a connected determinant in molecular systems.

        use determinants, only: det_info
        use excitations, only: excit
        use excitations, only: find_excitation_permutation2
        use hamiltonian_ueg, only: slater_condon2_ueg_excit
        use system, only: ndim

        use dSFMT_interface, only: dSFMT_t, get_rand_close_open

        type(det_info), intent(in) :: cdet
        type(dSFMT_t), intent(inout) :: rng
        real(p), intent(out) :: pgen, hmatel
        type(excit), intent(out) :: connection
        logical :: allowed_excitation

        integer :: ij_k(ndim), ij_spin, max_na
        real(dp) :: r

        ! 1. Must have a double excitation.
        connection%nexcit = 2

        ! 2. Select orbitals to excite from and into.
        call choose_ij_k(rng, cdet%occ_list, connection%from_orb(1), connection%from_orb(2), ij_k, ij_spin)
        call choose_ab_ueg(rng, cdet%f, connection%from_orb(1), connection%from_orb(2), ij_k, ij_spin, max_na, &
                           connection%to_orb(1), connection%to_orb(2), allowed_excitation)
        
        if (allowed_excitation) then

            ! 3. Return generation probability and conecting matrix element.

            ! TODO.
!            pgen = calc_pgen_ueg()

            call find_excitation_permutation2(cdet%f, connection)
            hmatel = slater_condon2_ueg_excit(connection%from_orb(1), connection%from_orb(2), &
                                              connection%to_orb(1), connection%to_orb(2), connection%perm)

        else

            ! Carelessly selected ij with no available excitations.  Not worth
            ! renormalising.  Discard: return a null excitation.
            hmatel = 0.0_p
            pgen = 1.0_p

        end if

    end subroutine gen_excit_ueg

    subroutine gen_excit_ueg_no_renorm(rng, cdet, pgen, connection, hmatel)

        ! Create a random excitation from cdet and calculate both the probability
        ! of selecting that excitation and the Hamiltonian matrix element.

        ! This doesn't exclude the case where, having selected all orbitals
        ! involved in the excitation, the final orbital selected is already
        ! occupied and so cannot be excited into.  Whilst this is somewhat
        ! wasteful (generating excitations which can't be performed), there is
        ! a balance between the cost of generating forbidden excitations and the
        ! O(N) cost of renormalising the generation probabilities.

        ! In:
        !    cdet: info on the current determinant (cdet) that we will gen
        !        from.
        !    parent_sign: sign of the population on the parent determinant (i.e.
        !        either a positive or negative integer).
        ! In/Out:
        !    rng: random number generator.
        ! Out:
        !    pgen: probability of generating the excited determinant from cdet.
        !    connection: excitation connection between the current determinant
        !        and the child determinant, on which progeny are gened.
        !    hmatel: < D | H | D' >, the Hamiltonian matrix element between a
        !    determinant and a connected determinant in molecular systems.

        use determinants, only: det_info
        use excitations, only: excit
        use excitations, only: find_excitation_permutation2
        use hamiltonian_ueg, only: slater_condon2_ueg_excit
        use system, only: ndim

        use dSFMT_interface, only: dSFMT_t, get_rand_close_open
        use bit_utils

        type(det_info), intent(in) :: cdet
        type(dSFMT_t), intent(inout) :: rng
        real(p), intent(out) :: pgen, hmatel
        type(excit), intent(out) :: connection

        logical :: allowed_excitation
        integer :: ij_k(ndim), ij_spin, max_na
        real(dp) :: r

        if (sum(count_set_bits(cdet%f)) /= 4) then
            write (6,*) 'HUH?'
            write (6,'(2(B32.32,2X))') cdet%f
            stop
        end if

        ! 1. Must have a double excitation.
        connection%nexcit = 2

        ! 2. Select orbitals to excite from and into.
        call choose_ij_k(rng, cdet%occ_list, connection%from_orb(1), connection%from_orb(2), ij_k, ij_spin)
        call find_ab_ueg(rng, cdet%f, connection%from_orb(1), connection%from_orb(2), ij_k, ij_spin, max_na, &
                           connection%to_orb(1), connection%to_orb(2), allowed_excitation)
        
        if (allowed_excitation) then

            ! 3. Return generation probability and conecting matrix element.

            pgen = calc_pgen_ueg_no_renorm(max_na)

            call find_excitation_permutation2(cdet%f, connection)
            hmatel = slater_condon2_ueg_excit(connection%from_orb(1), connection%from_orb(2), &
                                              connection%to_orb(1), connection%to_orb(2), connection%perm)


        else

            ! Carelessly selected ij with no available excitations.  Not worth
            ! renormalising.  Discard: return a null excitation.
            hmatel = 0.0_p
            pgen = 1.0_p

        end if

    end subroutine gen_excit_ueg_no_renorm

!--- Select random orbitals involved in a valid double excitation ---

    subroutine choose_ab_ueg(rng, f, i, j, ij_k, ij_spin, max_na, a, b, allowed_excitation)

        use basis, only: basis_fns, basis_length, nbasis, basis_lookup
        use const, only: i0_end
        use bit_utils, only: count_set_bits
        use system, only: nel, ndim
        use dSFMT_interface, only: dSFMT_t, get_rand_close_open
        use ueg_system, only: ueg_ternary_conserve, ueg_basis_index
        use utils, only: tri_ind

        type(dSFMT_t), intent(inout) :: rng
        integer(i0), intent(in) :: f(basis_length)
        integer, intent(in) :: i, j, ij_k(ndim), ij_spin
        integer, intent(out) :: a, b, max_na
        logical, intent(out) :: allowed_excitation

        integer, parameter :: max_attempts = 200
        integer :: iattempt

        call find_ab_ueg(rng, f, i, j, ij_k, ij_spin, max_na, a, b, allowed_excitation)
        if (max_na == 0) then
            ! No available excitations.
        else if (.not.allowed_excitation) then
            ! Attempt max_attempts times to find an excitation.
            do iattempt = 2, max_attempts
                call find_ab_ueg(rng, f, i, j, ij_k, ij_spin, max_na, a, b, allowed_excitation)
                if (allowed_excitation) exit
            end do
        end if

    end subroutine choose_ab_ueg

    subroutine choose_ij_k(rng, occ_list, i, j, ij_k, ij_spin)

        ! Randomly choose a pair of spin-orbitals for 1-band systems with
        ! Bloch orbitals.
        !
        ! See choose_ij_hub_k for a specific procedure for the momentum space
        ! formulation of the hubbard model.
        !
        ! In:
        !    occ_list: integer list of occupied spin-orbitals in a determinant.
        ! In/Out:
        !    rng: random number generator.
        ! Out:
        !    i, j: randomly selected spin-orbitals.
        !    ij_spin: -2 if (i,j) are both beta, +2 if (i,j) are both alpha and
        !        0 if it's a "mixed" excitation, i.e. alpha, beta or beta,
        !        alpha.

        use basis, only: basis_fns
        use system, only: nel, ndim
        use dSFMT_interface, only: dSFMT_t, get_rand_close_open

        integer, intent(in) :: occ_list(nel)
        type(dSFMT_t), intent(inout) :: rng
        integer, intent(out) :: i,j, ij_k(ndim), ij_spin
        integer :: ind
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

        r = get_rand_close_open(rng)
        ind = int(r*nel*(nel-1)/2) + 1

        ! i,j initially refer to the indices in the lists of occupied spin-orbitals
        ! rather than the spin-orbitals.
        ! Note that the indexing scheme for the strictly lower triangular array
        ! assumes j>i.  As occ_list is ordered, this means that we will return
        ! i,j (referring to spin-orbitals) where j>i.  This ordering is
        ! convenient subsequently, e.g. is assumed in the
        ! find_excitation_permutation2 routine.
        j = int(1.50_p + sqrt(2*ind-1.750_p))
        i = ind - ((j-1)*(j-2))/2

        ! i,j are the electrons we're exciting.  Find the occupied corresponding
        ! spin-orbitals.
        i = occ_list(i)
        j = occ_list(j)

        ij_spin = basis_fns(i)%Ms + basis_fns(j)%Ms
        ij_k = basis_fns(i)%l + basis_fns(j)%l 

    end subroutine choose_ij_k

!--- Select random orbitals in double excitations ---

    subroutine find_ab_ueg(rng, f, i, j, ij_k, ij_spin, max_na, a, b, allowed_excitation)

        use basis, only: basis_fns, basis_length, nbasis, basis_lookup, bit_lookup
        use const, only: i0_end
        use bit_utils, only: count_set_bits
        use system, only: nel, ndim
        use dSFMT_interface, only: dSFMT_t, get_rand_close_open
        use ueg_system, only: ueg_ternary_conserve, ueg_basis_index
        use utils, only: tri_ind

        type(dSFMT_t), intent(inout) :: rng
        integer(i0), intent(in) :: f(basis_length)
        integer, intent(in) :: i, j, ij_k(ndim), ij_spin
        integer, intent(out) :: a, b, max_na
        logical, intent(out) :: allowed_excitation

        integer :: fac, shift, ij_ind, ibp, ibe, n, ind, kb(ndim)
        integer(i0) :: poss_a(basis_length)

        ! Let's just check there are possible a,b first!
        if (j>i) then
            ij_ind = tri_ind((j+1)/2, (i+1)/2)
        else
            ij_ind = tri_ind((i+1)/2, (j+1)/2)
        end if

        ! Adjust for spin.  Allow a to be up (without bias) if possible.
        if (ij_spin == -2) then
            poss_a = iand(not(f), ishft(ueg_ternary_conserve(1:,ij_ind),1))
        else
            poss_a = iand(not(f), ueg_ternary_conserve(1:,ij_ind))
        end if
        max_na = sum(count_set_bits(poss_a))

        if (max_na > 0) then

            ! Select a random a.
            a = int(max_na*get_rand_close_open(rng)) + 1

            n = 0
            finda: do ibe = 1, basis_length
                if (poss_a(ibe) /= 0_i0) then
                    do ibp = 0, i0_end
                        if (btest(poss_a(ibe), ibp)) n = n + 1
                        if (n == a) then
                            ! Found our basis funtion!
                            a = basis_lookup(ibp, ibe)
                            exit finda
                        end if
                    end do
                end if
            end do finda

            ! Ok, b is now completely defined.
            kb = ij_k - basis_fns(a)%l

            ! Now, we have to be careful!
            ! max_na is the maximum number of ways of choosing a.
            ! However, we then use it to calculate the number of ways of
            ! generating this unique excitation.  As the UEG is a one-band
            ! model (and hence b is uniquely defined), then if a and b are of
            ! the same spin, they were both possible choices for a.  We hence
            ! need to half max_na in order to give the number of unique possible
            ! excitations given the choice of i and j.
            select case(ij_spin)
            case(2)
                b = ueg_basis_index(kb, 1)
                max_na = max_na / 2
            case(0)
                b = ueg_basis_index(kb, -1)
            case(-2)
                b = ueg_basis_index(kb, -1)
                max_na = max_na / 2
            end select

            ! Excitation is forbidden if b is already occupied!
            ibp = bit_lookup(1,b)
            ibe = bit_lookup(2,b)
            allowed_excitation = .not.btest(f(ibe), ibp)

            ! It is useful to return a,b ordered (e.g. for the find_excitation_permutation2 routine).
            if (a > b) then
                ind = a
                a = b
                b = ind
            end if

        else

            allowed_excitation = .false.

        end if

    end subroutine find_ab_ueg

!--- Excitation generation probabilities ---

    pure function calc_pgen_ueg_no_renorm(max_na) result(pgen)

        use system, only: nel

        real(p) :: pgen
        integer, intent(in) :: max_na

        ! We explicitly reject excitations i,j->a,b where b is already
        ! occupied.

        ! **WARNING**
        ! As an implementation detail, we have required a to be up spin if possible.
        ! Hence, given we require a to be unoccupied (but not b), the generation
        ! probability is asymmetric:
        ! pgen = p(i,j) p(a|i,j) p(b|i,j,a)

        ! p(i,j) = 1/binom(nel,2) = 2/(nel*(nel-1))
        ! p(a|i,j) is the number of unoccupied a orbitals which i,j can excite
        ! into.
        ! UEG is a one-band system, so p(b|i,j,a) = 1.

        pgen = 2.0_p / ( nel*(nel-1) * max_na )

    end function calc_pgen_ueg_no_renorm

end module excit_gen_ueg

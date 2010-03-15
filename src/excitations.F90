module excitations

! Module for dealing with excitations.

use const

implicit none

! A handy type for containing the excitation information needed to connect one
! determinant to another.
type excit
    ! Excitation level.
    integer :: nexcit
    ! Orbitals which are excited from and to.
    ! Only used for single and double excitations; undefined otherwise.
    integer :: from_orb(2), to_orb(2)
    ! True if a total odd number of permutations is required to align
    ! the determinants.  Only used for single and double excitations.
    ! Undefined otherwise.
    logical :: perm
end type excit

contains

    pure function get_excitation(f1,f2) result(excitation)

        ! In: 
        !    f1(basis_length): bit string representation of the Slater
        !        determinant.
        !    f2(basis_length): bit string representation of the Slater
        !        determinant.
        ! Returns:
        !    excitation: excit type containing the following information---
        !        excitation%nexcit: excitation level.
        !
        !    If the excitation is a single or double excitation then it also
        !    includes:
        ! 
        !        excitation%from_orbs(2): orbitals excited from in f1.
        !        excitation%to_orbs(2): orbitals excited to in f2.
        !        excitation%perm: true if an odd number of permutations are
        !            reqiured to align the determinants.
        !        The second element of from_orbs and to_orbs is zero for single
        !        excitations.

        use bit_utils
        use basis, only: basis_length, basis_lookup
        use system, only: nel

        type(excit) :: excitation
        integer(i0), intent(in) :: f1(basis_length), f2(basis_length)
        integer :: i, j, iexcit1, iexcit2, perm, iel1, iel2, shift
        logical :: test_f1, test_f2

        excitation = excit(0, 0, 0, .false.)

        if (any(f1/=f2)) then

            iexcit1 = 0
            iexcit2 = 0
            iel1 = 0
            iel2 = 0
            perm = 0

            ! Excitation level...
#ifdef _PGI
            ! Work round an *insane* bug in PGI where intrinsic bit operations
            ! return an integer(4) if the arguments are of a kind smaller than
            ! 4.  PGI gets it right if the kind is larger than 4, but that
            ! doesn't help us always in this case...
            excitation%nexcit = sum(count_set_bits(int(ieor(f1,f2),i0)))/2
#else
            excitation%nexcit = sum(count_set_bits(ieor(f1,f2)))/2
#endif

            ! Finding the permutation to align the determinants is non-trivial.
            ! It turns out to be quite easy with bit operations.
            ! The idea is to do a "dumb" permutation where the determinants are 
            ! expressed in two sections: orbitals not involved in the excitation
            ! and those that are.  Each section is stored in ascending index
            ! order.
            ! To obtain such ordering requires (for each orbital that is
            ! involved in the excitation) a total of
            ! nel - iel - nexcit + iexcit
            ! where nel is the number of electrons, iel is the position of the 
            ! orbital within the list of occupied states in the determinant,
            ! nexcit is the total number of excitations and iexcit is the number
            ! of the "current" orbital involved in excitations.
            ! e.g. Consider (1, 2, 3, 4, 5) -> (1, 3, 5, 6, 7).
            ! (1, 2, 3, 4) goes to (1, 3, 2, 4).
            ! 2 is the first (iexcit=1) orbital found involved in the excitation
            ! and so requires 5 - 2 - 2 + 1 = 2 permutation to shift it to the
            ! first "slot" in the excitation "block" in the list of states.
            ! 4 is the second orbital found and requires 5 - 4 - 2 + 2 = 1
            ! permutation to shift it the end (last "slot" in the excitation
            ! block).
            ! Whilst the resultant number of permutations isn't necessarily the
            ! minimal number for the determinants to align, this is irrelevant
            ! as the Slater--Condon rules only care about whether the number of
            ! permutations is even or odd.
            shift = nel - excitation%nexcit 

            if (excitation%nexcit <= 2) then

                do i = 1, basis_length
                    ! Bonus optimisation: skip bit strings which aren't changed.
                    ! This modifies the algorithm above by skipping basis
                    ! functions which are already lined up.  Doing so is
                    ! guaranteed to introduce an additional even number of
                    ! permuations (as we underestimate both iel1 and iel2 by the
                    ! same number of electrons).
                    if (f1(i) == f2(i)) cycle
                    do j = 0, i0_end

                        test_f1 = btest(f1(i),j)
                        test_f2 = btest(f2(i),j)

                        if (test_f2) iel2 = iel2 + 1

                        if (test_f1) then
                            iel1 = iel1 + 1
                            if (.not.test_f2) then
                                ! occupied in f1 but not in f2
                                iexcit1 = iexcit1 + 1
                                excitation%from_orb(iexcit1) = basis_lookup(j,i)
                                perm = perm + (shift - iel1 + iexcit1)
                            end if
                        else
                            if (test_f2) then
                                ! occupied in f1 but not in f2
                                iexcit2 = iexcit2 + 1
                                excitation%to_orb(iexcit2) = basis_lookup(j,i)
                                perm = perm + (shift - iel2 + iexcit2)
                            end if
                        end if

                    end do
                end do

                ! It seems that this test is faster than btest(perm,0)!
                excitation%perm = mod(perm,2) == 1

            end if
        end if

    end function get_excitation

    function get_excitation_level(f1, f2) result(level)

        ! In: 
        !    f1(basis_length): bit string representation of the Slater
        !        determinant.
        !    f2(basis_length): bit string representation of the Slater
        !        determinant.
        ! Returns:
        !    Excitation level connecting determinants f1 and f2.     

        use basis, only: basis_length
        use bit_utils, only: count_set_bits

        integer :: level
        integer(i0), intent(in) :: f1(basis_length), f2(basis_length)

#ifdef _PGI
        ! Work round an *insane* bug in PGI where intrinsic bit operations
        ! return an integer(4) if the arguments are of a kind smaller than
        ! 4.  PGI gets it right if the kind is larger than 4, but that
        ! doesn't help us always in this case...
        level = sum(count_set_bits(int(ieor(f1,f2),i0)))/2
#else
        level = sum(count_set_bits(ieor(f1,f2)))/2
#endif

    end function get_excitation_level

    subroutine find_excitation_permutation1(occ_list, excitation)

        ! Find the parity of the permutation required to maximally line up
        ! a determinant with an excitation of it, as needed for use with the
        ! Slater--Condon rules.
        !
        ! This version is for single excitations of a determinant.
        !
        ! In:
        !    occ_list: integer list of occupied orbitals in the Slater determinant.
        !    excitation: excit type specifying how the excited determinant is
        !        connected to the rdeterminant given in occ_list.
        ! Out:
        !    excitation: excit type with the parity of the permutation also
        !        specified.
        use system, only: nel
    
        integer, intent(in) :: occ_list(nel)
        type(excit), intent(inout) :: excitation

        integer :: i, shift, perm, j, j1, nholes

        ! Adapt algorithm from get_excitation and find_excitation_permutation2.

        ! As we only consider single excitations, this is much easier than
        ! find_excitation_permutation2.

        ! The comments in find_excitation_permutation2 apply here, except that
        ! there only one slot is needed at the end of the occ_list to hold
        ! orbitals involved in the excitation.
        
        shift = nel - excitation%nexcit 

        perm = 0
        nholes = 0
        do i = 1, nel
            j = occ_list(i)

            ! Check for inserting excited orbital.
            if (j > excitation%to_orb(1) .and. j1 < excitation%to_orb(1)) then
                ! Number of permutations to get to the end of the list.
                perm = perm + shift - (i - nholes) + 1
            end if

            ! Check for orbital exciting from.
            if (j == excitation%from_orb(1)) then
                ! Number of permutations to get to the end of the list.
                perm = perm + shift - i + 1
                nholes = 1
            end if

            j1 = j
        end do

        excitation%perm = mod(perm,2) == 1

    end subroutine find_excitation_permutation1

    subroutine find_excitation_permutation2(occ_list, excitation)

        ! Find the parity of the permutation required to maximally line up
        ! a determinant with an excitation of it, as needed for use with the
        ! Slater--Condon rules.
        !
        ! This version is for double excitations of a determinant.
        !
        ! In:
        !    occ_list: integer list of occupied orbitals in the Slater determinant.
        !    excitation: excit type specifying how the excited determinant is
        !        connected to the rdeterminant given in occ_list.
        !        Note that we require the lists of orbitals excited from/into to
        !        be ordered.
        ! Out:
        !    excitation: excit type with the parity of the permutation also
        !        specified.

        use system, only: nel
    
        integer, intent(in) :: occ_list(nel)
        type(excit), intent(inout) :: excitation

        integer :: i, j, j1, shift, perm, nholes

        shift = nel - excitation%nexcit 

        ! Adapt algorithm from get_excitation.
        ! The idea is to count the number of permutations involved in shifting
        ! orbitals involved in the excitation to the end of the list (thus
        ! giving the 2 determinants with maximum coincidence).

        ! The last 2 positions in the list are for orbitals involved in the
        ! excitation.
        ! We count the number of permutations required to shift i and a into the
        ! first slot and j and b into the second slot.

        ! The number of permutations required is simply the number of orbitals
        ! between a given orbital and the slot it's going into.
        ! This is given by:
        !   nel - 2 - position + slot_number
        ! where 2 comes from the fact we're considering double excitations and 
        ! slot_number is 1 or 2.

        ! The position of the i and j orbitals is easy, as occ_list is ordered.
        ! The position of the a and b orbitals is harder, as occ_list refers to
        ! the determinant we're exciting from.  However, this is still more
        ! efficient than generating the list for the excited determinant or
        ! testing all spin orbitals (as in get_excitation).

        ! The positions of a and b are given by:
        !   k - nholes
        ! where k is the position at which a/b is inserted and nholes is the
        ! number of orbitals we've excited from *before* k (i.e. not including
        ! k).

        ! This is all a counting exercise: occ_list isn't actually altered.

        perm = 0
        nholes = 0
        j1 = 0
        do i = 1, nel
            j = occ_list(i)

            ! First check if we insert a new electron.
            ! This works round the problem if we have to insert a
            ! before j and then b replaces j.  This order keeps the number of
            ! holes correct.
            if (j > excitation%to_orb(1) .and. j1 < excitation%to_orb(1)) then
                ! Number of permutations to get to the end of the list.
                perm = perm + shift - (i - nholes) + 1
                nholes = nholes - 1
            end if
            if (j > excitation%to_orb(2) .and. j1 < excitation%to_orb(2)) then
                ! Number of permutations to get to the end of the list.
                perm = perm + shift - (i - nholes) + 2
            end if

            if (j == excitation%from_orb(1)) then
                ! Number of permutations to get to the end of the list.
                perm = perm + shift - i + 1
                nholes = nholes + 1
            else if (j == excitation%from_orb(2)) then
                ! Number of permutations to get to the end of the list.
                perm = perm + shift - i + 2
                nholes = nholes + 1
            end if

            j1 = j
        end do

        excitation%perm = mod(perm,2) == 1

    end subroutine find_excitation_permutation2

    subroutine create_excited_det(f_in, connection, f_out)

        ! Generate a determinant from another determinant and the excitation
        ! information connecting the two determinants.
        ! In: 
        !    f_in(basis_length): bit string representation of the reference
        !        Slater determinant.
        !    connection: excitation connecting f_in to f_out.  Note that
        !        the perm field is not used.
        ! Out:
        !    f_out(basis_length): bit string representation of the excited
        !        Slater determinant.

        use basis, only: basis_length, bit_lookup

        integer(i0), intent(in) :: f_in(basis_length)
        type(excit), intent(in) :: connection
        integer(i0), intent(out) :: f_out(basis_length)

        integer :: i, orb, bit_pos, bit_element

        ! Unset the orbitals which are excited from and set the orbitals which
        ! are excited into.
        f_out = f_in
        do i = 1, connection%nexcit
            ! Clear i/j orbital.
            orb = connection%from_orb(i)
            bit_pos = bit_lookup(1,orb)
            bit_element = bit_lookup(2,orb)
            f_out(bit_element) = ibclr(f_out(bit_element), bit_pos)
            ! Set a/b orbital.
            orb = connection%to_orb(i)
            bit_pos = bit_lookup(1,orb)
            bit_element = bit_lookup(2,orb)
            f_out(bit_element) = ibset(f_out(bit_element), bit_pos)
        end do

    end subroutine create_excited_det

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
        use system, only: nvirt_alpha, nvirt_beta, nalpha, nbeta, nel
        use symmetry, only: sym_table, inv_sym

        real(dp) :: pgen
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

        pgen = 2.0_dp/(nalpha*nbeta*(nbasis-nel-forbidden_excitations))

    end function calc_pgen_hub_k

    function calc_pgen_hub_real(occ_list, f, nvirt_avail) result(pgen)

        use basis, only: basis_length
        use system, only: nel
        use hubbard_real, only: connected_orbs

        real(dp) :: pgen
        integer, intent(in) :: occ_list(nel)
        integer(i0), intent(in) :: f(basis_length)
        integer, intent(in) :: nvirt_avail
        integer :: i, no_excit

        ! For single excitations
        !   pgen = p(i) p(a|i) \chi_r
        ! where
        !   p(i) is the probability of choosing the i-th electron to excite
        !   from.
        !   p(i) = 1/nel
        !
        !   p(a|i) is the probability of choosing to excite into orbital a given
        !   that the i-th electron has been chosen.
        !   p(a|i) is the number of virtual orbitals connected to i and is
        !   calculated when the random excitation is chosen.

        !   \chi_r is a renormalisation to take into account the fact that not
        !   all electrons may be excited from (e.g. no connected orbitals are
        !   vacant).
        !   \chi_r = nel/(nel - no_excit)
        !   where no_excit is the number of occupied orbitals which have no
        !   connected excitations.

        ! \chi_r is a constant for a given determinant, so an optimisation is to
        ! calculate this once per determinant rather than for each walker on the
        ! same determinant.

        no_excit = 0
        do i = 1, nel
            ! See if there are any allowed excitations from this electron.
            ! (see notes in choose_ia_hub_real for how this works)
            if (all(ieor(f, connected_orbs(:,occ_list(i))) == 0)) then
                ! none allowed from this orbial
                no_excit = no_excit + 1
            end if
        end do

        pgen = 1.0_dp/(nvirt_avail * (nel - no_excit))

    end function calc_pgen_hub_real

end module excitations

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
        use determinants, only: basis_length, basis_lookup
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
        use symmetry, only: sym_table

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
        !  nalpha*nbeta + nbeta*nalpha (i.e. i can be alpha or beta).
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
        ! pgen = 1/(2*nalpha*nbeta) [ 1/(nbasis-nel-delta_d) + 1/(basis-nel-delta_d) ]
        !                       1
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
            b = 2*sym_table(ab_sym, ka)
            b_pos = bit_lookup(1,b)
            b_el = bit_lookup(2,b)
            ! Are (a,b) both unoccupied?
            if (btest(f(b_el), b_pos)) forbidden_excitations = forbidden_excitations + 1
        end do
        do b = 1, nvirt_beta
            kb = unocc_beta(b)/2
            a = 2*sym_table(ab_sym, kb) - 1
            a_pos = bit_lookup(1,a)
            a_el = bit_lookup(2,a)
            ! Are (a,b) both unoccupied?
            if (btest(f(a_el), a_pos)) forbidden_excitations = forbidden_excitations + 1
        end do

        pgen = 1.0_dp/(nalpha*nbeta*(nbasis-nel-forbidden_excitations))

    end function calc_pgen_hub_k

end module excitations

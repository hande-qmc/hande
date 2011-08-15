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
    ! In the case of folded spectrum code, there can be up to a quadruple excitation
    ! *NB* the only subroutines which work sensibly when handed a triple or quadruple
    !      excitation are those pointed to by create_spawned_particle_ptr.
    integer :: from_orb(4), to_orb(4)
    ! True if a total odd number of permutations is required to align
    ! the determinants.  Only used for single and double excitations.
    ! Undefined otherwise.
    logical :: perm
end type excit

! excit_mask(:,i) is a bit field with bits corresponding to all orbitals with
! a higher index than i set.
integer(i0), allocatable :: excit_mask(:,:) ! (basis_length, nbasis)

contains

    subroutine init_excitations()

        ! Allocate and initialise data in excit_mask.

        use basis, only: bit_lookup, nbasis, basis_length
        use checking, only: check_allocate

        integer :: ibasis, jbasis, ipos, iel, jpos, jel, ierr

        allocate(excit_mask(basis_length, nbasis), stat=ierr)
        call check_allocate('excit_mask', basis_length*nbasis, ierr)

        excit_mask = 0

        do ibasis = 1, nbasis
            ipos = bit_lookup(1, ibasis)
            iel = bit_lookup(2, ibasis)
            ! Set bits corresponding to all orbitals above ibasis.
            ! Sure, there are quicker (and probably more elegant) ways of doing
            ! this, but it's a one-off...
            ! Loop from jbasis=1 as separate_strings means that even if
            ! jbasis<ibasis, it can still come after ibasis in the bit string.
            do jbasis = 1, nbasis
                jpos = bit_lookup(1, jbasis)
                jel = bit_lookup(2, jbasis)
                if ( (jel==iel .and. jpos > ipos) .or. jel>iel) &
                    excit_mask(jel, ibasis) = ibset(excit_mask(jel, ibasis), jpos)
            end do
        end do

    end subroutine init_excitations


    subroutine end_excitations()

        ! Deallocate excit_mask.

        use checking, only: check_deallocate

        integer :: ierr

        deallocate(excit_mask, stat=ierr)
        call check_deallocate('excit_mask', ierr)

    end subroutine end_excitations


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
        !        excitation%from_orb(2): orbitals excited from in f1.
        !        excitation%to_orb(2): orbitals excited to in f2.
        !        excitation%perm: true if an odd number of permutations are
        !            reqiured to align the determinants.
        !        The second element of from_orb and to_orb is zero for single
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
#ifdef PGI
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

    pure function get_excitation_level(f1, f2) result(level)

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

#ifdef PGI
        ! Work round an *insane* bug in PGI where intrinsic bit operations
        ! return an integer(4) if the arguments are of a kind smaller than
        ! 4.  PGI gets it right if the kind is larger than 4, but that
        ! doesn't help us always in this case...
        level = sum(count_set_bits(int(ieor(f1,f2),i0)))/2
#else
        level = sum(count_set_bits(ieor(f1,f2)))/2
#endif

    end function get_excitation_level

    pure subroutine find_excitation_permutation1(f, excitation)

        ! Find the parity of the permutation required to maximally line up
        ! a determinant with an excitation of it, as needed for use with the
        ! Slater--Condon rules.
        !
        ! This version is for single excitations of a determinant.
        !
        ! In:
        !    f: bit string representation of the determinant.
        !    excitation: excit type specifying how the excited determinant is
        !        connected to the determinant given in occ_list.
        ! Out:
        !    excitation: excit type with the parity of the permutation also
        !        specified.

        use basis, only: basis_length
        use bit_utils, only: count_set_bits
    
        integer(i0), intent(in) :: f(basis_length)
        type(excit), intent(inout) :: excitation

        integer :: perm
        integer(i0) :: ia(basis_length)

        ! This is just a simplification of find_excitation_permutation2.  See
        ! the comments there (and ignore any that refer to j and b...).

        ia = ieor(excit_mask(:,excitation%from_orb(1)),excit_mask(:,excitation%to_orb(1)))
        perm = sum(count_set_bits(iand(f,ia)))
        if (excitation%from_orb(1) > excitation%to_orb(1)) perm = perm - 1
        excitation%perm = mod(perm,2) == 1

    end subroutine find_excitation_permutation1

    pure subroutine find_excitation_permutation2(f, excitation)

        ! Find the parity of the permutation required to maximally line up
        ! a determinant with an excitation of it, as needed for use with the
        ! Slater--Condon rules.
        !
        ! This version is for double excitations of a determinant.
        !
        ! In:
        !    f: bit string representation of the determinant.
        !    excitation: excit type specifying how the excited determinant is
        !        connected to the determinant described by f.
        !        Note that we require the lists of orbitals excited from/into to
        !        be ordered.
        ! Out:
        !    excitation: excit type with the parity of the permutation also
        !        specified.

        use basis, only: basis_length
        use bit_utils, only: count_set_bits
    
        integer(i0), intent(in) :: f(basis_length)
        type(excit), intent(inout) :: excitation

        integer :: perm
        integer(i0) :: ia(basis_length), jb(basis_length)

        ! Fast way of getting the parity of the permutation required to align
        ! two determinants given one determinant and the connecting exctitation.
        ! This is hard to generalise to all cases, but we actually only care
        ! about single and double excitations.  The idea is quite different from
        ! that used in get_excitation (where we also need to find the orbitals
        ! involved in the excitation).

        ! In the following & represents the bitwise and operation; ^ the bitwise exclusive or
        ! operation; xmask is a mask with all bits representing orbitals above
        ! x set; f is the string representing the determinant from which we
        ! excite and the excitation is defined by (i,j)->(a,b), where i<j and
        ! a<b.

        ! imask ^ amask returns a bit string with bits corresponding to all
        ! orbitals between i and a set, with max(i,a) set and min(i,a) cleared.
        ! Thus f & (imask ^ amask) returns a bit string with only bits set for
        ! the occupied orbitals which are between i and a (possibly including i)
        ! and so the popcount of this gives the number of orbitals between i and
        ! a (possibly one larger than the actual answer) number of permutations needed to
        ! align i and a in the same 'slot' in the determinant string.  We need
        ! to subtract one if i>a to correct for the overcounting.

        ! An analagous approach counts the number of permutations required so
        ! j and b are coincident.

        ! Finally, we need to account for some more overcounting/undercounting.
        ! If j is between i and a, then it is counted yet j can either be moved
        ! before i (resulting in the actual number of permutations being one
        ! less than that counted) or after i (resulting in moving j taking one
        ! more permutation than counted).  It doesn't matter which we do, as we
        ! are only interested in whether the number of permutations is odd or
        ! even.  We similarly need to take into account the case where i is
        ! between j and b.

        ia = ieor(excit_mask(:,excitation%from_orb(1)),excit_mask(:,excitation%to_orb(1)))
        jb = ieor(excit_mask(:,excitation%from_orb(2)),excit_mask(:,excitation%to_orb(2)))

        perm = sum(count_set_bits(iand(f,ia))) + sum(count_set_bits(iand(f,jb)))

        if (excitation%from_orb(1) > excitation%to_orb(1)) perm = perm - 1
        if (excitation%from_orb(1) > excitation%to_orb(2)) perm = perm - 1
        if (excitation%from_orb(2) > excitation%to_orb(2) .or. &
            excitation%from_orb(2) < excitation%to_orb(1)) perm = perm - 1

        excitation%perm = mod(perm,2) == 1

    end subroutine find_excitation_permutation2

    pure subroutine create_excited_det(f_in, connection, f_out)

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

    subroutine create_excited_det_complete(cdet_in, connection, cdet_out)
    
        ! Generate a complete excited determinant from another determinant and 
        !the excitation information connecting the two determinants.
        ! In: 
        !    cdet_in: info on the current determinant that we will excite
        !        from.  The f field must be set.
        !    connection: excitation connecting cdet_in to cdet_out.  Note that
        !        the perm field is not used.
        ! Out:
        !    cdet_out info: on the determinant that we will excite to
        use determinants, only : det_info, decode_det_spinocc_spinunocc

        type(det_info), intent(in)  :: cdet_in
        type(excit), intent(in)     :: connection
        type(det_info), intent(out) :: cdet_out

        ! Create the excited determinant bit string representation
        call create_excited_det(cdet_in%f, connection, cdet_out%f)

        ! Decode the excited determinant bit string representation
        call decode_det_spinocc_spinunocc(cdet_out%f,cdet_out)

        

    end subroutine create_excited_det_complete









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
        use symmetry, only: sym_table, inv_sym

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

    pure function calc_pgen_hub_real(occ_list, f, nvirt_avail) result(pgen)

        ! Calculate the generation probability of a given excitation for the
        ! Hubbard model in real space.
        !
        ! Note that all the information required for input should be available
        ! during the FCIQMC algorithm and should not be needed to be calculated.
        !
        ! Further, we assume only allowed excitations are generated.
        !
        ! In:
        !    occ_list: integer list of occupied orbitals in the Slater determinant.
        !    f: bit string representation of the determinant we're exciting
        !        from.
        !    nvirt_avail: the number of available orbitals into which an
        !        electron can be excited, given the choice of the orbital which 
        !        is being excited from (i.e. having chosen i, how many
        !        possibilities are there for a, where i is occupied and
        !        a occupied and D_i^a is connected to D.
        ! Returns:
        !    pgen: the generation probability of the excitation.  See notes in
        !        spawning.

        use basis, only: basis_length
        use system, only: nel
        use hubbard_real, only: connected_orbs

        use errors

        real(p) :: pgen
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
        !   nvirt_avail is the number of virtual orbitals connected to i and is
        !   calculated when the random excitation is chosen.
        !   p(a|i) = 1/nvirt_avail

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
            if (all(iand(not(f), connected_orbs(:,occ_list(i))) == 0)) then
                ! none allowed from this orbial
                no_excit = no_excit + 1
            end if
        end do

        pgen = 1.0_p/(nvirt_avail * (nel - no_excit))

    end function calc_pgen_hub_real

    pure subroutine enumerate_all_excitations_hub_real(cdet, max_excit, excitations)

        ! Find all excitations connected to a determinant constructed from the
        ! real-space (atomic) spin-orbitals.

        ! In:
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.  The f and occ_list fields must be set.
        ! Out:
        !    max_excit: the number of possible excitations from the determinant.
        !    excitations: array of excit variables containing the excitation
        !        information.  Note that only single excitations are allowed, so
        !        the nexcit field is not set, the second element in the from_orb
        !        and to_orb filed is not set, and the permutation field is also
        !        not set, as it's quite expensive to evaluate and not necessary
        !        for most of the excitations.  The array must be at least of the
        !        size of the maximum number of excitations: 2*ndim*nel.

        use basis, only: basis_length, bit_lookup
        use determinants, only: det_info
        use hubbard_real, only: connected_sites
        use system, only: ndim, nel

        type(det_info), intent(in) :: cdet
        integer, intent(out) :: max_excit
        type(excit), intent(out) :: excitations(:)

        integer :: ii, i, ia, a, a_pos, a_el

        max_excit = 0

        do ii = 1, nel
            i = cdet%occ_list(i)
            ! Each orbital is connected to (at most) 2 orbitals in each
            ! direction.
            do ia = 1, 2*ndim
                a = connected_sites(ia,a)
                if (a == 0) then
                    ! run out of possible excitations from this electron
                    cycle
                else
                    a_pos = bit_lookup(1,a)
                    a_el = bit_lookup(2,a)
                    if (.not.btest(cdet%f(a_el), a_pos)) then
                        ! a is unoccupied.  Have an excitation.
                        max_excit = max_excit+1
                        excitations(max_excit)%from_orb(1) = i
                        excitations(max_excit)%to_orb(1) = a
                    end if
                end if
            end do
        end do

    end subroutine enumerate_all_excitations_hub_real

    pure subroutine enumerate_all_excitations_hub_k(cdet, max_excit, excitations)

        ! Find all excitations connected to a determinant constructed from the
        ! momentum-space (Bloch) spin-orbitals.

        ! In:
        !    cdet: info on the current determinant (cdet) that we will spawn
        !        from.  The f, occ_list_alpha, occ_list_beta and
        !        occ_list_unocc_alpha fields must be set.
        ! Out:
        !    max_excit: the number of possible excitations from the determinant.
        !    excitations: array of excit variables containing the excitation
        !        information.  Note that only double excitations are allowed, so
        !        the nexcit field is not set and the permutation field is also
        !        not set, as it's quite expensive to evaluate and not necessary
        !        for most of the excitations.  The array must be at least of the
        !        size of the maximum number of excitations:
        !        nalpha*nbeta*min(nsites-nalpha,nsites-nbeta).
        !        WARNING: the from_orb and to_orb are not ordered and must be
        !        ordered before (e.g.) find_excitation_permutation2 is called.

        use basis, only: bit_lookup
        use determinants, only: det_info
        use symmetry, only: sym_table, inv_sym
        use system, only: ndim, nel, nalpha, nvirt_alpha, nbeta

        type(det_info), intent(in) :: cdet
        integer, intent(out) :: max_excit
        type(excit), intent(out) :: excitations(:)

        integer :: ii, i, ij, j, ij_sym, ia, a, b, b_pos, b_el

        max_excit = 0

        do ii = 1, nalpha
            i = cdet%occ_list_alpha(ii)
            do ij = 1, nbeta
                j = cdet%occ_list_beta(ij)
                ij_sym = sym_table((i+1)/2,(j+1)/2)
                ! Either a must be alpha and b beta or vice versa.  Without loss of
                ! generality, choose a to always be alpha.
                do ia = 1, nvirt_alpha
                    a = cdet%unocc_list_alpha(ia)
                    ! The choice of i,j,a fixes the wavevector and spin of b.
                    b = 2*sym_table(ij_sym, inv_sym((a+1)/2))
                    b_pos = bit_lookup(1,b)
                    b_el = bit_lookup(2,b)
                    if (.not.btest(cdet%f(b_el), b_pos)) then
                        ! If b is unoccupied then have found the excitation.
                        max_excit = max_excit + 1
                        excitations(max_excit)%from_orb(1:2) = (/ i, j /)
                        excitations(max_excit)%to_orb(1:2) = (/ a, b /)
                    end if
                end do
            end do
        end do

    end subroutine enumerate_all_excitations_hub_k

end module excitations

module determinants

! Generation, inspection and manipulation of Slater determinants. 

use const
use system
use basis
use parallel

implicit none

type det
    ! Bit string representation of the Slater determinant.
    ! f is used throughout to indicate a Slater determinant
    ! represented as a bit string.
    integer(i0), pointer :: f(:) ! (basis_length)
    ! Total spin of the determinant in units of electron spin (1/2).   
    integer(i0) :: Ms
end type det

! Number of possible determinants in the system.
integer :: ndets

! Store of all determinants.
! This will quickly become a memory issue, but for dealing with the FCI of small
! systems it is ok.
type(det), allocatable :: dets(:) ! (ndets)

! A handy type for containing the excitation information needed to connect one
! determinant to another.
type excit
    ! Excitation level.
    integer :: nexcit
    ! Orbitals which are excited from and to.
    ! Only used for single and double excitations.
    integer :: from_orb(2), to_orb(2)
    ! True if a total odd number of permutations is required to align
    ! the determinants.  Only used for single and double excitations.
    logical :: perm
end type excit

! If true then the determinant list is written to determinant_file.
logical :: write_determinants
character(255) :: determinant_file = 'DETS'

contains

    subroutine init_determinants()

        ! Initialise determinant information: number of determinants, how
        ! they are stored as bit strings, lookup arrays for converting from
        ! integer list of orbitals to bit strings and vice versa.

        use comb_m, only: binom

        integer :: i, bit_pos, bit_element

        ndets = binom(nbasis, nel)

        ! See note in basis.
        basis_length = nbasis/8
        if (mod(nbasis,8) /= 0) basis_length = basis_length + 1

        if (iproc == parent) then
            write (6,'(1X,a20,i4)') 'Number of electrons:', nel
            write (6,'(1X,a26,i4)') 'Number of basis functions:', nbasis
            write (6,'(1X,a26,i8,/)') 'Size of determinant space:', ndets
        end if

        ! Lookup arrays.
        allocate(bit_lookup(2,nbasis))
        allocate(basis_lookup(0:7,basis_length))
        basis_lookup = 0

        do i = 1, nbasis
            bit_pos = mod(i, 8) - 1
            if (bit_pos == -1) bit_pos = 7
            bit_element = (i+7)/8
            bit_lookup(:,i) = (/ bit_pos, bit_element /)
            basis_lookup(bit_pos, bit_element) = i
        end do

    end subroutine init_determinants

    subroutine find_all_determinants()
    
        ! Find all possible Slater determinants that can be formed from the
        ! basis functions.

        use comb_m, only: comb
        use utils, only: get_free_unit

        integer :: i, c(nel), ierr, iunit

        allocate(dets(ndets), stat=ierr)
        
        if (write_determinants) then
            iunit = get_free_unit()
            open(iunit, file=determinant_file, status='unknown')
        end if

        do i = 1, ndets
            c = comb(nbasis, nel, i)
            call init_det(c, dets(i))
            write (iunit,'(1X,i8,4X)',advance='no') i
            call write_det(dets(i)%f, iunit, new_line=.true.)
        end do

        close(iunit, status='keep')

    end subroutine find_all_determinants

    subroutine init_det(occ_list, d)

        ! Initialise a Slater determinant.
        ! In:
        !    occ_list(nel): list of occupied orbitals in the Slater determinant.
        ! Out:
        !    d: variable of type det containing a bit repesentation of the 
        !       occupied orbitals and the total spin.

        integer, intent(in) :: occ_list(nel)
        type(det), intent(out) :: d

        integer :: ierr

        allocate(d%f(basis_length), stat=ierr)

        d%f = encode_det(occ_list)
        d%Ms = det_spin(d%f)
    
    end subroutine init_det

    pure function encode_det(occ_list) result(bit_list)

        ! In:
        !    occ_list(nel): integer list of occupied orbitals in the Slater determinant.
        ! Returns:
        !    bit_list(basis_length): a bit string representation of the occupied
        !        orbitals.   The first element contains the first 8 basis
        !        functions, the second element the next 8 and so on.  A basis
        !        function is occupied if the relevant bit is set.

        integer(i0) :: bit_list(basis_length)
        integer, intent(in) :: occ_list(nel)
        integer :: i, orb, bit_pos, bit_element

        bit_list = 0
        do i = 1, nel
            orb = occ_list(i)
            bit_pos = bit_lookup(1,orb)
            bit_element = bit_lookup(2,orb)
            bit_list(bit_element) = ibset(bit_list(bit_element), bit_pos)
        end do
        
    end function encode_det

    pure function decode_det(f) result(occ_list)

        ! In:
        !    f(basis_length): bit string representation of the Slater
        !        determinant.
        ! Returns:
        !    occ_list(nel): integer list of occupied orbitals in the Slater determinant.

        integer :: occ_list(nel)
        integer(i0), intent(in) :: f(basis_length)
        integer :: i, j, iorb

        iorb = 1
        outer: do i = 1, basis_length
            do j = 0, 7
                if (btest(f(i), j)) then
                    occ_list(iorb) = basis_lookup(j, i)
                    if (iorb == nel) exit outer
                    iorb = iorb + 1
                end if
            end do
        end do outer

    end function decode_det

    subroutine write_det(f, iunit, new_line)

        ! Write out a determinant as a list of occupied orbitals in the 
        ! Slater determinant.
        ! In:
        !    f(basis_length): bit string representation of the Slater
        !        determinant.
        !    iunit (optional): io unit to which the output is written.
        !        Default: 6 (stdout).
        !    new_line (optional): if true, then a new line is written at
        !        the end of the list of occupied orbitals.  Default: no
        !        new line.

        integer(i0), intent(in) :: f(basis_length)
        integer, intent(in), optional :: iunit
        logical, intent(in), optional :: new_line
        integer :: occ_list(nel), io, i

        if (present(iunit)) then
            io = iunit
        else
            io = 6
        end if

        occ_list = decode_det(f)

        write (io,'("(")', advance='no')
        do i = 1, nel-1
            write (io,'(i4,",")', advance='no') occ_list(i)
        end do
        write (io,'(i4,")")', advance='no') occ_list(i)
        if (present(new_line)) then
            if (new_line) write (io,'()')
        end if

    end subroutine write_det

    pure function det_spin(f) result(Ms)

        ! In:
        !    f(basis_length): bit string representation of the Slater
        !        determinant.
        ! Returns:
        !    Ms: total spin of the determinant in units of electron spin (1/2).   

        use bit_utils, only: count_set_bits

        integer :: Ms
        integer(i0), intent(in) :: f(basis_length)
        integer(i0) :: alpha_mask, beta_mask, a, b
        integer :: i

        ! Alpha basis functions are in the even bits.
        alpha_mask = 85 !01010101
        ! Beta basis functions are in the odd bits.
        beta_mask = -86 !10101010
        Ms = 0
        do i = 1, basis_length
            ! Find bit string of all alpha orbitals.
            a = iand(f(i), alpha_mask)
            ! Find bit string of all beta orbitals.
            b = iand(f(i), beta_mask)
            Ms = Ms + count_set_bits(a) - count_set_bits(b)
        end do

    end function det_spin

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
            excitation%nexcit = sum(count_set_bits(ieor(f1,f2)))/2

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
            ! permutations are odd or even.
            shift = nel - excitation%nexcit 

            if (excitation%nexcit <= 2) then

                do i = 1, basis_length
                    if (f1(i) == f2(i)) cycle
                    do j = 0, 7

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

end module determinants

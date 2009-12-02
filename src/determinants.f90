module determinants

! Generation, inspection and manipulation of Slater determinants. 

use const
use system
use basis

implicit none

type det
    ! f(basis_length): bit string representation of the Slater determinant.
    integer(i0), pointer :: f(:)
    ! Ms: total spin of the determinant in units of electron spin (1/2).   
    integer(i0) :: Ms
end type det

type(det), allocatable :: dets(:)

type excit
    integer :: nexcit
    integer :: from_orb(2), to_orb(2)
    ! perm is true if a total odd number of permutations is required to align
    ! the determinants.
    logical :: perm
end type excit

integer :: basis_length
integer :: ndets

integer, allocatable :: bit_lookup(:,:) ! (2, nbasis)
integer, allocatable :: basis_lookup(:,:) ! (8, basis_length)

contains

    subroutine init_determinants()

        use comb_m, only: binom

        integer :: i, bit_pos, bit_element

        ndets = binom(nbasis, nel)

        basis_length = nbasis/8
        if (mod(nbasis,8) /= 0) basis_length = basis_length + 1

        write (6,*) 'Number of electrons: ', nel
        write (6,*) 'Number of basis functions: ', nbasis
        write (6,*) 'Number of determinants: ', ndets
        write (6,*) 'basis_length: ', basis_length

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

        integer :: i, c(nel), ierr
        type(excit) :: excitation

        allocate(dets(ndets), stat=ierr)

        do i = 1, ndets
            c = comb(nbasis, nel, i)
            call init_det(c, dets(i))
        end do

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

        alpha_mask = 85 !01010101
        beta_mask = -86 !10101010
        Ms = 0
        do i = 1, basis_length
            a = iand(f(i), alpha_mask)
            b = iand(f(i), beta_mask)
            Ms = Ms + count_set_bits(a) - count_set_bits(b)
        end do

    end function det_spin

    pure function get_excitation(f1,f2) result(excitation)

        ! 

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
            excitation%nexcit = sum(count_set_bits(ieor(f1,f2)))/2
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

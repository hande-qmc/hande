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
end type

type(det), allocatable :: dets(:)

integer :: basis_length

contains

    subroutine find_all_determinants()
    
        ! Find all possible Slater determinants that can be formed from the
        ! basis functions.

        use comb_m

        integer :: ndets, i, c(nel), ierr

        ndets = binom(nbasis, nel)

        basis_length = nbasis/8
        if (mod(nbasis,8) /= 0) basis_length = basis_length + 1

        write (6,*) 'Number of electrons: ', nel
        write (6,*) 'Number of basis functions: ', nbasis
        write (6,*) 'Number of determinants: ', ndets
        write (6,*) 'basis_length: ', basis_length

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
            bit_pos = mod(orb, 8) - 1
            if (bit_pos == -1) bit_pos = 7
            bit_element = (orb+7)/8
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
        integer :: bit_element, i, bit_pos, iorb

        bit_element = 1
        iorb = 1
        occ_list = 0
        do i = 1, nbasis
            bit_pos = mod(i, 8) - 1
            if (bit_pos == -1) bit_pos = 7
            if (btest(f(bit_element), bit_pos)) then
                occ_list(iorb) = i
                iorb = iorb + 1
            end if
            if (bit_pos == 7) bit_element = bit_element + 1
        end do

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

end module determinants

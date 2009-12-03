module basis

! Basis function information.

use kpoints

implicit none

! Store of information about the (spin) basis functions of the system.
type(kpoint), allocatable :: basis_fns(:) ! (nbasis)

! number of basis functions.  Equal to 2*number of sites as there are
! 2 spin orbitals per site.
integer :: nbasis

! The determinants are stored as a bit string.  Each element of an array is
! a byte (containing 8 bytes) as this was the smallest possible integer type
! that standard fortran can currently handle. (The bit type has just been
! deleted from the forthcoming F2008 standard, so we won't hold out breath...)
! basis_length is the length of the byte array necessary to contain a bit for
! each basis function, i.e. ceiling(nbasis/8).
integer :: basis_length

! A determinant is stored in the array f(nbasis).  A basis function is occupied
! in the determinant if the relevant bit is set.  The relevant bit is given by
! bit_element, the element of the array which contains the bit corresponding to
! the basis function, and bit_position, which contains the position of the bit
! within the given element.  bit_lookup(:,i) gives the (/ bit_position,
! bit_element /) of the i-th basis function.
! Note fortran numbers bits starting from 0.
integer, allocatable :: bit_lookup(:,:) ! (2, nbasis)

! The reverse lookup to bit_lookup.
! basis_lookup(i,j) gives the basis function corresponding to
! the i-th bit in the j-th element of a determinant array.
integer, allocatable :: basis_lookup(:,:) ! (8, basis_length)

contains

    pure function spin_symmetry(i, j) result(spin_match)

        ! In:
        !    i: index of a basis function
        !    j: index of a basis function
        ! Returns:
        !    true if i and j refer to basis functions of the same spin.

        logical :: spin_match
        integer, intent(in) :: i, j

        spin_match = basis_fns(i)%ms == basis_fns(j)%ms

    end function spin_symmetry

end module basis

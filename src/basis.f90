module basis

! Basis function information.

use kpoints

implicit none

type(kpoint), allocatable :: basis_fns(:)
integer :: nbasis
integer :: basis_length

integer, allocatable :: bit_lookup(:,:) ! (2, nbasis)
integer, allocatable :: basis_lookup(:,:) ! (8, basis_length)

contains

    pure function spin_symmetry(i, j) result(spin_match)

        logical :: spin_match
        integer, intent(in) :: i, j

        spin_match = basis_fns(i)%ms == basis_fns(j)%ms

    end function spin_symmetry

end module basis

module basis

! Basis function information.

use kpoints

implicit none

type(kpoint), allocatable :: basis_fns(:)
integer :: nbasis

end module basis

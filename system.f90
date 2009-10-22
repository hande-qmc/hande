module system

use const

implicit none

integer :: ndim ! 1, 2 or 3 dimensions.
integer :: n_sites ! Number of sites in crystal cell.
integer, allocatable :: lattice(:,:)  ! ndim, ndim.  Gives lattice vectors of crystal cell.
real(dp), allocatable :: box_length(:) ! ndim.  Gives lengths of lattice vectors.

integer :: nel = 0 ! # of electrons

end module system

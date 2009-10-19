module hubbard

use kpoints
use system

implicit none

type(kpoint), allocatable :: basis_fns(:)

contains

    subroutine init_basis_fns()
    ! Produce the basis functions.  The number of basis functions is
    ! equal to the number of sites in the crystal cell (ie the number
    ! of k-points used to sample the FBZ of the primitive cell).
    ! From the cell parameters and the "tilt" used (if any) generate
    ! the list of wavevectors and hence the kinetic energy associated
    ! with each basis function (one per wavevector).
    end subroutine init_basis_fns()

end module hubbard

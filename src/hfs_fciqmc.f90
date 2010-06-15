module hellmann_feynmann_sampling

! Module for performing Hellmann--Feynmann sampling in FCIQMC.

use const

implicit none

!--- Input options. ---

! Specifiy the magnitude squared of the l quantum vector which specifies a set
! of symmetry-related orbitals.  The occupation number of this set is then
! sampled using HF-FCIQMC.
integer :: lmag2

!--- Operator parameters. ---

! Bit string mask corresponding to be orbitals selected by lmag2.
integer(i0), allocatable :: lmask(:)

contains

    subroutine init_hellmann_feynmann_sampling()

        ! Initialisation of HF sampling: setup parameters for the operators being
        ! sampled.

        use basis, only: basis_length, set_orb_mask

        integer :: ierr

        allocate(lmask(basis_length), stat=ierr)

        call set_orb_mask(lmag2, lmask)

    end subroutine init_hellmann_feynmann_sampling

end module hellmann_feynmann_sampling

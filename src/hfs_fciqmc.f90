module hellmann_feynmann_sampling

! Module for performing Hellmann--Feynmann sampling in FCIQMC.

use const

use hfs_data

implicit none

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

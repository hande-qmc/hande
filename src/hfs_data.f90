module hfs_data

! Module for storing data specific to Hellmann--Feynmann sampling.

use const, only: p, i0

implicit none

!--- Input options. ---

! Specifiy the magnitude squared of the l quantum vector which specifies a set
! of symmetry-related orbitals.  The occupation number of this set is then
! sampled using HF-FCIQMC.
integer :: lmag2

!--- Operator parameters. ---

! Bit string mask corresponding to be orbitals selected by lmag2.
integer(i0), allocatable :: lmask(:)

!--- HFS-specific variables. ---

real(p) :: hf_shift

end module hfs_data

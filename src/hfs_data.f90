module hfs_data

! Module for storing data specific to Hellmann--Feynman sampling.

use const, only: p, i0, lint

implicit none

!--- Input options. ---

! Specifiy the magnitude squared of the l quantum vector which specifies a set
! of symmetry-related orbitals.  The occupation number of this set is then
! sampled using HF-FCIQMC.
integer :: lmag2 = -1

! Which operator are we sampling?
integer :: hf_operator

!--- Avaiable operators. ---

! Note that not all operators are implemented for all systems.

! Hamiltonian operator.
! Of course, we can obtain this via standard FCIQMC, but it's useful for
! debugging.
integer, parameter :: hamiltonian_operator = 2**0

! Kinetic operator, T.
integer, parameter :: kinetic_operator = 2**1

!--- Operator parameters. ---

! Bit string mask corresponding to be orbitals selected by lmag2.
integer(i0), allocatable :: lmask(:)

!--- HFS-specific variables. ---

! TODO: comment variables.

real(p) :: hf_shift = 0.0_p
real(p) :: proj_hf_O_hpsip
real(p) :: proj_hf_H_hfpsip

real(p) :: D0_hf_population
real(p) :: O00

integer(lint) :: hf_signed_pop

end module hfs_data

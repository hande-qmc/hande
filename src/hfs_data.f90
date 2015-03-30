module hfs_data

! Module for storing data specific to Hellmann--Feynman sampling.

use const, only: p

implicit none

!--- Input options. ---

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

! Double occupancy operator, D.
integer, parameter :: double_occ_operator = 2**2

! Dipole occupancy operator, \mu_i, where the component i depends upon the
! supplied integrals.
integer, parameter :: dipole_operator = 2**3

!--- Operator parameters. ---

!--- HFS-specific variables. ---

! Hellmann--Feynman shift.
real(p) :: hf_shift = 0.0_p

end module hfs_data

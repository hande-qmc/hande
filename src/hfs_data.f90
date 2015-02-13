module hfs_data

! Module for storing data specific to Hellmann--Feynman sampling.

use const, only: p, i0, int_64, dp

implicit none

!--- Input options. ---

! Which operator are we sampling?
integer :: hf_operator

! alpha -> 0 from + or from -?
integer :: alpha0

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

! The Hellmann--Feynmann equivalent of the projected estimator contains several
! terms which must be averaged separately.
! See documentation/theory/hellmann_feynman/hf.tex for details.
! \sum_i <D_0|O|D_i> N_i
real(p) :: proj_hf_O_hpsip
! \sum_i <D_0|H|D_i> \tilde{N}_i
real(p) :: proj_hf_H_hfpsip

! Hellmann-Feynman psip population on the reference determinant, D0.
real(p) :: D0_hf_population

! Value of <D0|O|D0>, where O is the operator we are sampling.
real(p) :: O00

! Signed population of Hellmann--Feynman particles
!     \sum_i sign(N_i) \tilde{N}_i,
! where
!     N_i is the Hamiltonian population on i 
!     \tilde{N}_i is the Hellmann--Feynman population on i 
real(p) :: hf_signed_pop

end module hfs_data

module fciqmc_data

! Data for fciqmc calculations and procedures which manipulate fciqmc and only
! fciqmc data.

use const
implicit none

!--- Input data ---

! number of monte carlo cycles/report loop
integer :: ncycles
! number of report cycles
integer :: nreport

! timestep
real(dp) :: tau = 0.10_dp

! shift
real(dp) :: shift = 0.0_dp

! Array sizes
integer :: walker_length
integer :: spawned_walker_length

! Current number of walkers stored in the main list.
integer :: tot_walkers

!--- Walker data ---

! Walker information: main list.
! a) determinants
integer(i0), allocatable :: walker_dets(:,:) ! (basis_length, walker_length)
! b) walker population
integer, allocatable :: walker_population(:) ! (walker_length)
! c) Diagonal matrix elements, K_ii.  Storing them avoids recalculation.
! K_ii = < D_i | H | D_i > - E_0, where E_0 = <D_0 | H | D_0> and |D_0> is the
! reference determinant.
real(dp), allocatable :: walker_energies(:)

! Walker information: spawned list.
! a) determinants.
integer(i0), allocatable :: spawned_walker_dets(:,:) ! (basis_length, spawned_walker_length)
! b) walker population.
integer, allocatable :: spawned_walker_population(:) ! (spawned_walker_length)
! c) next empty slot in the spawning array.
integer :: spawning_head

end module fciqmc_data

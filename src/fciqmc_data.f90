module fciqmc_data

use const
implicit none

! timestep
real(dp) :: tau

! shift
real(dp) :: shift

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
! a) determinants
integer(i0), allocatable :: spawned_walker_dets(:,:) ! (basis_length, walker_length)
! b) walker population
integer, allocatable :: spawned_walker_population(:) ! (walker_length)

end module fciqmc_data

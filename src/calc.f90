module calc

use const
use parallel, only: blacs_info

implicit none

! Flags for doing exact and/or Lanczos diagonalisation.
logical :: t_exact = .false., t_lanczos = .false.

! Hamiltonian matrix.  Clearly the scaling of the memory demands with system
! size is horrendous.  We only store one symmetry block at a time though.
! Best contained within a module to allow easy access from the Lanczos
! matrix-vector multiplication routines.
real(dp), allocatable :: hamil(:,:) ! (ndets, ndets)

! If true, then the eigenvectors are found during exact diagonalisation as well
! as the eigenvalues.  Doing this is substantially more expensive.
logical :: find_eigenvectors = .false.

! If true then the non-zero elements of the Hamiltonian matrix are written to hamiltonian_file.
logical :: write_hamiltonian = .false.
character(255) :: hamiltonian_file = 'HAMIL'

! BLACS info for diagonalisation
type(blacs_info) :: proc_blacs_info

! Distribution of Hamiltonian matrix across the processors. 
! No distribution.
integer, parameter :: distribute_off = 0
! Block cyclic distribution (see comments in parallel.F90 and the blacs and
! scalapack documentation).  Used for parallel exact diagonalisation.
integer, parameter :: distribute_blocks = 1
! Distribute matrix by columns.  Used for parallel Lanczos.
integer, parameter :: distribute_cols = 2

! Flag which stores which distribution mode is in use.
integer :: distribute = distribute_off

end module calc

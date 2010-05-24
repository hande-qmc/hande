module calc

use const
use parallel, only: blacs_info

implicit none

!--- Calculation info ---

! A calculation type is performed if the relevant bit (defined by the subsequent
! parameters) is set in calc_type.
! This can be easily tested using the doing_calc function.
integer :: calc_type = 0
! Flags for doing exact and/or Lanczos diagonalisation.
integer, parameter :: exact_diag = 2**0
integer, parameter :: lanczos_diag = 2**1
! Use the incredibly simple and naive FCIQMC or the optimised implementation?
integer, parameter :: fciqmc_calc = 2**2
integer, parameter :: simple_fciqmc_calc = 2**3
! Estimate the size of the Hilbert space using Monte Carlo?
integer, parameter :: mc_hilbert_space = 2**4

! Ms of determinants.  If not set, then all possible values of Ms are considered
! in FCI.  FCIQMC assumes ms = 0 if not given in input.
integer :: ms_in = huge(1)

! Symmetry block of determinants.  Ignored for real space formulation.  Refers
! to a wavevector in momentum space formulation.  If not set, then determinants
! of all possible momenta are considered in FCI.  FCIQMC assumes determinants
! with 0 momentum are to be considered if not specified in input.
integer :: sym_in = huge(1)

!--- Info for FCI calculations ---

! Hamiltonian matrix.  Clearly the scaling of the memory demands with system
! size is horrendous.  We only store one symmetry block at a time though.
! Best contained within a module to allow easy access from the Lanczos
! matrix-vector multiplication routines.
real(p), allocatable :: hamil(:,:) ! (ndets, ndets)

! If true, then the eigenvectors are found during exact diagonalisation as well
! as the eigenvalues.  Doing this is substantially more expensive.
logical :: find_eigenvectors = .false.

! If true then the non-zero elements of the Hamiltonian matrix are written to hamiltonian_file.
logical :: write_hamiltonian = .false.
character(255) :: hamiltonian_file = 'HAMIL'

! BLACS info for diagonalisation
type(blacs_info) :: proc_blacs_info

!--- Parallel info for FCI calculations ---

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

!--- Info for stocastic calculations ---

! Seed used to initialise the dSFMT random number generator.
integer :: seed = 7

contains

    function doing_calc(calc_param) result(doing)

        ! In:
        !    calc_param: integer corresponding to a type of calculation, e.g.
        !      lanczos diagonalisation or FCIQMC.
        !      It is possible to test to see if one or more out of a group of
        !      calculation types are being performed by setting calc_param to be
        !      the sum of the group of calculation types.
        ! Returns:
        !    true if the supplied calculation type is specifed in calc_type.
        
        logical :: doing
        integer, intent(in) :: calc_param

        doing = iand(calc_param, calc_type) /= 0

    end function doing_calc

end module calc

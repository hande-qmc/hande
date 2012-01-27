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
! Doing initiator-FCIQMC?
integer, parameter :: initiator_fciqmc = 2**4
! Doing continuous-time FCIQMC?
integer, parameter :: ct_fciqmc_calc = 2**5
! Doing Hellmann--Feynman sampling?
integer, parameter :: hfs_fciqmc_calc = 2**6
! Estimate the size of the Hilbert space using Monte Carlo?
integer, parameter :: mc_hilbert_space = 2**7
! Doing a folded spectrum calculation?
integer, parameter :: folded_spectrum = 2**8
! Doing Density Matrix Monte Carlo?
integer, parameter :: dmqmc_calc = 2**9

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

! If either of the following options are true, then the eigenvectors are found
! during exact diagonalisation as well as the eigenvalues.  Doing this is
! substantially more expensive.

! Print out the ground state wavefunction.
logical :: print_ground_state = .false.

! Analyse the ground state wavefunction.
logical :: analyse_ground_state = .false.

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

!--- Input data: Hilbert space truncation ---

! CI/CIQMC:
! If true, truncate the Slater determinant space such that it contains
! determinants which differ from the reference determinant (e.g. Hartree--Fock
! determinant) by at most truncation_level excitations.
! truncation_level excitations.  A truncation level equal to the number of
! electrons corresponds to the full space
! DMQMC:
! If true, truncate the density matrix space such that it only contains matrix
! elements corresponding to two determinants which differ by at most
! truncation_level excitations.  A truncation level equal to the number of
! electrons corresponds to the full space.
logical :: truncate_space = .false.
integer :: truncation_level

!--- Info for stocastic calculations ---

! Seed used to initialise the dSFMT random number generator.
integer :: seed = 7

!--- Info for DMQMC calculations ---

! For DMQMC, the user may want to calculate many different combinations
! of estimators. The above variable, calc-type, does a similar thing
! for the types of calculation used. dmqmc_calc_type is a variable
! which works in the same way, to allow the user to choose any combination
! of estimators in a general way. The new function doing_dmqmc_calc works
! in exactly the same way to doing_calc.
integer :: dmqmc_calc_type = 0
integer, parameter :: dmqmc_energy = 2**0
integer, parameter :: dmqmc_staggered_magnetisation = 2**1
integer, parameter :: dmqmc_energy_squared = 2**2
integer, parameter :: dmqmc_correlation = 2**3

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

    function doing_dmqmc_calc(calc_param) result(doing)

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

        doing = iand(calc_param, dmqmc_calc_type) /= 0

    end function doing_dmqmc_calc

end module calc

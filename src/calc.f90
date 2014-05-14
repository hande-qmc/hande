module calc

use const
use csr, only: csrpsy_t
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
! Doing continuous-time FCIQMC?
integer, parameter :: ct_fciqmc_calc = 2**4
! Doing Hellmann--Feynman sampling?
integer, parameter :: hfs_fciqmc_calc = 2**5
! Estimate the size of the Hilbert space using Monte Carlo?
integer, parameter :: mc_hilbert_space = 2**6
! Doing a folded spectrum calculation?
integer, parameter :: folded_spectrum = 2**7
! Doing Density Matrix Monte Carlo?
integer, parameter :: dmqmc_calc = 2**8
! Doing Coupled Cluster Monte Carlo?
integer, parameter :: ccmc_calc = 2**9

! Using the initiator approximation in FCIQMC or CCMC?
logical :: initiator_approximation = .false.

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
! For single-node calculations, we can also store the hamiltonian matrix in
! sparse format.  This is typically *far* smaller than hamil and hence it is
! usually better to run on a single node than distribute the Hamiltonian matrix.
! Of course, perhaps one day a kind soul will implement a parallel sparse matrix
! storage format.
type(csrpsy_t) :: hamil_csr
! Use sparse matrix rather than dense matrix?
logical :: use_sparse_hamil = .false.

! If either of the following options are true, then the eigenvectors are found
! during exact diagonalisation as well as the eigenvalues.  Doing this is
! substantially more expensive.

! Number of FCI wavefunctions to print out.
integer :: print_fci_wfn = 0
! ...and file to write them to.
character(255) :: print_fci_wfn_file = 'FCI_WFN'

! Number of FCI wavefunctions to compute properties of.
integer :: analyse_fci_wfn = 0

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
! Flag for using non-blocking communications.
! Default: False.
logical :: non_blocking_comm = .false.

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

! RAS-CI (only in QMC currently)
! The ras1 space can have at most truncation_level holes.
! The ras3 space can have at most ras3_max electrons.
! If all orbitals are either in the frozen core, ras1 or ras3 spaces, then
! ras3_max must equal truncation_level.
! The RAS2 space is analogous to the complete active space (CAS).
! Number of spatial orbitals in RAS1 and RAS3 spaces (counting from lowest
! orbitals up and highest orbitals down, respectively).
integer :: RAS(2) = (/ -1, -1 /)
! Min number of electrons in RAS1 space.
integer :: ras1_min
! Max number of electrons in RAS3 space.
integer :: ras3_max
! Bit masks for showing only orbitals in RAS1 and RAS3 spaces.
integer(i0), allocatable :: ras1(:), ras3(:) ! (basis_length)

!--- Info for stocastic calculations ---

! Seed used to initialise the dSFMT random number generator.
! Default: hash of global UUID and time.
integer :: seed

! --- QMC reference state and trial (importance-sampling) functions ---

! For the Heisenberg model, several different trial functions can be used in the
! energy estimator. Only a single determinant can be used for the Hubbard model.
integer, parameter :: single_basis = 0
integer, parameter :: neel_singlet = 1

! Which trial function are we using? Only relevant to the Heisneberg model.
! trial_function will always be 0 for other models to represent a single determinant.
integer :: trial_function = 0

! For the Heisenberd model, a guiding function may be used,
! |psi_G> = \sum_{i} a_i |psi_i>, so that the new Hamiltonian matrix elements are
! H_ij^new = (a_i*H_ij)/a_j. This is just importance sampling. These functions
! represent the different types of functions whihc may be used.
integer, parameter :: no_guiding = 0
! Note that when we use the Neel singlet state as a guiding function, it must also
! be used as the trial function in calculating the projected energy.
integer, parameter :: neel_singlet_guiding = 1

! If we are not using importance sampling, this is set to 0. Else it is set to one
! of the above values to specify the corresponding guiding function being used.
integer :: guiding_function = 0

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
integer, parameter :: dmqmc_rdm_r2 = 2**4
integer, parameter :: dmqmc_full_r2 = 2**5

! Combine information required for non-blocking report loop quantities
! into one type for convenience.
type nb_rep_t
    ! Array to store report loop estimators such as projected energy
    ! etc.
    ! This array must not be deallocated, copied or inspected in any
    ! way in between report loop communication.
    real(dp), allocatable :: rep_info(:)
    ! Array whose entries will contain:
    ! 1. The total number of spawned walkers in a given report loop
    !    which is to be used for calculating the spawning rate.
    ! 2. The number of walkers spawned from a given processor
    !    to all other processors except the current one, which
    !    is used for calculating the total number of walkers for a give
    !    report loop.
    integer :: nb_spawn(2)
    ! Array of requests used for non blocking communications.
    ! This array must not be deallocated, copied or inspected in any
    ! way in between report loop communication. request(nprocs)
    integer, allocatable :: request(:)
end type nb_rep_t

contains

    subroutine init_calc_defaults()

        ! Initialise calculation defaults which cannot be set at compile-time.

        use iso_c_binding, only: c_loc, c_ptr, c_char, c_int
        use report, only: GLOBAL_UUID
        use hashing, only: MurmurHash2
        use utils, only: fstring_to_carray

        character(len=len(GLOBAL_UUID)+10) :: seed_data
        character(kind=c_char), target :: cseed_data(len(seed_data)+1)
        type(c_ptr) :: cseed_ptr
        integer(c_int) :: n

        call date_and_time(time=seed_data(:10))
        seed_data(11:) = GLOBAL_UUID

        cseed_data = fstring_to_carray(seed_data)
        cseed_ptr = c_loc(cseed_data)
        n = size(cseed_data)-1

        seed = int(MurmurHash2(cseed_ptr, n, 12345_c_int))

    end subroutine init_calc_defaults

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

    subroutine alloc_nb_rep_t(ntypes, report_loop)

        ! Allocate nb_rep_t object.

        ! In:
        !    ntypes:
        !       number of types of walkers sampled (see sampling_size).
        ! Out:
        !    report_loop: initialies object for storing and communicating
        !        report loop quantities for non-blocking communications.

        use parallel, only: nprocs
        use checking, only: check_allocate

        integer, intent(in) :: ntypes
        type(nb_rep_t), intent(out) :: report_loop

        integer :: ierr

        allocate(report_loop%rep_info(ntypes*nprocs+7), stat=ierr)
        call check_allocate('report_loop%rep_info', size(report_loop%rep_info), ierr)

        report_loop%nb_spawn = 0

        allocate(report_loop%request(0:nprocs-1), stat=ierr)
        call check_allocate('report_loop%request', size(report_loop%request), ierr)

    end subroutine alloc_nb_rep_t

    subroutine dealloc_nb_rep_t(rep_comm)

        ! Deallocate nb_rep_t object.

        ! In/Out:
        !    rep_comm: nb_rep_t object to be deallocated.

        use checking, only: check_deallocate

        type(nb_rep_t), intent(inout) :: rep_comm

        integer :: ierr

        deallocate(rep_comm%rep_info, stat=ierr)
        call check_deallocate('rep_comm%rep_info', ierr)

        deallocate(rep_comm%request, stat=ierr)
        call check_deallocate('rep_comm%request', ierr)

    end subroutine dealloc_nb_rep_t

end module calc

module calc

use const
use csr, only: csrp_t
use parallel, only: blacs_info

implicit none

!--- Calculation metadata ---

! This is a copy, initialised in init_calc_defaults, of settings in
! environment_report.  By storing them here we avoid a cascade compilation on
! every compilation when the git hash is updated at compile-time.

type metadata_t
    character(40) :: git_sha1
    character(36) :: uuid
end type metadata_t

type(metadata_t) :: GLOBAL_META

!--- Calculation info ---

! A calculation type is performed if the relevant bit (defined by the subsequent
! parameters) is set in calc_type.
! This can be easily tested using the doing_calc function.

! [review] - AJWT: This global is a bit worrying as it appears that only one type of calc can happen at a time.  Possibly for long-term thought.
! [comment] - Verena: Note that when doing simple FCIQMC, calc_type = simple_fciqmc_calc + fciqmc_calc.
integer :: calc_type = 0
! Flag for doing exact diagonalisation.
integer, parameter :: exact_diag = 2**0
! 2**1 used to be lanczos diagonalisation which has been removed from HANDE.
! Use the incredibly simple and naive FCIQMC or the optimised implementation?
integer, parameter :: fciqmc_calc = 2**2
integer, parameter :: simple_fciqmc_calc = 2**3
! Doing continuous-time FCIQMC?
integer, parameter :: ct_fciqmc_calc = 2**4
! Doing Hellmann--Feynman sampling?
integer, parameter :: hfs_fciqmc_calc = 2**5
! Estimate the size of the Hilbert space using Monte Carlo?
integer, parameter :: mc_hilbert_space = 2**6
! Doing Density Matrix Monte Carlo?
integer, parameter :: dmqmc_calc = 2**8
! Doing Coupled Cluster Monte Carlo?
integer, parameter :: ccmc_calc = 2**9
! Monte Carlo estimate of thermal kinetic energy?
integer, parameter :: mc_canonical_estimates = 2**10
! Redistributing restart info.
integer, parameter :: restart_redistribute = 2**11
! Doing Unitary Coupled Cluster Monte Carlo?
integer, parameter :: uccmc_calc = 2**12
! Doing Trotterized Unitary Coupled Cluster Monte Carlo?
integer, parameter :: trot_uccmc_calc = 2**13
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
integer, parameter :: dmqmc_kinetic_energy = 2**6
integer, parameter :: dmqmc_H0_energy = 2**7
integer, parameter :: dmqmc_potential_energy = 2**8
integer, parameter :: dmqmc_HI_energy = 2**9

!--- global data (to deal with in HANDE 1.1)

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
integer(i0), allocatable :: ras1(:), ras3(:) ! (tot_string_len)

contains

    subroutine init_calc_defaults(git_sha1, uuid)

        ! Initialise calculation defaults which cannot be set at compile-time.

        ! In:
        !    git_sha1: git SHA1 hash of the calculation.  Note that only the first 40
        !       characters (the length of the SHA1 hash) are actually used.
        !    uuid: UUID of the calculation.

        character(*), intent(in) :: git_sha1
        character(36), intent(in) :: uuid

        GLOBAL_META = metadata_t(git_sha1, uuid)

    end subroutine init_calc_defaults

    function gen_seed(uuid) result(randish_seed)

        ! In:
        !    uuid: UUID of the calculation.
        ! Returns:
        !    A random(ish) seed based upon the hash of the time and the calculation UUID.

        use iso_c_binding, only: c_loc, c_ptr, c_char, c_int
        use hashing, only: MurmurHash2
        use utils, only: fstring_to_carray
        use parallel

        integer :: randish_seed
        character(36), intent(in) :: uuid

        character(len=len(uuid)+10) :: seed_data
        character(kind=c_char), target :: cseed_data(len(seed_data)+1)
        type(c_ptr) :: cseed_data_ptr
        integer(c_int) :: n
#ifdef PARALLEL
        integer :: ierr
#endif

        if (parent) then
            call date_and_time(time=seed_data(:10))
            seed_data(11:) = uuid

            cseed_data = fstring_to_carray(seed_data)
            cseed_data_ptr = c_loc(cseed_data)
            n = size(cseed_data)-1 ! Don't hash terminating null character.

            randish_seed = int(MurmurHash2(cseed_data_ptr, n, 12345_c_int))
        end if

#ifdef PARALLEL
        call mpi_bcast(randish_seed, 1, mpi_integer, root, MPI_COMM_WORLD, ierr)
#endif

    end function gen_seed

    function doing_calc(calc_param) result(doing)

        ! In:
        !    calc_param: integer corresponding to a type of calculation, e.g.
        !      exact diagonalisation or FCIQMC.
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
        !      exact diagonalisation or FCIQMC.
        !      It is possible to test to see if one or more out of a group of
        !      calculation types are being performed by setting calc_param to be
        !      the sum of the group of calculation types.
        ! Returns:
        !    true if the supplied calculation type is specifed in calc_type.

        logical :: doing
        integer, intent(in) :: calc_param

        doing = iand(calc_param, dmqmc_calc_type) /= 0

    end function doing_dmqmc_calc

    function get_calculation_string(calc_num) result(calc_name)

        ! Returns string corresponding to name of integer given.

        use const, only: int_32

        integer(int_32), intent(in) :: calc_num
        character(255) :: calc_name

        select case(calc_num)
        case(exact_diag)
            calc_name = "exact diagonalisation"
        case(fciqmc_calc)
            calc_name = "FCIQMC"
        case(simple_fciqmc_calc)
            calc_name = "simple FCIQMC"
        case(ct_fciqmc_calc)
            calc_name = "continuous time FCIQMC"
        case(hfs_fciqmc_calc)
            calc_name = "Hellmann-Feynman FCIQMC"
        case(mc_hilbert_space)
            calc_name = "Monte Carlo Hilbert space estimation"
        case(dmqmc_calc)
            calc_name = "DMQMC"
        case(ccmc_calc)
            calc_name = "CCMC"
        case(mc_canonical_estimates)
            calc_name = "Monte Carlo canonical estimates"
        case(restart_redistribute)
            calc_name = "Restart Redistribution"
        case(uccmc_calc)
            calc_name = "UCCMC"
        case(trot_uccmc_calc)
            calc_name = "Trotterized UCCMC"
        end select

    end function get_calculation_string

end module calc

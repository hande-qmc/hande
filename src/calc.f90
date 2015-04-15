module calc

use const
use csr, only: csrp_t
use dmqmc_data, only: rdm_t
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
! Doing Density Matrix Monte Carlo?
integer, parameter :: dmqmc_calc = 2**8
! Doing Coupled Cluster Monte Carlo?
integer, parameter :: ccmc_calc = 2**9
! Monte Carlo estimate of thermal kinetic energy?
integer, parameter :: mc_canonical_kinetic_energy = 2**10

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

!--- global data (to deal with)

! Symmetry block of determinants.  Ignored for real space formulation.  Refers
! to a wavevector in momentum space formulation.  If not set, then determinants
! of all possible momenta are considered in FCI.  FCIQMC assumes determinants
! with 0 momentum are to be considered if not specified in input.
integer :: sym_in = huge(1)

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
integer(i0), allocatable :: ras1(:), ras3(:) ! (string_len)

contains

    subroutine init_calc_defaults(git_sha1, uuid, seed)

        ! Initialise calculation defaults which cannot be set at compile-time.

        ! In:
        !    git_sha1: git SHA1 hash of the calculation.  Note that only the first 40
        !       characters (the length of the SHA1 hash) are actually used.
        !    uuid: UUID of the calculation.
        ! Out:
        !    seed: seed used to initialise the dSFMT random number generator.

        character(*), intent(in) :: git_sha1
        character(36), intent(in) :: uuid
        integer, intent(out) :: seed

        GLOBAL_META = metadata_t(git_sha1, uuid)
        seed = gen_seed(uuid)

    end subroutine init_calc_defaults

    function gen_seed(uuid) result(randish_seed)

        ! In:
        !    uuid: UUID of the calculation.
        ! Returns:
        !    A random(ish) seed based upon the hash of the time and the calculation UUID.

        use iso_c_binding, only: c_loc, c_ptr, c_char, c_int
        use hashing, only: MurmurHash2
        use utils, only: fstring_to_carray

        integer :: randish_seed
        character(36), intent(in) :: uuid

        character(len=len(uuid)+10) :: seed_data
        character(kind=c_char), target :: cseed_data(len(seed_data)+1)
        type(c_ptr) :: cseed_data_ptr
        integer(c_int) :: n

        call date_and_time(time=seed_data(:10))
        seed_data(11:) = uuid

        cseed_data = fstring_to_carray(seed_data)
        cseed_data_ptr = c_loc(cseed_data)
        n = size(cseed_data)-1 ! Don't hash terminating null character.

        randish_seed = int(MurmurHash2(cseed_data_ptr, n, 12345_c_int))

    end function gen_seed

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

    subroutine init_parallel_t(ntypes, ndata, non_blocking_comm, par_calc, nslots)

        ! Allocate and initialise a parallel_t object.

        ! In:
        !    ntypes: number of types of walkers sampled (see sampling_size).
        !    ndata: the number of additional data elements accumulated over all
        !        processors (nparticles_start_ind-1).
        !    non_blocking_comm: true if using non-blocking communications
        !    nslots: number of slots we divide proc_map into
        ! In/Out:
        !    par_calc: type containing parallel information for calculation
        !        see definitions above.

        use parallel, only: nprocs
        use checking, only: check_allocate
        use qmc_data, only: parallel_t

        integer, intent(in) :: ntypes, ndata
        logical, intent(in) :: non_blocking_comm
        integer, intent(in) :: nslots
        type(parallel_t), intent(inout) :: par_calc

        integer :: i, ierr

        associate(nb=>par_calc%report_comm)
            if (non_blocking_comm) then
                allocate(nb%rep_info(ntypes*nprocs+ndata), stat=ierr)
                call check_allocate('nb%rep_info', size(nb%rep_info), ierr)
                allocate(nb%request(0:nprocs-1), stat=ierr)
                call check_allocate('nb%request', size(nb%request), ierr)
            end if
        end associate

        call init_proc_map_t(nslots, par_calc%load%proc_map)

    end subroutine init_parallel_t

    subroutine init_proc_map_t(nslots, pm)

        ! In:
        !    nslots: number of slots (per processor) we divide proc_map_t into.
        ! In/Out:
        !    pm: proc_map_t object containing a mapping of slot to processor index.

        use checking, only: check_allocate

        use parallel, only: nprocs
        use spawn_data, only: proc_map_t

        integer, intent(in) :: nslots
        type(proc_map_t), intent(out) :: pm

        integer :: i, ierr

        pm%nslots = nslots
        allocate(pm%map(0:nslots*nprocs-1), stat=ierr)
        call check_allocate('pm%map', size(pm%map), ierr)
        forall (i=0:nslots*nprocs-1) pm%map(i) = modulo(i,nprocs)

    end subroutine init_proc_map_t

    subroutine dealloc_parallel_t(non_blocking_comm, par_calc)

        ! Deallocate parallel_t object.

        ! In:
        !    non_blocking_comm: true if using non-blocking communications
        ! In/Out:
        !    par_calc: type containing parallel information for calculation
        !        see definitions above.

        use checking, only: check_deallocate
        use qmc_data, only: parallel_t

        logical, intent(in) :: non_blocking_comm
        type(parallel_t), intent(inout) :: par_calc

        integer :: ierr

        associate(nb=>par_calc%report_comm)
            if (non_blocking_comm) then
                if (allocated(nb%rep_info)) then
                    deallocate(nb%rep_info, stat=ierr)
                    call check_deallocate('nb%rep_info', ierr)
                end if
                if (allocated(nb%request)) then
                    deallocate(nb%request, stat=ierr)
                    call check_deallocate('nb%request', ierr)
                end if
            end if
        end associate

        associate(proc_map=>par_calc%load%proc_map)
            if (allocated(proc_map%map)) then
                deallocate(proc_map%map, stat=ierr)
                call check_deallocate('proc_map%map', ierr)
            end if
        end associate

    end subroutine dealloc_parallel_t

end module calc

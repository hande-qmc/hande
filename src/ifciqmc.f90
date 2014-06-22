module ifciqmc

! Module for methods required for the initiator approximation to FCIQMC.
! This contains only core functionality that is *not* in the standard FCIQMC
! algorithm.  Other specialisations (e.g. annihilation, spawning) are found in
! the relevant source files for those parts of the algorithm.

use const

implicit none

! The complete active space (CAS) is given as (N_cas,N_active), where
! N_cas is the number of electrons in the N_active orbitals.
! The N-N_cas electrons occupy the lowest energy orbitals ("core"
! orbitals) for all determinants within the CAS.
! The 2M-N_core-N_active highest energy orbitals are inactive and are
! not occupied in any determinants within the CAS.
! Thus ANDing a determinant with cas_mask gives the electrons in the
! core or inactive orbitals.  The determinant is only in the CAS if the
! result is identical to the cas_core mask (i.e. all the core orbitals
! are filled and no electrons are in the inactive orbitals).
integer(i0), allocatable :: cas_mask(:), cas_core(:)

contains

    subroutine init_ifciqmc(nel)

        ! Allocate and initialise data required for i-FCIQMC.

        ! In:
        !    nel: number of electrons in system.

        use basis, only: basis_length, bit_lookup, nbasis
        use checking, only: check_allocate
        use fciqmc_data, only: initiator_CAS

        integer, intent(in) :: nel
        integer :: ierr, i, bit_pos, bit_element

        allocate(cas_mask(basis_length), stat=ierr)
        call check_allocate('cas_mask', basis_length, ierr)
        allocate(cas_core(basis_length), stat=ierr)
        call check_allocate('cas_core', basis_length, ierr)

        ! Create a mask which has bits set for all core electrons and a mask
        ! which has bits set for all inactive orbitals.
        cas_mask = 0
        cas_core = 0
        ! Set core obitals.
        do i = 1, nel - initiator_CAS(1)
            bit_pos = bit_lookup(1,i)
            bit_element = bit_lookup(2,i)
            cas_mask = ibset(cas_mask(bit_element), bit_pos)
            cas_core = ibset(cas_core(bit_element), bit_pos)
        end do
        ! Set inactive obitals.
        do i = nel - initiator_CAS(1) + 2*initiator_CAS(2) + 1, nbasis
            bit_pos = bit_lookup(1,i)
            bit_element = bit_lookup(2,i)
            cas_mask = ibset(cas_mask(bit_element), bit_pos)
        end do

    end subroutine init_ifciqmc

    subroutine set_parent_flag(parent_population, f, determ_flag, parent_flag)

        ! Test whether the parent determinant is an initiator.
        !
        ! In:
        !    parent_population: current population of walkers on the parent
        !                       determinant.
        !    f: bit string representation of the parent determinant.
        !    determ_flag: 0 if f is deterministic and 1 otherwise.
        ! Out:
        !    parent_flag: set to 0 if the determinant is an initiator and 1 otherwise.

        use basis, only: basis_length
        use fciqmc_data, only: initiator_population

        real(p), intent(in) :: parent_population
        integer(i0), intent(in) :: f(basis_length)
        integer, intent(in) :: determ_flag
        integer, intent(out) :: parent_flag

        if (abs(parent_population) > initiator_population) then
            ! Has a high enough population to be an initiator.
            parent_flag = 0
        else if (all(iand(f,cas_mask) == cas_core)) then
            ! Is in the complete active space.
            parent_flag = 0
        else if (determ_flag == 0) then
            ! Is a deterministic state.
            parent_flag = 0
        else
            ! Isn't an initiator.
            parent_flag = 1
        end if

    end subroutine set_parent_flag

    subroutine set_parent_flag_dummy(parent_population, f, determ_flag, parent_flag)

        ! A deliberately empty routine.

        ! Use the same interface as set_parent_flag, but don't actually do
        ! anything.  This call *should* then get optimised away during standard
        ! FCIQMC calculations.

        use basis, only: basis_length

        real(p), intent(in) :: parent_population
        integer(i0), intent(in) :: f(basis_length)
        integer, intent(in) :: determ_flag
        integer, intent(out) :: parent_flag
        
        ! Simple statements (optimised away?) which remove any unused/unset
        ! variable warnings for compilation with -Werror.
        parent_flag = f(1)
        parent_flag = parent_population
        parent_flag = determ_flag
        parent_flag = 0

    end subroutine set_parent_flag_dummy

    subroutine end_ifciqmc()

        ! Deallocate data required for i-FCIQMC.

        use checking, only: check_deallocate

        integer :: ierr

        if (allocated(cas_mask)) then
            deallocate(cas_mask, stat=ierr)
            call check_deallocate('cas_mask', ierr)
        end if
        if (allocated(cas_core)) then
            deallocate(cas_core, stat=ierr)
            call check_deallocate('cas_core', ierr)
        end if

    end subroutine end_ifciqmc

end module ifciqmc

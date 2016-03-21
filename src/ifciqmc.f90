module ifciqmc

! Module for methods required for the initiator approximation to FCIQMC.
! This contains only core functionality that is *not* in the standard FCIQMC
! algorithm.  Other specialisations (e.g. annihilation, spawning) are found in
! the relevant source files for those parts of the algorithm.

use const

implicit none

contains

    pure subroutine set_parent_flag_real(parent_population, initiator_pop, determ_flag, parent_flag)

        ! Test whether the parent determinant is an initiator.
        !
        ! In:
        !    parent_population: current population of walkers on the parent
        !                       determinant in each space.
        !    initiator_pop: the population above which a determinant is an initiator.
        !    determ_flag: 0 if determinant is deterministic and 1 otherwise.
        ! Out:
        !    parent_flag: set to 0 if the determinant is an initiator and 1 otherwise.

        real(p), intent(in) :: parent_population(:), initiator_pop
        integer, intent(in) :: determ_flag
        integer, intent(out) :: parent_flag
        integer :: i

        parent_flag = 0_int_p
        do i = 1, size(parent_population)
            if (.not.(abs(parent_population(i)) > initiator_pop .or. determ_flag == 0)) then
                ! Isn't an initiator.
                parent_flag = ibset(parent_flag, i - 1)
            ! Otherwise has a high enough population to be an initiator in this space,
            ! or is a deterministic state (which must always be an initiator).
            end if
        end do
    end subroutine set_parent_flag_real

    pure subroutine set_parent_flag_complex(parent_population, initiator_pop, determ_flag, parent_flag)

        ! Test whether the parent determinant is an initiator.
        !
        ! In:
        !    parent_population: current population of walkers on the parent
        !                       determinant in each space.
        !    initiator_pop: the population above which a determinant is an initiator.
        !    determ_flag: 0 if determinant is deterministic and 1 otherwise.
        ! Out:
        !    parent_flag: set to 0 if the determinant is an initiator and 1 otherwise.

        use const, only: p

        real(p), intent(in) :: parent_population(:), initiator_pop
        integer, intent(in) :: determ_flag
        integer, intent(out) :: parent_flag
        integer :: i

        parent_flag = 0_int_p
        do i = 1, size(parent_population), 2
            if ((.not. ((parent_population(i) ** 2) + (parent_population(i+1) ** 2)) > &
                    initiator_pop ** 2) .or. (.not.(determ_flag == 0))) then
                ! Isn't an initiator.
                parent_flag = ibset(parent_flag, i - 1)
                parent_flag = ibset(parent_flag, i)
            end if
            ! Otherwise has a high enough population to be an initiator in this space,
            ! or is a deterministic state (which must always be an initiator).
        end do
    end subroutine set_parent_flag_complex
end module ifciqmc

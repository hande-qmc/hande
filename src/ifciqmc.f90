module ifciqmc

! Module for methods required for the initiator approximation to FCIQMC.
! This contains only core functionality that is *not* in the standard FCIQMC
! algorithm.  Other specialisations (e.g. annihilation, spawning) are found in
! the relevant source files for those parts of the algorithm.

use const

implicit none

contains

    subroutine set_parent_flag(parent_population, initiator_pop, determ_flag, parent_flag)

        ! Test whether the parent determinant is an initiator.
        !
        ! In:
        !    parent_population: current population of walkers on the parent
        !                       determinant in each space.
        !    initiator_pop: the population above which a determinant is an initiator.
        !    determ_flag: 0 if determinant is deterministic and 1 otherwise.
        ! Out:
        !    parent_flag: set to 0 if the determinant is an initiator and 1 otherwise.

        real(p), intent(in) :: parent_population, initiator_pop
        integer, intent(in) :: determ_flag
        integer, intent(out) :: parent_flag
        integer :: i

        if (abs(parent_population) > initiator_pop) then
            ! Has a high enough population to be an initiator in this space.
            parent_flag = 0
        else if (determ_flag == 0) then
            ! Is a deterministic state (which must always be an initiator).
            parent_flag = 0
        else
            ! Isn't an initiator.
            parent_flag = 1
        end if

    end subroutine set_parent_flag

end module ifciqmc

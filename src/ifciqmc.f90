module ifciqmc

! Module for methods required for the initiator approximation to FCIQMC.
! This contains only core functionality that is *not* in the standard FCIQMC
! algorithm.  Other specialisations (e.g. annihilation, spawning) are found in
! the relevant source files for those parts of the algorithm.

use const

implicit none

contains

    subroutine set_parent_flag(parent_population, initiator_pop, f, determ_flag, parent_flag)

        ! Test whether the parent determinant is an initiator.
        !
        ! In:
        !    parent_population: current population of walkers on the parent
        !                       determinant.
        !    initiator_pop: the population above which a determinant is an initiator.
        !    f: bit string representation of the parent determinant.
        !    determ_flag: 0 if f is deterministic and 1 otherwise.
        ! Out:
        !    parent_flag: set to 0 if the determinant is an initiator and 1 otherwise.

        real(p), intent(in) :: parent_population, initiator_pop
        integer(i0), intent(in) :: f(:)
        integer, intent(in) :: determ_flag
        integer, intent(out) :: parent_flag

        if (abs(parent_population) > initiator_pop) then
            ! Has a high enough population to be an initiator.
            parent_flag = 0
        else if (determ_flag == 0) then
            ! Is a deterministic state (which must always be an initiator).
            parent_flag = 0
        else
            ! Isn't an initiator.
            parent_flag = 1
        end if

    end subroutine set_parent_flag

    subroutine set_parent_flag_dummy(parent_population, initiator_pop, f, determ_flag, parent_flag)

        ! A deliberately empty routine.

        ! Use the same interface as set_parent_flag, but don't actually do
        ! anything.  This call *should* then get optimised away during standard
        ! FCIQMC calculations.

        real(p), intent(in) :: parent_population, initiator_pop
        integer(i0), intent(in) :: f(:)
        integer, intent(in) :: determ_flag
        integer, intent(out) :: parent_flag
        
        ! Simple statements (optimised away?) which remove any unused/unset
        ! variable warnings for compilation with -Werror.
        parent_flag = f(1)
        parent_flag = parent_population
        parent_flag = determ_flag
        parent_flag = 0

    end subroutine set_parent_flag_dummy

end module ifciqmc

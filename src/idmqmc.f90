module idmqmc

! Module for methods required for the initiator approximation to DMQMC.
! This contains only core functionality that is *not* in the standard DMQMC
! algorithm.  Other specialisations (e.g. annihilation, spawning) are found in
! the relevant source files for those parts of the algorithm.

use const

implicit none

contains

    subroutine set_parent_flag_dmqmc(parent_population, initiator_pop, f1, f2, level, parent_flag)

        ! Test whether the parent determinant is an initiator.
        !
        ! In:
        !    parent_population: current population of walkers on the parent
        !                       determinant.
        !    initiator_pop: the population above which a determinant is an initiator.
        !    f1: bit string representation of the parent determinant bra/ket.
        !    f2: bit string representation of the parent determinant bra/ket.
        !    level: excitation level at which to set determinant to be an initiator.
        ! Out:
        !    parent_flag: set to 0 if the determinant is an initiator and 1 otherwise.

        use excitations, only: get_excitation_level
        use determinants, only: det_info_t
        use system, only: sys_t

        real(p), intent(in) :: parent_population, initiator_pop
        integer(i0), intent(in) :: f1(:), f2(:)
        integer, intent(in) :: level
        integer, intent(out) :: parent_flag
        integer :: excitation

        excitation = get_excitation_level(f1, f2)

        if (abs(parent_population) > initiator_pop) then
            ! Has a high enough population to be an initiator.
            parent_flag = 0
        else if (level >= excitation) then
            ! Is on excitation level "level" or below.
            parent_flag = 0
        else
            ! Isn't an initiator.
            parent_flag = 1
        end if

    end subroutine set_parent_flag_dmqmc

end module idmqmc

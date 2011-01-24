module folded_spectrum_utils

use const

use proc_pointers

implicit none

contains

    subroutine fs_spawner(cdet, parent_sign, nspawn, connection)

        type(det_info), intent(in) :: cdet
        integer, intent(in) :: parent_sign
        integer, intent(out) :: nspawn(:)
        type(excit), intent(out) :: connection(:)

    end subroutine fs_spawner

    subroutine fs_stochastic_death(Kii, population, tot_population, ndeath)

        ! Particles will attempt to die with probability
        !  p_d = tau*M_ii
        ! where tau is the timestep and M_ii is the appropriate diagonal
        ! matrix element.
        ! For folded spectrum FCIQMC, M_ii = K_ii - S where S is the shift and K_ii is
        !  K_ii =  (< D_i | H | D_i > - E_0 - E_offset)^2.
        ! We store < D_i | H | D_i > - E_0 - E_offset, so just need to square it
        ! before passing it to the standard death routine.

        ! In:
        !    Kii: < D_i | H | D_i > - E_0, where D_i is the determinant on
        !         which the particles reside.
        ! In/Out:
        !    population: number of particles on determinant D_i.
        !    tot_population: total number of particles.
        ! Out:
        !    ndeath: running total of number of particles died/cloned.

        use death, only: stochastic_death

        real(p), intent(in) :: Kii
        integer, intent(inout) :: population, tot_population
        integer, intent(out) :: ndeath

        call stochastic_death(Kii**2, population, tot_population, ndeath)

    end subroutine fs_stochastic_death

end module folded_spectrum_utils

module death

! Module for performing the death part of the FCIQMC algorithm.

use const
use fciqmc_data
implicit none

contains

    subroutine stochastic_death(Mii, population)

        ! Particles will attempt to die with probability
        !  p_d = tau*M_ii
        ! where tau is the timestep and M_ii is the appropriate diagonal
        ! matrix element.
        ! For FCIQMC M_ii = K_ii - S where S is the shift and K_ii is
        !  K_ii =  < D_i | H | D_i > - E_0.
        ! For Hellmann--Feynmann sampling M_ii = < D_i | O | D_i >, where O is
        ! the operator being sampled.

        ! In:

        use dSFMT_interface, only: genrand_real2

        real(p), intent(in) :: Mii
        integer, intent(inout) :: population

        real(p) :: pd
        real(dp) :: r
        integer :: kill

        ! Optimisation: the number of particles on a given determinant can die
        ! stochastically...
        ! This amounts to multplying p_d by the population.  int(p_d) is thus
        ! the number that definitely die and the fractional part of p_d is the
        ! probability of an additional death.

        pd = tau*Mii

        ! This will be the same for all particles on the determinant, so we can
        ! attempt all deaths in one shot.
        pd = pd*abs(population)

        ! Number that definitely die...
        kill = int(pd)

        ! In addition, stochastic death (bad luck! ;-))
        pd = pd - kill ! Remaining chance...
        r = genrand_real2()
        if (abs(pd) > r) then
            if (pd > 0.0_p) then
                ! die die die!
                kill = kill + 1
            else
                ! clone clone clone! doesn't quite have the same ring to it.
                kill = kill - 1
            end if
        end if

        ! Find the new number of particles.
        ! If kill is positive then particles are killed (i.e. the population
        ! should move towards zero).
        ! Also update the total number of particles on the processor.
        if (population < 0) then
            population = population + kill
        else
            population = population - kill
        end if
        nparticles = nparticles - kill

    end subroutine stochastic_death

end module death

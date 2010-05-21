module death

! Module for performing the death part of the FCIQMC algorithm.

use const
use fciqmc_data
implicit none

contains

    subroutine stochastic_death(cpos)

        ! Particles will attempt to die with probability
        !  p_d = tau(K_ii - S)
        ! where tau is the timestep, S is the shift and K_ii is
        !  K_ii =  < D_i | H | D_i > - E_0
        ! In:
        !    cpos: current position within the main walkers array.

        use dSFMT_interface, only: genrand_real2

        integer, intent(in) :: cpos
        real(p) :: pd
        real(dp) :: r
        integer :: current_pop, new_pop, kill

        ! Optimisation: the number of particles on a given determinant can die
        ! stochastically...

        pd = tau*(walker_energies(cpos)-shift)

        current_pop = walker_population(cpos)

        ! This will be the same for all particles on the determinant, so we can
        ! attempt all deaths in one shot.
        pd = pd*abs(current_pop)

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

        ! Find the new total of particles.
        ! If kill is positive then particles are killed (i.e. the population
        ! should move towards zero).
        ! Also update the number of particles on the processor.
        if (current_pop < 0) then
            new_pop = current_pop + kill
        else
            new_pop = current_pop - kill
        end if
        nparticles = nparticles - abs(current_pop) + abs(new_pop)

        walker_population(cpos) = new_pop

    end subroutine stochastic_death

end module death

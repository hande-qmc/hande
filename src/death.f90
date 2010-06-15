module death

! Module for performing the death part of the FCIQMC algorithm.

use const
use fciqmc_data
implicit none

contains

    subroutine stochastic_death(Kii, population, tot_population)

        ! Particles will attempt to die with probability
        !  p_d = tau*M_ii
        ! where tau is the timestep and M_ii is the appropriate diagonal
        ! matrix element.
        ! For FCIQMC M_ii = K_ii - S where S is the shift and K_ii is
        !  K_ii =  < D_i | H | D_i > - E_0.

        ! In:
        !    Kii: < D_i | H | D_i > - E_0, where D_i is the determinant on
        !         which the particles reside.
        ! In/Out:
        !    population: number of particles on determinant D_i.
        !    tot_population: total number of particles.
        
        ! Note that population and tot_population refer to a single 'type' of
        ! population, i.e. either a set of Hamiltonian walkers or a set of
        ! Hellmann--Feynmann walkers.

        use dSFMT_interface, only: genrand_real2

        real(p), intent(in) :: Kii
        integer, intent(inout) :: population, tot_population

        real(p) :: pd
        real(dp) :: r
        integer :: kill, old_population

        ! Optimisation: the number of particles on a given determinant can die
        ! stochastically...
        ! This amounts to multplying p_d by the population.  int(p_d) is thus
        ! the number that definitely die and the fractional part of p_d is the
        ! probability of an additional death.

        pd = tau*(Kii-shift)

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
        old_population = population
        if (population < 0) then
            population = population + kill
        else
            population = population - kill
        end if
        tot_population = tot_population - abs(old_population) + abs(population)

    end subroutine stochastic_death

end module death

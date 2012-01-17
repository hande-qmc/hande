module death

! Module for performing the death part of the FCIQMC algorithm.

use const
use fciqmc_data
implicit none

contains

    subroutine stochastic_death(Kii, population, tot_population, ndeath)

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
        ! Out:
        !    ndeath: running total of number of particles died/cloned.
        
        ! Note that population and tot_population refer to a single 'type' of
        ! population, i.e. either a set of Hamiltonian walkers or a set of
        ! Hellmann--Feynman walkers.

        use dSFMT_interface, only: genrand_real2

        real(p), intent(in) :: Kii
        integer, intent(inout) :: population, ndeath
        integer(lint), intent(inout) :: tot_population

        real(p) :: pd
        real(dp) :: r
        integer :: kill, old_population

        ! Optimisation: the number of particles on a given determinant can die
        ! stochastically...
        ! This amounts to multplying p_d by the population.  int(p_d) is thus
        ! the number that definitely die and the fractional part of p_d is the
        ! probability of an additional death.
        
        ! dmqmc_factor below is set to 1.0 when not performing a DMQMC calculation, and so
        ! can be ignored in these cases.
        ! When performing dmqmc calculations, dmqmc_factor = 2.0. This factor is included
        ! here because in DMQMC calculations, instead of spawning from one end with
        ! the full probability, we spawn from two different ends with half probability each.
        ! Hence, tau is set to tau/2 in DMQMC calculations, so that an extra factor is not
        ! required in every spawning routine. In this death step however, we use kii, which
        ! has a factor of 1/2 included for convenience already, for conveniece elsewhere.
        ! Hence we have to multiply by an extra factor of 2 to account for the extra 1/2 in tau.

        pd = tau*(Kii-shift)*dmqmc_factor

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
        ndeath = ndeath + abs(kill)

    end subroutine stochastic_death

    subroutine stochastic_hf_cloning(Oii, hamiltonian_pop, hf_pop, tot_hf_pop)

        ! Clone Hellmann--Feynman particles from Hamiltonian particles.
        ! HF particles are created from Hamiltonian particles on the same
        ! determinant with probability 
        !   tau(O_ii - \tilde{S})
        ! where
        !   O_ii = < D_i | O | D_i >
        !   \tilde{S} is the HF shift.

        ! In:
        !    Oii: < D_i | O | D_i > (stored in the appropriate element of
        !        walker_enegies).
        !    hamiltonian_pop: number of Hamiltonian particles on determinant
        !        D_i.
        ! In/Out:
        !    hf_pop: number of Hellmann--Feynman particles on determinant D_i.
        !    tot_hf_pop: total number of Hellmann--Feynman particles.

        use dSFMT_interface, only: genrand_real2

        use hfs_data, only: hf_shift

        real(p), intent(in) :: Oii
        integer, intent(in) :: hamiltonian_pop
        integer, intent(inout) :: hf_pop
        integer(lint), intent(inout) :: tot_hf_pop

        real(p) :: pd, matel
        real(dp) :: r
        integer :: clone, old_pop

        ! Attempt to clone with a Hellmann--Feynman walker with probability
        !  p = tau*(Oii-shift)
        ! This will be the same for all particles on the determinant, so we can
        ! attempt all cloning in one shot.
        matel = Oii-hf_shift
        pd = tau*abs(matel*hamiltonian_pop)

        ! Number that are definitely cloned...
        clone = int(pd)

        ! In addition, stochastic cloning.
        pd = pd - clone
        r = genrand_real2()
        if (pd > r) clone = clone + 1

        ! Hellmann--Feynman offsping have the same sign as the Hamiltonian
        ! parents if Oii-hf_shift is negative, otherwise they have the opposite sign.
        old_pop = hf_pop
        if (matel > 0.0_p) then
            hf_pop = hf_pop - sign(clone, hamiltonian_pop)
        else
            hf_pop = hf_pop + sign(clone, hamiltonian_pop)
        end if
        tot_hf_pop = tot_hf_pop - abs(old_pop) + abs(hf_pop)
        if (abs(clone)>1) write (17,*) 'CLONING', clone, old_pop, hf_pop, hamiltonian_pop, Oii, Oii-hf_shift

    end subroutine stochastic_hf_cloning

end module death

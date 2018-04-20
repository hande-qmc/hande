module death

! Module for performing the death part of the FCIQMC algorithm.

use const
use qmc_io
implicit none

contains

    subroutine stochastic_death(rng, sys, qs, dfock, Kii, loc_shift, proj_energy, logging_info, population, &
                                tot_population, ndeath)

        ! Particles will attempt to die with probability
        !  p_d = tau*M_ii
        ! where tau is the timestep and M_ii is the appropriate diagonal
        ! matrix element.
        ! For FCIQMC M_ii = K_ii - S where S is the shift and K_ii is
        !  K_ii =  < D_i | H | D_i > - E_0.

        ! In:
        !    sys: the system
        !    qs: qmc_state_t object. tau, propagator and dmqmc_factor are used.
        !    dfock: \sum_i (f_i - f^0_i), where f_i (f^0_i) is the Fock eigenvalue of the i-th orbital occupied in D_i (D_0).
        !    Kii: < D_i | H | D_i > - E_0, where D_i is the determinant on
        !         which the particles reside.
        !    loc_shift: The value of the shift to be used in the death step.
        !    proj_energy: The current estimate of the energy.
        ! In/Out:
        !    rng: random number generator.
        !    population: (unscaled) number of particles on determinant D_i.
        !    tot_population: total number of particles.
        ! Out:
        !    ndeath: running total of (unscaled) number of particles died/cloned.

        ! Note:

        ! * population and tot_population refer to a single 'type' of
        !   population, e.g. either a set of Hamiltonian walkers or a set of
        !   Hellmann--Feynman walkers.
        ! * the population and ndeath should be unscaled (ie not divided by
        !   pop_real_factor) so to avoid a scaling and unscaling step.

        use dSFMT_interface, only: dSFMT_t, get_rand_close_open
        use qmc_data, only: qmc_state_t
        use system, only: sys_t
        use spawning, only: calc_qn_weighting
        use const, only: debug
        use logging, only: write_logging_death, logging_t

        type(sys_t), intent(in) :: sys
        real(p), intent(in) :: Kii, dfock
        type(qmc_state_t), intent(in) :: qs
        type(dSFMT_t), intent(inout) :: rng
        real(p), intent(in) :: loc_shift, proj_energy
        integer(int_p), intent(inout) :: population, ndeath
        real(dp), intent(inout) :: tot_population
        type(logging_t), intent(in) :: logging_info

        real(p) :: pd, pd_saved
        real(dp) :: r
        integer(int_p) :: kill, old_population
        real(p) :: weight

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

        weight = calc_qn_weighting(qs%propagator, dfock)
        pd = qs%tau*((Kii-proj_energy)*weight+(proj_energy-loc_shift))*qs%dmqmc_factor

        pd_saved = pd

        ! This will be the same for all particles on the determinant, so we can
        ! attempt all deaths in one shot.
        pd = pd*abs(population)

        ! Number that definitely die...
        kill = int(pd, int_p)

        ! In addition, stochastic death (bad luck! ;-))
        pd = pd - kill ! Remaining chance...
        r = get_rand_close_open(rng)
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
        if (debug) call write_logging_death(logging_info, Kii, proj_energy, loc_shift, weight, kill, pd_saved, &
                                        real(old_population,p)/qs%psip_list%pop_real_factor, &
                                        real(population,p)/qs%psip_list%pop_real_factor)
        tot_population = tot_population + &
            real(abs(population) - abs(old_population),p)/qs%psip_list%pop_real_factor
        ndeath = ndeath + abs(kill)

    end subroutine stochastic_death

    subroutine stochastic_hf_cloning(rng, tau, Oii, hf_shift, hamiltonian_pop, ncloned)

        ! Clone Hellmann--Feynman particles from Hamiltonian particles.
        ! HF particles are created from Hamiltonian particles on the same
        ! determinant with probability
        !   tau(O_ii - \tilde{S})
        ! where
        !   O_ii = < D_i | O | D_i >
        !   \tilde{S} is the HF shift.

        ! In/Out:
        !    rng: random number generator.
        ! In:
        !    tau: timestep being used.
        !    Oii: < D_i | O | D_i > (stored in the appropriate element of
        !        walker_enegies).
        !    hf_shift: the Hellmann--Feynman shift, \tilde{S}.
        !    hamiltonian_pop: number of Hamiltonian particles on determinant
        !        D_i.
        ! Out:
        !    ncloned: number of Hellmann--Feynman particles created on determinant D_i.

        use dSFMT_interface, only: dSFMT_t, get_rand_close_open

        type(dSFMT_t), intent(inout) :: rng
        real(p), intent(in) :: tau, Oii, hf_shift
        integer(int_p), intent(in) :: hamiltonian_pop
        integer(int_p), intent(out) :: ncloned

        real(p) :: pd, matel
        real(dp) :: r

        ! Attempt to clone with a Hellmann--Feynman walker with probability
        !  p = tau*(Oii-shift)
        ! This will be the same for all particles on the determinant, so we can
        ! attempt all cloning in one shot.
        matel = Oii-hf_shift
        pd = tau*abs(matel*hamiltonian_pop)

        ! Number that are definitely cloned...
        ncloned = int(pd, int_p)

        ! In addition, stochastic cloning.
        pd = pd - ncloned
        r = get_rand_close_open(rng)
        if (pd > r) ncloned = ncloned + 1

        ! Hellmann--Feynman offsping have the same sign as the Hamiltonian
        ! parents if Oii-hf_shift is negative, otherwise they have the opposite sign.
        if (matel > 0.0_p) then
            ncloned = -sign(ncloned, hamiltonian_pop)
        else
            ncloned = +sign(ncloned, hamiltonian_pop)
        end if

    end subroutine stochastic_hf_cloning

end module death

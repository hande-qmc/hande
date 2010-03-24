module energy_evaluation

! This module contains procedure for evaluating and estimating the energy of
! a system based upon the population dynamics of an FCIQMC calculation.

use const

implicit none

contains

    subroutine update_shift(nparticles_old, nparticles,nupdate_steps)

        ! Update the shift according to:
        !  shift(beta) = shift(beta-A*tau) - xi*log(N_w(tau)/N_w(beta-A*tau))/(A*tau)
        ! where
        !  * shift(beta) is the shift at imaginary time beta;
        !  * A*tau is the amount of imaginary time between shift-updates (=# of
        !    Monte Carlo cycles between updating the shift);
        !  * xi is a damping factor (0.05-0.10 is appropriate) to prevent large fluctations;
        !  * N_w(beta) is the total number of particles at imaginary time beta.
        ! The running average of the shift is also updated.
        ! In:
        !    nparticles_old: N_w(beta-A*tau).
        !    nparticles: N_w(beta).

        use fciqmc_data, only: shift, tau, shift_damping, av_shift

        integer, intent(in) :: nparticles_old, nparticles, nupdate_steps

        shift = shift - log(real(nparticles,8)/nparticles_old)*shift_damping/(tau*nupdate_steps)
        av_shift = av_shift + shift

    end subroutine update_shift

    subroutine update_proj_energy_hub_k(idet, inst_proj_energy)

        ! Add the contribution of the current determinant to the projected
        ! energy.
        ! The correlation energy given by the projected energy is:
        !   \sum_{i \neq 0} <D_i|H|D_0> N_i/N_0
        ! where N_i is the population on the i-th determinant, D_i,
        ! and 0 refers to the reference determinant.
        ! During a MC cycle we store
        !   \sum_{i \neq 0} <D_i|H|D_0> N_i
        ! If the current determinant is the reference determinant, then
        ! N_0 is stored as D0_population.  This makes normalisation very
        ! efficient.
        ! This procedure is for the Hubbard model in momentum space only.
        ! In:
        !    idet: index of current determinant in the main walker list.
        ! In/Out:
        !    inst_proj_energy: running total of the \sum_{i \neq 0} <D_i|H|D_0> N_i.
        !    This is updated if D_i is connected to D_0 (and isn't D_0).

        use fciqmc_data, only: walker_dets, walker_population, f0, D0_population, proj_energy
        use excitations, only: excit, get_excitation
        use hamiltonian, only: slater_condon2_hub_k

        integer, intent(in) :: idet
        real(p), intent(inout) :: inst_proj_energy
        type(excit) :: excitation
        real(p) :: hmatel

        excitation = get_excitation(walker_dets(:,idet), f0)

        if (excitation%nexcit == 0) then
            ! Have reference determinant.
            D0_population = walker_population(idet)
        else if (excitation%nexcit == 2) then
            ! Have a determinant connected to the reference determinant: add to 
            ! projected energy.
            hmatel = slater_condon2_hub_k(excitation%from_orb(1), excitation%from_orb(2), &
                                       & excitation%to_orb(1), excitation%to_orb(2),excitation%perm)
            inst_proj_energy = inst_proj_energy + hmatel*walker_population(idet)
        end if

    end subroutine update_proj_energy_hub_k

    subroutine update_proj_energy_hub_real(idet, inst_proj_energy)

        ! Add the contribution of the current determinant to the projected
        ! energy.
        ! The correlation energy given by the projected energy is:
        !   \sum_{i \neq 0} <D_i|H|D_0> N_i/N_0
        ! where N_i is the population on the i-th determinant, D_i,
        ! and 0 refers to the reference determinant.
        ! During a MC cycle we store
        !   \sum_{i \neq 0} <D_i|H|D_0> N_i
        ! If the current determinant is the reference determinant, then
        ! N_0 is stored as D0_population.  This makes normalisation very
        ! efficient.
        ! This procedure is for the Hubbard model in real space only.
        ! In:
        !    idet: index of current determinant in the main walker list.
        ! In/Out:
        !    inst_proj_energy: running total of the \sum_{i \neq 0} <D_i|H|D_0> N_i.
        !    This is updated if D_i is connected to D_0 (and isn't D_0).

        use fciqmc_data, only: walker_dets, walker_population, f0, D0_population, proj_energy
        use excitations, only: excit, get_excitation
        use hamiltonian, only: slater_condon1_hub_real

        integer, intent(in) :: idet
        real(p), intent(inout) :: inst_proj_energy
        type(excit) :: excitation
        real(p) :: hmatel

        excitation = get_excitation(walker_dets(:,idet), f0)

        if (excitation%nexcit == 0) then
            ! Have reference determinant.
            D0_population = walker_population(idet)
        else if (excitation%nexcit == 1) then
            ! Have a determinant connected to the reference determinant: add to 
            ! projected energy.
            hmatel = slater_condon1_hub_real(excitation%from_orb(1), excitation%to_orb(1), excitation%perm)
            inst_proj_energy = inst_proj_energy + hmatel*walker_population(idet)
        end if

    end subroutine update_proj_energy_hub_real

end module energy_evaluation

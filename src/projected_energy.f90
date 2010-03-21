module projected_energy

use const

use fciqmc_data, only: walker_dets, walker_population, f0, D0_population, proj_energy
use excitations, only: excit, get_excitation
implicit none

contains

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

        use hamiltonian, only: slater_condon2_hub_k

        integer, intent(in) :: idet
        real(dp), intent(inout) :: inst_proj_energy
        type(excit) :: excitation
        real(dp) :: hmatel

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

        use hamiltonian, only: slater_condon1_hub_real

        integer, intent(in) :: idet
        real(dp), intent(inout) :: inst_proj_energy
        type(excit) :: excitation
        real(dp) :: hmatel

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

end module projected_energy

module projected_energy

use const

use fciqmc_data, only: walker_dets, walker_population, f0, D0_population, proj_energy
use excitations, only: excit, get_excitation
implicit none

contains

    subroutine update_proj_energy_hub_k(idet)

        use hamiltonian, only: slater_condon2_hub_k
        integer, intent(in) :: idet
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
            proj_energy = proj_energy + hmatel*walker_population(idet)
        end if

    end subroutine update_proj_energy_hub_k

end module projected_energy

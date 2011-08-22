module heisenberg_estimators

use const

implicit none

contains

subroutine update_proj_energy_heisenberg_basic(idet)

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
        ! This procedure is for the Heisenberg model only
        ! In:
        !    idet: index of current determinant in the main walker list.

        use fciqmc_data, only: walker_dets, walker_population, f0, D0_population, &
                               proj_energy, calculate_magnetisation
        use excitations, only: excit, get_excitation
        use basis, only: bit_lookup
        use system, only: J_coupling
        use hubbard_real, only: connected_orbs

        integer, intent(in) :: idet
        type(excit) :: excitation
        integer :: bit_position, bit_element

        excitation = get_excitation(walker_dets(:,idet), f0)

        if (excitation%nexcit == 0) then
            ! Have reference determinant.
            D0_population = D0_population + walker_population(1,idet)
        else if (excitation%nexcit == 1) then
            ! Have a determinant connected to the reference determinant: add to 
            ! projected energy.
            
            bit_position = bit_lookup(1,excitation%from_orb(1))
            bit_element = bit_lookup(2,excitation%from_orb(1))
            
            if (btest(connected_orbs(bit_element, excitation%to_orb(1)), bit_position)) &
                     proj_energy = proj_energy - 2.0_p*J_coupling*walker_population(1,idet)
        end if
        
        if (calculate_magnetisation) call update_magnetisation_heisenberg(idet)

    end subroutine update_proj_energy_heisenberg_basic
    
    subroutine update_proj_energy_heisenberg_positive(idet)

        use fciqmc_data, only: walker_dets, walker_population, walker_energies, &
                               proj_energy, calculate_magnetisation, D0_population
        use system, only: J_coupling, ndim, nsites

        integer, intent(in) :: idet
        
        proj_energy = proj_energy + (J_coupling*ndim*nsites+2*walker_energies(1,idet))* &
                                     walker_population(1,idet)
                                     
        D0_population = D0_population + walker_population(1,idet)
        
        if (calculate_magnetisation) call update_magnetisation_heisenberg(idet)

    end subroutine update_proj_energy_heisenberg_positive
    
    subroutine update_proj_energy_heisenberg_neel_singlet(idet)
    
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
        ! This procedure is for the Heisenberg model only
        ! In:
        !    idet: index of current determinant in the main walker list.

        use fciqmc_data, only: walker_dets, walker_population, walker_energies, &
                               walker_reference_data, calculate_magnetisation, &
                               proj_energy, neel_singlet_amp, D0_population
        use system, only: nsites, ndim, J_coupling

        integer, intent(in) :: idet
        integer :: i, n, ipos, lattice_1_up, lattice_2_up
        
        n = walker_reference_data(1,idet)
        lattice_1_up = walker_reference_data(2,idet)
        
        ! Deduce the number of 0-1 bonds, where the 1's are on the
        ! second sublattice
        lattice_2_up = ((ndim*nsites) + nint(walker_energies(1,idet)/J_coupling))/2 - lattice_1_up
        
        proj_energy = proj_energy + (neel_singlet_amp(n+1) * walker_energies(1,idet) * walker_population(1,idet))
        if (n > 0) proj_energy = proj_energy - (2 * J_coupling * lattice_1_up * &
                                 walker_population(1,idet) * neel_singlet_amp(n))
        if (n < (nsites/2)) proj_energy = proj_energy - (2 * J_coupling * lattice_2_up * &
                            walker_population(1,idet) * neel_singlet_amp(n+2))
        
        D0_population = D0_population + walker_population(1,idet)*neel_singlet_amp(n+1)
        
        if (calculate_magnetisation) call update_magnetisation_heisenberg(idet)

    end subroutine update_proj_energy_heisenberg_neel_singlet
    
    pure function neel_singlet_data(idet) result(spin_config_data)
    
        ! This subroutine calculates the number of up spins on one of the
        ! sublattices so that it can be stored, along with the total number
        ! of spins which are up on the sublattice. These are later used in the
        ! update_proj_energy_heisenberg_neel_singlet subroutine.
    
        use basis, only: basis_length, basis_lookup
        use bit_utils, only: count_set_bits
        use determinants, only: lattice_mask
        use fciqmc_data, only: walker_dets
        use hubbard_real, only: connected_orbs
        use system, only: nsites
    
        integer, intent(in) :: idet
        integer(i0) :: f_mask(basis_length), f_not(basis_length), g(basis_length)
        integer :: lattice_1_up, n, i, ipos, basis_find
        integer :: spin_config_data(2)
    
        f_mask = iand(walker_dets(:,idet), lattice_mask)
        n = sum(count_set_bits(f_mask))
        if (n > nsites/2) n = nsites/2 - n
    
        lattice_1_up = 0
        f_not = not(walker_dets(:,idet))
        do i = 1, basis_length
            do ipos = 0, i0_end
                if (btest(f_mask(i), ipos)) then
                    basis_find = basis_lookup(ipos, i)
                    g = iand(f_not, connected_orbs(:,basis_find))
                    lattice_1_up = lattice_1_up + sum(count_set_bits(g))
                end if
            end do
        end do
        
        spin_config_data(1) = n
        spin_config_data(2) = lattice_1_up
    
    end function neel_singlet_data
    
    subroutine update_magnetisation_heisenberg(idet)

        ! Add the contribution of the current determinant to the projected
        ! magnetisation. This only gets run when calculate_magnetisation
        ! is true, and only in the Heisenberg model.
        ! The staggered magnetisation squared in the z direction is:
        !   \sum_{i \neq 0} <D_i|M_z^2|D_i> N_i^2/Magnitude^2
        ! where N_i is the population on the i-th basis function,
        ! and Magnitude^2 is the sum of the squares of the components
        ! of the wavefunction, ie, the sum of the squares of the psips.
        ! During a MC cycle we store
        !   \sum_{i \neq 0} <D_i|M_z^2|D_i> N_i^2 and Magnitude^2
        ! This procedure is for the Heisenberg model only
        ! In:
        !    idet: index of current determinant in the main walker list.
        
        use fciqmc_data, only: walker_dets, population_squared, &
                               walker_population, average_magnetisation
        use basis, only: basis_length
        use determinants, only: lattice_mask
        use system, only: nel
        use bit_utils, only: count_set_bits

        integer, intent(in) :: idet
        integer(i0) :: f_mask(basis_length)
        integer :: sublattice1_up_spins

        f_mask = iand(walker_dets(:,idet), lattice_mask)
        sublattice1_up_spins = sum(count_set_bits(f_mask))
        average_magnetisation = average_magnetisation + &
             (walker_population(1,idet)**2)*((4*sublattice1_up_spins - 2*nel)**2)
        population_squared = population_squared + walker_population(1,idet)**2

    end subroutine update_magnetisation_heisenberg

end module heisenberg_estimators

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
    
        ! Add the contribution of the current basis fucntion to the
        ! projected energy.
        ! This uses the trial function
        ! |psi> = \sum_{i} |D_i>
        ! Meaning that every single basis function has a weight of one in
        ! the sum. A unitary transformation is applied to h when using this
        ! estimator, so that all components of the true ground state
        ! are positive, and hence we get a good overlap.
        ! This procedure is for the Heisenberg model only.
        ! In:
        !    idet: index of current determinant in the main walker list.

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
    
        ! Add the contribution of the current basis fucntion to the
        ! projected energy.
        ! This uses the Neel singlet state as a trial function.
        ! This function is an integral over the Neel states in all directions
        ! (giving equal weight to each). Hence it is rotationally invariant
        ! and hence is an S=0 eigenstate. The ground state is known to be an
        ! S = 0 eigenstate also (ie of the total spin squared operator). This
        ! makes this state a suitable trial function. For details on the state, see
        ! K. Runge, Phys. Rev. B 45, 7229 (1992).
        ! This procedure is for the Heisenberg model only.
        ! In:
        !    idet: index of current determinant in the main walker list.

        use fciqmc_data, only: walker_dets, walker_population, walker_energies, &
                               walker_reference_data, calculate_magnetisation, &
                               proj_energy, neel_singlet_amp, D0_population, &
                               importance_sampling
        use system, only: nsites, ndim, J_coupling

        integer, intent(in) :: idet
        integer :: i, n, ipos, lattice_1_up, lattice_2_up
        real(dp) :: importance_sampling_factor = 1.0
        
        n = walker_reference_data(1,idet)
        lattice_1_up = walker_reference_data(2,idet)
        
        ! If importance sampling is applied then the psip amplitudes, n_i,
        ! will represent the quantities
        ! f_i = c_i^T * c_i
        ! where c_i the amplitude of |D_i> in the true ground state and
        ! c_i^T is the amplitude of |D_i> in the trial ground state. Hence, we
        ! must remove this extra factor if we wish to calculate the projected eneergy
        ! in the same way. This is done with the factor, importance_sampling_factor.
        if (importance_sampling) importance_sampling_factor = 1.0/neel_singlet_amp(n)
        
        ! Deduce the number of 0-1 bonds, where the 1's are on the
        ! second sublattice:
        ! The total number of 0-1 bonds, n(0-1) can be found from the diagonal 
        ! element of the current basis function:
        ! hmatel = -J_coupling*(ndim*nsites - 2*n(0-1))
        ! This means we can avoid calculating n(0-1) again, which is expensive.
        ! We know the number of 0-1 bonds where the 1 (the spin up) is on sublattice 1,
        ! so can then deduce the number where the 1 is on sublattice 2.
        lattice_2_up = ((ndim*nsites) + nint(walker_energies(1,idet)/J_coupling))/2 - lattice_1_up
        
        ! There are three contributions to add to the projected energy from
        ! the current basis function. Consider the Neel singlet state:
        ! |NS> = \sum_{i} a_i|D_i>
        ! The amplitude a_i only depend on the number of spins up on sublattice 1.
        ! We want to calculate \sum_{i} (a_i * <D_i|H|D_j> * n_j) where |D_j> is the
        ! current basis function, and then add this to the current projected energy.
        
        ! Firstly, consider the diagonal term:
        ! We have <D_j|H|D_j> stored, so this is simple:
        proj_energy = proj_energy + (neel_singlet_amp(n) * walker_energies(1,idet) * &
                                          walker_population(1,idet) * importance_sampling_factor)
        
        ! Now, to find all other basis functions connected to |D_j>, we find 0-1 bonds
        ! and then flip both of these spins. The resulting basis function, |D_i> will be
        ! connected. The amplitude a_i only depends on the number of spins up on
        ! sublattice 1. This will depend on whether, for the 0-1 bond flipped,
        ! the 1 (up spin) was one sublattice 1 or 2. If the up spin was on lattice
        ! 1, there will be one less up spin on sublattice 1 after the flip. If it
        ! was on sublattice 2, there will be one extra spin. So we can have either of
        ! the two amplitudes, neel_singlet_amp(n-1) or neel_singlet_amp(n+1) 
        ! respectively. We just need to know how many of each type of 0-1 bonds there
        ! are. But we have these already - they are stored as lattice_1_up and lattice_2_up.
        ! Finally note that the matrix element is -2*J_coupling, and we can put this together...
        
        ! From 0-1 bonds where the 1 is on sublattice 1, we have:
        proj_energy = proj_energy - (2 * J_coupling * lattice_1_up * &
                    walker_population(1,idet) * neel_singlet_amp(n-1) * importance_sampling_factor)
                                 
        ! And from 1-0 bond where the 1 is on sublattice 2, we have:
        proj_energy = proj_energy - (2 * J_coupling * lattice_2_up * &
                    walker_population(1,idet) * neel_singlet_amp(n+1) * importance_sampling_factor)
        
        ! Now we just need to find the contribution to the denominator. The total
        ! denominator is
        ! \sum_{i} (a_i * n_i)
        ! Hence from this paritcular basis function, |D_j>, we just add (a_j * n_j)
        
        D0_population = D0_population + &
                          walker_population(1,idet)*neel_singlet_amp(n)*importance_sampling_factor
        
        if (calculate_magnetisation) call update_magnetisation_heisenberg(idet)

    end subroutine update_proj_energy_heisenberg_neel_singlet
    
    pure function neel_singlet_data(f) result(spin_config_data)
    
        ! This subroutine calculates the number of up spins on one of the
        ! sublattices so that it can be stored, and also the total number
        ! of 0-1 bonds where the up spins lie on this same sublattice. These are
        ! later used in the update_proj_energy_heisenberg_neel_singlet subroutine.
        ! This procedure is for the Heisenberg model only.
        ! In:
        !    f: bit string representation of the basis function.
        ! Returns: 
        !    spin_config_data: A two component integer vector, where the first
        !       component stores the total number of spins up on the first sublattice 
        !       for this basis function, and the second component stores the numbers 
        !       of 0-1 bonds, where the 1 (the up spin) is on the first sublattice.
    
        use basis, only: basis_length, basis_lookup
        use bit_utils, only: count_set_bits
        use determinants, only: lattice_mask
        use hubbard_real, only: connected_orbs
        use system, only: nsites
    
        integer(i0), intent(in) :: f(basis_length)
        integer(i0) :: f_mask(basis_length), f_not(basis_length), g(basis_length)
        integer :: lattice_1_up, n, i, ipos, basis_find
        integer :: spin_config_data(2)
    
        ! Calculate the number of up spins on the first sublattice.
        f_mask = iand(f, lattice_mask)
        n = sum(count_set_bits(f_mask))
        if (n > nsites/2) n = nsites/2 - n
    
        ! Find the number of 0-1 bonds where the 1 lies on the first sublattice.
        lattice_1_up = 0
        f_not = not(f)
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
        ! magnetisation. This is only run when calculate_magnetisation
        ! is true, and only in the Heisenberg model.
        ! The staggered magnetisation squared in the z direction is:
        !   \sum_{i} <D_i|M_z^2|D_i> N_i^2/Magnitude^2
        ! where N_i is the population on the i-th basis function,
        ! and Magnitude^2 is the sum of the squares of the components
        ! of the wavefunction, ie, the sum of the squares of the psips.
        ! During a MC cycle we store
        !   \sum_{i} <D_i|M_z^2|D_i> N_i^2 and Magnitude^2
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

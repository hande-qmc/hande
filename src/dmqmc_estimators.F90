module dmqmc_estimators

use const

implicit none

contains

   subroutine update_dmqmc_estimators(ireport, ntot_particles_old)

        ! Update the shift and average the shift and estimators

        ! Should be called every report loop in an DMQMC calculation.

        ! In:
        !    ireport: index of the current report loop.
        ! Inout:
        !    ntot_particles_old: total number (across all processors) of
        !        particles in the simulation at end of the previous report loop.
        !        Returns the current total number of particles for use in the
        !        next report loop.

        use calc, only: doing_dmqmc_calc, dmqmc_energy, dmqmc_staggered_magnetisation
        use energy_evaluation, only: update_shift
        use fciqmc_data, only: nparticles, sampling_size, target_particles, rspawn
        use fciqmc_data, only: shift, av_shift, vary_shift, start_vary_shift
        use fciqmc_data, only: nreport, ncycles, trace, thermal_energy, thermal_staggered_mag
        use fciqmc_data, only: total_trace, total_estimator_numerators, number_dmqmc_estimators
        use parallel

        integer, intent(in) :: ireport
        integer, intent(inout) :: ntot_particles_old(sampling_size)
        integer :: counter

#ifdef PARALLEL
        real(dp) :: ir(sampling_size+1+((number_dmqmc_estimators+1)*ncycles))
        real(dp) :: ir_sum(sampling_size+1+((number_dmqmc_estimators+1)*ncycles))
        integer :: ntot_particles(sampling_size), ierr

        ! Need to sum the number of particles and other quantites over all processors.
        ir(1:sampling_size) = nparticles
        ir(sampling_size+1) = rspawn
        ! The trace is stored for each iteration within the report loop, so will occupy
        ! components from sampling_size+2 to sampling_size+1+ncycles, since it will
        ! have a total of ncycle components.
        ir(sampling_size+2:sampling_size+1+ncycles) = trace
        counter = 0
        if (doing_dmqmc_calc(dmqmc_energy)) then
            counter = counter + 1
            ! Similarly for all the other estimators, as with trace, they will occupy
            ! a total of ncycle components, one for each iteration, so the components
            ! which need to be filled for each estimator, is as below.
            ir( sampling_size+2+(counter*ncycles):sampling_size+1+((counter+1)*ncycles) )=thermal_energy
        end if
        if (doing_dmqmc_calc(dmqmc_staggered_magnetisation)) then
            counter = counter + 1
            ir( sampling_size+2+(counter*ncycles):sampling_size+1+((counter+1)*ncycles) )=thermal_staggered_mag
        end if
        ! Merge the lists from each processor together.
        call mpi_allreduce(ir, ir_sum, size(ir), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
        ntot_particles = nint(ir_sum(1:sampling_size))
        rspawn = ir_sum(sampling_size+1)
        if (parent) then
            ! Let the parent processor store the merged data lists, so that these
            ! can be output afterwards by the parent.
            total_trace = ir_sum(sampling_size+2:sampling_size+1+ncycles)
            do counter = 1, number_dmqmc_estimators
               total_estimator_numerators(:,counter) = &
                  + ir_sum( sampling_size+2+(counter*ncycles):sampling_size+1+((counter+1)*ncycles) )
            end do
        end if

        if (vary_shift) then
            call update_shift(ntot_particles_old(1), ntot_particles(1), ncycles)
        end if
        ntot_particles_old = ntot_particles
        if (ntot_particles(1) > target_particles .and. .not.vary_shift) then
            vary_shift = .true.
            start_vary_shift = ireport
        end if

#else
        ! If only one a single core, don't need to merge the data
        ! from several processors, simply let the total values be equal
        ! to the standard values, ready to be output.

        total_trace = trace
        counter = 0
        if (doing_dmqmc_calc(dmqmc_energy)) then
            counter = counter+1
            total_estimator_numerators(:,counter) = thermal_energy
        end if
        if (doing_dmqmc_calc(dmqmc_staggered_magnetisation)) then
            counter = counter+1
            total_estimator_numerators(:,counter) = thermal_staggered_mag
        end if

        if (vary_shift) call update_shift(ntot_particles_old(1), nparticles(1), ncycles)
        ntot_particles_old = nparticles
        if (nparticles(1) > target_particles .and. .not.vary_shift) then
            vary_shift = .true.
            start_vary_shift = ireport
        end if        
#endif

            rspawn = rspawn/(ncycles*nprocs)

   end subroutine update_dmqmc_estimators

   subroutine call_dmqmc_estimators(idet, beta_index)
       
       ! This function calls the processes to update the estimators which
       ! have been requested by the user to be calculated.
       ! First, calculate the excitation level between the two bitsrtings
       ! corresponding to the the two ends. Then add the contribution from
       ! the current density matrix element to the trace, which is always
       ! calculated. Then call other estimators, as required.

       ! In:
       !    idet: Current position in the main bitsrting list.
       !    beta_index: Current iteration within the report loop.

       use basis, only: basis_length, total_basis_length
       use calc, only: doing_dmqmc_calc, dmqmc_energy, dmqmc_staggered_magnetisation
       use excitations, only: get_excitation, excit
       use fciqmc_data, only: walker_dets, walker_population, trace, doing_reduced_dm
       use proc_pointers, only: update_dmqmc_energy_ptr, update_dmqmc_stag_mag_ptr

       integer, intent(in) :: idet, beta_index
       type(excit) :: excitation     
   
       ! Get excitation.
       excitation = get_excitation(walker_dets(:basis_length,idet), &
                        walker_dets((1+basis_length):total_basis_length,idet))

       ! If diagonal element, add to the trace.
       if (excitation%nexcit == 0) trace(beta_index)=trace(beta_index) + walker_population(1,idet)
       ! See which estimators are to be calculated, and call the corresponding procedures.
       if (doing_dmqmc_calc(dmqmc_energy)) call update_dmqmc_energy_ptr(idet, beta_index, excitation)
       if (doing_dmqmc_calc(dmqmc_staggered_magnetisation)) &
                                     call update_dmqmc_stag_mag_ptr(idet, beta_index, excitation)
       if (doing_reduced_dm) call update_reduced_density_matrix(idet)

   end subroutine call_dmqmc_estimators

   subroutine dmqmc_energy_heisenberg(idet, beta_index, excitation)

       ! For the Heisenberg model only.
       ! Add the contribution from the current density matrix element
       ! to the thermal energy estimate.

       ! In:
       !    idet: Current position in the main bitsrting list.
       !    beta_index: Current iteration within the report loop.
       !    excitation: excit type variable which stores information on
       !    the excitation between the two bitstring ends, corresponding
       !    to the two labels for the density matrix element.
 
       use basis, only: basis_length, total_basis_length
       use basis, only: bit_lookup
       use excitations, only: excit
       use fciqmc_data, only: walker_dets, walker_population
       use fciqmc_data, only: walker_energies, thermal_energy, H00
       use hubbard_real, only: connected_orbs
       use system, only: J_coupling

       integer, intent(in) :: idet, beta_index
       type(excit), intent(in) :: excitation
       integer :: bit_element, bit_position

       ! If no excitation, we have a diagonal element, so add elements
       ! which involve the diagonal element of the Hamiltonian.
       if (excitation%nexcit == 0) then
           thermal_energy(beta_index)=thermal_energy(beta_index) + &
                             (walker_energies(1,idet)+H00)*walker_population(1,idet)
       else if (excitation%nexcit == 1) then
       ! If not a diagonal element, but only a single excitation, then the corresponding
       ! Hamiltonian element may be non-zero. Calculate if the flipped spins are 
       ! neighbours on the lattice, and if so, add the contirbution from this site.
           bit_position = bit_lookup(1,excitation%from_orb(1))
           bit_element = bit_lookup(2,excitation%from_orb(1))
           if (btest(connected_orbs(bit_element, excitation%to_orb(1)), bit_position)) &
                 thermal_energy(beta_index)=thermal_energy(beta_index) - &
                             (2.0*J_coupling*walker_population(1,idet))
       end if

   end subroutine dmqmc_energy_heisenberg

   subroutine dmqmc_energy_squared_heisenberg(idet, beta_index,excitation)

       ! For the Heisenberg model only.
       ! Add the contribution from the current density matrix element
       ! to the thermal energy squared estimate.

       ! In:
       !    idet: Current position in the main bitsrting list.
       !    beta_index: Current iteration within the report loop.
       !    excitation: excit type variable which stores information on
       !    the excitation between the two bitstring ends, corresponding
       !    to the two labels for the density matrix element.

       use basis, only: basis_length, total_basis_length
       use basis, only: bit_lookup
       use excitations, only: excit
       use fciqmc_data, only: walker_dets, walker_population
       use fciqmc_data, only: walker_energies, thermal_energy_squared, H00
       use hubbard_real, only: connected_orbs, next_nearest_orbs
       use hubbard_real, only: double_next_nearest_orbs
       use system, only: J_coupling

       integer, intent(in) :: idet, beta_index
       type(excit), intent(in) :: excitation
       integer :: bit_element1, bit_position1, bit_element2, bit_position_2
       integer :: sum_H1_H2

       sum_H1_H2 = 0

       if (excitation%nexcit == 0) then
           ! If there are 0 excitations then either nothing happens twice, or we
           ! flip the same pair of spins twice. The Hamiltonian element for doing nothing
           ! is just the diagonal element. For each possible pairs of spins which can be
           ! flipped, there is a mtarix element of -2*J_coupling, so we just need to count
           ! the number of such pairs, which can be found simply from the diagonal element.

           sum_H1_H2 = (walker_energies(1,idet)+H00)**2
           sum_H1_H2 = sum_H1_H2 - J_coupling*nbonds - (walker_energies(1,idet)+H00)

       else if (excitation%nexcit == 1) then
           ! If there is only one excitation (2 spins flipped) then the contribution to H^2
           ! depend on the positions of the spins relative to one another.
           ! If the the spins are nearest neighbors then we could either do nothing and then
           ! flip the pair, or flip the pair and then do nothing.
           ! If next nearest neighbors and there is only one two-bond path to get from one
           ! spin to the other, we first flip the pair on the first bond, then flip the pair
           ! on the second bond. This flipping can only be done in exactly one order, not both -
           ! the two spins which change are opposite, so the middle spin will initially only
           ! be the same as one or the other spin. This is nice, because we don't have check
           ! which way up the intermeddiate spin is - there will always be one order which contributes.
           ! If there are two such paths, then this could happen by either paths, but again, the two
           ! intermeddiate spins will only allow one order of spin flipping for each path, no
           ! matter which way up they are, so we only need to check if there are two possible paths.

           if (next_nearest_orbs(excitation%from_orb(1),excitation%to_orb(1)) /= 0) then
               ! Contribution for next nearest neighbors
               sum_H1_H2 = -2.0*J_coupling*next_nearest_orbs(excitation%from_orb(1),excitation%to_orb(1))
           else
               ! Contributions for nearest neighbors
               bit_position1 = bit_lookup(1,excitation%from_orb(1))
               bit_element1 = bit_lookup(2,excitation%from_orb(1))
               if (btest(connected_orbs(bit_element1, excitation%to_orb(1)), bit_position1)) &
                   sum_H1_H2 = -2.0*J_coupling*(walker_energies(1,idet)+walker_energies(2,idet)+(2*H00))
           end if

       else if (excitation%nexcit == 2) then
           ! If there are two excitations (4 spins flipped) then, once again, the contribution
           ! to the thermal energy squared will depend on the positions of the spins. If there
           ! are two pairs of spins flipped which are separated then there is one way
           ! for this to happen - by flipping one pair, and then the other (this also requires
           ! that the two spins within each neighboring pair are opposite, as ever for the
           ! Heisenberg model). These two flips can happen in either order.
           ! If the you can draw a path from the first spin which goes through every
           ! spin and comes back to the original spin, then there are two ways to pair the
           ! four spins - however, once again, only one way to pair them such that the spins
           ! within a pair are opposite. This is because the four spins must have a total spin of
           ! ms=0 (else flipping them all will give a configuration with a different total ms)
           ! so there must be two up and two down. So for a square lattice, they can only be
           ! paired in one way. Again, the flipping can be done in either order, giving a factor of 2.

           bit_position1 = bit_lookup(1,excitation%from_orb(1))
           bit_element1 = bit_lookup(2,excitation%from_orb(1))
           bit_position2 = bit_lookup(1,excitation%from_orb(2))
           bit_element2 = bit_lookup(2,excitation%from_orb(2))
           if (btest(connected_orbs(bit_element1, excitation%to_orb(1)), bit_position1) .and. &
           btest(connected_orbs(bit_element2, excitation%to_orb(2)), bit_position2) then
               sum_H1_H2 = -4.0*J_coupling
           else if (btest(connected_orbs(bit_element1, excitation%to_orb(2)), bit_position1) .and. &
           btest(connected_orbs(bit_element2, excitation%to_orb(1)), bit_position2)) then
               sum_H1_H2 = -4.0*J_coupling
           end if

           thermal_energy_squared(beta_index) = thermal_energy_squared(beta_index) + &
               sum_H1_H2*walker_population(1,idet)

       end if

   end subroutine dmqmc_energy_squared_heisenberg

   subroutine dmqmc_energy_hub_real(idet, beta_index, excitation)

       ! For the Heisenberg model only.
       ! Add the contribution from the current density matrix element
       ! to the thermal energy estimate.

       ! In:
       !    idet: Current position in the main bitsrting list.
       !    beta_index: Current iteration within the report loop.
       !    excitation: excit type variable which stores information on
       !    the excitation between the two bitstring ends, corresponding
       !    to the two labels for the density matrix element.

       use basis, only: basis_length, total_basis_length
       use excitations, only: excit
       use fciqmc_data, only: walker_dets, walker_population
       use fciqmc_data, only: walker_energies, thermal_energy, H00
       use hamiltonian, only: slater_condon1_hub_real        

       integer, intent(in) :: idet, beta_index
       type(excit), intent(in) :: excitation
       integer :: bit_element, bit_position
       real(p) :: hmatel

       ! If no excitation, we have a diagonal element, so add elements
       ! which involve the diagonal element of the Hamiltonian.
       if (excitation%nexcit == 0) then
           thermal_energy(beta_index)=thermal_energy(beta_index) + &
                             (walker_energies(1,idet)+H00)*walker_population(1,idet)
       else if (excitation%nexcit == 1) then
       ! If not a diagonal element, but only a single excitation, then the corresponding
       ! Hamiltonian element may be non-zero. Calculate if the flipped spins are
       ! neighbours on the lattice, and if so, add the contirbution from this site.
           hmatel = slater_condon1_hub_real(excitation%from_orb(1), excitation%to_orb(1), excitation%perm)
           thermal_energy(beta_index)=thermal_energy(beta_index) + (hmatel*walker_population(1,idet))
       end if

   end subroutine dmqmc_energy_hub_real

   subroutine dmqmc_stag_mag_heisenberg(idet, beta_index, excitation)

       ! For the Heisenberg model only.
       ! Add the contribution from the current density matrix element
       ! to the thermal staggered magnetisation estimate.

       ! In:
       !    idet: Current position in the main bitsrting list.
       !    beta_index: Current iteration within the report loop.
       !    excitation: excit type variable which stores information on
       !    the excitation between the two bitstring ends, corresponding
       !    to the two labels for the density matrix element.

       use basis, only: basis_length, total_basis_length
       use bit_utils, only: count_set_bits
       use determinants, only: lattice_mask
       use excitations, only: excit
       use fciqmc_data, only: walker_dets, walker_population, thermal_staggered_mag
       use system, only: nel

       integer, intent(in) :: idet, beta_index
       type(excit), intent(in) :: excitation
       integer :: bit_element, bit_position, stag_mag_el, sublattice1_up_spins
       integer(i0) :: f_mask(basis_length)

       ! The staggered magnetisation operator in the z-direction is diagonal
       ! in a basis of spins which are eigenstates of the z-direction spin
       ! operator. Hence, we only have contributions from the diagonal
       ! elements, where there is no excitation.
       if (excitation%nexcit == 0) then
           ! Calculate the number of up spins on sublattice 1. From this, and
           ! the total number of spins up, nel, we may calculate the matrix
           ! element for the staggered magnetisation operator for this spin
           ! configuartion.
           ! See the function 'diagonal_element_heisenberg_staggered' in
           ! hamiltonian.f90 for an explnation of why the matrix element is
           ! 4*sublattice1_ip_spins - 2*nel.
           f_mask = iand(walker_dets(:basis_length,idet), lattice_mask)
           sublattice1_up_spins = sum(count_set_bits(f_mask))
           ! Matrix element of staggered magnetisation.
           stag_mag_el = (4*sublattice1_up_spins - 2*nel)
           thermal_staggered_mag(beta_index)=thermal_staggered_mag(beta_index) + &
                             stag_mag_el*walker_population(1,idet)

       end if

   end subroutine dmqmc_stag_mag_heisenberg

   subroutine update_reduced_density_matrix(idet)

       use basis, only: basis_length, total_basis_length
       use dmqmc_procedures, only: decode_dm_bitstring
       use fciqmc_data, only: subsystem_B_mask, reduced_density_matrix
       use fciqmc_data, only: walker_dets, walker_population

       integer, intent(in) :: idet
       integer(i0) :: f1(basis_length), f2(basis_length)
       integer(i0) :: end1, end2

       f1 = iand(subsystem_B_mask,walker_dets(:basis_length,idet))
       f2 = iand(subsystem_B_mask,walker_dets(basis_length+1:total_basis_length,idet))

       if (sum(abs(f1-f2)) == 0) then
           call decode_dm_bitstring(walker_dets(:,idet),end1,end2)
           reduced_density_matrix(end1,end2) = reduced_density_matrix(end1,end2) + walker_population(1,idet)
       end if

    end subroutine update_reduced_density_matrix


end module dmqmc_estimators

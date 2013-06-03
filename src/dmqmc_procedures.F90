module dmqmc_procedures

use const
implicit none

contains

    subroutine init_dmqmc()

         use basis, only: basis_length, total_basis_length, bit_lookup, basis_lookup
         use calc, only: doing_dmqmc_calc, dmqmc_calc_type, dmqmc_energy, dmqmc_energy_squared
         use calc, only: dmqmc_staggered_magnetisation, dmqmc_correlation
         use checking, only: check_allocate, check_deallocate
         use fciqmc_data, only: trace, energy_index, energy_squared_index, correlation_index
         use fciqmc_data, only: staggered_mag_index, estimator_numerators, doing_reduced_dm
         use fciqmc_data, only: dmqmc_factor, number_dmqmc_estimators, ncycles, tau, dmqmc_weighted_sampling
         use fciqmc_data, only: correlation_mask, correlation_sites, half_density_matrix
         use fciqmc_data, only: dmqmc_sampling_probs, dmqmc_accumulated_probs, flip_spin_matrix
         use fciqmc_data, only: doing_concurrence, calculate_excit_distribution, excit_distribution
         use fciqmc_data, only: nreport, average_shift_until, shift_profile, dmqmc_vary_weights
         use fciqmc_data, only: finish_varying_weights, weight_altering_factors, dmqmc_find_weights
         use parallel, only: parent
         use system, only: max_number_excitations

         integer :: ierr, i, bit_position, bit_element

         number_dmqmc_estimators = 0
         trace = 0

         if (doing_dmqmc_calc(dmqmc_energy)) then
             number_dmqmc_estimators = number_dmqmc_estimators + 1
             energy_index = number_dmqmc_estimators
         end if
         if (doing_dmqmc_calc(dmqmc_energy_squared)) then
             number_dmqmc_estimators = number_dmqmc_estimators + 1
             energy_squared_index = number_dmqmc_estimators
         end if
         if (doing_dmqmc_calc(dmqmc_correlation)) then
             number_dmqmc_estimators = number_dmqmc_estimators + 1
             correlation_index = number_dmqmc_estimators
             allocate(correlation_mask(1:basis_length), stat=ierr)
             call check_allocate('correlation_mask',basis_length,ierr)
             correlation_mask = 0
             do i = 1, 2
                 bit_position = bit_lookup(1,correlation_sites(i))
                 bit_element = bit_lookup(2,correlation_sites(i))
                 correlation_mask(bit_element) = ibset(correlation_mask(bit_element), bit_position)
             end do
         end if
         if (doing_dmqmc_calc(dmqmc_staggered_magnetisation)) then
             number_dmqmc_estimators = number_dmqmc_estimators + 1
             staggered_mag_index = number_dmqmc_estimators
         end if

         allocate(estimator_numerators(1:number_dmqmc_estimators), stat=ierr)
         call check_allocate('estimator_numerators',number_dmqmc_estimators,ierr)
         estimator_numerators = 0

         if (calculate_excit_distribution .or. dmqmc_find_weights) then
             allocate(excit_distribution(0:max_number_excitations), stat=ierr)
             call check_allocate('excit_distribution',max_number_excitations+1,ierr)             
             excit_distribution = 0.0_p
         end if

         if (average_shift_until > 0) then
             allocate(shift_profile(1:nreport+1), stat=ierr)
             call check_allocate('shift_profile',nreport+1,ierr)
             shift_profile = 0.0_p
         end if

         ! In DMQMC we want the spawning probabilities to have an extra factor of a half,
         ! because we spawn from two different ends with half probability. To avoid having
         ! to multiply by an extra variable in every spawning routine to account for this, we
         ! multiply the time step by 0.5 instead, then correct this in the death step (see below).
         tau = tau*0.5_p
         ! Set dmqmc_factor to 2 so that when probabilities in death.f90 are multiplied
         ! by this factor it cancels the factor of 0.5 introduced into the timestep in DMQMC.
         ! Every system uses the same death routine, so this factor only needs to be added once.
         ! This factor is also used in updated the shift, where the true tau is needed.
         dmqmc_factor = 2.0_p

         if (dmqmc_weighted_sampling) then
             ! dmqmc_sampling_probs stores the factors by which probabilities are to
             ! be reduced when spawning away from the diagonal. The trial function required
             ! from these probabilities, for use in importance sampling, is actually that of
             ! the accumulated factors, ie, if dmqmc_sampling_probs = (a, b, c, ...) then
             ! dmqmc_accumulated_factors = (1, a, ab, abc, ...). This is the array which we
             ! need to create and store. dmqmc_sampling_probs is no longer needed and so can
             ! be deallocated. Also, the user may have only input factors for the first few
             ! excitation levels, but we need to store factors for all levels, as done below.
             if (.not.allocated(dmqmc_sampling_probs)) then
                 allocate(dmqmc_sampling_probs(1:max_number_excitations), stat=ierr)
                 call check_allocate('dmqmc_sampling_probs',max_number_excitations,ierr)
                 dmqmc_sampling_probs = 1.0_p
             end if
             if (half_density_matrix) dmqmc_sampling_probs(1) = dmqmc_sampling_probs(1)*2.0_p
             allocate(dmqmc_accumulated_probs(0:max_number_excitations), stat=ierr)
             call check_allocate('dmqmc_accumulated_probs',max_number_excitations+1,ierr)
             dmqmc_accumulated_probs(0) = 1.0_p
             do i = 1, size(dmqmc_sampling_probs)
                 dmqmc_accumulated_probs(i) = dmqmc_accumulated_probs(i-1)*dmqmc_sampling_probs(i)
             end do
             dmqmc_accumulated_probs(size(dmqmc_sampling_probs)+1:max_number_excitations) = &
                                    dmqmc_accumulated_probs(size(dmqmc_sampling_probs))
             if (dmqmc_vary_weights) then
                 ! Allocate an array to store the factors by which the weights will change each
                 ! iteration.
                 allocate(weight_altering_factors(0:max_number_excitations), stat=ierr)
                 call check_allocate('weight_altering_factors',max_number_excitations+1,ierr) 
                 weight_altering_factors = dble(dmqmc_accumulated_probs)**(1/dble(finish_varying_weights))
                 ! If varying the weights, start the accumulated probabilties as all 1.0
                 ! initially, and then alter them gradually later.
                 dmqmc_accumulated_probs = 1.0_p
             end if
         else
             allocate(dmqmc_accumulated_probs(0:max_number_excitations), stat=ierr)
             call check_allocate('dmqmc_accumulated_probs',max_number_excitations+1,ierr)
             dmqmc_accumulated_probs = 1.0_p
             if (half_density_matrix) dmqmc_accumulated_probs(1:max_number_excitations) &
                                      = 2.0_p*dmqmc_accumulated_probs(1:max_number_excitations)
         end if

         ! If doing a reduced density matrix calculation, allocate and define the bit masks that
         ! have 1's at the positions referring to either subsystems A or B.
         if (doing_reduced_dm) call setup_rdm_arrays()

         ! If doing concurrence calculation then construct and store the 4x4 flip spin matrix i.e.
         ! \sigma_y \otimes \sigma_y
 
         if (doing_concurrence) then
             allocate(flip_spin_matrix(4,4), stat=ierr)
             call check_allocate('flip_spin_matrix', 16,ierr)
             flip_spin_matrix = 0._p
             flip_spin_matrix(1,4) = -1._p
             flip_spin_matrix(4,1) = -1._p
             flip_spin_matrix(3,2) = 1._p
             flip_spin_matrix(2,3) = 1._p    
         end if

    end subroutine init_dmqmc

    subroutine setup_rdm_arrays()

        ! Setup the bit masks needed for RDM calculations. These are masks for the bits referring
        ! to either subsystem A or B. Also calculate the positions and elements of the sites
        ! in subsyetsm A, and finally allocate the RDM itself.

        use basis, only: basis_length, bit_lookup, basis_lookup
        use checking, only: check_allocate, check_deallocate
        use errors
        use fciqmc_data, only: subsystem_A_mask, subsystem_B_mask, subsystem_A_bit_positions
        use fciqmc_data, only: subsystem_A_list, reduced_density_matrix, subsystem_A_size
        use fciqmc_data, only: output_rdm, rdm_unit
        use system, only: system_type, heisenberg, nsites
        use utils, only: get_free_unit

        integer :: i, ierr, ipos, basis_find, bit_position, bit_element

        subsystem_A_size = ubound(subsystem_A_list,1)
        allocate(subsystem_A_mask(1:basis_length), stat=ierr)
        call check_allocate('subsystem_A_mask',basis_length,ierr)
        allocate(subsystem_B_mask(1:basis_length), stat=ierr)
        call check_allocate('subsystem_B_mask',basis_length,ierr)
        allocate(subsystem_A_bit_positions(subsystem_A_size,2), stat=ierr)
        call check_allocate('subsystem_A_bit_positions',2*subsystem_A_size,ierr)
        subsystem_A_mask = 0
        subsystem_B_mask = 0
        subsystem_A_bit_positions = 0

        ! For the Heisenberg model only currently.
        if (system_type==heisenberg) then
            do i = 1, subsystem_A_size
                bit_position = bit_lookup(1,subsystem_A_list(i))
                bit_element = bit_lookup(2,subsystem_A_list(i))
                subsystem_A_mask(bit_element) = ibset(subsystem_A_mask(bit_element), bit_position)
                subsystem_A_bit_positions(i,1) = bit_position
                subsystem_A_bit_positions(i,2) = bit_element
            end do
            subsystem_B_mask = subsystem_A_mask
            ! We cannot just flip the mask for system A to get that for system B, because
            ! there the trailing bits on the end don't refer to anything and should be
            ! set to 0. So, first set these to 1 and then flip all the bits.
            do ipos = 0, i0_end
                basis_find = basis_lookup(ipos, basis_length)
                if (basis_find == 0) then
                    subsystem_B_mask(basis_length) = ibset(subsystem_B_mask(basis_length),ipos)
                end if
            end do
            subsystem_B_mask = not(subsystem_B_mask)
        else
            call stop_all("setup_rdm_arrays","The use of RDMs is currently only implemented for the &
                           &Heisenberg model.")
        end if

        if (subsystem_A_size <= int(nsites/2)) then
            ! In this case, for an ms = 0 subspace (as the ground state of the Heisenberg model
            ! will be) then any combination of spins can occur in the subsystem, from all spins
            ! down to all spins up. Hence the total size of the reduced density matrix will be
            ! 2**(number of spins in subsystem A).
            allocate(reduced_density_matrix(2**subsystem_A_size,2**subsystem_A_size), stat=ierr)
            call check_allocate('reduced_density_matrix', 2**(2*subsystem_A_size),ierr)
            reduced_density_matrix = 0
        else if (subsystem_A_size == nsites) then
            allocate(reduced_density_matrix(2**subsystem_A_size,2**subsystem_A_size), stat=ierr)
            call check_allocate('reduced_density_matrix', 2**(2*subsystem_A_size),ierr)
            reduced_density_matrix = 0
        end if
        
        if (output_rdm) then
            rdm_unit = get_free_unit()
            open(rdm_unit, file='reduced_dm', status='replace')
            write(rdm_unit,'(a37,i6)') "# Number of elements in each RDM row:", &
                                       ubound(reduced_density_matrix,1)
            close(rdm_unit)
        end if

    end subroutine setup_rdm_arrays

    subroutine random_distribution_heisenberg()

        ! For the Heisenberg model only. Distribute the initial number of psips
        ! along the main diagonal. Each diagonal element should be chosen
        ! with the same probability.

        ! Currently this creates psips with Ms = ms_in only.

        ! If we have number of sites = nsites,
        ! and total spin value = ms_in,
        ! then number of up spins is equal to up_spins = (ms_in + nsites)/2.

        ! Start from state with all spins down, then choose the above number of
        ! spins to flip up with equal probability.

        use basis, only: nbasis, basis_length, bit_lookup
        use calc, only: ms_in
        use dSFMT_interface, only:  genrand_real2
        use fciqmc_data, only: D0_population
        use parallel
        use system, only: nsites

        integer :: i, up_spins, rand_basis, bits_set
        integer :: bit_element, bit_position, npsips
        integer(i0) :: f(basis_length)
        real(dp) :: rand_num

        up_spins = (ms_in+nsites)/2
        npsips = int(D0_population/nprocs)
        ! If initial number of psips does not split evenly between all processors,
        ! add the leftover psips to the first processors in order.
        if (D0_population-(nprocs*int(D0_population/nprocs)) > iproc) npsips = npsips+1

        do i = 1, npsips

            ! Start with all spins down.
            f = 0
            bits_set = 0

            do
                ! If half the spins are now flipped up, we have our basis
                ! function fully created, so exit the loop.
                if (bits_set==up_spins) exit
                ! Choose a random spin to flip.
                rand_num = genrand_real2()
                rand_basis = ceiling(rand_num*nbasis)
                ! Find the corresponding positions for this spin.
                bit_position = bit_lookup(1,rand_basis)
                bit_element = bit_lookup(2,rand_basis)
                if (.not. btest(f(bit_element),bit_position)) then
                    ! If not flipped up, flip the spin up.
                    f(bit_element) = ibset(f(bit_element),bit_position)
                    bits_set = bits_set + 1
                end if
            end do

            ! Now call a routine to add the corresponding diagonal element to
            ! the spawned walkers list.
            call create_particle(f,f,1)

        end do

    end subroutine random_distribution_heisenberg

    subroutine create_particle(f1,f2,nspawn)

        ! Create a psip on a diagonal element of the density
        ! matrix by adding it to the spawned walkers list. This
        ! list can then be sorted correctly by the direct_annihilation
        ! routine

        ! In:
        !    f_new: Bit string representation of index of the diagonal
        !           element upon which a new psip shall be placed.

        use hashing
        use basis, only: basis_length, total_basis_length
        use fciqmc_data, only: spawned_walkers, spawning_head, spawned_pop
        use parallel

        integer(i0), intent(in) :: f1(basis_length), f2(basis_length)
        integer, intent(in) :: nspawn
        integer(i0) :: f_new(total_basis_length)
#ifndef PARALLEL
        integer, parameter :: iproc_spawn = 0
#else
        integer :: iproc_spawn
#endif

        ! Create the bitstring of the psip.
        f_new = 0
        f_new(:basis_length) = f1
        f_new((basis_length+1):(total_basis_length)) = f2

#ifdef PARALLEL
        ! Need to determine which processor the spawned walker should be sent to.
        iproc_spawn = modulo(murmurhash_bit_string(f_new, &
                                (total_basis_length)), nprocs)
#endif

        ! Move to the next position in the spawning array.
        spawning_head(iproc_spawn) = spawning_head(iproc_spawn) + 1

        ! Set info in spawning array.
        ! Zero it as not all fields are set.
        spawned_walkers(:,spawning_head(iproc_spawn)) = 0
        ! indices 1 to total_basis_length store the bitstring.
        spawned_walkers(:(2*basis_length),spawning_head(iproc_spawn)) = f_new
        ! The final index stores the number of psips created.
        spawned_walkers((2*basis_length)+1,spawning_head(iproc_spawn)) = nspawn

    end subroutine create_particle

    subroutine decode_dm_bitstring(f, index1, index2)

        ! This function maps an input DMQMC bitstring to two indices
        ! giving the corresponding position of the bitstring in the reduced
        ! density matrix.

        use basis, only: basis_length
        use fciqmc_data, only: subsystem_A_bit_positions, subsystem_A_bit_positions
        use fciqmc_data, only: subsystem_A_size

        integer(i0), intent(in) :: f(basis_length*2)
        integer(i0), intent(out) :: index1, index2
        integer :: i

        ! Start from all bits down, so that we can flip bits up one by one.
        index1 = 0
        index2 = 0

        ! Loop over all the sites in the sublattice considered for the reduced density matrix.
        do i = 1, subsystem_A_size
            ! If the spin is up, flip the corresponding bit in the first index up.
            if (btest(f(subsystem_A_bit_positions(i,2)),subsystem_A_bit_positions(i,1))) &
                index1 = ibset(index1,i-1)
            ! Similarly for the second index, by looking at the second end of the bitstring.
            if (btest(f(subsystem_A_bit_positions(i,2)+basis_length),subsystem_A_bit_positions(i,1))) &
                index2 = ibset(index2,i-1)
        end do

        ! The process above maps to numbers between 0 and 2^subsystem_A_size-1, but the smallest
        ! and largest of the reduced density matrix are one more than these, so add one...
        index1 = index1+1
        index2 = index2+1

    end subroutine decode_dm_bitstring
 
    subroutine update_sampling_weights()
        
        ! This routine updates the values of the weights used in importance sampling. It also
        ! removes or adds psips from the various levels accordingly.

        ! In:
        !    iteration: The current iteration in the beta loop.

        use annihilation, only: remove_unoccupied_dets
        use basis, only: basis_length, total_basis_length
        use excitations, only: get_excitation_level
        use fciqmc_data, only: dmqmc_accumulated_probs, finish_varying_weights
        use fciqmc_data, only: weight_altering_factors, tot_walkers, walker_dets, walker_population
        use fciqmc_data, only: nparticles
        use dSFMT_interface, only:  genrand_real2

        integer :: idet, excit_level, nspawn, sign_factor, old_population
        real(p) :: new_factor
        real(dp) :: rand_num, prob

        ! Alter weights for the next iteration.
        dmqmc_accumulated_probs = dble(dmqmc_accumulated_probs)*weight_altering_factors

        ! When the weights for an excitation level are increased by a factor, the number
        ! of psips on that level has to decrease by the same factor, else the wavefunction
        ! which the psips represent will not be the correct importance sampled wavefunction
        ! for the new weights. The code below loops over every psips and destorys (or creates)
        ! it with the appropriate probability.
        do idet = 1, tot_walkers
            excit_level = get_excitation_level(walker_dets(1:basis_length,idet),&
                    walker_dets(basis_length+1:total_basis_length,idet))
            old_population = abs(walker_population(1,idet))
            rand_num = genrand_real2()
            ! If weight_altering_factors(excit_level) > 1, need to kill psips.
            ! If weight_altering_factors(excit_level) < 1, need to create psips.
            prob = abs(1.0_dp - weight_altering_factors(excit_level)**(-1))*old_population
            nspawn = int(prob)
            prob = prob - nspawn
            if (rand_num < prob) nspawn = nspawn + 1
            if (weight_altering_factors(excit_level) > 1.0_dp) then
                sign_factor = -1
            else
                sign_factor = +1
            end if
            nspawn = sign(nspawn,walker_population(1,idet)*sign_factor)
            ! Update the population on this determinant.
            walker_population(1,idet) = walker_population(1,idet) + nspawn
            ! Update the total number of walkers
            nparticles(1) = nparticles(1) - old_population + abs(walker_population(1,idet))
        end do

        ! Call the annihilation routine to update the main walker list, as some
        ! sites will now have no psips on and so need removing from the simulation.
        call remove_unoccupied_dets()

    end subroutine update_sampling_weights

    subroutine output_and_alter_weights()

        ! This routine will alter and output the sampling weights used in importance 
        ! sampling. It uses the excitation distribution, calculated on the beta loop
        ! which has just finished, and finds the weights needed so that each excitation
        ! level will have roughly equal numbers of psips in the next loop. For example,
        ! to find the weights of psips on the 1st excitation level, divide the number of
        ! psips on the 1st excitation level by the number on the 0th level, then multiply
        ! the old sampling weight by this number to give the new weight. This can be used
        ! when the weights are being introduced gradually each beta loop, too. The weights
        ! are output and can then be used in future DMQMC runs.

        use fciqmc_data, only: dmqmc_sampling_probs, dmqmc_accumulated_probs
        use fciqmc_data, only: excit_distribution, finish_varying_weights
        use fciqmc_data, only: dmqmc_vary_weights, weight_altering_factors
        use parallel
        use system, only: max_number_excitations

        integer :: i, ierr
#ifdef PARALLEL
        real(p) :: merged_excit_dist(max_number_excitations) 
        call mpi_allreduce(excit_distribution, merged_excit_dist, max_number_excitations, &
            MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
        
        excit_distribution = merged_excit_dist        
#endif

        ! It is assumed that there is an even maximum number of excitations.
        do i = 1, (max_number_excitations/2)
            ! Don't include levels where there are very few psips accumulated.
            if (excit_distribution(i-1) > 10.0_p .and. excit_distribution(i) > 10.0_p) then
                ! Alter the sampling weights using the relevant excitation distribution.
                dmqmc_sampling_probs(i) = dmqmc_sampling_probs(i)*&
                    (excit_distribution(i)/excit_distribution(i-1))
                dmqmc_sampling_probs(max_number_excitations+1-i) = dmqmc_sampling_probs(i)**(-1)
            end if
        end do
        
        ! Recalculate dmqmc_accumulated_probs with the new weights.
        do i = 1, max_number_excitations
            dmqmc_accumulated_probs(i) = dmqmc_accumulated_probs(i-1)*dmqmc_sampling_probs(i)
        end do

        ! If dmqmc_vary_weights is true then the weights are to be introduced gradually at the
        ! start of each beta loop. This required redefining weight_altering_factors to coincide
        ! with the new sampling weights.
        if (dmqmc_vary_weights) then
            weight_altering_factors = dble(dmqmc_accumulated_probs)**(1/dble(finish_varying_weights))
            ! Reset the weights for the next loop.
            dmqmc_accumulated_probs = 1.0_p
        end if

        if (parent) then
            ! Print out weights in a form which can be copied into an input file.
            write(6, '(a31,2X)', advance = 'no') ' # Importance sampling weights:'
            do i = 1, max_number_excitations
                write (6, '(es12.4,2X)', advance = 'no') dmqmc_sampling_probs(i)
            end do
            write (6, '()', advance = 'yes')
        end if

    end subroutine output_and_alter_weights

end module dmqmc_procedures

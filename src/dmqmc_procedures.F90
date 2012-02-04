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
         use fciqmc_data, only: staggered_mag_index, estimator_numerators, subsystem_A_size
         use fciqmc_data, only: subsystem_A_mask, subsystem_B_mask, subsystem_A_bit_positions
         use fciqmc_data, only: subsystem_A_list, dmqmc_factor, number_dmqmc_estimators, ncycles
         use fciqmc_data, only: reduced_density_matrix, doing_reduced_dm, tau
         use fciqmc_data, only: correlation_mask, correlation_sites, half_density_matrix
         use fciqmc_data, only: dmqmc_sampling_probs, dmqmc_accumulated_probs
         use parallel, only: parent
         use system, only: system_type, heisenberg, nsites, nel

         integer :: ierr
         integer :: i, ipos, basis_find, bit_position, bit_element
         integer :: max_number_excitations

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

         if (allocated(dmqmc_sampling_probs)) then
             if (half_density_matrix) dmqmc_sampling_probs(1) = dmqmc_sampling_probs(1)*2.0_p
             ! dmqmc_sampling_probs stores the factors by which probabilities are to
             ! be reduced when spawning away from the diagonal. The trial function required
             ! from these probabilities, for use in importance sampling, is actually that of
             ! the accumulated factors, ie, if dmqmc_sampling_probs = (a, b, c, ...) then
             ! dmqmc_accumulated_factors = (1, a, ab, abc, ...). This is the array which we
             ! need to create and store. dmqmc_sampling_probs is no longer needed and so can
             ! be deallocated. Also, the user may have only input factors for the first few
             ! excitation levels, but we need to store factors for all levels, as done below.
             max_number_excitations = min(nel, (nsites-nel))
             allocate(dmqmc_accumulated_probs(0:max_number_excitations), stat=ierr)
             call check_allocate('dmqmc_accumulated_probs',max_number_excitations+1,ierr)
             dmqmc_accumulated_probs(0) = 1.0_p
             do i = 1, size(dmqmc_sampling_probs)
                 dmqmc_accumulated_probs(i) = dmqmc_accumulated_probs(i-1)*dmqmc_sampling_probs(i)
             end do
             dmqmc_accumulated_probs(size(dmqmc_sampling_probs):max_number_excitations) = &
                                        dmqmc_sampling_probs(size(dmqmc_sampling_probs))
             deallocate(dmqmc_sampling_probs, stat=ierr)
            call check_deallocate('dmqmc_sampling_probs',ierr)
         else
             max_number_excitations = min(nel, (nsites-nel))
             allocate(dmqmc_accumulated_probs(0:max_number_excitations), stat=ierr)
             call check_allocate('dmqmc_accumulated_probs',max_number_excitations+1,ierr)
             dmqmc_accumulated_probs = 1.0_p
             if (half_density_matrix) dmqmc_accumulated_probs(1:max_number_excitations) &
                                      = 2.0_p*dmqmc_accumulated_probs(1:max_number_excitations)
         end if

         ! If doing a reduced density matrix calculation, then allocate and define the
         ! bit masks that have 1's at the positions referring to either subsystems A or B.
         if (doing_reduced_dm) then
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
                 subsystem_A_mask = 0
                 do i = 1, subsystem_A_size
                     bit_position = bit_lookup(1,subsystem_A_list(i))
                     bit_element = bit_lookup(2,subsystem_A_list(i))
                     subsystem_A_mask(bit_element) = ibset(subsystem_A_mask(bit_element), bit_position)
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
             end if
             if (subsystem_A_size <= int(nsites/2)) then
                 ! In this case, for an ms = 0 subspace (as the ground state of the Heisenberg model
                 ! will be) then any combination of spins can occur in the subsystem, from all spins
                 ! down to all spins up. Hence the total size of the reduced density matrix will be
                 ! 2**(number of spins in subsystem A).
                 allocate(reduced_density_matrix(2**subsystem_A_size,2**subsystem_A_size), stat=ierr)
                 call check_allocate('reduced_density_matrix', 2**(2*subsystem_A_size),ierr)
                 reduced_density_matrix = 0
             end if
             do i = 1, subsystem_A_size
                 bit_position = bit_lookup(1,subsystem_A_list(i))
                 bit_element = bit_lookup(2,subsystem_A_list(i))
                 subsystem_A_bit_positions(i,1) = bit_position
                 subsystem_A_bit_positions(i,2) = bit_element
             end do
         end if

    end subroutine init_dmqmc

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
                if (btest(f(bit_element),bit_position)) then
                    ! If already flipped up, flip back down.
                    f(bit_element) = ibclr(f(bit_element),bit_position)
                    bits_set = bits_set - 1
                else
                    ! If not flipped up, flip the spin up.
                    f(bit_element) = ibset(f(bit_element),bit_position)
                    bits_set = bits_set + 1
                end if
            end do

            ! Now call a routine to add the corresponding diagonal element to
            ! the spawned walkers list.
            call create_diagonal_particle(f)

        end do

    end subroutine random_distribution_heisenberg

    subroutine create_diagonal_particle(f_new)

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

        integer(i0), intent(in) :: f_new(basis_length)
        integer(i0) :: f_new_diagonal(total_basis_length)
#ifndef PARALLEL
        integer, parameter :: iproc_spawn = 0
#else
        integer :: iproc_spawn
#endif

        ! Create the bitstring of a psip on a diagonal element.
        f_new_diagonal = 0
        f_new_diagonal(:basis_length) = f_new
        f_new_diagonal((basis_length+1):(total_basis_length)) = f_new

#ifdef PARALLEL
        ! Need to determine which processor the spawned walker should be sent to.
        iproc_spawn = modulo(murmurhash_bit_string(f_new_diagonal, &
                                (total_basis_length)), nprocs)
#endif

        ! Move to the next position in the spawning array.
        spawning_head(iproc_spawn) = spawning_head(iproc_spawn) + 1

        ! Set info in spawning array.
        ! Zero it as not all fields are set.
        spawned_walkers(:,spawning_head(iproc_spawn)) = 0
        ! indices 1 to total_basis_length store the bitstring.
        spawned_walkers(:(2*basis_length),spawning_head(iproc_spawn)) = f_new_diagonal
        ! The final index stores the number of psips created, always 1 in this situation.
        spawned_walkers((2*basis_length)+1,spawning_head(iproc_spawn)) = 1

    end subroutine create_diagonal_particle

    subroutine decode_dm_bitstring(f, index1, index2)

        use basis, only: total_basis_length, basis_length
        use fciqmc_data, only: subsystem_A_bit_positions, subsystem_A_bit_positions
        use fciqmc_data, only: subsystem_A_size

        integer(i0), intent(in) :: f(total_basis_length)
        integer(i0), intent(out) :: index1, index2
        integer :: i

        index1 = 0
        index2 = 0

        do i = 1, subsystem_A_size
            if (btest(f(subsystem_A_bit_positions(i,2)),subsystem_A_bit_positions(i,1))) &
                index1 = ibset(index1,i-1)
            if (btest(f(subsystem_A_bit_positions(i,2)+basis_length),subsystem_A_bit_positions(i,1))) &
                index2 = ibset(index2,i-1)
        end do

        index1 = index1+1
        index2 = index2+1

    end subroutine decode_dm_bitstring
    
end module dmqmc_procedures

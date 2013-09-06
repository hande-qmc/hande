module dmqmc_estimators

use const

implicit none

contains

   subroutine update_dmqmc_estimators(ntot_particles_old, ireport, nreplica)

        ! Update the shift and average the shift and estimators.

        ! This is called every report loop in an DMQMC calculation.

        ! Inout:
        !    ntot_particles_old: total number (across all processors) of
        !        particles in the simulation at end of the previous report loop.
        !        Returns the current total number of particles for use in the
        !        next report loop.
        ! In:
        !    ireport: The number of the report loop currently being performed.
        !    nreplica: The total number of replica simulations being performed.

        use spawn_data, only: annihilate_wrapper_spawn_t
        use calc, only: doing_dmqmc_calc, dmqmc_energy, dmqmc_staggered_magnetisation
        use calc, only: dmqmc_energy_squared, dmqmc_renyi_2
        use checking, only: check_allocate
        use dmqmc_procedures, only: rdms
        use energy_evaluation, only: update_shift
        use fciqmc_data, only: nparticles, sampling_size, target_particles, rspawn
        use fciqmc_data, only: shift, vary_shift, nreport, replica_tricks
        use fciqmc_data, only: estimator_numerators, number_dmqmc_estimators
        use fciqmc_data, only: nreport, ncycles, trace, calculate_excit_distribution
        use fciqmc_data, only: excit_distribution, average_shift_until, shift_profile
        use fciqmc_data, only: calc_inst_rdm, rdm_spawn, nrdms, rdm_traces, renyi_2
        use hash_table, only: reset_hash_table
        use parallel

        integer(lint), intent(inout) :: ntot_particles_old(sampling_size)
        integer, intent(in) :: ireport, nreplica
        integer(lint) :: ntot_particles(sampling_size)
        integer :: irdm, ireplica, i

#ifdef PARALLEL
        real(dp), allocatable :: ir(:)
        real(dp), allocatable :: ir_sum(:)
        integer :: ierr, array_size, min_ind, max_ind
#endif

        ! rdm array, and calculate any desired estimators.
        if (calc_inst_rdm) then
            ! WARNING: cannot pass rdm_spawn%spawn to procedures expecting an
            ! array of type spawn_t due to a bug in gfortran which results in
            ! memory deallocations!
            ! See https://groups.google.com/forum/#!topic/comp.lang.fortran/VuFvOsLs6hE
            ! and http://gcc.gnu.org/bugzilla/show_bug.cgi?id=58310.
            ! The explicit loop is also meant to be more efficient anyway, as it
            ! prevents any chance of copy-in/copy-out...
            do irdm = 1, nrdms
                call annihilate_wrapper_spawn_t(rdm_spawn(irdm)%spawn, .false.)
                ! Now is also a good time to reset the hash table (otherwise we
                ! attempt to lookup non-existent data in the next cycle!).
                call reset_hash_table(rdm_spawn(irdm)%ht)
                ! spawn_t comms changes the memory used by spawn%sdata.  Make
                ! sure the hash table always uses the currently 'active'
                ! spawning memory.
                rdm_spawn(irdm)%ht%data_label => rdm_spawn(irdm)%spawn%sdata
            end do
            call calculate_rdm_traces(rdms, rdm_spawn%spawn, rdm_traces)
            if (doing_dmqmc_calc(dmqmc_renyi_2)) call calculate_renyi_2(rdms, rdm_spawn%spawn, renyi_2)
            do irdm = 1, nrdms
                rdm_spawn(irdm)%spawn%head = rdm_spawn(irdm)%spawn%head_start
            end do
        end if


#ifdef PARALLEL
        array_size = 2*sampling_size+1+number_dmqmc_estimators
        if (calculate_excit_distribution) array_size = array_size + size(excit_distribution)
        if (calc_inst_rdm) array_size = array_size + size(rdm_traces(:,1))
        if (doing_dmqmc_calc(dmqmc_renyi_2)) array_size = array_size + size(renyi_2)

        allocate(ir(1:array_size), stat=ierr)
        call check_allocate('ir',array_size,ierr)
        allocate(ir_sum(1:array_size), stat=ierr)
        call check_allocate('ir_sum',array_size,ierr)

        ! Need to sum the number of particles and other quantites over all processors.
        min_ind = 1; max_ind = sampling_size
        ir(min_ind:max_ind) = nparticles
        min_ind = max_ind + 1; max_ind = min_ind
        ir(min_ind) = rspawn
        min_ind = max_ind + 1; max_ind = min_ind + sampling_size - 1
        ir(min_ind:max_ind) = trace
        min_ind = max_ind + 1; max_ind = min_ind + number_dmqmc_estimators - 1
        ir(min_ind:max_ind) = estimator_numerators
        if (calculate_excit_distribution) then
            min_ind = max_ind + 1; max_ind = min_ind + size(excit_distribution) - 1
            ir(min_ind:max_ind) = excit_distribution
        end if
        if (calc_inst_rdm) then
            min_ind = max_ind + 1; max_ind = min_ind + size(rdm_traces) - 1
            ir(min_ind:max_ind) = rdm_traces(:,1)
        end if
        if (calc_inst_rdm) then
            min_ind = max_ind + 1; max_ind = min_ind + size(renyi_2) - 1
            ir(min_ind:max_ind) = renyi_2
        end if

        ! Merge the lists from each processor together.
        call mpi_allreduce(ir, ir_sum, size(ir), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

        min_ind = 1; max_ind = sampling_size
        ntot_particles = nint(ir_sum(min_ind:max_ind))
        min_ind = max_ind + 1; max_ind = min_ind
        rspawn = ir_sum(min_ind)
        min_ind = max_ind + 1; max_ind = min_ind + sampling_size - 1
        trace = nint(ir_sum(min_ind:max_ind))
        min_ind = max_ind + 1; max_ind = min_ind + number_dmqmc_estimators - 1
        estimator_numerators = ir_sum(min_ind:max_ind)
        if (calculate_excit_distribution) then
            min_ind = max_ind + 1; max_ind = min_ind + size(excit_distribution) - 1
            excit_distribution = ir_sum(min_ind:max_ind)
        end if
        if (calc_inst_rdm) then
            min_ind = max_ind + 1; max_ind = min_ind + size(rdm_traces) - 1
            rdm_traces(:,1) = real(ir_sum(min_ind:max_ind),p)
        end if
        if (calc_inst_rdm) then
            min_ind = max_ind + 1; max_ind = min_ind + size(renyi_2) - 1
            renyi_2 = real(ir_sum(min_ind:max_ind),p)
        end if
#else
        ntot_particles = nparticles
#endif

        ! If average_shift_until = -1 then it means that the shift should be updated to
        ! use the values of shift stored in shift_profile. Otherwise, use the standard update
        ! routine.
        if (average_shift_until == -1) then
            if (ireport < nreport) shift = shift_profile(ireport+1)
        else
            if (vary_shift) then
                do ireplica = 1, nreplica
                    call update_shift(shift(ireplica), ntot_particles_old(ireplica), &
                        ntot_particles(ireplica), ncycles)
                end do
            end if
            if (ntot_particles(1) > target_particles .and. .not.vary_shift) vary_shift = .true.
        end if

        ntot_particles_old = ntot_particles
        rspawn = rspawn/(ncycles*nprocs)

   end subroutine update_dmqmc_estimators

   subroutine call_dmqmc_estimators(idet, iteration)

       ! This function calls the processes to update the estimators which
       ! have been requested by the user to be calculated.
       ! First, calculate the excitation level between the two bitsrtings
       ! corresponding to the the two ends. Then add the contribution from
       ! the current density matrix element to the trace, which is always
       ! calculated. Then call other estimators, as required.

       ! In:
       !    idet: Current position in the main bitstring list.
       !    iteration: current Monte Carlo cycle.

       use basis, only: basis_length, total_basis_length
       use calc, only: doing_dmqmc_calc, dmqmc_energy, dmqmc_staggered_magnetisation
       use calc, only: dmqmc_energy_squared, dmqmc_correlation
       use excitations, only: get_excitation, excit
       use fciqmc_data, only: walker_dets, walker_population, trace, doing_reduced_dm
       use fciqmc_data, only: dmqmc_accumulated_probs, start_averaging, dmqmc_find_weights
       use fciqmc_data, only: calculate_excit_distribution, excit_distribution
       use fciqmc_data, only: sampling_size, dmqmc_accumulated_probs_old
       use proc_pointers, only: update_dmqmc_energy_ptr, update_dmqmc_stag_mag_ptr
       use proc_pointers, only: update_dmqmc_energy_squared_ptr, update_dmqmc_correlation_ptr

       integer, intent(in) :: idet, iteration
       type(excit) :: excitation
       real(p) :: unweighted_walker_pop(sampling_size)

       ! Get excitation.
       excitation = get_excitation(walker_dets(:basis_length,idet), &
                        walker_dets((1+basis_length):total_basis_length,idet))

       ! When performing importance sampling the result is that certain excitation
       ! levels have smaller psips populations than the true density matrix by some
       ! factor. In these cases, we want to multiply the psip population by this factor
       ! to calculate the contribution from these excitation levels correctly.

       ! In the case of no importance sampling, unweighted_walker_pop = walker_population(1,idet).
       unweighted_walker_pop = walker_population(:,idet)*dmqmc_accumulated_probs(excitation%nexcit)

       ! The following only use the populations with ireplica = 1, so don't waste time calculating
       ! them if there won't be a contribution.
       if (abs(unweighted_walker_pop(1)) > 0) then
           ! If diagonal element, add to the trace.
           if (excitation%nexcit == 0) trace = trace + walker_population(:,idet)
           ! See which estimators are to be calculated, and call the corresponding procedures.
           ! Energy
           If (doing_dmqmc_calc(dmqmc_energy)) call update_dmqmc_energy_ptr&
                   &(idet, excitation, unweighted_walker_pop(1))
           ! Energy squared
           if (doing_dmqmc_calc(dmqmc_energy_squared)) call update_dmqmc_energy_squared_ptr&
                   &(idet, excitation, unweighted_walker_pop(1))
           ! Spin-spin correlation function
           if (doing_dmqmc_calc(dmqmc_correlation)) call update_dmqmc_correlation_ptr&
                   &(idet, excitation, unweighted_walker_pop(1))
           ! Staggered magnetisation
           if (doing_dmqmc_calc(dmqmc_staggered_magnetisation)) call update_dmqmc_stag_mag_ptr&
                   &(idet, excitation, unweighted_walker_pop(1))
           ! Excitation distribution
           if (calculate_excit_distribution) excit_distribution(excitation%nexcit) = &
                   excit_distribution(excitation%nexcit) + abs(walker_population(1,idet))
           ! Excitation_distribtuion for calculating importance sampling weights
           if (dmqmc_find_weights .and. iteration > start_averaging) excit_distribution(excitation%nexcit) = &
                   excit_distribution(excitation%nexcit) + abs(walker_population(1,idet))
       end if

       ! Reduced density matrix
       if (doing_reduced_dm) call update_reduced_density_matrix_heisenberg&
               &(idet, excitation, walker_population(:,idet), iteration)

       dmqmc_accumulated_probs_old = dmqmc_accumulated_probs

   end subroutine call_dmqmc_estimators

   subroutine dmqmc_energy_heisenberg(idet, excitation, walker_pop)

       ! For the Heisenberg model only.
       ! Add the contribution from the current density matrix element
       ! to the thermal energy estimate.

       ! In:
       !    idet: Current position in the main bitstring (density matrix) list.
       !    excitation: excit type variable which stores information on
       !        the excitation between the two bitstring ends, corresponding
       !        to the two labels for the density matrix element.
       !    walker_pop: number of particles on the current density matrix
       !        element.

       use basis, only: basis_length, total_basis_length
       use basis, only: bit_lookup
       use excitations, only: excit
       use fciqmc_data, only: walker_dets
       use fciqmc_data, only: walker_data, H00
       use fciqmc_data, only: estimator_numerators, energy_index
       use hubbard_real, only: connected_orbs
       use system, only: J_coupling

       integer, intent(in) :: idet
       type(excit), intent(in) :: excitation
       real(p) , intent(in) :: walker_pop
       integer :: bit_element, bit_position

       ! If no excitation, we have a diagonal element, so add elements
       ! which involve the diagonal element of the Hamiltonian.
       if (excitation%nexcit == 0) then
           estimator_numerators(energy_index) = estimator_numerators(energy_index) + &
                                    (walker_data(1,idet)+H00)*walker_pop
       else if (excitation%nexcit == 1) then
       ! If not a diagonal element, but only a single excitation, then the corresponding
       ! Hamiltonian element may be non-zero. Calculate if the flipped spins are
       ! neighbours on the lattice, and if so, add the contirbution from this site.
           bit_position = bit_lookup(1,excitation%from_orb(1))
           bit_element = bit_lookup(2,excitation%from_orb(1))
           if (btest(connected_orbs(bit_element, excitation%to_orb(1)), bit_position)) &
                 estimator_numerators(energy_index) = estimator_numerators(energy_index) - &
                                   (2.0*J_coupling*walker_pop)
       end if

   end subroutine dmqmc_energy_heisenberg

   subroutine dmqmc_energy_squared_heisenberg(idet, excitation, walker_pop)

       ! For the Heisenberg model only.
       ! Add the contribution from the current density matrix element
       ! to the thermal energy squared estimate.

       ! In:
       !    idet: Current position in the main bitstring (density matrix) list.
       !    excitation: excit type variable which stores information on
       !        the excitation between the two bitstring ends, corresponding
       !        to the two labels for the density matrix element.
       !    walker_pop: number of particles on the current density matrix
       !        element.

       use basis, only: basis_length, total_basis_length
       use basis, only: bit_lookup
       use excitations, only: excit
       use fciqmc_data, only: walker_dets
       use fciqmc_data, only: walker_data, H00
       use fciqmc_data, only: estimator_numerators, energy_squared_index
       use hubbard_real, only: connected_orbs, next_nearest_orbs
       use system, only: J_coupling, nbonds

       integer, intent(in) :: idet
       type(excit), intent(in) :: excitation
       real(p) , intent(in) :: walker_pop
       integer :: bit_element1, bit_position1, bit_element2, bit_position2
       real(p) :: sum_H1_H2, J_coupling_squared

       sum_H1_H2 = 0
       J_coupling_squared = J_coupling**2

       if (excitation%nexcit == 0) then
           ! If there are 0 excitations then either nothing happens twice, or we
           ! flip the same pair of spins twice. The Hamiltonian element for doing nothing
           ! is just the diagonal element. For each possible pairs of spins which can be
           ! flipped, there is a mtarix element of -2*J_coupling, so we just need to count
           ! the number of such pairs, which can be found simply from the diagonal element.

           sum_H1_H2 = (walker_data(1,idet)+H00)**2
           sum_H1_H2 = sum_H1_H2 + 2.0*J_coupling_squared*nbonds + 2.0*J_coupling*(walker_data(1,idet)+H00)

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
           ! which way up the intermediate spin is - there will always be one order which contributes.
           ! If there are two such paths, then this could happen by either paths, but again, the two
           ! intermediate spins will only allow one order of spin flipping for each path, no
           ! matter which way up they are, so we only need to check if there are two possible paths.

           if (next_nearest_orbs(excitation%from_orb(1),excitation%to_orb(1)) /= 0) then
               ! Contribution for next-nearest neighbors.
               sum_H1_H2 = 4.0*J_coupling_squared*next_nearest_orbs(excitation%from_orb(1),excitation%to_orb(1))
           end if
           ! Contributions for nearest neighbors.
           ! Note, for certain lattices, such as the triangular lattice, two spins can be both
           ! nearest neighbors *and* next-nearest neighbors. Therefore, it is necessary in general
           ! to check for both situations.
           bit_position1 = bit_lookup(1,excitation%from_orb(1))
           bit_element1 = bit_lookup(2,excitation%from_orb(1))
           if (btest(connected_orbs(bit_element1, excitation%to_orb(1)), bit_position1)) &
                   sum_H1_H2 = sum_H1_H2 - 4.0*J_coupling*(walker_data(1,idet)+H00)

       else if (excitation%nexcit == 2) then
           ! If there are two excitations (4 spins flipped) then, once again, the contribution
           ! to the thermal energy squared will depend on the positions of the spins. If there
           ! are two pairs of spins flipped which are separated then there is one way
           ! for this to happen - by flipping one pair, and then the other (this also requires
           ! that the two spins within each neighboring pair are opposite, as ever for the
           ! Heisenberg model). These two flips can happen in either order.
           ! In some cases the spins may be such that we may pair the spins in more than one way.
           ! For example, if the four spins are in a square shape, or for a 4-by-4 Heisenberg
           ! model, the spins could be connected across the whole lattice, forming a ring due
           ! to the periodic boundaries. In these cases it may be possible to perform the spin
           ! flips by pairing them in either of two ways. To account for this possibility we
           ! have to try and pair the spins in both ways, so we always check both if statements
           ! below. Again, once these pairings have been chosen, the flips can be performed in
           ! either order.

           bit_position1 = bit_lookup(1,excitation%from_orb(1))
           bit_element1 = bit_lookup(2,excitation%from_orb(1))
           bit_position2 = bit_lookup(1,excitation%from_orb(2))
           bit_element2 = bit_lookup(2,excitation%from_orb(2))
           if (btest(connected_orbs(bit_element1, excitation%to_orb(1)), bit_position1) .and. &
           btest(connected_orbs(bit_element2, excitation%to_orb(2)), bit_position2)) &
               sum_H1_H2 = 8.0*J_coupling_squared
           if (btest(connected_orbs(bit_element1, excitation%to_orb(2)), bit_position1) .and. &
           btest(connected_orbs(bit_element2, excitation%to_orb(1)), bit_position2)) &
               sum_H1_H2 = sum_H1_H2 + 8.0*J_coupling_squared

       end if

       estimator_numerators(energy_squared_index) = estimator_numerators(energy_squared_index) + &
               sum_H1_H2*walker_pop

   end subroutine dmqmc_energy_squared_heisenberg

   subroutine dmqmc_energy_hub_real(idet, excitation, walker_pop)

       ! For the Heisenberg model only.
       ! Add the contribution from the current density matrix element
       ! to the thermal energy estimate.

       ! In:
       !    idet: Current position in the main bitstring (density matrix) list.
       !    excitation: excit type variable which stores information on
       !        the excitation between the two bitstring ends, corresponding
       !        to the two labels for the density matrix element.
       !    walker_pop: number of particles on the current density matrix
       !        element.

       use basis, only: basis_length, total_basis_length
       use excitations, only: excit
       use fciqmc_data, only: walker_dets
       use fciqmc_data, only: walker_data, H00
       use fciqmc_data, only: estimator_numerators, energy_index
       use hamiltonian_hub_real, only: slater_condon1_hub_real

       integer, intent(in) :: idet
       type(excit), intent(in) :: excitation
       real(p) , intent(in) :: walker_pop
       real(p) :: hmatel

       ! If no excitation, we have a diagonal element, so add elements
       ! which involve the diagonal element of the Hamiltonian.
       if (excitation%nexcit == 0) then
           estimator_numerators(energy_index) = estimator_numerators(energy_index) + &
                                 (walker_data(1,idet)+H00)*walker_pop
       else if (excitation%nexcit == 1) then
       ! If not a diagonal element, but only a single excitation, then the corresponding
       ! Hamiltonian element may be non-zero. Calculate if the flipped spins are
       ! neighbours on the lattice, and if so, add the contirbution from this site.
           hmatel = slater_condon1_hub_real(excitation%from_orb(1), excitation%to_orb(1), excitation%perm)
           estimator_numerators(energy_index) = estimator_numerators(energy_index) + &
                                 (hmatel*walker_pop)
       end if

   end subroutine dmqmc_energy_hub_real

   subroutine dmqmc_correlation_function_heisenberg(idet, excitation, walker_pop)

       ! For the Heisenberg model only.
       ! Add the contribution from the current density matrix element
       ! to the thermal spin correlation function estimator.

       ! In:
       !    idet: Current position in the main bitstring (density matrix) list.
       !    excitation: excit type variable which stores information on
       !        the excitation between the two bitstring ends, corresponding
       !        to the two labels for the density matrix element.
       !    walker_pop: number of particles on the current density matrix
       !        element.

       use basis, only: basis_length, total_basis_length
       use basis, only: bit_lookup
       use bit_utils, only: count_set_bits
       use excitations, only: excit
       use fciqmc_data, only: walker_dets
       use fciqmc_data, only: walker_data, H00, correlation_mask
       use fciqmc_data, only: estimator_numerators, correlation_index
       use hubbard_real, only: connected_orbs
       use system, only: J_coupling

       integer, intent(in) :: idet
       type(excit), intent(in) :: excitation
       real(p) , intent(in) :: walker_pop
       integer(i0) :: f(basis_length)
       integer :: bit_element1, bit_position1, bit_element2, bit_position2
       integer :: sign_factor

       ! If no excitation, we want the diagonal element of the correlation function operator.
       if (excitation%nexcit == 0) then
           ! If the two spis i and j are the same, the matrix element is +1/4.
           ! If they are different, the matrix element is -1/4.
           ! So we want sign_factor to be +1 in the former case, and -1 in the latter case.
           ! f as calculated below will have 0's at sites other than i and j, and the same values
           ! as walker_dets at i and j. Hence, if f has two 1's or no 1's, we want
           ! sign_factor = +1. Else if we have one 1, we want sign_factor = -1.
           f = iand(walker_dets(:basis_length,idet), correlation_mask)
           ! Count if we have zero, one or two 1's.
           sign_factor = sum(count_set_bits(f))
           ! The operation below will map 0 and 2 to +1, and will map 1 to -1, as is easily checked.
           sign_factor = (mod(sign_factor+1,2)*2)-1
           ! Hence sign_factor can be used to find the matrix element, as used below.
           estimator_numerators(correlation_index) = estimator_numerators(correlation_index) + &
                                    (sign_factor*(walker_pop/4))
       else if (excitation%nexcit == 1) then
       ! If not a diagonal element, but only a single excitation, then the corresponding
       ! matrix element will be 1/2 if and only if the two sites which are flipped are
       ! sites i and j, else it will be 0. We assume that excitations will only be set if
       ! i and j are opposite (else they could not be flipped, for ms=0).
           bit_position1 = bit_lookup(1,excitation%from_orb(1))
           bit_element1 = bit_lookup(2,excitation%from_orb(1))
           bit_position2 = bit_lookup(1,excitation%to_orb(1))
           bit_element2 = bit_lookup(2,excitation%to_orb(1))
           if (btest(correlation_mask(bit_element1), bit_position1).and.&
               btest(correlation_mask(bit_element2), bit_position2)) &
                 estimator_numerators(correlation_index) = &
                     estimator_numerators(correlation_index) + (walker_pop/2)
       end if

   end subroutine dmqmc_correlation_function_heisenberg

   subroutine dmqmc_stag_mag_heisenberg(idet, excitation, walker_pop)

       ! For the Heisenberg model only.
       ! Add the contribution from the current density matrix element
       ! to the thermal staggered magnetisation estimate.

       ! In:
       !    idet: Current position in the main bitstring (density matrix) list.
       !    excitation: excit type variable which stores information on
       !        the excitation between the two bitstring ends, corresponding
       !        to the two labels for the density matrix element.
       !    walker_pop: number of particles on the current density matrix
       !        element.

       use basis, only: basis_length, total_basis_length, bit_lookup
       use bit_utils, only: count_set_bits
       use determinants, only: lattice_mask
       use excitations, only: excit
       use fciqmc_data, only: walker_dets
       use fciqmc_data, only: estimator_numerators, staggered_mag_index
       use system, only: nel, nsites

       integer, intent(in) :: idet
       type(excit), intent(in) :: excitation
       real(p) , intent(in) :: walker_pop
       integer :: bit_element1, bit_position1, bit_element2, bit_position2
       integer(i0) :: f(basis_length)
       integer :: n_up_plus
       integer :: total_sum

       total_sum = 0

       ! This is for the staggered magnetisation squared. This is given by
       ! M^2 = M_x^2 + M_y^2 + M_z^2
       ! where M_x = \sum{i}(-1)^i S_i^x, etc... and S_i^x is the spin operator
       ! at site i, in the x direction.
       if (excitation%nexcit == 0) then
           ! Need to calculate the number of spins up on sublattice 1:
           ! N_u(+) - Number up on + sublattice (where (-1)^i above is +1)
           ! Note the number of down spins on a sublattice is easily obtained from
           ! N_u(+) since there are N/2 spins on each - and the number of spins up on
           ! a different sublattice is easily obtained since there are nel spins
           ! up in total. Hence the matrix element will be written only in terms
           ! of the number of up spins on sublattice 1, to save computation.
           f = iand(walker_dets(:basis_length,idet), lattice_mask)
           n_up_plus = sum(count_set_bits(f))
           ! Below, the term in brackets and middle term come from the z component (the
           ! z operator is diagonal) and one nsites/4 factor comes from the x operator,
           ! the other nsites/4 factor from the y operator.
           total_sum = (2*n_up_plus-nel)**2 + (nsites/2)
       else if (excitation%nexcit == 1) then
           ! Off-diagonal elements from the y and z operators. For the pair of spins
           ! that are flipped, if they are on the same sublattice, we get a factor of
           ! 1, or if on different sublattices, a factor of -1.
           bit_position1 = bit_lookup(1,excitation%from_orb(1))
           bit_element1 = bit_lookup(2,excitation%from_orb(1))
           bit_position2 = bit_lookup(1,excitation%to_orb(1))
           bit_element2 = bit_lookup(2,excitation%to_orb(1))
           if (btest(lattice_mask(bit_element1), bit_position1)) total_sum = total_sum+1
           if (btest(lattice_mask(bit_element2), bit_position2)) total_sum = total_sum+1
           ! The operation below will map 0 and 2 to +1, and will map 1 to -1, as is easily checked.
           ! We want this - if both or no spins on this sublattice, then both on same sublattice
           ! either way, so plus one. Else they are on different sublattices, so we want a factor
           ! of -1, as we get.
           total_sum = (mod(total_sum+1,2)*2)-1
       end if

       estimator_numerators(staggered_mag_index) = estimator_numerators(staggered_mag_index) + &
                                  (real(total_sum)/real(nsites**2))*walker_pop

   end subroutine dmqmc_stag_mag_heisenberg

   subroutine update_reduced_density_matrix_heisenberg(idet, excitation, walker_pop, iteration)

       ! Add a contribution from the current walker to the reduced density
       ! matrix estimator, which is produced by 'tracing out' a given set of
       ! spins.

       ! Applicable only to the Heisenberg model.

       ! This procedure takes the two determinants bitstrings for the current
       ! walker and, if the two bitstrings of the B subsystem are identical,
       ! adds the walker population to the corresponding reduced density matrix
       ! element.

       ! In:
       !    idet: Current position in the main bitstring (density matrix) list.
       !    walker_pop: number of particles on the current density matrix
       !        element. Note that this walker population is still weighted
       !        by the importance sampling factors. These factors must be
       !        removed before any estimates can be calculated.

       use basis, only: basis_length, total_basis_length
       use dmqmc_procedures, only: decode_dm_bitstring, rdms
       use excitations, only: excit
       use fciqmc_data, only: reduced_density_matrix, walker_dets, walker_population
       use fciqmc_data, only: sampling_size, calc_inst_rdm, calc_ground_rdm, nrdms
       use fciqmc_data, only: nsym_vec, start_averaging, rdm_spawn, dmqmc_accumulated_probs
       use spawning, only: create_spawned_particle_rdm

       integer, intent(in) :: idet, iteration
       integer, intent(in) :: walker_pop(sampling_size)
       type(excit), intent(in) :: excitation
       real(p) :: unweighted_walker_pop(sampling_size)
       integer :: irdm, isym, ireplica
       integer(i0) :: f1(basis_length), f2(basis_length)

       if (.not. (iteration > start_averaging .or. calc_inst_rdm)) return

       ! Loop over all RDMs to be calculated.
       do irdm = 1, nrdms
           ! Loop over every symmetry-equivalent subsystem for this RDM.
           do isym = 1, nsym_vec

               ! Apply the mask for the B subsystem to set all sites in the A subsystem to 0.
               f1 = iand(rdms(irdm)%B_masks(isym,:),walker_dets(:basis_length,idet))
               f2 = iand(rdms(irdm)%B_masks(isym,:),walker_dets(basis_length+1:total_basis_length,idet))

               ! Once this is done, check if the resulting bitstring (which can only possibly
               ! have 1's in the B subsystem) are identical. If they are, then this psip gives
               ! a contibution to the reduced density matrix for subsystem A. This is because we
               ! get the reduced density matrix for A by 'tracing out' over B, which in practice
               ! means only keeping matrix elements that are on the diagonal for subsystem B.
               if (sum(abs(f1-f2)) == 0) then
                   ! Call a function which maps the subsystem A state to two RDM bitstrings.
                   call decode_dm_bitstring(walker_dets(:,idet),irdm,isym)

                   if (calc_ground_rdm) then
                       ! The above routine actually maps to numbers between 0 and 2^rdms(1)%A_nsites-1,
                       ! but the smallest and largest reduced density matrix indices are one more than
                       ! these, so add one.
                       rdms(irdm)%end1 = rdms(irdm)%end1 + 1
                       rdms(irdm)%end2 = rdms(irdm)%end2 + 1
                       unweighted_walker_pop = walker_population(:,idet)*dmqmc_accumulated_probs(excitation%nexcit)
                       ! Note, when storing the entire RDM (as done here), the maximum value of
                       ! rdms(i)%rdm_basis_length is 1, so we only consider this one element here.
                       reduced_density_matrix(rdms(irdm)%end1(1),rdms(irdm)%end2(1)) = &
                           reduced_density_matrix(rdms(irdm)%end1(1),rdms(irdm)%end2(1)) + unweighted_walker_pop(1)
                   end if

                   if (calc_inst_rdm) then
                      do ireplica = 1, sampling_size
                          if (abs(walker_pop(ireplica)) > 0) then
                              call create_spawned_particle_rdm(irdm, walker_pop(ireplica), &
                                      ireplica, rdm_spawn(irdm))
                          end if
                      end do
                   end if

               end if
           end do
       end do

    end subroutine update_reduced_density_matrix_heisenberg

    subroutine call_ground_rdm_procedures(beta_cycle)

        ! Wrapper for calling relevant reduced density matrix procedures.

        use checking, only: check_allocate, check_deallocate
        use dmqmc_procedures, only: rdms
        use fciqmc_data, only: reduced_density_matrix
        use fciqmc_data, only: doing_von_neumann_entropy, doing_concurrence
        use fciqmc_data, only: output_rdm, rdm_unit, rdm_traces
        use parallel
        use utils, only: get_free_unit, append_ext, int_fmt

        integer, intent(in) :: beta_cycle
        real(p), allocatable :: old_rdm_elements(:)
        integer :: i, j, k, ierr, new_unit
        character(10) :: rdm_filename
        ! If in paralell then merge the reduced density matrix onto one processor.
#ifdef PARALLEL

        real(dp), allocatable :: dm(:,:)
        real(dp), allocatable :: dm_sum(:,:)
        integer :: num_eigv

        num_eigv = 2**rdms(1)%A_nsites

        allocate(dm(num_eigv,num_eigv), stat=ierr)
        call check_allocate('dm',num_eigv**2,ierr)
        allocate(dm_sum(num_eigv,num_eigv), stat=ierr)
        call check_allocate('dm_sum',num_eigv**2,ierr)

        dm = reduced_density_matrix

        call mpi_allreduce(dm, dm_sum, size(dm), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

        reduced_density_matrix = dm_sum

        deallocate(dm)
        call check_deallocate('dm',ierr)
        deallocate(dm_sum)
        call check_deallocate('dm_sum',ierr)

#endif

        rdm_traces = 0.0_p

        if (parent) then
            ! Force the reduced desnity matrix to be symmetric by averaging the upper and
            ! lower triangles.
            do i = 1, ubound(reduced_density_matrix,1)
                do j = 1, i-1
                    reduced_density_matrix(i,j) = 0.5_p*(reduced_density_matrix(i,j) +&
                            reduced_density_matrix(j,i))
                    reduced_density_matrix(j,i) = reduced_density_matrix(i,j)
                end do
                ! Add current contirbution to the trace.
                rdm_traces(1,1) = rdm_traces(1,1) + reduced_density_matrix(i,i)
            end do

            ! Call the routines to calculate the desired quantities.
            if (doing_von_neumann_entropy) call calculate_vn_entropy(rdm_traces(1,1))
            if (doing_concurrence) call calculate_concurrence()

            write (6,'(1x,a12,1X,f22.12)') "# RDM trace=", rdm_traces(1,1)

            if (output_rdm) then
                new_unit = get_free_unit()
                call append_ext('rdm', beta_cycle, rdm_filename)
                open(new_unit, file=trim(rdm_filename), status='replace')
                write(new_unit,'(a5,1x,es15.8)') "Trace", rdm_traces(1,1)
                do i = 1, ubound(reduced_density_matrix,1)
                    do j = i, ubound(reduced_density_matrix,1)
                        write(new_unit,'(a1,'//int_fmt(i,0)//',a1,'//int_fmt(j,0)//',a1,1x,es15.8)') &
                            "(",i,",",j,")", reduced_density_matrix(i,j)
                    end do
                end do
                close(new_unit)
            end if

        end if

    end subroutine call_ground_rdm_procedures

    subroutine calculate_vn_entropy(trace_rdm)

        ! Calculate the Von Neumann Entropy. Use lapack to calculate the
        ! eigenvalues {\lambda_j} of the reduced density matrix.
        ! Then VN Entropy S = -\sum_j\lambda_j\log_2{\lambda_j}.

        ! Need to paralellise for large subsystems and introduce test to
        ! check whether diagonalisation should be performed in serial or
        ! paralell.

        use checking, only: check_allocate, check_deallocate
        use dmqmc_procedures, only: rdms
        use fciqmc_data, only: reduced_density_matrix

        real(p), intent(in) :: trace_rdm
        integer :: i, rdm_size
        integer :: info, ierr, lwork
        real(p), allocatable :: work(:)
        real(p), allocatable :: dm_tmp(:,:)
        real(p) :: eigv(2**rdms(1)%A_nsites)
        real(p) :: vn_entropy
        
        rdm_size = 2**rdms(1)%A_nsites
        vn_entropy = 0._p

        ! Find the optimal size of the workspace.
        allocate(work(1), stat=ierr)
        call check_allocate('work',1,ierr)
#ifdef SINGLE_PRECISION
        call ssyev('N', 'U', rdm_size, reduced_density_matrix, rdm_size, eigv, work, -1, info)
#else
        call dsyev('N', 'U', rdm_size, reduced_density_matrix, rdm_size, eigv, work, -1, info)
#endif

        lwork = nint(work(1))
        deallocate(work)
        call check_deallocate('work',ierr)

        ! Now perform the diagonalisation.
        allocate(work(lwork), stat=ierr)
        call check_allocate('work',lwork,ierr)

        ! The matrix input into the following diagonalisation routines will have their upper half
        ! (including the diagonal) destroyed. We might want reduced_desntiy_matrix later, so
        ! use some temporary space:
        allocate(dm_tmp(rdm_size,rdm_size), stat=ierr)
        call check_allocate('dm_tmp',rdm_size**2,ierr)
        dm_tmp = reduced_density_matrix

#ifdef SINGLE_PRECISION
        call ssyev('N', 'U', rdm_size, dm_tmp, rdm_size, eigv, work, lwork, info)
#else
        call dsyev('N', 'U', rdm_size, dm_tmp, rdm_size, eigv, work, lwork, info)
#endif
        write(6,'(a24)',advance='no') "Eigenvalues thrown away:"
        do i = 1, ubound(eigv,1)
            if (eigv(i) < 0.0_p) then
                write(6,'(es15.8,2x)',advance='no') eigv(i)/trace_rdm
                cycle
            end if
            vn_entropy = vn_entropy - eigv(i)*(log(eigv(i))/log(2.0_p))
        end do
        write(6,'()',advance='yes')
        write (6,'(1x,a36,1X,f22.9)') "# Unnormalised von Neumann entropy= ", vn_entropy

        deallocate(dm_tmp)
        call check_deallocate('dm_tmp',ierr)

    end subroutine calculate_vn_entropy
    
    subroutine calculate_concurrence()

        ! Calculate the concurrence of a qubit. For a reduced density matrix \rho,
        ! the concurrence, C =  max(0, \lamda_1 - \lambda_2 - \lambda_3 -\lambda_4) where
        ! \lambda_i are the eigenvalues of the matrix, R = \sqrt{\sqrt{\rho}\~{\rho}\sqrt{\rho}}.
        ! Where \~\rho = {\sigma_y \otimes \sigma_y} \rho^{\ast} {\sigma_y \otimes \sigma_y}. 
        ! \lambda_1 > ... > \lambda_4.

        ! This can be simplified to finding the square root of the eigenvalues \{\lambda_i\} of \rho\~{\rho} 
        ! and in the case where \rho is a real, symmetric matrix then we can further simplify the problem
        ! to finding the eigenvalues of R = \rho \sigma_y \otimes \sigma_y.

        ! Below we have named {\sigma_y \otimes \sigma_y} flip_spin_matrix as in the literature.

        use checking, only: check_allocate, check_deallocate
        use fciqmc_data, only: reduced_density_matrix, flip_spin_matrix
        integer :: info, ierr, lwork
        real(p), allocatable :: work(:)
        real(p) :: reigv(4), ieigv(4)
        real(p) :: concurrence
        real(p) :: rdm_spin_flip(4,4), rdm_spin_flip_tmp(4,4), VL(4,4), VR(4,4)
        
        ! Make rdm_spin_flip_tmp because sgeev and dgeev delete input matrix.
        rdm_spin_flip_tmp = matmul(reduced_density_matrix, flip_spin_matrix)
        rdm_spin_flip = rdm_spin_flip_tmp
        
        ! Find the optimal size of the workspace.
        allocate(work(1), stat=ierr)
        call check_allocate('work',1,ierr)
#ifdef SINGLE_PRECISION
        call sgeev('N', 'N', 4, rdm_spin_flip_tmp, 4, reigv, ieigv, VL, 1, VR, 1, work, -1, info)
#else
        call dgeev('N', 'N', 4, rdm_spin_flip_tmp, 4, reigv, ieigv, VL, 1, VR, 1,  work, -1, info)
#endif
        lwork = nint(work(1))
        deallocate(work)
        call check_deallocate('work',ierr)

        ! Now perform the diagonalisation.
        allocate(work(lwork), stat=ierr)
        call check_allocate('work',lwork,ierr)
#ifdef SINGLE_PRECISION
        call sgeev('N', 'N', 4, rdm_spin_flip, 4, reigv, ieigv, VL, 1, VR, 1, work, lwork, info)
#else
        call dgeev('N', 'N', 4, rdm_spin_flip, 4, reigv, ieigv, VL, 1, VR, 1, work, lwork, info)
#endif
        ! Calculate the concurrence. Take abs of eigenvalues so that this is equivelant to sqauring
        ! and then square-rooting.
        concurrence = 2._p*maxval(abs(reigv)) - sum(abs(reigv)) 
        concurrence = max(0._p, concurrence)
        write (6,'(1x,a28,1X,f22.12)') "# Unnormalised concurrence= ", concurrence

    end subroutine calculate_concurrence
    
    subroutine calculate_rdm_traces(rdm_data, rdm_lists, traces)

        use dmqmc_procedures, only: rdm
        use excitations, only: get_excitation_level
        use spawn_data, only: spawn_t

        type(rdm), intent(in) :: rdm_data(:)
        type(spawn_t), intent(in) :: rdm_lists(:)
        real(p) :: traces(:,:)

        integer :: irdm, i, rdm_bl
        integer, parameter :: thread_id = 0

        traces = 0.0_p 

        do irdm = 1, size(rdm_data)
            rdm_bl = rdm_data(irdm)%rdm_basis_length
            do i = 1, rdm_lists(irdm)%head(thread_id,0)

                if (all( rdm_lists(irdm)%sdata(1:rdm_bl,i) == rdm_lists(irdm)%sdata(rdm_bl+1:2*rdm_bl,i))) then
                    traces(:,irdm) = traces(:,irdm) + &
                        real(rdm_lists(irdm)%sdata(rdm_lists(irdm)%bit_str_len+1:rdm_lists(irdm)%element_len,i),p)
                end if
            end do
        end do

    end subroutine calculate_rdm_traces

    subroutine calculate_renyi_2(rdm_data, rdm_lists, r2)

        use dmqmc_procedures, only: rdm
        use excitations, only: get_excitation_level
        use fciqmc_data, only: dmqmc_accumulated_probs_old
        use spawn_data, only: spawn_t

        type(rdm), intent(in) :: rdm_data(:)
        type(spawn_t), intent(in) :: rdm_lists(:)
        real(p) :: r2(:)
        integer :: i, irdm, excit_level, rdm_bl
        real(p) :: unweighted_pop_1, unweighted_pop_2
        integer, parameter :: thread_id = 0

        r2 = 0.0_p 

        do irdm = 1, size(rdm_data)
            rdm_bl = rdm_data(irdm)%rdm_basis_length
            do i = 1, rdm_lists(irdm)%head(thread_id,0)

                excit_level = get_excitation_level(rdm_lists(irdm)%sdata(1:rdm_bl,i),&
                    rdm_lists(irdm)%sdata(rdm_bl+1:2*rdm_bl,i))

                unweighted_pop_1 = rdm_lists(irdm)%sdata(rdm_lists(irdm)%bit_str_len+1,i)*dmqmc_accumulated_probs_old(excit_level)
                unweighted_pop_2 = rdm_lists(irdm)%sdata(rdm_lists(irdm)%bit_str_len+2,i)*dmqmc_accumulated_probs_old(excit_level)

                if (excit_level == 0) then
                    r2(irdm) = r2(irdm) + unweighted_pop_1*unweighted_pop_2
                else
                    ! As we only hold RDM elements above the diagonal, off-diagonal elements must
                    ! be counted twice.
                    r2(irdm) = r2(irdm) + 2*unweighted_pop_1*unweighted_pop_2
                end if

            end do
        end do

    end subroutine calculate_renyi_2

end module dmqmc_estimators

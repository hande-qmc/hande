module dmqmc_estimators

use const

implicit none

contains

   subroutine communicate_dmqmc_estimates()

        ! Sum together the contributions to the various DMQMC estimators (and
        ! some other non-physical quantities such as the rate of spawning and
        ! total number of walkers) across all MPI processes.

        ! This is called every report loop in a DMQMC calculation.

        use spawn_data, only: annihilate_wrapper_spawn_t
        use calc, only: doing_dmqmc_calc, dmqmc_energy, dmqmc_staggered_magnetisation
        use calc, only: dmqmc_energy_squared, dmqmc_rdm_r2
        use checking, only: check_allocate
        use dmqmc_procedures, only: rdms
        use fciqmc_data, only: nparticles, sampling_size, rspawn, shift, replica_tricks
        use fciqmc_data, only: estimator_numerators, number_dmqmc_estimators
        use fciqmc_data, only: ncycles, trace, calculate_excit_distribution
        use fciqmc_data, only: excit_distribution, tot_nparticles
        use fciqmc_data, only: calc_inst_rdm, rdm_spawn, nrdms, rdm_traces, renyi_2
        use hash_table, only: reset_hash_table
        use parallel

        integer :: irdm

#ifdef PARALLEL
        real(dp), allocatable :: ir(:)
        real(dp), allocatable :: ir_sum(:)
        integer :: ierr, array_size, min_ind, max_ind
#endif

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
            if (doing_dmqmc_calc(dmqmc_rdm_r2)) call calculate_rdm_renyi_2(rdms, rdm_spawn%spawn, renyi_2)
            do irdm = 1, nrdms
                rdm_spawn(irdm)%spawn%head = rdm_spawn(irdm)%spawn%head_start
            end do
        end if

#ifdef PARALLEL
        ! Put all the quantities to be communicated together in one array.

        array_size = 2*sampling_size+1+number_dmqmc_estimators
        if (calculate_excit_distribution) array_size = array_size + size(excit_distribution)
        if (calc_inst_rdm) array_size = array_size + size(rdm_traces)
        if (doing_dmqmc_calc(dmqmc_rdm_r2)) array_size = array_size + size(renyi_2)

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
            do irdm = 1, nrdms
                min_ind = max_ind + 1; max_ind = min_ind + size(rdm_traces(:,irdm)) - 1
                ir(min_ind:max_ind) = rdm_traces(:,irdm)
            end do
        end if
        if (doing_dmqmc_calc(dmqmc_rdm_r2)) then
            min_ind = max_ind + 1; max_ind = min_ind + size(renyi_2) - 1
            ir(min_ind:max_ind) = renyi_2
        end if

        ! Sum the data from each processor together.
        call mpi_allreduce(ir, ir_sum, size(ir), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

        ! Extract the summed data from the combined array.
        min_ind = 1; max_ind = sampling_size
        tot_nparticles = nint(ir_sum(min_ind:max_ind))
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
            do irdm = 1, nrdms
                min_ind = max_ind + 1; max_ind = min_ind + size(rdm_traces(:,irdm)) - 1
                rdm_traces(:,irdm) = real(ir_sum(min_ind:max_ind),p)
            end do
        end if
        if (doing_dmqmc_calc(dmqmc_rdm_r2)) then
            min_ind = max_ind + 1; max_ind = min_ind + size(renyi_2) - 1
            renyi_2 = real(ir_sum(min_ind:max_ind),p)
        end if
#else
        tot_nparticles = nparticles
#endif

        rspawn = rspawn/(ncycles*nprocs)

   end subroutine communicate_dmqmc_estimates

   subroutine update_shift_dmqmc(loc_tot_nparticles, loc_tot_nparticles_old, ireport)

        ! In/Out:
        !    loc_tot_nparticles: total number (across all processors) of
        !        particles in the simulation at end of the previous report loop.
        !    loc_tot_nparticles_old: total number (across all processors) of
        !        particles in the simulation currently.
        ! In:
        !    ireport: The number of the report loop currently being performed.

        use energy_evaluation, only: update_shift
        use fciqmc_data, only: shift, shift_profile, average_shift_until, ncycles
        use fciqmc_data, only: target_particles, vary_shift, nreport

        real(dp), intent(in) :: loc_tot_nparticles(:)
        real(dp), intent(in) :: loc_tot_nparticles_old(:)
        integer, intent(in) :: ireport
        integer :: ireplica

        ! If average_shift_until = -1 then it means that the shift should be
        ! updated to use the values of shift stored in shift_profile. Otherwise,
        ! use the standard update routine.
        if (average_shift_until == -1) then
            if (ireport < nreport) shift = shift_profile(ireport+1)
        else
            if (vary_shift) then
                do ireplica = 1, size(loc_tot_nparticles)
                    call update_shift(shift(ireplica), loc_tot_nparticles_old(ireplica), &
                        loc_tot_nparticles(ireplica), ncycles)
                end do
            end if
            if (loc_tot_nparticles(1) > target_particles .and. (.not. vary_shift)) &
                vary_shift = .true.
        end if

   end subroutine update_shift_dmqmc

   subroutine update_dmqmc_estimators(sys, idet, iteration)

       ! This function calls the processes to update the estimators which have
       ! been requested by the user to be calculated. First, calculate the
       ! excitation level between the two bitstrings corresponding to the the
       ! two ends. Then add the contribution from the current density matrix
       ! element to the trace, which is always calculated. Then call other
       ! estimators, as required.

       ! In:
       !    sys: system being studied.
       !    idet: Current position in the main bitstring list.
       !    iteration: current Monte Carlo cycle.

       use calc, only: doing_dmqmc_calc, dmqmc_energy, dmqmc_staggered_magnetisation
       use calc, only: dmqmc_energy_squared, dmqmc_correlation, dmqmc_full_r2
       use excitations, only: get_excitation, excit
       use fciqmc_data, only: walker_dets, walker_population, trace, doing_reduced_dm
       use fciqmc_data, only: dmqmc_accumulated_probs, start_averaging, dmqmc_find_weights
       use fciqmc_data, only: calculate_excit_distribution, excit_distribution
       use fciqmc_data, only: sampling_size, dmqmc_accumulated_probs_old, real_factor
       use proc_pointers, only: update_dmqmc_energy_ptr, update_dmqmc_stag_mag_ptr
       use proc_pointers, only: update_dmqmc_energy_squared_ptr, update_dmqmc_correlation_ptr
       use system, only: sys_t

       type(sys_t), intent(in) :: sys
       integer, intent(in) :: idet, iteration
       type(excit) :: excitation
       real(p) :: unweighted_walker_pop(sampling_size)

       ! Get excitation.
       excitation = get_excitation(sys%nel, walker_dets(:sys%basis%basis_length,idet), &
                        walker_dets((1+sys%basis%basis_length):sys%basis%total_basis_length,idet))

       ! When performing importance sampling the result is that certain
       ! excitation levels have smaller psips populations than the true density
       ! matrix by some factor. In these cases, we want to multiply the psip
       ! population by this factor to calculate the contribution from these
       ! excitation levels correctly.

       ! In the case of no importance sampling, unweighted_walker_pop = walker_population(1,idet).
       unweighted_walker_pop = real(walker_population(:,idet),dp)*dmqmc_accumulated_probs(excitation%nexcit)/&
                                real_factor

       ! If diagonal element, add to the trace.
       if (excitation%nexcit == 0) trace = trace + real(walker_population(:,idet),p)/real_factor

       ! The following only use the populations with ireplica = 1, so only call
       ! them if the determinant is occupied in the first replica.
       if (abs(unweighted_walker_pop(1)) > 0) then
           ! See which estimators are to be calculated, and call the
           ! corresponding procedures.
           ! Energy
           If (doing_dmqmc_calc(dmqmc_energy)) call update_dmqmc_energy_ptr&
                   &(sys, idet, excitation, unweighted_walker_pop(1))
           ! Energy squared.
           if (doing_dmqmc_calc(dmqmc_energy_squared)) call update_dmqmc_energy_squared_ptr&
                   &(sys, idet, excitation, unweighted_walker_pop(1))
           ! Spin-spin correlation function.
           if (doing_dmqmc_calc(dmqmc_correlation)) call update_dmqmc_correlation_ptr&
                   &(sys, idet, excitation, unweighted_walker_pop(1))
           ! Staggered magnetisation.
           if (doing_dmqmc_calc(dmqmc_staggered_magnetisation)) call update_dmqmc_stag_mag_ptr&
                   &(sys, idet, excitation, unweighted_walker_pop(1))
           ! Excitation distribution.
           if (calculate_excit_distribution) excit_distribution(excitation%nexcit) = &
                   excit_distribution(excitation%nexcit) + real(abs(walker_population(1,idet)),p)/real_factor
           ! Excitation distribtuion for calculating importance sampling weights.
           if (dmqmc_find_weights .and. iteration > start_averaging) excit_distribution(excitation%nexcit) = &
                   excit_distribution(excitation%nexcit) + real(abs(walker_population(1,idet)),p)/real_factor
       end if

       ! Full Renyi entropy (S_2).
       if (doing_dmqmc_calc(dmqmc_full_r2)) call update_full_renyi_2(unweighted_walker_pop)

       ! Reduced density matrices.
       if (doing_reduced_dm) call update_reduced_density_matrix_heisenberg&
               &(idet, excitation, walker_population(:,idet), iteration)

       dmqmc_accumulated_probs_old = dmqmc_accumulated_probs

   end subroutine update_dmqmc_estimators

   subroutine dmqmc_energy_heisenberg(sys, idet, excitation, walker_pop)

       ! For the Heisenberg model only.
       ! Add the contribution from the current density matrix element to the 
       ! thermal energy estimate.

       ! In:
       !    sys: system being studied.
       !    idet: Current position in the main bitstring (density matrix) list.
       !    excitation: excit type variable which stores information on
       !        the excitation between the two bitstring ends, corresponding
       !        to the two labels for the density matrix element.
       !    walker_pop: number of particles on the current density matrix
       !        element.

       use excitations, only: excit
       use fciqmc_data, only: walker_dets
       use fciqmc_data, only: walker_data, H00
       use fciqmc_data, only: estimator_numerators, energy_index
       use real_lattice, only: connected_orbs
       use system, only: sys_t

       type(sys_t), intent(in) :: sys
       integer, intent(in) :: idet
       type(excit), intent(in) :: excitation
       real(p), intent(in) :: walker_pop
       integer :: bit_element, bit_position

       ! If no excitation, we have a diagonal element, so add elements which
       ! involve the diagonal element of the Hamiltonian.
       if (excitation%nexcit == 0) then
           estimator_numerators(energy_index) = estimator_numerators(energy_index) + &
                                    (walker_data(1,idet)+H00)*walker_pop
       else if (excitation%nexcit == 1) then
       ! If not a diagonal element, but only a single excitation, then the
       ! corresponding Hamiltonian element may be non-zero. Calculate if the
       ! flipped spins are neighbours on the lattice, and if so, add the
       ! contribution from this site.
           bit_position = sys%basis%bit_lookup(1,excitation%from_orb(1))
           bit_element = sys%basis%bit_lookup(2,excitation%from_orb(1))
           if (btest(connected_orbs(bit_element, excitation%to_orb(1)), bit_position)) &
                 estimator_numerators(energy_index) = estimator_numerators(energy_index) - &
                                   (sys%heisenberg%J*walker_pop/2)
       end if

   end subroutine dmqmc_energy_heisenberg

   subroutine dmqmc_energy_squared_heisenberg(sys, idet, excitation, walker_pop)

       ! For the Heisenberg model only.
       ! Add the contribution from the current density matrix element to the
       ! thermal energy squared estimate.

       ! In:
       !    sys: system being studied.
       !    idet: Current position in the main bitstring (density matrix) list.
       !    excitation: excit type variable which stores information on the
       !        excitation between the two bitstring ends, corresponding to the
       !        two labels for the density matrix element.
       !    walker_pop: number of particles on the current density matrix
       !        element.

       use excitations, only: excit
       use fciqmc_data, only: walker_dets
       use fciqmc_data, only: walker_data, H00
       use fciqmc_data, only: estimator_numerators, energy_squared_index
       use real_lattice, only: connected_orbs, next_nearest_orbs
       use system, only: sys_t

       type(sys_t), intent(in) :: sys
       integer, intent(in) :: idet
       type(excit), intent(in) :: excitation
       real(p), intent(in) :: walker_pop
       integer :: bit_element1, bit_position1, bit_element2, bit_position2
       real(p) :: sum_H1_H2, J_coupling_squared

       sum_H1_H2 = 0.0_p
       J_coupling_squared = sys%heisenberg%J**2/16

       if (excitation%nexcit == 0) then
           ! If there are 0 excitations then either nothing happens twice, or we
           ! flip the same pair of spins twice. The Hamiltonian element for doing
           ! nothing is just the diagonal element. For each possible pairs of
           ! spins which can be flipped, there is a mtarix element of -J/2, so
           ! we just need to count the number of such pairs, which can be found 
           ! simply from the diagonal element.

           sum_H1_H2 = (walker_data(1,idet)+H00)**2
           associate(sh=>sys%heisenberg)
               sum_H1_H2 = sum_H1_H2 + 2*J_coupling_squared*sh%nbonds + sh%J*(walker_data(1,idet)+H00)/2
           end associate

       else if (excitation%nexcit == 1) then
           ! If there is only one excitation (2 spins flipped) then the
           ! contribution to H^2 depend on the positions of the spins relative
           ! to one another. If the the spins are nearest neighbors then we
           ! could either do nothing and then flip the pair, or flip the pair
           ! and then do nothing. If next nearest neighbors and there is only
           ! one two-bond path to get from one spin to the other, we first flip
           ! the pair on the first bond, then flip the pair on the second bond.
           ! This flipping can only be done in exactly one order, not both - the
           ! two spins which change are opposite, so the middle spin will
           ! initially only be the same as one or the other spin. This is nice, 
           ! because we don't have check which way up the intermediate spin is -
           ! there will always be one order which contributes. If there are two
           ! such paths, then this could happen by either paths, but again, the
           ! two intermediate spins will only allow one order of spin flipping
           ! for each path, no matter which way up they are, so we only need to
           ! check if there are two possible paths.

           if (next_nearest_orbs(excitation%from_orb(1),excitation%to_orb(1)) /= 0_i0) then
               ! Contribution for next-nearest neighbors.
               sum_H1_H2 = 4.0_p*J_coupling_squared*next_nearest_orbs(excitation%from_orb(1),excitation%to_orb(1))
           end if
           ! Contributions for nearest neighbors.
           ! Note, for certain lattices, such as the triangular lattice, two
           ! spins can be both nearest neighbors *and* next-nearest neighbors.
           ! Therefore, it is necessary in general to check for both situations.
           bit_position1 = sys%basis%bit_lookup(1,excitation%from_orb(1))
           bit_element1 = sys%basis%bit_lookup(2,excitation%from_orb(1))
           if (btest(connected_orbs(bit_element1, excitation%to_orb(1)), bit_position1)) &
                   sum_H1_H2 = sum_H1_H2 - sys%heisenberg%J*(walker_data(1,idet)+H00)

       else if (excitation%nexcit == 2) then
           ! If there are two excitations (4 spins flipped) then, once again,
           ! the contribution to the thermal energy squared will depend on the
           ! positions of the spins. If there are two pairs of spins flipped
           ! which are separated then there is one way for this to happen - by
           ! flipping one pair, and then the other (this also requires that the
           ! two spins within each neighboring pair are opposite, as ever for
           ! the Heisenberg model). These two flips can happen in either order.
           ! In some cases the spins may be such that we may pair the spins in
           ! more than one way. For example, if the four spins are in a square
           ! shape, or for a 4-by-4 Heisenberg model, the spins could be
           ! connected across the whole lattice, forming a ring due to the
           ! periodic boundaries. In these cases it may be possible to perform
           ! the spin flips by pairing them in either of two ways. To account
           ! for this possibility we have to try and pair the spins in both
           ! ways, so we always check both if statements below. Again, once
           ! these pairings have been chosen, the flips can be performed in
           ! either order.

           bit_position1 = sys%basis%bit_lookup(1,excitation%from_orb(1))
           bit_element1 = sys%basis%bit_lookup(2,excitation%from_orb(1))
           bit_position2 = sys%basis%bit_lookup(1,excitation%from_orb(2))
           bit_element2 = sys%basis%bit_lookup(2,excitation%from_orb(2))
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

   subroutine dmqmc_energy_hub_real(sys, idet, excitation, walker_pop)

       ! Add the contribution from the current density matrix element to the
       ! thermal energy estimate.

       ! In:
       !    sys: system being studied.
       !    idet: Current position in the main bitstring (density matrix) list.
       !    excitation: excit type variable which stores information on
       !        the excitation between the two bitstring ends, corresponding to
       !        the two labels for the density matrix element.
       !    walker_pop: number of particles on the current density matrix
       !        element.

       use excitations, only: excit
       use fciqmc_data, only: walker_dets
       use fciqmc_data, only: walker_data, H00
       use fciqmc_data, only: estimator_numerators, energy_index
       use hamiltonian_hub_real, only: slater_condon1_hub_real
       use system, only: sys_t

       type(sys_t), intent(in) :: sys
       integer, intent(in) :: idet
       type(excit), intent(in) :: excitation
       real(p), intent(in) :: walker_pop
       real(p) :: hmatel

       ! If no excitation, we have a diagonal element, so add elements which
       ! involve the diagonal element of the Hamiltonian.
       if (excitation%nexcit == 0) then
           estimator_numerators(energy_index) = estimator_numerators(energy_index) + &
                                 (walker_data(1,idet)+H00)*walker_pop
       else if (excitation%nexcit == 1) then
       ! If not a diagonal element, but only a single excitation, then the
       ! corresponding Hamiltonian element may be non-zero. Calculate if the
       ! flipped spins are neighbours on the lattice, and if so, add the
       ! contribution from this site.
           hmatel = slater_condon1_hub_real(sys, excitation%from_orb(1), excitation%to_orb(1), excitation%perm)
           estimator_numerators(energy_index) = estimator_numerators(energy_index) + &
                                 (hmatel*walker_pop)
       end if

   end subroutine dmqmc_energy_hub_real

   subroutine dmqmc_correlation_function_heisenberg(sys, idet, excitation, walker_pop)

       ! For the Heisenberg model only.
       ! Add the contribution from the current density matrix element to the
       ! thermal spin correlation function estimator.

       ! In:
       !    sys: system being studied.
       !    idet: Current position in the main bitstring (density matrix) list.
       !    excitation: excit type variable which stores information on
       !        the excitation between the two bitstring ends, corresponding to
       !        the two labels for the density matrix element.
       !    walker_pop: number of particles on the current density matrix
       !        element.

       use bit_utils, only: count_set_bits
       use excitations, only: excit
       use fciqmc_data, only: walker_dets
       use fciqmc_data, only: walker_data, H00, correlation_mask
       use fciqmc_data, only: estimator_numerators, correlation_index
       use real_lattice, only: connected_orbs
       use system, only: sys_t

       type(sys_t), intent(in) :: sys
       integer, intent(in) :: idet
       type(excit), intent(in) :: excitation
       real(p), intent(in) :: walker_pop
       integer(i0) :: f(sys%basis%basis_length)
       integer :: bit_element1, bit_position1, bit_element2, bit_position2
       integer :: sign_factor

       ! If no excitation, we want the diagonal element of the correlation
       ! function operator.
       if (excitation%nexcit == 0) then
           ! If the two spins i and j are the same, the matrix element is +1/4.
           ! If they are different, the matrix element is -1/4. So we want
           ! sign_factor to be +1 in the former case, and -1 in the latter case.
           ! f as calculated below will have 0's at sites other than i and j,
           ! and the same values as walker_dets at i and j. Hence, if f has
           ! two 1's or no 1's, we want sign_factor = +1. Else if we have one 1,
           ! we want sign_factor = -1.
           f = iand(walker_dets(:sys%basis%basis_length,idet), correlation_mask)
           ! Count if we have zero, one or two 1's.
           sign_factor = sum(count_set_bits(f))
           ! The operation below will map 0 and 2 to +1, and will map 1 to -1,
           ! as is easily checked.
           sign_factor = (mod(sign_factor+1,2)*2)-1
           ! Hence sign_factor can be used to find the matrix element, as used
           ! below.
           estimator_numerators(correlation_index) = estimator_numerators(correlation_index) + &
                                    (sign_factor*(walker_pop/4))
       else if (excitation%nexcit == 1) then
           ! If not a diagonal element, but only a single excitation, then the
           ! corresponding matrix element will be 1/2 if and only if the two
           ! sites which are flipped are sites i and j, else it will be 0. We
           ! assume that excitations will only be set if i and j are opposite
           ! (else they could not be flipped, for ms=0).
           bit_position1 = sys%basis%bit_lookup(1,excitation%from_orb(1))
           bit_element1 = sys%basis%bit_lookup(2,excitation%from_orb(1))
           bit_position2 = sys%basis%bit_lookup(1,excitation%to_orb(1))
           bit_element2 = sys%basis%bit_lookup(2,excitation%to_orb(1))
           if (btest(correlation_mask(bit_element1), bit_position1) .and. btest(correlation_mask(bit_element2), bit_position2)) &
                 estimator_numerators(correlation_index) = estimator_numerators(correlation_index) + (walker_pop/2)
       end if

   end subroutine dmqmc_correlation_function_heisenberg

   subroutine dmqmc_stag_mag_heisenberg(sys, idet, excitation, walker_pop)

       ! For the Heisenberg model only.
       ! Add the contribution from the current density matrix element to the
       ! thermal staggered magnetisation estimate.

       ! In:
       !    sys: system being studied.
       !    idet: Current position in the main bitstring (density matrix) list.
       !    excitation: excit type variable which stores information on
       !        the excitation between the two bitstring ends, corresponding
       !        to the two labels for the density matrix element.
       !    walker_pop: number of particles on the current density matrix
       !        element.

       use bit_utils, only: count_set_bits
       use determinants, only: lattice_mask
       use excitations, only: excit
       use fciqmc_data, only: walker_dets
       use fciqmc_data, only: estimator_numerators, staggered_mag_index
       use system, only: sys_t

       type(sys_t), intent(in) :: sys
       integer, intent(in) :: idet
       type(excit), intent(in) :: excitation
       real(p), intent(in) :: walker_pop
       integer :: bit_element1, bit_position1, bit_element2, bit_position2
       integer(i0) :: f(sys%basis%basis_length)
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
           ! Note the number of down spins on a sublattice is easily obtained
           ! from N_u(+) since there are N/2 spins on each - and the number of
           ! spins up on a different sublattice is easily obtained since there
           ! are nel spins up in total. Hence the matrix element will be written
           ! only in terms of the number of up spins on sublattice 1, to save
           ! computation.
           f = iand(walker_dets(:sys%basis%basis_length,idet), lattice_mask)
           n_up_plus = sum(count_set_bits(f))
           ! Below, the term in brackets and middle term come from the z
           ! component (the z operator is diagonal) and one nsites/4 factor
           ! comes from the x operator, the other nsites/4 factor from the y
           ! operator.
           total_sum = (2*n_up_plus-sys%nel)**2 + (sys%lattice%nsites/2)
       else if (excitation%nexcit == 1) then
           ! Off-diagonal elements from the y and z operators. For the pair of
           ! spins that are flipped, if they are on the same sublattice, we get
           ! a factor of 1, or if on different sublattices, a factor of -1.
           bit_position1 = sys%basis%bit_lookup(1,excitation%from_orb(1))
           bit_element1 = sys%basis%bit_lookup(2,excitation%from_orb(1))
           bit_position2 = sys%basis%bit_lookup(1,excitation%to_orb(1))
           bit_element2 = sys%basis%bit_lookup(2,excitation%to_orb(1))
           if (btest(lattice_mask(bit_element1), bit_position1)) total_sum = total_sum+1
           if (btest(lattice_mask(bit_element2), bit_position2)) total_sum = total_sum+1
           ! The operation below will map 0 and 2 to +1, and will map 1 to -1,
           ! as is easily checked. We want this - if both or no spins on this
           ! sublattice, then both on same sublattice either way, so plus one.
           ! Else they are on different sublattices, so we want a factor of -1,
           ! as we get.
           total_sum = (mod(total_sum+1,2)*2)-1
       end if

       estimator_numerators(staggered_mag_index) = estimator_numerators(staggered_mag_index) + &
                                  (real(total_sum)/real(sys%lattice%nsites**2))*walker_pop

   end subroutine dmqmc_stag_mag_heisenberg

   subroutine update_full_renyi_2(walker_pop)

       ! Add the contribution from the current density matrix element to the
       ! Renyi entropy (S_2) of the full density matrix.

       ! In:
       !    walker_pop: number of particles on the current density matrix
       !        element, for both replicas.

       use fciqmc_data, only: estimator_numerators, full_r2_index

       real(p), intent(in) :: walker_pop(:)

       estimator_numerators(full_r2_index) = estimator_numerators(full_r2_index) + &
                                  walker_pop(1)*walker_pop(2)

   end subroutine update_full_renyi_2

   subroutine update_reduced_density_matrix_heisenberg(idet, excitation, walker_pop, iteration)

       ! Add the contribution from the current walker to the reduced density
       ! matrices being sampled. This is performed by 'tracing out' the
       ! subsystem B spins.

       ! Applicable only to the Heisenberg model.

       ! This procedure takes the two determinants bitstrings for the current
       ! psips and, if the two bitstrings of the B subsystem are identical,
       ! adds the walker population to the corresponding reduced density matrix
       ! element.

       ! In:
       !    idet: current position in the main bitstring (density matrix) list.
       !    excitation: excit type variable which stores information on
       !        the excitation between the two bitstring ends, corresponding to
       !        the two labels for the density matrix element.
       !    walker_pop: number of particles on the current density matrix
       !        element. Note that this walker population is still weighted
       !        by the importance sampling factors. These factors must be
       !        removed before any estimates can be calculated.
       !    iteration: interation number.  No accumulation of the RDM is
       !        performed if iteration <= start_averaging.

       use basis, only: basis_global
       use dmqmc_procedures, only: decode_dm_bitstring, rdms
       use excitations, only: excit
       use fciqmc_data, only: reduced_density_matrix, walker_dets, walker_population
       use fciqmc_data, only: sampling_size, calc_inst_rdm, calc_ground_rdm, nrdms
       use fciqmc_data, only: start_averaging, rdm_spawn, dmqmc_accumulated_probs
       use fciqmc_data, only: nsym_vec, real_factor
       use spawning, only: create_spawned_particle_rdm

       integer, intent(in) :: idet, iteration
       integer(int_p), intent(in) :: walker_pop(sampling_size)
       type(excit), intent(in) :: excitation
       real(p) :: unweighted_walker_pop(sampling_size)
       integer :: irdm, isym, ireplica
       integer(i0) :: f1(basis_global%basis_length), f2(basis_global%basis_length)

       if (.not. (iteration > start_averaging .or. calc_inst_rdm)) return

       ! Loop over all RDMs to be calculated.
       do irdm = 1, nrdms
           ! Loop over every symmetry-equivalent subsystem for this RDM.
           do isym = 1, nsym_vec

               ! Apply the mask for the B subsystem to set all sites in the A
               ! subsystem to 0.
               f1 = iand(rdms(irdm)%B_masks(:,isym),walker_dets(:basis_global%basis_length,idet))
               f2 = iand(rdms(irdm)%B_masks(:,isym),walker_dets(basis_global%basis_length+1:basis_global%total_basis_length,idet))

               ! Once this is done, check if the resulting bitstrings (which can
               ! only possibly have 1's in the B subsystem) are identical. If
               ! they are, then this psip contributes to the reduced density
               ! matrix for subsystem A. This is because we get the reduced
               ! density matrix for A by 'tracing out' over B, which in practice
               ! means only keeping matrix elements that are on the diagonal for
               ! subsystem B.
               if (sum(abs(f1-f2)) == 0_i0) then
                   ! Call a function which maps the subsystem A state to two RDM
                   ! bitstrings.
                   call decode_dm_bitstring(walker_dets(:,idet),irdm,isym)

                   if (calc_ground_rdm) then
                       ! The above routine actually maps to numbers between 0
                       ! and 2^rdms(1)%A_nsites-1, but the smallest and largest
                       ! reduced density matrix indices are one more than these,
                       ! so add one.
                       rdms(irdm)%end1 = rdms(irdm)%end1 + 1
                       rdms(irdm)%end2 = rdms(irdm)%end2 + 1
                       unweighted_walker_pop = real(walker_population(:,idet),p)*&
                                                dmqmc_accumulated_probs(excitation%nexcit)/real_factor
                       ! Note, when storing the entire RDM (as done here), the
                       ! maximum value of rdms(i)%rdm_basis_length is 1, so we
                       ! only consider this one element here.
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

        ! Wrapper for calling ground-state RDM procedures (*not*
        ! beta-dependent RDMs).

       ! In:
       !    beta_cycle: index of the beta loop being performed.

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
        ! If in parallel then merge the reduced density matrix onto one
        ! processor.
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
            ! Force the reduced density matrix to be symmetric by averaging the
            ! upper and lower triangles.
            do i = 1, ubound(reduced_density_matrix,1)
                do j = 1, i-1
                    reduced_density_matrix(i,j) = 0.5_p*(reduced_density_matrix(i,j) + reduced_density_matrix(j,i))
                    reduced_density_matrix(j,i) = reduced_density_matrix(i,j)
                end do
                ! Add current contribution to the trace.
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

        ! Need to paralellise for large subsystems and introduce test to check
        ! whether diagonalisation should be performed in serial or paralell.

        ! In:
        !    trace_rdm: The trace of the RDM being considered.

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
        vn_entropy = 0.0_p

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

        ! The matrix input into the following diagonalisation routines will have
        ! their upper half (including the diagonal) destroyed. We might want
        ! reduced_desntiy_matrix later, so use some temporary space:
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

        ! Calculate the concurrence of a qubit. For a reduced density matrix
        ! \rho, the concurrence,
        ! C =  max(0, \lamda_1 - \lambda_2 - \lambda_3 -\lambda_4)
        ! where \lambda_i are the eigenvalues of the matrix,
        ! R = \sqrt{\sqrt{\rho}\~{\rho}\sqrt{\rho}},
        ! and
        ! \~\rho = {\sigma_y \otimes \sigma_y} \rho^{\ast} {\sigma_y \otimes \sigma_y}. 
        ! \lambda_1 > ... > \lambda_4.

        ! This can be simplified to finding the square root of the eigenvalues
        ! \{\lambda_i\} of \rho\~{\rho} and in the case where \rho is a real,
        ! symmetric matrix then we can further simplify the problem to finding
        ! the eigenvalues of R = \rho \sigma_y \otimes \sigma_y.

        ! Below we have named {\sigma_y \otimes \sigma_y} flip_spin_matrix as in
        ! the literature.

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
        ! Calculate the concurrence. Take abs of eigenvalues so that this is
        ! equivalant to sqauring and then square-rooting.
        concurrence = 2.0_p*maxval(abs(reigv)) - sum(abs(reigv)) 
        concurrence = max(0.0_p, concurrence)
        write (6,'(1x,a28,1X,f22.12)') "# Unnormalised concurrence= ", concurrence

    end subroutine calculate_concurrence
    
    subroutine calculate_rdm_traces(rdm_data, rdm_lists, traces)

        ! In:
        !    rdm_data: Array of rdm derived types, holding information about
        !        the various subsystems for which RDMs are being estimated.
        !    rdm_lists: Array of rdm_spawn_t derived types, which hold all of
        !        the RDM psips which belong to this processor.
        ! Out:
        !    r2: The calculated RDM traces.

        use dmqmc_procedures, only: rdm
        use excitations, only: get_excitation_level
        use spawn_data, only: spawn_t

        type(rdm), intent(in) :: rdm_data(:)
        type(spawn_t), intent(in) :: rdm_lists(:)
        real(p), intent(out) :: traces(:,:)
        integer :: irdm, i, rdm_bl
        integer, parameter :: thread_id = 0

        traces = 0.0_p 

        ! Loop over all RDMs being calculated.
        do irdm = 1, size(rdm_data)
            rdm_bl = rdm_data(irdm)%rdm_basis_length
            ! Loop over the total population of RDM psips on this processor.
            do i = 1, rdm_lists(irdm)%head(thread_id,0)
                ! If on the diagonal of the RDM...
                if (all( rdm_lists(irdm)%sdata(1:rdm_bl,i) == rdm_lists(irdm)%sdata(rdm_bl+1:2*rdm_bl,i))) then
                    traces(:,irdm) = traces(:,irdm) + &
                        real(rdm_lists(irdm)%sdata(rdm_lists(irdm)%bit_str_len+1:rdm_lists(irdm)%element_len,i),p)
                end if
            end do
        end do

    end subroutine calculate_rdm_traces

    subroutine calculate_rdm_renyi_2(rdm_data, rdm_lists, r2)

        ! Calculate the Renyi entropy (S_2) for all instantaneous RDMs being
        ! calculated.

        ! In:
        !    rdm_data: Array of rdm derived types, holding information about the
        !        various subsystems for which RDMs are being estimated.
        !    rdm_lists: Array of rdm_spawn_t derived types, which hold all of
        !        the RDM psips which belong to this processor.
        ! Out:
        !    r2: The calculated Renyi entropies (S_2).

        use dmqmc_procedures, only: rdm
        use excitations, only: get_excitation_level
        use fciqmc_data, only: dmqmc_accumulated_probs_old
        use spawn_data, only: spawn_t

        type(rdm), intent(in) :: rdm_data(:)
        type(spawn_t), intent(in) :: rdm_lists(:)
        real(p), intent(out) :: r2(:)
        integer :: i, irdm, excit_level, rdm_bl
        real(p) :: unweighted_pop_1, unweighted_pop_2
        integer, parameter :: thread_id = 0

        r2 = 0.0_p 

        ! Loop over all RDMs being calculated.
        do irdm = 1, size(rdm_data)
            rdm_bl = rdm_data(irdm)%rdm_basis_length
            ! Loop over the total population of RDM psips on this processor.
            do i = 1, rdm_lists(irdm)%head(thread_id,0)

                excit_level = get_excitation_level(int(rdm_lists(irdm)%sdata(1:rdm_bl,i), i0), &
                    int(rdm_lists(irdm)%sdata(rdm_bl+1:2*rdm_bl,i), i0) )

                ! Renormalise the psip populations to correct for the importance
                ! sampling procedure used.
                unweighted_pop_1 = rdm_lists(irdm)%sdata(rdm_lists(irdm)%bit_str_len+1,i)*dmqmc_accumulated_probs_old(excit_level)
                unweighted_pop_2 = rdm_lists(irdm)%sdata(rdm_lists(irdm)%bit_str_len+2,i)*dmqmc_accumulated_probs_old(excit_level)

                ! As we only hold RDM elements above the diagonal, off-diagonal
                ! elements must be counted twice.
                if (excit_level == 0) then
                    r2(irdm) = r2(irdm) + unweighted_pop_1*unweighted_pop_2
                else
                    r2(irdm) = r2(irdm) + 2*unweighted_pop_1*unweighted_pop_2
                end if

            end do
        end do

    end subroutine calculate_rdm_renyi_2

end module dmqmc_estimators

module dmqmc_estimators

use const

implicit none

contains

   subroutine update_dmqmc_estimators(ntot_particles_old)

        ! Update the shift and average the shift and estimators

        ! Should be called every report loop in an DMQMC calculation.

        ! Inout:
        !    ntot_particles_old: total number (across all processors) of
        !        particles in the simulation at end of the previous report loop.
        !        Returns the current total number of particles for use in the
        !        next report loop.

        use calc, only: doing_dmqmc_calc, dmqmc_energy, dmqmc_staggered_magnetisation
        use calc, only: dmqmc_energy_squared
        use energy_evaluation, only: update_shift
        use fciqmc_data, only: nparticles, sampling_size, target_particles, rspawn
        use fciqmc_data, only: shift, vary_shift
        use fciqmc_data, only: estimator_numerators, number_dmqmc_estimators
        use fciqmc_data, only: nreport, ncycles, trace
        use parallel

        integer(lint), intent(inout) :: ntot_particles_old(sampling_size)

#ifdef PARALLEL
        real(dp) :: ir(sampling_size+1+((number_dmqmc_estimators+1)*ncycles))
        real(dp) :: ir_sum(sampling_size+1+((number_dmqmc_estimators+1)*ncycles))
        integer :: ierr
        integer(lint) :: ntot_particles(sampling_size)

        ! Need to sum the number of particles and other quantites over all processors.
        ir(1:sampling_size) = nparticles
        ir(sampling_size+1) = rspawn
        ir(sampling_size+2) = trace
        ir(sampling_size+3:sampling_size+2+number_dmqmc_estimators) = estimator_numerators
        ! Merge the lists from each processor together.
        call mpi_allreduce(ir, ir_sum, size(ir), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
        ntot_particles = nint(ir_sum(1:sampling_size),lint)
        rspawn = ir_sum(sampling_size+1)
        if (parent) then
            ! Let the parent processor store the merged data lists, so that these
            ! can be output afterwards by the parent.
            trace = nint(ir_sum(sampling_size+2))
            estimator_numerators = ir_sum(sampling_size+3:sampling_size+2+number_dmqmc_estimators)
        end if

        if (vary_shift) then
            call update_shift(ntot_particles_old(1), ntot_particles(1), ncycles)
        end if
        ntot_particles_old = ntot_particles
        if (ntot_particles(1) > target_particles .and. .not.vary_shift) then
            vary_shift = .true.
        end if

#else
        ! If only one a single core, don't need to merge the data
        ! from several processors. Can simply output everything as it is.

        if (vary_shift) call update_shift(ntot_particles_old(1), nparticles(1), ncycles)
        ntot_particles_old = nparticles
        if (nparticles(1) > target_particles .and. .not.vary_shift) then
            vary_shift = .true.
        end if
#endif

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
       !    idet: Current position in the main bitsrting list.

       use basis, only: basis_length, total_basis_length
       use calc, only: doing_dmqmc_calc, dmqmc_energy, dmqmc_staggered_magnetisation
       use calc, only: dmqmc_energy_squared, dmqmc_correlation
       use excitations, only: get_excitation, excit
       use fciqmc_data, only: walker_dets, walker_population, trace, doing_reduced_dm
       use fciqmc_data, only: dmqmc_accumulated_probs, reduced_dm_start_averaging
       use proc_pointers, only: update_dmqmc_energy_ptr, update_dmqmc_stag_mag_ptr
       use proc_pointers, only: update_dmqmc_energy_squared_ptr, update_dmqmc_correlation_ptr

       integer, intent(in) :: idet, iteration
       type(excit) :: excitation
       real(p) :: unweighted_walker_pop

       ! Get excitation.
       excitation = get_excitation(walker_dets(:basis_length,idet), &
                        walker_dets((1+basis_length):total_basis_length,idet))

       ! When performing importance sampling the result is that certain excitation
       ! levels have smaller psips populations than the true density matrix by some
       ! factor. In these cases, we want to multiply the psip population by this factor
       ! to calculate the contribution from these excitation levels correctly.
       
       ! In the case of no importance sampling, unweighted_walker_pop = walker_population(1,idet)
       ! and so this change can be ignored.
       unweighted_walker_pop = walker_population(1,idet)*dmqmc_accumulated_probs(excitation%nexcit)

       ! If diagonal element, add to the trace.
       if (excitation%nexcit == 0) trace = trace + walker_population(1,idet)
       ! See which estimators are to be calculated, and call the corresponding procedures.
       ! Energy
       if (doing_dmqmc_calc(dmqmc_energy)) call update_dmqmc_energy_ptr&
               &(idet, excitation, unweighted_walker_pop)
       ! Energy squared
       if (doing_dmqmc_calc(dmqmc_energy_squared)) call update_dmqmc_energy_squared_ptr&
               &(idet, excitation, unweighted_walker_pop)
       ! Spin-spin correlation function
       if (doing_dmqmc_calc(dmqmc_correlation)) call update_dmqmc_correlation_ptr&
               &(idet, excitation, unweighted_walker_pop)
       ! Staggered magnetisation
       if (doing_dmqmc_calc(dmqmc_staggered_magnetisation)) call update_dmqmc_stag_mag_ptr&
               &(idet, excitation, unweighted_walker_pop)
       ! Reduced density matrix
       if (doing_reduced_dm .and. iteration > reduced_dm_start_averaging) &
               call update_reduced_density_matrix(idet, unweighted_walker_pop)

   end subroutine call_dmqmc_estimators

   subroutine dmqmc_energy_heisenberg(idet, excitation, walker_pop)

       ! For the Heisenberg model only.
       ! Add the contribution from the current density matrix element
       ! to the thermal energy estimate.

       ! In:
       !    idet: Current position in the main bitsrting list.
       !    excitation: excit type variable which stores information on
       !    the excitation between the two bitstring ends, corresponding
       !    to the two labels for the density matrix element.

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
       !    idet: Current position in the main bitsrting list.
       !    excitation: excit type variable which stores information on
       !    the excitation between the two bitstring ends, corresponding
       !    to the two labels for the density matrix element.

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
           ! which way up the intermeddiate spin is - there will always be one order which contributes.
           ! If there are two such paths, then this could happen by either paths, but again, the two
           ! intermeddiate spins will only allow one order of spin flipping for each path, no
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
       !    idet: Current position in the main bitsrting list.
       !    excitation: excit type variable which stores information on
       !    the excitation between the two bitstring ends, corresponding
       !    to the two labels for the density matrix element.

       use basis, only: basis_length, total_basis_length
       use excitations, only: excit
       use fciqmc_data, only: walker_dets
       use fciqmc_data, only: walker_data, H00
       use fciqmc_data, only: estimator_numerators, energy_index
       use hamiltonian, only: slater_condon1_hub_real

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
       !    idet: Current position in the main bitsrting list.
       !    excitation: excit type variable which stores information on
       !    the excitation between the two bitstring ends, corresponding
       !    to the two labels for the density matrix element.

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
           if (btest(correlation_mask(bit_element1), bit_position1).and.btest(correlation_mask(bit_element2), bit_position2)) &
                 estimator_numerators(correlation_index) = estimator_numerators(correlation_index) + &
                                   (walker_pop/2)
       end if

   end subroutine dmqmc_correlation_function_heisenberg

   subroutine dmqmc_stag_mag_heisenberg(idet, excitation, walker_pop)

       ! For the Heisenberg model only.
       ! Add the contribution from the current density matrix element
       ! to the thermal staggered magnetisation estimate.

       ! In:
       !    idet: Current position in the main bitsrting list.
       !    excitation: excit type variable which stores information on
       !    the excitation between the two bitstring ends, corresponding
       !    to the two labels for the density matrix element.

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

   subroutine update_reduced_density_matrix(idet, walker_pop)

       ! Add a contirbution from the current walker to the reduced density
       ! matrix estimator. This function takes the two bitstrings for the current
       ! walker and, if the two bitstrings of the B subsystem are identical,
       ! adds the walker population to the corresponding reduced density matrix element.

       use basis, only: basis_length, total_basis_length
       use dmqmc_procedures, only: decode_dm_bitstring
       use fciqmc_data, only: subsystem_B_mask, reduced_density_matrix
       use fciqmc_data, only: walker_dets, walker_population

       integer, intent(in) :: idet
       real(p), intent(in) :: walker_pop
       integer(i0) :: f1(basis_length), f2(basis_length)
       integer(i0) :: end1, end2

       ! Apply the mask for the B subsystem to set all sites in the A subsystem to 0.
       f1 = iand(subsystem_B_mask,walker_dets(:basis_length,idet))
       f2 = iand(subsystem_B_mask,walker_dets(basis_length+1:total_basis_length,idet))

       ! Once this is done, check if the resulting bitstring (which can only possibly
       ! have 1's in the B subsystem) are identical. If they are, then this psip gives
       ! a contibution to the reduced density matrix for subsystem A. This is because we
       ! get the reduced density matrix for A by 'tracing out' over B, which in practice
       ! means only keeping matrix elements that are on the diagonal for subsystem B.
       if (sum(abs(f1-f2)) == 0) then
           ! We need to assign positions in the reduced density matrix for this walker,
           ! so call a function which maps the possible subsystem A spin configurations
           ! to indices.
           call decode_dm_bitstring(walker_dets(:,idet),end1,end2)
           reduced_density_matrix(end1,end2) = reduced_density_matrix(end1,end2) + walker_pop
       end if

    end subroutine update_reduced_density_matrix

    subroutine output_reduced_density_matrix()

        ! Normalise and ouput the reduced density matrix to the screen. This is
        ! called at the end of each beta loop in DMQMC calculations, when the
        ! reduced density matrix is requested.

        use fciqmc_data, only: reduced_density_matrix, subsystem_A_size
        use parallel
        use utils, only: print_matrix

        integer :: i, rdm_size
        real(p) :: trace_rdm
   
#ifdef PARALLEL
        real(dp) :: dm(2**subsystem_A_size,2**subsystem_A_size)
        real(dp) :: dm_sum(2**subsystem_A_size,2**subsystem_A_size)
        integer :: ierr

        dm = reduced_density_matrix

        call mpi_allreduce(dm, dm_sum, size(dm), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

        reduced_density_matrix = dm_sum

#endif

        trace_rdm = 0.0_p
        rdm_size = ubound(reduced_density_matrix,1)

        if (parent) then
            do i = 1, rdm_size
                trace_rdm = trace_rdm + reduced_density_matrix(i,i)
            end do
            reduced_density_matrix = reduced_density_matrix/trace_rdm
            write (6, '(a23)') 'Reduced density matrix:'
            call print_matrix(reduced_density_matrix)
        end if

    end subroutine output_reduced_density_matrix

    subroutine calculate_vn_entropy()

        ! Calculate the Von Neumann Entropy. Use lapack to calculate the 
        ! eigenvalues {\lambda_j} of the reduced density matrix. 
        ! Then VN Entropy S = -\sum_j\lambda_j\log_2{\lambda_j}

        ! Need to paralellise for large subsystems and introduce test to
        ! check whether diagonalisation should be performed in serial or
        ! paralell.
       
        use checking, only: check_allocate, check_deallocate
        use fciqmc_data, only: reduced_density_matrix, subsystem_A_size
        use parallel

        integer :: i, rdm_size
        integer :: info, ierr, lwork
        real(p), allocatable :: work(:)
        real(p) :: eigv(2**subsystem_A_size)
        real(p) :: trace_rdm, vn_entropy
#ifdef PARALLEL
        real(dp) :: dm(2**subsystem_A_size,2**subsystem_A_size)
        real(dp) :: dm_sum(2**subsystem_A_size,2**subsystem_A_size)

        vn_entropy = 0._p
        dm = reduced_density_matrix

        call mpi_allreduce(dm, dm_sum, size(dm), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

        reduced_density_matrix = dm_sum

#endif

        trace_rdm = 0.0_p
        rdm_size = ubound(reduced_density_matrix,1)

        if (parent) then
            do i = 1, rdm_size
                trace_rdm = trace_rdm + reduced_density_matrix(i,i)
            end do
            reduced_density_matrix = reduced_density_matrix/trace_rdm
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

#ifdef SINGLE_PRECISION
            call ssyev('N', 'U', rdm_size, reduced_density_matrix, rdm_size, eigv, work, lwork, info)
#else
            call dsyev('N', 'U', rdm_size, reduced_density_matrix, rdm_size, eigv, work, lwork, info)
#endif
            do i = 1, ubound(eigv,1)
                vn_entropy = vn_entropy - eigv(i)*(log(eigv(i))/log(2.))             
            end do
            write (6,'(1x,a23,1X,f22.12)') "# Von-Neumann Entropy= ", vn_entropy
        end if
        
    end subroutine calculate_vn_entropy

end module dmqmc_estimators

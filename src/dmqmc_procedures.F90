module dmqmc_procedures

use const
implicit none

! This type contains information for the RDM corresponding to a
! given subsystem. It takes translational symmetry into account by
! storing information for all subsystems which are equivalent by
! translational symmetry.
type rdm
    ! The total number of sites in subsystem A.
    integer :: A_nsites
    ! Similar to basis_length, rdm_basis_length is the length of the
    ! byte array necessary to contain a bit for each subsystem-A basis
    ! function. An array of twice this length is stored to hold both
    ! RDM indices.
    integer :: rdm_basis_length
    ! The sites in subsystem A, as entered by the user.
    integer, allocatable :: subsystem_A(:)
    ! B_masks(:,i) has bits set at all bit positions corresponding to
    ! sites in version i of subsystem B, where the different 'versions'
    ! correspond to subsystems which are equivalent by symmetry.
    integer, allocatable :: B_masks(:,:)
    ! bit_pos(i,j,1) contains the position of the bit corresponding to
    ! site i in 'version' j of subsystem A.
    ! bit_pos(i,j,2) contains the element of the bit corresponding to
    ! site i in 'version' j of subsystem A.
    ! Note that site i in a given version is the site that corresponds to
    ! site i in all other versions of subsystem A (and so bit_pos(i,:,1)
    ! and bit_pos(i,:,2) will not be sorted). This is very important
    ! so that equivalent psips will contribute to the same RDM element.
    integer, allocatable :: bit_pos(:,:,:)
    ! Two bitstrings of length rdm_basis_length. To be used as temporary
    ! bitstrings to prevent having to regularly allocate different
    ! length bitstrings for different RDMs.
    integer(i0), allocatable :: end1(:), end2(:)
end type rdm

! This stores all the information for the various RDMs that the user asks
! to be calculated. Each element of this array corresponds to one of these RDMs.
type(rdm), allocatable :: rdms(:)

contains

    subroutine init_dmqmc(sys)

         ! In:
         !    sys: system being studied.

         use basis, only: basis_length, total_basis_length, bit_lookup, basis_lookup
         use calc, only: doing_dmqmc_calc, dmqmc_calc_type, dmqmc_energy, dmqmc_energy_squared
         use calc, only: dmqmc_staggered_magnetisation, dmqmc_correlation, dmqmc_full_r2
         use checking, only: check_allocate
         use fciqmc_data, only: trace, energy_index, energy_squared_index, correlation_index
         use fciqmc_data, only: staggered_mag_index, estimator_numerators, doing_reduced_dm
         use fciqmc_data, only: dmqmc_factor, number_dmqmc_estimators, ncycles, tau, dmqmc_weighted_sampling
         use fciqmc_data, only: correlation_mask, correlation_sites, half_density_matrix
         use fciqmc_data, only: dmqmc_sampling_probs, dmqmc_accumulated_probs, flip_spin_matrix
         use fciqmc_data, only: doing_concurrence, calculate_excit_distribution, excit_distribution
         use fciqmc_data, only: nreport, average_shift_until, shift_profile, dmqmc_vary_weights
         use fciqmc_data, only: finish_varying_weights, weight_altering_factors, dmqmc_find_weights
         use fciqmc_data, only: sampling_size, rdm_traces, nrdms, dmqmc_accumulated_probs_old
         use fciqmc_data, only: full_r2_index
         use system, only: sys_t

        type(sys_t), intent(in) :: sys

         integer :: ierr, i, bit_position, bit_element

         number_dmqmc_estimators = 0

         allocate(trace(sampling_size), stat=ierr)
         call check_allocate('trace',sampling_size,ierr)
         trace = 0.0_p

         allocate(rdm_traces(sampling_size,nrdms), stat=ierr)
         call check_allocate('rdm_traces',sampling_size*nrdms,ierr)
         rdm_traces = 0.0_p

         if (doing_dmqmc_calc(dmqmc_full_r2)) then
             number_dmqmc_estimators = number_dmqmc_estimators + 1
             full_r2_index = number_dmqmc_estimators
         end if
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
             allocate(excit_distribution(0:sys%max_number_excitations), stat=ierr)
             call check_allocate('excit_distribution',sys%max_number_excitations+1,ierr)
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
                 allocate(dmqmc_sampling_probs(1:sys%max_number_excitations), stat=ierr)
                 call check_allocate('dmqmc_sampling_probs',sys%max_number_excitations,ierr)
                 dmqmc_sampling_probs = 1.0_p
             end if
             if (half_density_matrix) dmqmc_sampling_probs(1) = dmqmc_sampling_probs(1)*2.0_p
             allocate(dmqmc_accumulated_probs(0:sys%max_number_excitations), stat=ierr)
             call check_allocate('dmqmc_accumulated_probs',sys%max_number_excitations+1,ierr)
             allocate(dmqmc_accumulated_probs_old(0:sys%max_number_excitations), stat=ierr)
             call check_allocate('dmqmc_accumulated_probs_old',sys%max_number_excitations+1,ierr)
             dmqmc_accumulated_probs(0) = 1.0_p
             dmqmc_accumulated_probs_old = 1.0_p
             do i = 1, size(dmqmc_sampling_probs)
                 dmqmc_accumulated_probs(i) = dmqmc_accumulated_probs(i-1)*dmqmc_sampling_probs(i)
             end do
             dmqmc_accumulated_probs(size(dmqmc_sampling_probs)+1:sys%max_number_excitations) = &
                                    dmqmc_accumulated_probs(size(dmqmc_sampling_probs))
             if (dmqmc_vary_weights) then
                 ! Allocate an array to store the factors by which the weights will change each
                 ! iteration.
                 allocate(weight_altering_factors(0:sys%max_number_excitations), stat=ierr)
                 call check_allocate('weight_altering_factors',sys%max_number_excitations+1,ierr) 
                 weight_altering_factors = dble(dmqmc_accumulated_probs)**(1/dble(finish_varying_weights))
                 ! If varying the weights, start the accumulated probabilties as all 1.0
                 ! initially, and then alter them gradually later.
                 dmqmc_accumulated_probs = 1.0_p
             end if
         else
             allocate(dmqmc_accumulated_probs(0:sys%max_number_excitations), stat=ierr)
             call check_allocate('dmqmc_accumulated_probs',sys%max_number_excitations+1,ierr)
             allocate(dmqmc_accumulated_probs_old(0:sys%max_number_excitations), stat=ierr)
             call check_allocate('dmqmc_accumulated_probs_old',sys%max_number_excitations+1,ierr)
             dmqmc_accumulated_probs = 1.0_p
             dmqmc_accumulated_probs_old = 1.0_p
             if (half_density_matrix) dmqmc_accumulated_probs(1:sys%max_number_excitations) = 2.0_p
         end if

         ! If doing a reduced density matrix calculation, allocate and define the bit masks that
         ! have 1's at the positions referring to either subsystems A or B.
         if (doing_reduced_dm) call setup_rdm_arrays(sys)

         ! If doing concurrence calculation then construct and store the 4x4 flip spin matrix i.e.
         ! \sigma_y \otimes \sigma_y
         if (doing_concurrence) then
             allocate(flip_spin_matrix(4,4), stat=ierr)
             call check_allocate('flip_spin_matrix',16,ierr)
             flip_spin_matrix = 0.0_p
             flip_spin_matrix(1,4) = -1.0_p
             flip_spin_matrix(4,1) = -1.0_p
             flip_spin_matrix(3,2) = 1.0_p
             flip_spin_matrix(2,3) = 1.0_p    
         end if

    end subroutine init_dmqmc

    subroutine setup_rdm_arrays(sys)

        ! Setup the bit masks needed for RDM calculations. These are masks for the bits referring
        ! to either subsystem A or B. Also calculate the positions and elements of the sites
        ! in subsyetsm A, and finally allocate the RDM itself.

        ! In:
        !    sys: system being studied.

        use calc, only: ms_in, doing_dmqmc_calc, dmqmc_rdm_r2
        use checking, only: check_allocate
        use errors
        use fciqmc_data, only: reduced_density_matrix, nrdms, calc_ground_rdm, calc_inst_rdm
        use fciqmc_data, only: replica_tricks, renyi_2, sampling_size, bit_shift, spawn_cutoff
        use fciqmc_data, only: spawned_rdm_length, rdm_spawn
        use hash_table, only: alloc_hash_table
        use parallel, only: parent
        use spawn_data, only: alloc_spawn_t
        use system, only: sys_t, heisenberg
        use utils, only: int_fmt

        type(sys_t), intent(in) :: sys

        integer :: i, ierr, ipos, basis_find, size_spawned_rdm, total_size_spawned_rdm
        integer :: bit_position, bit_element, nbytes_int

        ! For the Heisenberg model only currently.
        if (sys%system==heisenberg) then
            call find_rdm_masks(sys)
        else
            call stop_all("setup_rdm_arrays","The use of RDMs is currently only implemented for the &
                           &Heisenberg model.")
        end if

        total_size_spawned_rdm = 0
        nbytes_int = bit_size(i)/8

        if (calc_inst_rdm) then
            allocate(rdm_spawn(nrdms), stat=ierr)
            call check_allocate('rdm_spawn', nrdms, ierr)
        end if
        if (doing_dmqmc_calc(dmqmc_rdm_r2)) then
            allocate(renyi_2(nrdms), stat=ierr)
            call check_allocate('renyi_2', nrdms, ierr)
            renyi_2 = 0.0_p
        end if

        do i = 1, nrdms
            rdms(i)%rdm_basis_length = ceiling(real(rdms(i)%A_nsites)/i0_length)

            allocate(rdms(i)%end1(rdms(i)%rdm_basis_length), stat=ierr)
            call check_allocate('rdms(i)%end1', rdms(i)%rdm_basis_length, ierr)
            allocate(rdms(i)%end2(rdms(i)%rdm_basis_length), stat=ierr)
            call check_allocate('rdms(i)%end2', rdms(i)%rdm_basis_length, ierr)
            rdms(i)%end1 = 0_i0
            rdms(i)%end2 = 0_i0

            ! With the calc_ground_rdm option, the entire RDM is allocated. If the following condition
            ! is met then the number of rows is greater than the maximum integer accessible. This
            ! would clearly be too large, so abort in this case.
            if (calc_ground_rdm .and. rdms(i)%rdm_basis_length > 1) call stop_all("setup_rdm_arrays",&
                "A requested RDM is too large for all indices to be addressed by a single integer.")

            if (calc_inst_rdm) then
                size_spawned_rdm = (rdms(i)%rdm_basis_length*2+sampling_size)*i0_length/8
                total_size_spawned_rdm = total_size_spawned_rdm + size_spawned_rdm
                if (spawned_rdm_length < 0) then
                    ! Given in MB.  Convert.
                    ! Note that the factor of 2 is because two spawning arrays are stored, and
                    ! 21*nbytes_int is added because there are 21 integers in the hash table for each
                    ! spawned rdm slot.
                    spawned_rdm_length = int((-real(spawned_rdm_length,p)*10**6)/&
                                          (2*size_spawned_rdm + 21*nbytes_int))
                end if

                ! Note the initiator approximation is not implemented for density matrix calculations.
                call alloc_spawn_t(rdms(i)%rdm_basis_length*2, sampling_size, .false., spawned_rdm_length, &
                                     spawn_cutoff, bit_shift, 27, rdm_spawn(i)%spawn)
                ! Hard code hash table collision limit for now.  This should
                ! give an ok performance...
                ! We will only use the first 2*rdm_basis_length elements for the
                ! hash, even though rdm_spawn%spawn%sdata is larger than that in
                ! the first dimension...
                call alloc_hash_table(3*spawned_rdm_length, 7, rdms(i)%rdm_basis_length*2, &
                                     0, 0, 17, rdm_spawn(i)%ht, rdm_spawn(i)%spawn%sdata)
            end if
        end do

        if (parent) then
            write (6,'(1X,a58,f7.2)') 'Memory allocated per core for the spawned RDM lists (MB): ', &
                total_size_spawned_rdm*real(2*spawned_rdm_length,p)/10**6
            write (6,'(1X,a49,'//int_fmt(spawned_rdm_length,1)//',/)') &
                'Number of elements per core in spawned RDM lists:', spawned_rdm_length
        end if

        ! For an ms = 0 subspace, assuming less than or exactly half the spins in the subsystem are in
        ! the subsystem, then any combination of spins can occur in the subsystem, from all spins down
        ! to all spins up. Hence the total size of the reduced density matrix will be 2**(number of spins
        ! in subsystem A).
        if (calc_ground_rdm) then
            if (ms_in == 0 .and. rdms(1)%A_nsites <= floor(real(sys%lattice%nsites,p)/2.0_p)) then
                allocate(reduced_density_matrix(2**rdms(1)%A_nsites,2**rdms(1)%A_nsites), stat=ierr)
                call check_allocate('reduced_density_matrix', 2**(2*rdms(1)%A_nsites),ierr)
                reduced_density_matrix = 0.0_p
            else
                if (ms_in /= 0) then
                    call stop_all("setup_rdm_arrays","Reduced density matrices can only be used for Ms=0 &
                                   &calculations.")
                else if (rdms(1)%A_nsites > floor(real(sys%lattice%nsites,p)/2.0_p)) then
                    call stop_all("setup_rdm_arrays","Reduced density matrices can only be used for subsystems &
                                  &whose size is less than half the total system size.")
                end if
            end if
        end if
        
    end subroutine setup_rdm_arrays

    subroutine find_rdm_masks(sys)

        ! Initialise bit masks for converting a density matrix basis function
        ! into its corresponding reduced density matrix basis function. Bit
        ! masks will be calculated for the subsystems requested by the user, and
        ! also for all subsystems which are equivalent by translational symmetry.

        ! In:
        !    sys: system being studied.

        use basis, only: basis_length, bit_lookup, basis_lookup, basis_fns, nbasis
        use checking, only: check_allocate, check_deallocate
        use errors
        use fciqmc_data, only: nrdms, nsym_vec
        use hubbard_real, only: find_translational_symmetry_vecs, map_vec_to_cell
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        integer :: i, j, k, l, ipos, ierr
        integer :: basis_find, bit_position, bit_element
        integer :: A_mask(basis_length)
        real(p), allocatable :: sym_vecs(:,:)
        integer :: r(sys%lattice%ndim)

        ! Return all translational symmetry vectors in sym_vecs.
        call find_translational_symmetry_vecs(sys, sym_vecs, nsym_vec)

        ! Allocate the RDM arrays.
        do i = 1, nrdms
            if (rdms(i)%A_nsites == sys%lattice%nsites) then
                call stop_all('find_rdm_masks','You are attempting to use the full density matrix &
                              &as an RDM. This is not supported. You should use the &
                              &dmqmc_full_renyi_2 option to calculate the Renyi entropy of the &
                              &whole lattice.')
            else
                allocate(rdms(i)%B_masks(basis_length,nsym_vec), stat=ierr)
                call check_allocate('rdms(i)%B_masks', nsym_vec*basis_length,ierr)
                allocate(rdms(i)%bit_pos(rdms(i)%A_nsites,nsym_vec,2), stat=ierr)
                call check_allocate('rdms(i)%bit_pos', nsym_vec*rdms(i)%A_nsites*2,ierr)
            end if
            rdms(i)%B_masks = 0
            rdms(i)%bit_pos = 0
        end do

        ! Run through every site on every subsystem and add every translational symmetry vector.
        do i = 1, nrdms ! Over every subsystem.
            do j = 1, nsym_vec ! Over every symmetry vector.
                A_mask = 0
                do k = 1, rdms(i)%A_nsites ! Over every site in the subsystem.
                    r = basis_fns(rdms(i)%subsystem_A(k))%l
                    r = r + nint(sym_vecs(:,j))
                    ! If r is outside the cell considered in this simulation, shift it by the
                    ! appropriate lattice vector so that it is in this cell.
                    call map_vec_to_cell(sys%lattice%ndim, sys%lattice%lvecs, r)
                    ! Now need to find which basis function this site corresponds to. Simply loop
                    ! over all basis functions and check...
                    do l = 1, nbasis
                        if (all(basis_fns(l)%l == r)) then
                            bit_position = bit_lookup(1,l)
                            bit_element = bit_lookup(2,l)

                            A_mask(bit_element) = ibset(A_mask(bit_element), bit_position)

                            rdms(i)%bit_pos(k,j,1) = bit_position
                            rdms(i)%bit_pos(k,j,2) = bit_element
                        end if
                    end do
                end do

                rdms(i)%B_masks(:,j) = A_mask
                ! We cannot just flip the mask for system A to get that for system B,
                ! because the trailing bits on the end don't refer to anything and should
                ! be set to 0. So, first set these to 1 and then flip all the bits.
                do ipos = 0, i0_end
                    basis_find = basis_lookup(ipos, basis_length)
                    if (basis_find == 0) then
                        rdms(i)%B_masks(basis_length,j) = ibset(rdms(i)%B_masks(basis_length,j),ipos)
                    end if
                end do
                rdms(i)%B_masks(:,j) = not(rdms(i)%B_masks(:,j))

            end do
        end do

        deallocate(sym_vecs,stat=ierr)
        call check_deallocate('sym_vecs',ierr)

    end subroutine find_rdm_masks

    subroutine create_initial_density_matrix(rng, sys, nparticles_tot)

        ! Create a starting density matrix by sampling the elements of the (unnormalised)
        ! identity matrix. This is a sampling of the (unnormalised) infinite-temperature
        ! density matrix. This is done by picking determinants/spin configurations from
        ! with uniform probabilities in the space being considered.

        ! In/Out:
        !    rng: random number generator.
        ! In:
        !    sys: system being studied.
        ! Out:
        !    nparticles_tot: The total number of psips in the generated density matrix.

        use annihilation, only: direct_annihilation
        use calc, only: initiator_approximation
        use dSFMT_interface, only:  dSFMT_t, get_rand_close_open
        use errors
        use fciqmc_data, only: sampling_size, D0_population, all_sym_sectors
        use parallel
        use system, only: sys_t, heisenberg
        use utils, only: binom_r

        type(dSFMT_t), intent(inout) :: rng
        type(sys_t), intent(in) :: sys
        real(dp), intent(out) :: nparticles_tot(sampling_size)
        real(dp) :: nparticles_temp(sampling_size)
        integer :: nel, ireplica, ierr
        integer :: npsips, npsips_this_proc
        real(dp) :: total_size, sector_size
        real(dp) :: r, prob

        npsips_this_proc = nint(D0_population)/nprocs
        ! If the initial number of psips does not split evenly between all processors,
        ! add the leftover psips to the first processors in order.
        if (nint(D0_population)-(nprocs*int(D0_population/nprocs)) > iproc) &
              npsips_this_proc = npsips_this_proc + 1

        nparticles_temp = 0.0_dp

        do ireplica = 1, sampling_size
            select case(sys%system)
            case(heisenberg)
                if (all_sym_sectors) then
                    ! The size (number of configurations) of all symmetry sectors combined.
                    total_size = 2.0_dp**(real(sys%lattice%nsites,dp))

                    do nel = 0, sys%lattice%nsites
                        ! The size of this symmetry sector alone.
                        sector_size = binom_r(sys%lattice%nsites, nel)
                        prob = real(npsips_this_proc,dp)*sector_size/total_size
                        npsips = floor(prob)
                        ! If there are a non-integer number of psips to be spawned in this sector
                        ! then add an extra psip with the required probability.
                        prob = prob - npsips
                        r = get_rand_close_open(rng)
                        if (r < prob) npsips = npsips + 1

                        nparticles_temp(ireplica) = nparticles_temp(ireplica) + real(npsips, dp)
                        call random_distribution_heisenberg(rng, nel, npsips, ireplica)
                    end do
                else
                    ! This process will always create excatly D0_population psips.
                    nparticles_tot = nint(real(D0_population, dp))

                    call random_distribution_heisenberg(rng, sys%nel, npsips_this_proc, ireplica)
                end if
            case default
                call stop_all('create_initial_density_matrix','DMQMC not implemented for this system.')
            end select
        end do

#ifdef PARALLEL
        if (all_sym_sectors) then
            ! Finally, count the total number of particles across all processes.
#ifdef PARALLEL
            call mpi_allreduce(nparticles_temp, nparticles_tot, sampling_size, MPI_REAL8, MPI_SUM, &
                                MPI_COMM_WORLD, ierr)
#else
            nparticles_tot = nparticles_temp
#endif
        end if
#else
        nparticles_tot = nparticles_temp
#endif

        call direct_annihilation(sys, rng, initiator_approximation)

    end subroutine create_initial_density_matrix

    subroutine random_distribution_heisenberg(rng, spins_up, npsips, ireplica)

        ! For the Heisenberg model only. Distribute the initial number of psips
        ! along the main diagonal. Each diagonal element should be chosen
        ! with the same probability.

        ! Start from state with all spins down, then choose the above number of
        ! spins to flip up with equal probability.

        ! In/Out:
        !    rng: random number generator.
        ! In:
        !    spins_up: for the spin configurations generated, this number
        !       specifies how many of the spins shall be up.
        !    npsips: The total number of psips to be created.
        !    ireplica: index of replica (ie which of the possible concurrent
        !       DMQMC populations are we initialising)

        use basis, only: nbasis, basis_length, bit_lookup
        use calc, only: ms_in
        use dSFMT_interface, only: dSFMT_t, get_rand_close_open
        use fciqmc_data, only: encoding_factor
        use parallel
        use system

        type(dSFMT_t), intent(inout) :: rng
        integer, intent(in) :: spins_up, npsips, ireplica
        integer :: i, rand_basis, bits_set
        integer :: bit_element, bit_position
        integer(i0) :: f(basis_length)
        real(dp) :: rand_num

        do i = 1, npsips

            ! Start with all spins down.
            f = 0
            bits_set = 0

            do
                ! If half the spins are now flipped up, we have our basis
                ! function fully created, so exit the loop.
                if (bits_set==spins_up) exit
                ! Choose a random spin to flip.
                rand_num = get_rand_close_open(rng)
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
            call create_diagonal_density_matrix_particle(f,encoding_factor,ireplica)

        end do

    end subroutine random_distribution_heisenberg

    subroutine create_diagonal_density_matrix_particle(f, nspawn, particle_type)

        ! Create a psip on a diagonal element of the density matrix by adding
        ! it to the spawned walkers list. This list can then be sorted correctly
        ! by the direct_annihilation routine.

        ! In:
        !    f: bitstring representation of index of the diagonal element upon
        !        which a new psip shall be placed.
        !    nspawn: the number of particles to be added to this diagonal element.
        !    particle_type: the label of the replica to which this particle is
        !        to sample.

        use hashing
        use basis, only: basis_length, total_basis_length
        use fciqmc_data, only: qmc_spawn
        use parallel
        use errors, only: stop_all

        integer(i0), intent(in) :: f(basis_length)
        integer(int_p), intent(in) :: nspawn
        integer ::particle_type
        integer(i0) :: f_new(total_basis_length)
#ifndef PARALLEL
        integer, parameter :: iproc_spawn = 0
#else
        integer :: iproc_spawn
#endif

        ! Create the bitstring of the determinant.
        f_new = 0
        f_new(:basis_length) = f
        f_new((basis_length+1):(total_basis_length)) = f

#ifdef PARALLEL
        ! Need to determine which processor the spawned psip should be sent to.
        iproc_spawn = modulo(murmurhash_bit_string(f_new, &
                                total_basis_length, qmc_spawn%hash_seed), nprocs)
#endif

        ! Move to the next position in the spawning array.
        qmc_spawn%head(0,iproc_spawn) = qmc_spawn%head(0,iproc_spawn) + 1

        ! qmc_spawn%head_start(0,1) holds the number of slots in the spawning array per processor.
        if (qmc_spawn%head(0,iproc_spawn) - qmc_spawn%head_start(0,iproc_spawn) >= qmc_spawn%head_start(0,1)) &
            call stop_all('create_diagonal_density_matrix_particle', 'There is no space left in the spawning array.')

        ! Set info in spawning array.
        ! Zero it as not all fields are set.
        qmc_spawn%sdata(:,qmc_spawn%head(0,iproc_spawn)) = 0_int_s
        ! indices 1 to total_basis_length store the bitstring.
        qmc_spawn%sdata(:(total_basis_length),qmc_spawn%head(0,iproc_spawn)) = int(f_new, int_s)
        ! The final index stores the number of psips created.
        qmc_spawn%sdata((total_basis_length)+particle_type,qmc_spawn%head(0,iproc_spawn)) = int(nspawn, int_s)

    end subroutine create_diagonal_density_matrix_particle

    subroutine decode_dm_bitstring(f, irdm, isym)

        ! This function maps a full DMQMC bitstring to two bitstrings encoding the subsystem-A
        ! RDM bitstrings. These resulting bitstrings are stored in the end1 and end2 components
        ! of rdms(irdm).
        
        ! Crucially, the mapping is performed so that, if there are two subsystems which are
        ! equivalent by symmetry, then equivalent sites in those two subsystems will be mapped to
        ! the same RDM bitstrings. This is clearly a requirement to obtain a correct representation
        ! of the RDM when using translational symmetry.

        ! In:
        !    f: bitstring representation of the subsystem-A state.
        !    irdm: The label of the RDM being considered.
        !    isym: The label of the symmetry vector being considered.

        use basis, only: basis_length, bit_lookup

        integer(i0), intent(in) :: f(2*basis_length)
        integer, intent(in) :: irdm, isym
        integer :: i, bit_pos, bit_element

        ! Start from all bits down, so that we can flip bits up one by one.
        rdms(irdm)%end1 = 0
        rdms(irdm)%end2 = 0

        ! Loop over all the sites in the subsystem considered for the reduced density matrix.
        do i = 1, rdms(irdm)%A_nsites
            ! Find the final bit positions and elements.
            bit_pos = bit_lookup(1,i)
            bit_element = bit_lookup(2,i)

            ! If the spin is up, set the corresponding bit in the first bitstring.
            if (btest(f(rdms(irdm)%bit_pos(i,isym,2)),rdms(irdm)%bit_pos(i,isym,1))) &
                rdms(irdm)%end1(bit_element) = ibset(rdms(irdm)%end1(bit_element),bit_pos)
            ! Similarly for the second index, by looking at the second end of the bitstring.
            if (btest(f(rdms(irdm)%bit_pos(i,isym,2)+basis_length),rdms(irdm)%bit_pos(i,isym,1))) &
                rdms(irdm)%end2(bit_element) = ibset(rdms(irdm)%end2(bit_element),bit_pos)
        end do

    end subroutine decode_dm_bitstring
 
    subroutine update_sampling_weights(rng)
        
        ! This routine updates the values of the weights used in importance sampling. It also
        ! removes or adds psips from the various excitation levels accordingly.

        ! In/Out:
        !    rng: random number generator.

        use annihilation, only: remove_unoccupied_dets
        use basis, only: basis_length, total_basis_length
        use excitations, only: get_excitation_level
        use fciqmc_data, only: dmqmc_accumulated_probs, finish_varying_weights
        use fciqmc_data, only: weight_altering_factors, tot_walkers, walker_dets, walker_population
        use fciqmc_data, only: nparticles, sampling_size, encoding_factor
        use dSFMT_interface, only: dSFMT_t, get_rand_close_open

        type(dSFMT_t), intent(inout) :: rng
        integer :: idet, ireplica, excit_level, nspawn, sign_factor
        real(dp) :: new_population_target(sampling_size)
        integer(int_p) :: old_population(sampling_size), new_population(sampling_size)
        real(dp) :: r, pextra

        ! Alter weights for the next iteration.
        dmqmc_accumulated_probs = dble(dmqmc_accumulated_probs)*weight_altering_factors

        ! When the weights for an excitation level are increased by a factor, the number
        ! of psips on that level has to decrease by the same factor, else the wavefunction
        ! which the psips represent will not be the correct importance sampled wavefunction
        ! for the new weights. The code below loops over every psip and destroys (or creates)
        ! it with the appropriate probability.
        do idet = 1, tot_walkers

            excit_level = get_excitation_level(walker_dets(1:basis_length,idet), &
                    walker_dets(basis_length+1:total_basis_length,idet))

            old_population = abs(walker_population(:,idet))

            ! The new population that we are aiming for. If this is not an
            ! integer then we will have to round up or down to an integer with
            ! an unbiased probability.
            new_population_target = abs(real(walker_population(:,idet),dp))/weight_altering_factors(excit_level)
            new_population = int(new_population_target, int_p)

            ! If new_population_target is not an integer, round it up or down
            ! with an unbiased probability. Do this for each replica.
            do ireplica = 1, sampling_size

                pextra = new_population_target(ireplica) - new_population(ireplica)

                if (pextra > depsilon) then
                    r = get_rand_close_open(rng)
                    if (r < pextra) new_population(ireplica) = new_population(ireplica) + 1_int_p
                end if

                ! Finally, update the walker population.
                walker_population(ireplica,idet) = sign(new_population(ireplica), walker_population(ireplica,idet))

            end do

            ! Update the total number of walkers.
            nparticles = nparticles + real(new_population - old_population, dp)/encoding_factor

        end do

        ! Call the annihilation routine to update the main walker list, as some
        ! sites will have become unoccupied and so need removing from the simulation.
        call remove_unoccupied_dets(rng)

    end subroutine update_sampling_weights

    subroutine output_and_alter_weights(max_number_excitations)

        ! This routine will alter and output the sampling weights used in importance 
        ! sampling. It uses the excitation distribution, calculated on the beta loop
        ! which has just finished, and finds the weights needed so that each excitation
        ! level will have roughly equal numbers of psips in the next loop. For example,
        ! to find the weights of psips on the 1st excitation level, divide the number of
        ! psips on the 1st excitation level by the number on the 0th level, then multiply
        ! the old sampling weight by this number to give the new weight. This can be used
        ! when the weights are being introduced gradually each beta loop, too. The weights
        ! are output and can then be used in future DMQMC runs.

        ! In:
        !    max_number_excitations: maximum number of excitations possible (see
        !       sys_t type in system for details).

        use fciqmc_data, only: dmqmc_sampling_probs, dmqmc_accumulated_probs
        use fciqmc_data, only: excit_distribution, finish_varying_weights
        use fciqmc_data, only: dmqmc_vary_weights, weight_altering_factors
        use parallel

        integer, intent(in) :: max_number_excitations

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
        ! start of each beta loop. This requires redefining weight_altering_factors to coincide
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

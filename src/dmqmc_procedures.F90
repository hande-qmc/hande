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

    subroutine init_dmqmc()

         use basis, only: basis_length, total_basis_length, bit_lookup, basis_lookup
         use calc, only: doing_dmqmc_calc, dmqmc_calc_type, dmqmc_energy, dmqmc_energy_squared
         use calc, only: dmqmc_staggered_magnetisation, dmqmc_correlation
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
         use system

         integer :: ierr, i, bit_position, bit_element

         number_dmqmc_estimators = 0

         allocate(trace(sampling_size), stat=ierr)
         call check_allocate('trace',sampling_size,ierr)
         trace = 0.0_p

         allocate(rdm_traces(sampling_size,nrdms), stat=ierr)
         call check_allocate('rdm_traces',sampling_size*nrdms,ierr)
         rdm_traces = 0.0_p

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
             allocate(dmqmc_accumulated_probs_old(0:max_number_excitations), stat=ierr)
             call check_allocate('dmqmc_accumulated_probs_old',max_number_excitations+1,ierr)
             dmqmc_accumulated_probs(0) = 1.0_p
             dmqmc_accumulated_probs_old = 1.0_p
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
             allocate(dmqmc_accumulated_probs_old(0:max_number_excitations), stat=ierr)
             call check_allocate('dmqmc_accumulated_probs_old',max_number_excitations+1,ierr)
             dmqmc_accumulated_probs = 1.0_p
             dmqmc_accumulated_probs_old = 1.0_p
             if (half_density_matrix) dmqmc_accumulated_probs(1:max_number_excitations) = 2.0_p
         end if

         ! If doing a reduced density matrix calculation, allocate and define the bit masks that
         ! have 1's at the positions referring to either subsystems A or B.
         if (doing_reduced_dm) call setup_rdm_arrays()

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

    subroutine setup_rdm_arrays()

        ! Setup the bit masks needed for RDM calculations. These are masks for the bits referring
        ! to either subsystem A or B. Also calculate the positions and elements of the sites
        ! in subsyetsm A, and finally allocate the RDM itself.

        use calc, only: ms_in, doing_dmqmc_calc, dmqmc_renyi_2
        use checking, only: check_allocate
        use errors
        use fciqmc_data, only: reduced_density_matrix, nrdms, calc_ground_rdm, calc_inst_rdm
        use fciqmc_data, only: replica_tricks, renyi_2, sampling_size
        use fciqmc_data, only: spawned_rdm_length, rdm_spawn
        use hash_table, only: alloc_hash_table
        use parallel, only: parent
        use spawn_data, only: alloc_spawn_t
        use system

        integer :: i, ierr, ipos, basis_find, size_spawned_rdm, total_size_spawned_rdm
        integer :: bit_position, bit_element

        ! For the Heisenberg model only currently.
        if (sys_global%system==heisenberg) then
            call find_rdm_masks()
        else
            call stop_all("setup_rdm_arrays","The use of RDMs is currently only implemented for the &
                           &Heisenberg model.")
        end if

        total_size_spawned_rdm = 0

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

            if (doing_dmqmc_calc(dmqmc_renyi_2)) then
                allocate(renyi_2(nrdms), stat=ierr)
                call check_allocate('renyi_2', nrdms, ierr)
                renyi_2 = 0.0_p
            end if
            if (calc_inst_rdm) then
                allocate(rdm_spawn(nrdms), stat=ierr)
                call check_allocate('rdm_spawn', nrdms, ierr)

                size_spawned_rdm = (rdms(i)%rdm_basis_length*2+sampling_size)*i0_length/8
                total_size_spawned_rdm = total_size_spawned_rdm + size_spawned_rdm
                if (spawned_rdm_length < 0) then
                    ! Given in MB.  Convert.
                    ! Note that we store 2 arrays.
                    spawned_rdm_length = int((-real(spawned_rdm_length,p)*10**6)/(2*size_spawned_rdm))
                end if

                ! Note the initiator approximation is not implemented for density matrix calculations.
                call alloc_spawn_t(rdms(i)%rdm_basis_length*2, sampling_size, .false., &
                                 spawned_rdm_length, 7, rdm_spawn(i)%spawn)
                ! Hard code hash table collision limit for now.  This should
                ! give an ok performance...
                ! We will only use the first 2*rdm_basis_length elements for the
                ! hash, even though rdm_spawn%spawn%sdata is larger than that in
                ! the first dimension...
                call alloc_hash_table(nint(real(spawned_rdm_length)/3), 3, rdms(i)%rdm_basis_length*2, &
                                      0, 0, 17, rdm_spawn(i)%ht, rdm_spawn(i)%spawn%sdata)
            end if
        end do

        if (parent) write (6,'(1X,a58,f7.2)') 'Memory allocated per core for the spawned RDM lists (MB): ', &
                total_size_spawned_rdm*real(2*spawned_rdm_length,p)/10**6

        ! For an ms = 0 subspace, assuming less than or exactly half the spins in the subsystem are in
        ! the subsystem, then any combination of spins can occur in the subsystem, from all spins down
        ! to all spins up. Hence the total size of the reduced density matrix will be 2**(number of spins
        ! in subsystem A).
        if (calc_ground_rdm) then
            if (ms_in == 0 .and. rdms(1)%A_nsites <= floor(real(sys_global%lattice%nsites,p)/2.0_p)) then
                allocate(reduced_density_matrix(2**rdms(1)%A_nsites,2**rdms(1)%A_nsites), stat=ierr)
                call check_allocate('reduced_density_matrix', 2**(2*rdms(1)%A_nsites),ierr)
                reduced_density_matrix = 0.0_p
            else
                if (ms_in /= 0) then
                    call stop_all("setup_rdm_arrays","Reduced density matrices can only be used for Ms=0 &
                                   &calculations.")
                else if (rdms(1)%A_nsites > floor(real(sys_global%lattice%nsites,p)/2.0_p)) then
                    call stop_all("setup_rdm_arrays","Reduced density matrices can only be used for subsystems &
                                  &whose size is less than half the total system size.")
                end if
            end if
        end if
        
    end subroutine setup_rdm_arrays

    subroutine find_rdm_masks()

        use basis, only: basis_length, bit_lookup, basis_lookup, basis_fns, nbasis
        use checking, only: check_allocate, check_deallocate
        use errors
        use fciqmc_data, only: nrdms, nsym_vec
        use hubbard_real, only: map_vec_to_cell
        use system

        integer :: i, j, k, l, ipos, ierr
        integer :: basis_find, bit_position, bit_element
        integer :: r(sys_global%lattice%ndim), nvecs(3), A_mask(basis_length)
        real(p) :: v(sys_global%lattice%ndim), test_vec(sys_global%lattice%ndim), temp_vec(sys_global%lattice%ndim)
        real(p), allocatable :: trans_vecs(:,:)
        integer :: scale_fac

        ! The maximum number of translational symmetry vectors is nsites (for
        ! the case of a non-tilted lattice), so allocate this much storage.
        allocate(trans_vecs(sys_global%lattice%ndim,sys_global%lattice%nsites),stat=ierr)
        call check_allocate('trans_vecs',sys_global%lattice%ndim*sys_global%lattice%nsites,ierr)

        ! The number of symmetry vectors in each direction.
        nvecs = 0
        ! The total number of symmetry vectors.
        nsym_vec = 0

        do i = 1, sys_global%lattice%ndim
            scale_fac = maxval(abs(sys_global%lattice%lattice(:,i)))
            v = real(sys_global%lattice%lattice(:,i),p)/real(scale_fac,p)

            do j = 1, scale_fac-1
                test_vec = v*j
                if (all(.not. (abs(test_vec-real(nint(test_vec),p)) > 0.0_p) )) then
                    ! If test_vec has all integer components.
                    ! This is a symmetry vector, so store it.
                    nvecs(i) = nvecs(i) + 1
                    nsym_vec = nsym_vec + 1
                    trans_vecs(:,nsym_vec) = test_vec
                end if
            end do
        end do

        ! Next, add all combinations of the above generated vectors to form a closed group.

        ! Add all pairs of the above vectors.
        do i = 1, nvecs(1)
            do j = nvecs(1)+1, sum(nvecs)
                nsym_vec = nsym_vec + 1
                trans_vecs(:,nsym_vec) = trans_vecs(:,i)+trans_vecs(:,j)
            end do
        end do
        do i = nvecs(1)+1, nvecs(1)+nvecs(2)
            do j = nvecs(1)+nvecs(2)+1, sum(nvecs)
                nsym_vec = nsym_vec + 1
                trans_vecs(:,nsym_vec) = trans_vecs(:,i)+trans_vecs(:,j)
            end do
        end do

        ! Add all triples of the above vectors.
        do i = 1, nvecs(1)
            do j = nvecs(1)+1, nvecs(1)+nvecs(2)
                do k = nvecs(1)+nvecs(2)+1, sum(nvecs)
                    nsym_vec = nsym_vec + 1
                    trans_vecs(:,nsym_vec) = trans_vecs(:,i)+trans_vecs(:,j)+trans_vecs(:,k)
                end do
            end do
        end do

        ! Include the identity transformation vector in the first slot.
        trans_vecs(:,2:nsym_vec+1) = trans_vecs(:,1:nsym_vec)
        trans_vecs(:,1) = 0
        nsym_vec = nsym_vec + 1

        ! Allocate the RDM arrays.
        do i = 1, nrdms
            allocate(rdms(i)%B_masks(basis_length,nsym_vec), stat=ierr)
            call check_allocate('rdms(i)%B_masks', nsym_vec*basis_length,ierr)
            allocate(rdms(i)%bit_pos(rdms(i)%A_nsites,nsym_vec,2), stat=ierr)
            call check_allocate('rdms(i)%bit_pos', nsym_vec*rdms(i)%A_nsites*2,ierr)
            rdms(i)%B_masks = 0
            rdms(i)%bit_pos = 0
        end do

        ! Run through every site on every subsystem and add every translational symmetry vector.
        do i = 1, nrdms ! Over every subsystem.
            do j = 1, nsym_vec ! Over every symmetry vector.
                A_mask = 0
                do k = 1, rdms(i)%A_nsites ! Over every site in the subsystem.
                    r = basis_fns(rdms(i)%subsystem_A(k))%l
                    r = r + nint(trans_vecs(:,j))
                    ! If r is outside the cell considered in this simulation, shift it by the
                    ! appropriate lattice vector so that it is in this cell.
                    call map_vec_to_cell(r)
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

        deallocate(trans_vecs,stat=ierr)
        call check_deallocate('trans_vecs',ierr)

    end subroutine find_rdm_masks

    subroutine random_distribution_heisenberg(rng, ireplica)

        ! For the Heisenberg model only. Distribute the initial number of psips
        ! along the main diagonal. Each diagonal element should be chosen
        ! with the same probability.

        ! Currently this creates psips with Ms = ms_in only.

        ! If we have number of sites = nsites and total spin value = ms_in,
        ! then the number of up spins is equal to up_spins = (ms_in + nsites)/2.

        ! Start from state with all spins down, then choose the above number of
        ! spins to flip up with equal probability.

        ! In/Out:
        !    rng: random number generator.

        use basis, only: nbasis, basis_length, bit_lookup
        use calc, only: ms_in
        use dSFMT_interface, only:  dSFMT_t, get_rand_close_open
        use fciqmc_data, only: D0_population
        use parallel
        use system

        type(dSFMT_t), intent(inout) :: rng
        integer, intent(in) :: ireplica
        integer :: i, up_spins, rand_basis, bits_set
        integer :: bit_element, bit_position, npsips
        integer(i0) :: f(basis_length)
        real(dp) :: rand_num

        up_spins = (ms_in+sys_global%lattice%nsites)/2
        npsips = int(D0_population/nprocs)
        ! If the initial number of psips does not split evenly between all processors,
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
            call create_particle(f,f,1,ireplica)

        end do

    end subroutine random_distribution_heisenberg

    subroutine create_particle(f1, f2, nspawn, particle_type)

        ! Create a psip on a diagonal element of the density matrix by adding
        ! it to the spawned walkers list. This list can then be sorted correctly
        ! by the direct_annihilation routine.

        ! In:
        !    f_new: Bit string representation of index of the diagonal
        !           element upon which a new psip shall be placed.

        use hashing
        use basis, only: basis_length, total_basis_length
        use fciqmc_data, only: qmc_spawn
        use parallel

        integer(i0), intent(in) :: f1(basis_length), f2(basis_length)
        integer, intent(in) :: nspawn, particle_type
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
                                total_basis_length, qmc_spawn%hash_seed), nprocs)
#endif

        ! Move to the next position in the spawning array.
        qmc_spawn%head(0,iproc_spawn) = qmc_spawn%head(0,iproc_spawn) + 1

        ! Set info in spawning array.
        ! Zero it as not all fields are set.
        qmc_spawn%sdata(:,qmc_spawn%head(0,iproc_spawn)) = 0
        ! indices 1 to total_basis_length store the bitstring.
        qmc_spawn%sdata(:(total_basis_length),qmc_spawn%head(0,iproc_spawn)) = f_new
        ! The final index stores the number of psips created.
        qmc_spawn%sdata((total_basis_length)+particle_type,qmc_spawn%head(0,iproc_spawn)) = nspawn

    end subroutine create_particle

    subroutine decode_dm_bitstring(f, irdm, isym)

        ! This function maps a full DMQMC bitstring to two bitstrings encoding the subsystem-A
        ! RDM bitstrings. These resulting bitstrings are stored in the end1 and end2 components
        ! of rdms(irdm).
        
        ! Crucially, the mapping is performed so that, if there are two subsystems which are
        ! equivalent by symmetry, then equivalent sites in those two subsystems will be mapped to
        ! the same RDM bitstrings. This is clearly a requirement to obtain a correct representation
        ! of the RDM when using translational symmetry.

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
        use fciqmc_data, only: nparticles, sampling_size
        use dSFMT_interface, only: dSFMT_t, get_rand_close_open

        type(dSFMT_t), intent(inout) :: rng
        integer :: idet, ireplica, excit_level, nspawn, sign_factor, old_population
        real(p) :: new_factor
        real(dp) :: rand_num, prob

        ! Alter weights for the next iteration.
        dmqmc_accumulated_probs = dble(dmqmc_accumulated_probs)*weight_altering_factors

        ! When the weights for an excitation level are increased by a factor, the number
        ! of psips on that level has to decrease by the same factor, else the wavefunction
        ! which the psips represent will not be the correct importance sampled wavefunction
        ! for the new weights. The code below loops over every psips and destroys (or creates)
        ! it with the appropriate probability.
        do idet = 1, tot_walkers

            excit_level = get_excitation_level(walker_dets(1:basis_length,idet),&
                    walker_dets(basis_length+1:total_basis_length,idet))

            do ireplica = 1, sampling_size
                old_population = abs(walker_population(ireplica,idet))
                rand_num = get_rand_close_open(rng)
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
                nspawn = sign(nspawn,walker_population(ireplica,idet)*sign_factor)
                ! Update the population on this determinant.
                walker_population(ireplica,idet) = walker_population(ireplica,idet) + nspawn
                ! Update the total number of walkers.
                nparticles(ireplica) = nparticles(ireplica) - old_population + &
                        abs(walker_population(ireplica,idet))
            end do
        end do

        ! Call the annihilation routine to update the main walker list, as some
        ! sites will have become unoccupied and so need removing from the simulation.
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
        use system

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

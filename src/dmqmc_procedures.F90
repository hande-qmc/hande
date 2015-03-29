module dmqmc_procedures

use const
implicit none

contains

    subroutine init_dmqmc(sys, qmc_in, dmqmc_in, nreplicas, qs, dmqmc_estimates, weighted_sampling)

        ! In:
        !    sys: system being studied.
        !    nreplicas: number of replicas being used.
        ! In/Out:
        !    qmc_in: Input options relating to QMC methods.
        !    dmqmc_in: Input options relating to DMQMC.
        !    dmqmc_estimates: type containing estimates for observables.
        !    weighted_sampling: type containing information for weighted
        !        sampling.

        use calc, only: doing_dmqmc_calc, dmqmc_calc_type, dmqmc_energy, dmqmc_energy_squared
        use calc, only: dmqmc_staggered_magnetisation, dmqmc_correlation, dmqmc_full_r2
        use checking, only: check_allocate
        use fciqmc_data
        use system, only: sys_t

        use qmc_data, only: qmc_in_t, qmc_state_t
        use dmqmc_data, only: dmqmc_in_t, dmqmc_estimates_t, dmqmc_weighted_sampling_t

        type(sys_t), intent(in) :: sys
        integer, intent(in) :: nreplicas
        type(qmc_in_t), intent(in) :: qmc_in
        type(dmqmc_in_t), intent(inout) :: dmqmc_in
        type(qmc_state_t), intent(inout) :: qs
        type(dmqmc_estimates_t), intent(inout) :: dmqmc_estimates
        type(dmqmc_weighted_sampling_t), intent(inout) :: weighted_sampling

        integer :: ierr, i, bit_position, bit_element

        allocate(dmqmc_estimates%trace(nreplicas), stat=ierr)
        call check_allocate('dmqmc_estimates%trace',size(dmqmc_estimates%trace),ierr)
        dmqmc_estimates%trace = 0.0_p

        nrdms = dmqmc_in%rdm%nrdms
        allocate(rdm_traces(nreplicas,nrdms), stat=ierr)
        call check_allocate('rdm_traces',size(rdm_traces),ierr)
        rdm_traces = 0.0_p

        ! If calculating a correlaton function then set up the necessary bit
        ! mask. This has a bit set for each of the two sites/orbitals being
        ! considered in the correlation function.
        if (doing_dmqmc_calc(dmqmc_correlation)) then
            allocate(dmqmc_in%correlation_mask(1:sys%basis%string_len), stat=ierr)
            call check_allocate('dmqmc_in%correlation_mask',sys%basis%string_len,ierr)
            dmqmc_in%correlation_mask = 0_i0
            do i = 1, 2
                bit_position = sys%basis%bit_lookup(1,dmqmc_in%correlation_sites(i))
                bit_element = sys%basis%bit_lookup(2,dmqmc_in%correlation_sites(i))
                dmqmc_in%correlation_mask(bit_element) = ibset(dmqmc_in%correlation_mask(bit_element), bit_position)
            end do
        end if

        if (dmqmc_in%calc_excit_dist .or. dmqmc_in%find_weights) then
            allocate(dmqmc_estimates%excit_dist(0:sys%max_number_excitations), stat=ierr)
            call check_allocate('dmqmc_estimates%excit_dist',sys%max_number_excitations+1,ierr)
            dmqmc_estimates%excit_dist = 0.0_p
        end if

        ! When using an importance sampled initial density matrix we use then
        ! unsymmetrised version of Bloch's equation. This means we don't have
        ! to worry about these factors of 1/2.
        if (.not. dmqmc_in%propagate_to_beta) then
            ! In DMQMC we want the spawning probabilities to have an extra factor
            ! of a half, because we spawn from two different ends with half
            ! probability. To avoid having to multiply by an extra variable in
            ! every spawning routine to account for this, we multiply the time
            ! step by 0.5 instead, then correct this in the death step (see below).
            qs%tau = qs%tau*0.5_p
            ! Set dmqmc_factor to 2 so that when probabilities in death.f90 are
            ! multiplied by this factor it cancels the factor of 0.5 introduced
            ! into the timestep in DMQMC.cThis factor is also used in updated the
            ! shift, where the true tau is needed.
            qs%dmqmc_factor = 2.0_p
        end if
        ! Set the timestep to be the appropriate factor of ef so that results
        ! are at temperatures commensurate(ish) with the reduced (inverse) temperature
        ! Beta = 1\Theta = T/T_F, where T_F is the Fermi-Temperature. Also need
        ! to set the appropriate beta = Beta / T_F.
        if (dmqmc_in%fermi_temperature) then
            qs%tau = qs%tau / sys%ueg%ef
            dmqmc_in%init_beta = dmqmc_in%init_beta / sys%ueg%ef
        end if


        if (dmqmc_in%weighted_sampling) then
            ! sampling_probs stores the factors by which probabilities
            ! are to be reduced when spawning away from the diagonal. The trial
            ! function required from these probabilities, for use in importance
            ! sampling, is actually that of the accumulated factors, ie, if
            ! sampling_probs = (a, b, c, ...) then
            ! dmqmc_accumulated_factors = (1, a, ab, abc, ...).
            ! This is the array which we need to create and store.
            ! sampling_probs is no longer needed and so can be
            ! deallocated. Also, the user may have only input factors for the
            ! first few excitation levels, but we need to store factors for all
            ! levels, as done below.
            if (.not.allocated(dmqmc_in%sampling_probs)) then
                allocate(dmqmc_in%sampling_probs(1:sys%max_number_excitations), stat=ierr)
                call check_allocate('dmqmc_in%sampling_probs',sys%max_number_excitations,ierr)
                dmqmc_in%sampling_probs = 1.0_p
            end if
            allocate(weighted_sampling%probs(0:sys%max_number_excitations), stat=ierr)
            call check_allocate('weighted_sampling%probs',sys%max_number_excitations+1,ierr)
            allocate(weighted_sampling%probs_old(0:sys%max_number_excitations), stat=ierr)
            call check_allocate('weighted_sampling%probs_old',sys%max_number_excitations+1,ierr)
            weighted_sampling%probs(0) = 1.0_p
            weighted_sampling%probs_old = 1.0_p
            do i = 1, size(dmqmc_in%sampling_probs)
            weighted_sampling%probs(i) =  weighted_sampling%probs(i-1)*dmqmc_in%sampling_probs(i)
            end do
            weighted_sampling%probs(size(dmqmc_in%sampling_probs)+1:sys%max_number_excitations) = &
               weighted_sampling%probs(size(dmqmc_in%sampling_probs))
            if (dmqmc_in%vary_weights) then
                ! Allocate an array to store the factors by which the weights
                ! will change each iteration.
                allocate(weighted_sampling%altering_factors(0:sys%max_number_excitations), stat=ierr)
                call check_allocate('weighted_sampling%altering_factors',sys%max_number_excitations+1,ierr)
                weighted_sampling%altering_factors = real(weighted_sampling%probs,dp)**(1/real(dmqmc_in%finish_varying_weights,dp))
                ! If varying the weights, start the accumulated probabilties
                ! as all 1.0 initially, and then alter them gradually later.
                weighted_sampling%probs = 1.0_p
            end if
        else
            ! If not using the importance sampling procedure, turn it off by
            ! setting all amplitudes to 1.0 in the relevant arrays.
            allocate(weighted_sampling%probs(0:sys%max_number_excitations), stat=ierr)
            call check_allocate('weighted_sampling%probs',sys%max_number_excitations+1,ierr)
            allocate(weighted_sampling%probs_old(0:sys%max_number_excitations), stat=ierr)
            call check_allocate('weighted_sampling%probs_old',sys%max_number_excitations+1,ierr)
            weighted_sampling%probs = 1.0_p
            weighted_sampling%probs_old = 1.0_p
        end if

        ! If doing a reduced density matrix calculation, allocate and define
        ! the bit masks that have 1's at the positions referring to either
        ! subsystems A or B.
        if (dmqmc_in%rdm%doing_rdm) call setup_rdm_arrays(sys, .true., dmqmc_estimates%rdm_info, dmqmc_estimates%ground_rdm%rdm, &
                                                          qmc_in, dmqmc_in%rdm, nreplicas)

    end subroutine init_dmqmc

    subroutine setup_rdm_arrays(sys, called_from_dmqmc, rdm_info, ground_rdm, qmc_in, rdm_in, nreplicas)

        ! Setup the bit masks needed for RDM calculations. These are masks for
        ! the bits referring to either subsystem A or B. Also calculate the
        ! positions and elements of the sites in subsyetsm A, and finally
        ! allocate the RDM itself (including allocating the instances of the RDM
        ! spawning arrays and hash tables, for instantaneous RDM calculations).

        ! In:
        !    sys: system being studied.
        !    called_from_dmqmc: This variable should be true if this routine is
        !        called from DMQMC, false otherwise. This routine is also used
        !        by the FCI code, in which case qmc_in, rdm_in and nreplicas
        !        will not be passed in.
        !    qmc_in (optional): Input options relating to QMC methods.  Only
        !         needed for spawn_cutoff and if calc_inst_rdm is true.
        !    rdm_in (optional): Input options relating to reduced density matrices.
        !    nreplicas (optional): number of replicas being used.  Must be
        !        specified if qmc_in is.
        ! Out:
        !     ground_rdm: The array used to store the RDM in ground-state
        !         calculations.
        ! In/Out:
        !     rdm_info: information relating to RDM subsystems being studied.

        use calc, only: ms_in, doing_dmqmc_calc, dmqmc_rdm_r2, use_mpi_barriers
        use checking, only: check_allocate
        use dmqmc_data, only: rdm_t
        use errors
        use fciqmc_data, only: nrdms
        use fciqmc_data, only: renyi_2, real_bit_shift
        use fciqmc_data, only: rdm_spawn
        use hash_table, only: alloc_hash_table
        use parallel, only: parent
        use spawn_data, only: alloc_spawn_t
        use system, only: sys_t, heisenberg
        use utils, only: int_fmt

        use dmqmc_data, only: dmqmc_rdm_in_t
        use qmc_data, only: qmc_in_t

        type(sys_t), intent(in) :: sys
        logical, intent(in) :: called_from_dmqmc
        type(rdm_t), intent(inout) :: rdm_info(:)
        real(p), allocatable, intent(out) :: ground_rdm(:,:)
        type(qmc_in_t), intent(in), optional :: qmc_in
        type(dmqmc_rdm_in_t), intent(in), optional :: rdm_in
        integer, intent(in), optional :: nreplicas

        integer :: i, ierr, ipos, basis_find, size_spawned_rdm, total_size_spawned_rdm
        integer :: bit_position, bit_element, nbytes_int, spawn_length_loc
        logical :: calc_ground_rdm

        ! If this routine was not called from DMQMC then we must be doing a
        ! ground state RDM calculation.
        calc_ground_rdm = .not. called_from_dmqmc
        ! This should only be present if called_from_dmqmc is true.
        if (present(rdm_in)) calc_ground_rdm = rdm_in%calc_ground_rdm

        ! For the Heisenberg model only currently.
        if (sys%system == heisenberg) then
            call find_rdm_masks(sys, rdm_info)
        else
            call stop_all("setup_rdm_arrays","The use of RDMs is currently only implemented for &
                           &the Heisenberg model.")
        end if

        ! Loop over all subsystems for which we are calculating RDMs.
        ! Setup the rdms array.
        do i = 1, nrdms
            ! Initialise the instance of the rdm type for this subsystem.
            rdm_info(i)%string_len = ceiling(real(rdm_info(i)%A_nsites)/i0_length)
            allocate(rdm_info(i)%end1(rdm_info(i)%string_len), stat=ierr)
            call check_allocate('rdm_info(i)%end1', rdm_info(i)%string_len, ierr)
            allocate(rdm_info(i)%end2(rdm_info(i)%string_len), stat=ierr)
            call check_allocate('rdm_info(i)%end2', rdm_info(i)%string_len, ierr)
            rdm_info(i)%end1 = 0_i0
            rdm_info(i)%end2 = 0_i0

            ! With the calc_ground_rdm option, the entire RDM is allocated. If
            ! the following condition is met then the number of rows is greater
            ! than the maximum integer accessible. This would clearly be too
            ! large, so abort in this case.
            if (calc_ground_rdm .and. rdm_info(i)%string_len > 1) call stop_all("setup_rdm_arrays",&
                "A requested RDM is too large for all indices to be addressed by a single integer.")
        end do

        ! For an ms = 0 subspace, assuming less than or exactly half the spins
        ! in the subsystem are in the subsystem, then any combination of spins
        ! can occur in the subsystem, from all spins down to all spins up. Hence
        ! the total size of the reduced density matrix will be 2**(number of
        ! spins in subsystem A).
        if (calc_ground_rdm) then
            if (ms_in == 0 .and. rdm_info(1)%A_nsites <= floor(real(sys%lattice%nsites,p)/2.0_p)) then
                allocate(ground_rdm(2**rdm_info(1)%A_nsites,2**rdm_info(1)%A_nsites), stat=ierr)
                call check_allocate('ground_rdm', 2**(2*rdm_info(1)%A_nsites),ierr)
                ground_rdm = 0.0_p
            else
                if (ms_in /= 0) then
                    call stop_all("setup_rdm_arrays","Reduced density matrices can only be used for Ms=0 &
                                   &calculations.")
                else if (rdm_info(1)%A_nsites > floor(real(sys%lattice%nsites,p)/2.0_p)) then
                    call stop_all("setup_rdm_arrays","Reduced density matrices can only be used for subsystems &
                                  &whose size is less than half the total system size.")
                end if
            end if
        end if

        ! Setup code relating to instantaneous RDMs.
        if (present(rdm_in)) then
            ! Create the instances of the rdm_spawn_t type for instantaneous RDM
            ! calculatons.
            if (rdm_in%calc_inst_rdm) then
                allocate(rdm_spawn(nrdms), stat=ierr)
                call check_allocate('rdm_spawn', nrdms, ierr)

                ! If calculating Renyi entropy (S2).
                if (doing_dmqmc_calc(dmqmc_rdm_r2)) then
                    allocate(renyi_2(nrdms), stat=ierr)
                    call check_allocate('renyi_2', nrdms, ierr)
                    renyi_2 = 0.0_p
                end if

                total_size_spawned_rdm = 0
                nbytes_int = bit_size(i)/8

                do i = 1, nrdms
                    ! Allocate the spawn_t and hash table instances for this RDM.
                    if (.not.present(qmc_in)) call stop_all('setup_rdm_arrays', 'qmc_in not supplied.')
                    size_spawned_rdm = (rdm_info(i)%string_len*2+nreplicas)*int_s_length/8
                    total_size_spawned_rdm = total_size_spawned_rdm + size_spawned_rdm

                    spawn_length_loc = rdm_in%spawned_length

                    if (spawn_length_loc < 0) then
                        ! Given in MB.  Convert.
                        ! Note that the factor of 2 is because two spawning arrays
                        ! are stored, and 21*nbytes_int is added because there are
                        ! 21 integers in the hash table for each spawned rdm slot.
                        ! 21 was found to be appropriate after testing.
                        spawn_length_loc = int((-real(spawn_length_loc,p)*10**6)/&
                                              (2*size_spawned_rdm + 21*nbytes_int))
                    end if

                    ! Note the initiator approximation is not implemented for density matrix calculations.
                    call alloc_spawn_t(rdm_info(i)%string_len*2, nreplicas, .false., &
                                         spawn_length_loc, qmc_in%spawn_cutoff, real_bit_shift, &
                                         27, use_mpi_barriers, rdm_spawn(i)%spawn)
                    ! Hard code hash table collision limit for now.  The length of
                    ! the table is three times as large as the spawning arrays and
                    ! each hash value can have 7 clashes. This was found to give
                    ! reasonable performance.
                    call alloc_hash_table(3*spawn_length_loc, 7, rdm_info(i)%string_len*2, &
                                         0, 0, 17, rdm_spawn(i)%ht, rdm_spawn(i)%spawn%sdata)
                end do

                if (parent) then
                    write (6,'(1X,a58,f7.2)') 'Memory allocated per core for the spawned RDM lists (MB): ', &
                        total_size_spawned_rdm*real(2*spawn_length_loc,p)/10**6
                    write (6,'(1X,a49,'//int_fmt(spawn_length_loc,1)//',/)') &
                        'Number of elements per core in spawned RDM lists:', spawn_length_loc
                end if
            end if
        end if

    end subroutine setup_rdm_arrays

    subroutine find_rdm_masks(sys, rdm_info)

        ! Initialise bit masks for converting a density matrix basis function
        ! into its corresponding reduced density matrix basis function. Bit
        ! masks will be calculated for the subsystems requested by the user, and
        ! also for all subsystems which are equivalent by translational symmetry.

        ! In:
        !    sys: system being studied.

        use checking, only: check_allocate, check_deallocate
        use dmqmc_data, only: rdm_t
        use errors
        use fciqmc_data, only: nrdms, nsym_vec
        use real_lattice, only: find_translational_symmetry_vecs, map_vec_to_cell, enumerate_lattice_vectors
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        type(rdm_t), intent(inout) :: rdm_info(:)

        integer :: i, j, k, l, ipos, ierr
        integer :: basis_find, bit_position, bit_element
        integer(i0) :: A_mask(sys%basis%string_len)
        real(p), allocatable :: sym_vecs(:,:)
        integer :: r(sys%lattice%ndim)
        integer, allocatable :: lvecs(:,:)

        ! Return all translational symmetry vectors in sym_vecs.
        call find_translational_symmetry_vecs(sys, sym_vecs, nsym_vec)

        ! Allocate the RDM arrays.
        do i = 1, nrdms
            if (rdm_info(i)%A_nsites == sys%lattice%nsites) then
                call stop_all('find_rdm_masks','You are attempting to use the full density matrix &
                              &as an RDM. This is not supported. You should use the &
                              &dmqmc_full_renyi_2 option to calculate the Renyi entropy of the &
                              &whole lattice.')
            else
                allocate(rdm_info(i)%B_masks(sys%basis%string_len,nsym_vec), stat=ierr)
                call check_allocate('rdm_info(i)%B_masks', nsym_vec*sys%basis%string_len,ierr)
                allocate(rdm_info(i)%bit_pos(rdm_info(i)%A_nsites,nsym_vec,2), stat=ierr)
                call check_allocate('rdm_info(i)%bit_pos', nsym_vec*rdm_info(i)%A_nsites*2,ierr)
            end if
            rdm_info(i)%B_masks = 0_i0
            rdm_info(i)%bit_pos = 0
        end do

        ! Run through every site on every subsystem and add every translational
        ! symmetry vector.
        allocate(lvecs(sys%lattice%ndim,3**sys%lattice%ndim), stat=ierr)
        call check_allocate('lvecs', size(lvecs), ierr)
        call enumerate_lattice_vectors(sys%lattice, lvecs)
        do i = 1, nrdms ! Over every subsystem.
            do j = 1, nsym_vec ! Over every symmetry vector.
                A_mask = 0_i0
                do k = 1, rdm_info(i)%A_nsites ! Over every site in the subsystem.
                    r = sys%basis%basis_fns(rdm_info(i)%subsystem_A(k))%l
                    r = r + nint(sym_vecs(:,j))
                    ! If r is outside the cell considered in this simulation,
                    ! shift it by the appropriate lattice vector so that it is
                    ! in this cell.
                    call map_vec_to_cell(sys%basis%nbasis, sys%basis%basis_fns, sys%lattice%ndim, lvecs, r)
                    ! Now need to find which basis function this site 
                    ! corresponds to. Simply loopover all basis functions and
                    ! check...
                    do l = 1, sys%basis%nbasis
                        if (all(sys%basis%basis_fns(l)%l == r)) then
                            ! Found the correct basis function!
                            bit_position = sys%basis%bit_lookup(1,l)
                            bit_element = sys%basis%bit_lookup(2,l)
                            A_mask(bit_element) = ibset(A_mask(bit_element), bit_position)
                            rdm_info(i)%bit_pos(k,j,1) = bit_position
                            rdm_info(i)%bit_pos(k,j,2) = bit_element
                        end if
                    end do
                end do

                ! We cannot just flip the mask for system A to get that for
                ! system B, because the trailing bits on the end don't refer to
                ! anything and should be set to 0. So, first set these to 1 and
                ! then flip all the bits.
                rdm_info(i)%B_masks(:,j) = A_mask
                do ipos = 0, i0_end
                    basis_find = sys%basis%basis_lookup(ipos, sys%basis%string_len)
                    if (basis_find == 0) then
                        rdm_info(i)%B_masks(sys%basis%string_len,j) = ibset(rdm_info(i)%B_masks(sys%basis%string_len,j),ipos)
                    end if
                end do
                rdm_info(i)%B_masks(:,j) = not(rdm_info(i)%B_masks(:,j))

            end do
        end do

        deallocate(lvecs, stat=ierr)
        call check_deallocate('lvecs', ierr)
        deallocate(sym_vecs,stat=ierr)
        call check_deallocate('sym_vecs',ierr)

    end subroutine find_rdm_masks

    subroutine create_initial_density_matrix(rng, sys, qmc_in, dmqmc_in, reference, annihilation_flags, &
                                             target_nparticles_tot, psip_list, spawn, nload_slots)

        ! Create a starting density matrix by sampling the elements of the
        ! (unnormalised) identity matrix. This is a sampling of the
        ! (unnormalised) infinite-temperature density matrix. This is done by
        ! picking determinants/spin configurations with uniform probabilities in
        ! the space being considered.

        ! In/Out:
        !    rng: random number generator.
        !    psip_list: particle_t object conaining sample of initial density
        !       matrix on output.
        !    spawn: spawn_t object.  Reset on input and output.  Used to
        !       communicate the generated particles.
        ! In:
        !    sys: system being studied.
        !    qmc_in: input options relating to QMC methods.
        !    dmqmc_in: input options relating to DMQMC.
        !    reference: current reference determinant.
        !    annihilation_flags: calculation specific annihilation flags.
        !    target_nparticles_tot: The total number of psips to attempt to
        !        generate across all processes.
        !    nload_slots: number of load balancing slots (per processor).
        ! Out:
        !    nparticles_tot: The total number of psips in the generated density
        !        matrix across all processes, for all replicas.

        use annihilation, only: direct_annihilation
        use dSFMT_interface, only:  dSFMT_t, get_rand_close_open
        use errors
        use fciqmc_data, only: real_factor
        use parallel
        use system, only: sys_t, heisenberg, ueg, hub_k, hub_real
        use utils, only: binom_r
        use qmc_common, only: redistribute_particles
        use qmc_data, only: qmc_in_t, reference_t, particle_t, annihilation_flags_t
        use spawn_data, only:spawn_t
        use dmqmc_data, only: dmqmc_in_t
        use calc, only: sym_in

        type(dSFMT_t), intent(inout) :: rng
        type(sys_t), intent(in) :: sys
        type(qmc_in_t), intent(in) :: qmc_in
        type(dmqmc_in_t), intent(in) :: dmqmc_in
        type(reference_t), intent(in) :: reference
        type(annihilation_flags_t), intent(in) :: annihilation_flags
        integer(int_64), intent(in) :: target_nparticles_tot
        type(particle_t), intent(inout) :: psip_list
        type(spawn_t), intent(inout) :: spawn
        integer, intent(in) :: nload_slots

        real(p) :: nparticles_temp(psip_list%nspaces)
        integer :: nel, ireplica, ierr
        integer(int_64) :: npsips_this_proc, npsips
        real(dp) :: total_size, sector_size
        real(dp) :: r, prob

        npsips_this_proc = target_nparticles_tot/nprocs
        ! If the initial number of psips does not split evenly between all
        ! processors, add the leftover psips to the first processors in order.
        if (target_nparticles_tot-(nprocs*npsips_this_proc) > iproc) &
              npsips_this_proc = npsips_this_proc + 1_int_64

        nparticles_temp = 0.0_p

        do ireplica = 1, psip_list%nspaces
            select case(sys%system)
            case(heisenberg)
                if (dmqmc_in%all_spin_sectors) then
                    ! The size (number of configurations) of all symmetry
                    ! sectors combined.
                    total_size = 2.0_dp**(real(sys%lattice%nsites,dp))

                    do nel = 0, sys%lattice%nsites
                        ! The size of this symmetry sector alone.
                        sector_size = binom_r(sys%lattice%nsites, nel)
                        prob = real(npsips_this_proc,dp)*sector_size/total_size
                        npsips = floor(prob, int_64)
                        ! If there are a non-integer number of psips to be
                        ! spawned in this sector then add an extra psip with the
                        ! required probability.
                        prob = prob - npsips
                        r = get_rand_close_open(rng)
                        if (r < prob) npsips = npsips + 1_int_64

                        nparticles_temp(ireplica) = nparticles_temp(ireplica) + real(npsips, p)
                        call random_distribution_heisenberg(rng, sys%basis, nel, npsips, ireplica, spawn)
                    end do
                else
                    ! This process will always create excatly the target number
                    ! of psips.
                    call random_distribution_heisenberg(rng, sys%basis, sys%nel, npsips_this_proc, ireplica, spawn)
                end if
            case(ueg, hub_k)
                if (dmqmc_in%propagate_to_beta) then
                    ! Initially distribute psips along the diagonal according to
                    ! a guess.
                    if (dmqmc_in%grand_canonical_initialisation) then
                        call init_grand_canonical_ensemble(sys, dmqmc_in, sym_in, npsips_this_proc, spawn, rng)
                    else
                        call init_uniform_ensemble(sys, npsips_this_proc, sym_in, reference%f0, ireplica, &
                                                   dmqmc_in%all_sym_sectors, rng, spawn)
                    end if
                    ! Perform metropolis algorithm on initial distribution so
                    ! that we are sampling the trial density matrix.
                    if (dmqmc_in%metropolis_attempts > 0) call initialise_dm_metropolis(sys, rng, qmc_in, dmqmc_in, &
                                                                               npsips_this_proc, sym_in, ireplica, spawn)
                else
                    call random_distribution_electronic(rng, sys, sym_in, npsips_this_proc, ireplica, &
                                                        dmqmc_in%all_sym_sectors, spawn)
                end if
            case(hub_real)
                call random_distribution_electronic(rng, sys, sym_in, npsips_this_proc, ireplica, dmqmc_in%all_sym_sectors, spawn)
            case default
                call stop_all('create_initial_density_matrix','DMQMC not implemented for this system.')
            end select
        end do

        ! Finally, count the total number of particles across all processes.
        if (dmqmc_in%all_spin_sectors) then
#ifdef PARALLEL
            call mpi_allreduce(nparticles_temp, psip_list%tot_nparticles, psip_list%nspaces, MPI_PREAL, MPI_SUM, &
                                MPI_COMM_WORLD, ierr)
#else
            psip_list%tot_nparticles = nparticles_temp
#endif
        else
            psip_list%tot_nparticles = target_nparticles_tot
        end if

        call direct_annihilation(sys, rng, qmc_in, reference, annihilation_flags, psip_list, spawn)

        if (dmqmc_in%propagate_to_beta) then
            ! Reset the position of the first spawned particle in the spawning array
            spawn%head = spawn%head_start
            ! During the metropolis steps determinants originally in the correct
            ! portions of the spawned walker array are no longer there due to
            ! new determinants being accepted. So we need to reorganise the
            ! determinants appropriately.
            call redistribute_particles(psip_list%states, real_factor, psip_list%pops, &
                                        psip_list%nstates, psip_list%nparticles, spawn, nload_slots)
            call direct_annihilation(sys, rng, qmc_in, reference, annihilation_flags, psip_list, spawn)
        end if

    end subroutine create_initial_density_matrix

    subroutine random_distribution_heisenberg(rng, basis, spins_up, npsips, ireplica, spawn)

        ! For the Heisenberg model only. Distribute the initial number of psips
        ! along the main diagonal. Each diagonal element should be chosen
        ! with the same probability.

        ! Start from a state with all spins down, then choose the above number
        ! of spins to flip up with equal probability.

        ! In/Out:
        !    rng: random number generator.
        !    spawn: spawn_t object to hold spawned particles.
        ! In:
        !    basis: information about the single-particle basis.
        !    spins_up: for the spin configurations generated, this number
        !       specifies how many of the spins shall be up.
        !    npsips: The total number of psips to be created.
        !    ireplica: index of replica (ie which of the possible concurrent
        !       DMQMC populations are we initialising)

        use basis_types, only: basis_t
        use calc, only: ms_in
        use dSFMT_interface, only: dSFMT_t, get_rand_close_open
        use fciqmc_data, only: real_factor
        use spawn_data, only: spawn_t
        use parallel
        use system

        type(dSFMT_t), intent(inout) :: rng
        type(basis_t), intent(in) :: basis
        integer, intent(in) :: spins_up
        integer(int_64), intent(in) :: npsips
        integer, intent(in) :: ireplica
        type(spawn_t), intent(inout) :: spawn
        integer(int_64) :: i
        integer :: rand_basis, bits_set
        integer :: bit_element, bit_position
        integer(i0) :: f(basis%string_len)
        real(dp) :: rand_num

        do i = 1, npsips

            ! Start with all spins down.
            f = 0_i0
            bits_set = 0

            do
                ! If half the spins are now flipped up, we have our basis
                ! function fully created, so exit the loop.
                if (bits_set == spins_up) exit
                ! Choose a random spin to flip.
                rand_num = get_rand_close_open(rng)
                rand_basis = ceiling(rand_num*basis%nbasis)
                ! Find the corresponding positions for this spin.
                bit_position = basis%bit_lookup(1,rand_basis)
                bit_element = basis%bit_lookup(2,rand_basis)
                if (.not. btest(f(bit_element),bit_position)) then
                    ! If not flipped up, flip the spin up.
                    f(bit_element) = ibset(f(bit_element),bit_position)
                    bits_set = bits_set + 1
                end if
            end do

            ! Now call a routine to add the corresponding diagonal element to
            ! the spawned walkers list.
            call create_diagonal_density_matrix_particle(f,basis%string_len, &
                    basis%tensor_label_len, real_factor,ireplica, spawn)

        end do

    end subroutine random_distribution_heisenberg

    subroutine random_distribution_electronic(rng, sys, sym, npsips, ireplica, all_sym_sectors, spawn)

        ! For the electronic Hamiltonians only. Distribute the initial number of psips
        ! along the main diagonal. Each diagonal element should be chosen
        ! with the same probability.

        ! Determinants are generated uniformly in the Hilbert space associated
        ! to selected symmetry sector and spin polarisation.

        ! In:
        !    sys: system being studied.
        !    sym: only keep determinants in this symmetry sector.
        !    npsips: The total number of psips to be created.
        !    ireplica: index of replica (ie which of the possible concurrent
        !       DMQMC populations are we initialising)
        !    all_sym_sectors: create determinants in all symmetry sectors?
        ! In/Out:
        !    rng: random number generator
        !    spawn: spawn_t object to hold spawned particles.

        use dSFMT_interface, only: dSFMT_t, get_rand_close_open
        use symmetry, only: symmetry_orb_list
        use hilbert_space, only: gen_random_det_full_space
        use system, only: sys_t
        use fciqmc_data, only: real_factor
        use spawn_data, only: spawn_t

        type(dSFMT_t), intent(inout) :: rng
        type(sys_t), intent(in) :: sys
        integer, intent(in) :: sym
        integer(int_64), intent(in) :: npsips
        integer, intent(in) :: ireplica
        logical, intent(in) :: all_sym_sectors
        type(spawn_t), intent(inout) :: spawn

        integer(int_64) :: i
        integer(i0) :: f(sys%basis%string_len)
        integer :: occ_list(sys%nalpha+sys%nbeta)

        do i = 1, npsips
            do
                ! Generate a random determinant uniformly in this specific
                ! symmetry sector and spin polarisation.
                call gen_random_det_full_space(rng, sys, f, occ_list)
                if (all_sym_sectors .or. symmetry_orb_list(sys, occ_list) == sym) then
                    call create_diagonal_density_matrix_particle(f, sys%basis%string_len, &
                        sys%basis%tensor_label_len, real_factor, ireplica, spawn)
                    exit
                end if
            end do
        end do

    end subroutine random_distribution_electronic

    subroutine initialise_dm_metropolis(sys, rng, qmc_in, dmqmc_in, npsips, sym, ireplica, spawn)

        ! Attempt to initialise the temperature dependent trial density matrix
        ! using the metropolis algorithm. We either uniformly distribute psips
        ! on all excitation levels or use the grand canonical partition function
        ! as first guess and then use the Metropolis algorithm to distribute
        ! psips according to desired trial density matrix.
        ! The Metropolis algorithm works by generating a new determinant which
        ! is accepted / rejected based on the value of the total energy of the
        ! Slater determinant i.e. E_i = <D_i | H_T | D_i>, where H_T is the
        ! "trial" Hamiltonian.

        ! *** Warning ***: It is up to the user to decide whether enough
        !     metropolis steps have been carried out and that the trial density
        !     matrix is indeed being sampled correctly.

        ! In:
        !    sys: system being studied.
        !    qmc_in: input options relating to QMC methods.
        !    dmqmc_in: input options relating to DMQMC.
        !    beta: (inverse) temperature at which we're looking to sample the
        !        trial density matrix.
        !    sym: symmetry index of determinant space we wish to sample.
        !    npsips: number of psips to distribute in this sector.
        !    ireplica: replica index.
        ! In/Out:
        !    rng: random number generator.
        !    spawn: spawn_t object containing the initial distribution of
        !        psips on the diagonal.

        use dSFMT_interface, only: dSFMT_t, get_rand_close_open
        use system, only: sys_t
        use determinants, only: alloc_det_info_t, det_info_t, dealloc_det_info_t, decode_det_spinocc_spinunocc, &
                                encode_det
        use excitations, only: excit_t, create_excited_det
        use fciqmc_data, only: real_factor
        use parallel, only: nprocs, nthreads, parent
        use hilbert_space, only: gen_random_det_truncate_space
        use proc_pointers, only: trial_dm_ptr, gen_excit_ptr
        use qmc_data, only: qmc_in_t
        use utils, only: int_fmt
        use spawn_data, only: spawn_t
        use dmqmc_data, only: dmqmc_in_t

        type(sys_t), intent(in) :: sys
        type(qmc_in_t), intent(in) :: qmc_in
        type(dmqmc_in_t), intent(in) :: dmqmc_in
        integer, intent(in) :: sym
        integer(int_64), intent(in) :: npsips
        integer, intent(in) :: ireplica
        type(dSFMT_t), intent(inout) :: rng
        type(spawn_t), intent(inout) :: spawn

        integer :: occ_list(sys%nel), naccept
        integer :: idet, iattempt, nsuccess
        integer :: thread_id = 0, proc
        integer(i0) :: f_old(sys%basis%string_len), f_new(sys%basis%string_len)
        real(p), target :: tmp_data(1)
        real(p) :: pgen, hmatel, E_new, E_old, prob
        real(dp) :: r
        type(det_info_t) :: cdet
        type(excit_t) :: connection
        real(dp) :: move_prob(0:sys%nalpha, sys%max_number_excitations)

        naccept = 0 ! Number of metropolis moves which are accepted.
        nsuccess = 0 ! Number of successful proposal steps i.e. excluding null excitations.
        idet = 0

        call alloc_det_info_t(sys, cdet)

        ! The metropolis move is to excite [1:max_metropolis_move] electrons. This is achieved
        ! either by using the excitation generators when working in a
        ! symmetry contrained system or by uniformly exciting the electrons
        ! among the available levels. In the latter case we need to set the
        ! probabilities of a particular move (e.g. move two alpha spins),
        ! so do this here.
        if (dmqmc_in%all_sym_sectors) call set_level_probabilities(sys, move_prob, dmqmc_in%max_metropolis_move)

        ! Visit every psip metropolis_attempts times.
        do iattempt = 1, dmqmc_in%metropolis_attempts
            do proc = 0, nprocs-1
                do idet = spawn%head_start(nthreads-1,proc)+1, spawn%head(thread_id,proc)
                    cdet%f = spawn%sdata(:sys%basis%string_len,idet)
                    E_old = trial_dm_ptr(sys, cdet%f)
                    tmp_data(1) = E_old
                    cdet%data => tmp_data
                    call decode_det_spinocc_spinunocc(sys, cdet%f, cdet)
                    if (dmqmc_in%all_sym_sectors) then
                        call gen_random_det_truncate_space(rng, sys, dmqmc_in%max_metropolis_move, cdet, move_prob, occ_list)
                        nsuccess = nsuccess + 1
                        call encode_det(sys%basis, occ_list, f_new)
                    else
                        call gen_excit_ptr%full(rng, sys, qmc_in, cdet, pgen, connection, hmatel)
                        ! Check that we didn't generate a null excitation.
                        ! [todo] - Modify accordingly if pgen is ever calculated in for the ueg.
                        if (hmatel == 0) cycle
                        nsuccess = nsuccess + 1
                        call create_excited_det(sys%basis, cdet%f, connection, f_new)
                    end if
                    ! Accept new det with probability p = min[1,exp(-\beta(E_new-E_old))]
                    E_new = trial_dm_ptr(sys, f_new)
                    prob = exp(-1.0_p*dmqmc_in%init_beta*(E_new-E_old))
                    r = get_rand_close_open(rng)
                    if (prob > r) then
                        ! Accept the new determinant by modifying the entry
                        ! in spawned walker list.
                        naccept = naccept + 1
                        spawn%sdata(:sys%basis%string_len,idet) = f_new
                        spawn%sdata(sys%basis%string_len+1:sys%basis%tensor_label_len,idet) = f_new
                    end if
                end do
            end do
        end do

        if (parent) write (6,'(1X,"#",1X, "Average acceptance ratio: ",f8.7,1X," Average number of null excitations: ", f8.7)') &
                           real(naccept)/nsuccess, real(dmqmc_in%metropolis_attempts*npsips-nsuccess)/&
                                                   &(dmqmc_in%metropolis_attempts*npsips)

        call dealloc_det_info_t(cdet)

    end subroutine initialise_dm_metropolis

    subroutine init_uniform_ensemble(sys, npsips, sym, f0, ireplica, all_sym_sectors, rng, spawn)

        ! Create an initital distribution of psips along the diagonal.
        ! This subroutine will return a list of occupied determinants in
        ! spawn which can then be used for the metropolis initialisation.
        ! Psips are distributed equally among all excitation levels.

        ! In:
        !    sys: system being studied.
        !    npsips: number of psips to distribute on the diagonal.
        !    sym: symmetry of determinants being occupied.
        !    f0: bit string of reference determinant
        !    ireplica: replica index.
        !    all_sym_sectors: create determinants in all symmetry sectors?
        ! In/Out:
        !    rng: random number generator.
        !    spawn: spawn_t object containing list of occupied determinants.

        use spawn_data, only: spawn_t
        use determinants, only: encode_det, decode_det_spinocc_spinunocc, dealloc_det_info_t, &
                                det_info_t, alloc_det_info_t
        use fciqmc_data, only: real_factor
        use hilbert_space, only: gen_random_det_truncate_space
        use symmetry, only: symmetry_orb_list
        use system, only: sys_t
        use dSFMT_interface, only: dSFMT_t, get_rand_close_open

        type(sys_t), intent(in) :: sys
        integer(int_64), intent(in) :: npsips
        integer, intent(in) :: sym
        integer(i0), intent(in) :: f0(sys%basis%string_len)
        integer, intent(in) :: ireplica
        logical, intent(in) :: all_sym_sectors
        type(dSFMT_t), intent(inout) :: rng
        type(spawn_t), intent(inout) :: spawn

        integer :: occ_list(sys%nel), idet, ilevel
        integer(i0) :: f(sys%basis%string_len)
        type(det_info_t) :: det0
        integer(int_64) :: psips_per_level
        real(dp) :: ptrunc_level(0:sys%nalpha, sys%max_number_excitations)

        ! [todo] - Include all psips.
        psips_per_level = int(npsips/sys%max_number_excitations)
        call alloc_det_info_t(sys, det0)
        call decode_det_spinocc_spinunocc(sys, f0, det0)
        det0%f = f0
        ! gen_random_det_truncate_space does not produce determinants at
        ! excitation level zero, so take care of this explicitly.
        do idet = 1, psips_per_level
            call create_diagonal_density_matrix_particle(f0, sys%basis%string_len, &
                                                         sys%basis%tensor_label_len, real_factor, ireplica, spawn)
        end do

        ! Uniformly distribute the psips on all other levels.
        call set_level_probabilities(sys, ptrunc_level, sys%max_number_excitations)

        do idet = 1, npsips-psips_per_level
            ! Repeatedly attempt to create determinants on every other
            ! excitation level.
            do
                call gen_random_det_truncate_space(rng, sys, sys%max_number_excitations, det0, ptrunc_level(0:,:), occ_list)
                if (all_sym_sectors .or. symmetry_orb_list(sys, occ_list) == sym) then
                    call encode_det(sys%basis, occ_list, f)
                    call create_diagonal_density_matrix_particle(f, sys%basis%string_len, &
                                                                sys%basis%tensor_label_len, real_factor, ireplica, spawn)
                    exit
                end if
            end do
        end do

        call dealloc_det_info_t(det0)

    end subroutine init_uniform_ensemble

    subroutine init_grand_canonical_ensemble(sys, dmqmc_in, sym, npsips, spawn, rng)

        ! Initially distribute psips according to the grand canonical
        ! distribution function.

        ! In:
        !    sys: system being studied.
        !    dmqmc_in: input options for dmqmc.
        !    sym: symmetry sector under consideration.
        !    npsips: number of psips to create on the diagonal.
        ! In/Out:
        !    spawn: spawned list.
        !    rng: random number generator.

        use system, only: sys_t
        use spawn_data, only: spawn_t
        use fciqmc_data, only: real_factor
        use symmetry, only: symmetry_orb_list
        use dSFMT_interface, only: dSFMT_t, get_rand_close_open
        use determinants, only: encode_det
        use canonical_kinetic_energy, only: generate_allowed_orbital_list
        use dmqmc_data, only: dmqmc_in_t

        type(sys_t), intent(in) :: sys
        integer, intent(in) :: sym
        type(dmqmc_in_t), intent(in) :: dmqmc_in
        integer(int_64), intent(in) :: npsips
        type(spawn_t), intent(inout) :: spawn
        type(dSFMT_t), intent(inout) :: rng

        real(dp) :: p_single(sys%basis%nbasis/2)
        real(dp) :: r
        integer :: occ_list(sys%nel)
        integer(i0) :: f(sys%basis%string_len)
        integer :: ireplica, iorb, ipsip
        logical :: gen

        ireplica = 1

        ! Calculate orbital occupancies.
        ! * Warning *: We assume that we are dealing with a system without
        ! magnetic fields or other funny stuff, so the probabilty of occupying
        ! an alpha spin orbital is equal to that of occupying a beta spin
        ! orbital.
        forall(iorb=1:sys%basis%nbasis:2) p_single(iorb/2+1) = 1.0_p / &
                                          (1+exp(dmqmc_in%init_beta*(sys%basis%basis_fns(iorb)%sp_eigv-sys%ueg%chem_pot)))

        ! In the grand canoical ensemble the probability of occupying a
        ! determinant, |D_i>, is given by \prod_i^N p_i, where the p_i's are the
        ! Fermi factors calculated above in p_single(i). We normally work in the
        ! canonical ensemble, however, so we need to discard any determinant
        ! generated which does not contain nel electrons to obtain the correct
        ! normalisation and hence, the correct distribution. For small systems
        ! the number fluctuations are usually small, so this routine is quite
        ! fast.
        ipsip = 0
        do while (ipsip < npsips)
            occ_list = 0
            ! Select the alpha and beta spin orbitals and discard any
            ! determinant without the correct number of particles.
            if (sys%nalpha > 0) call generate_allowed_orbital_list(sys, rng, p_single, sys%nalpha, 1, occ_list(:sys%nalpha), gen)
            if (.not. gen) cycle
            if (sys%nbeta > 0) call generate_allowed_orbital_list(sys, rng, p_single, sys%nbeta, 0, occ_list(sys%nalpha+1:), gen)
            if (.not. gen) cycle
            ! Create the determinant.
            if (dmqmc_in%all_sym_sectors .or. symmetry_orb_list(sys, occ_list) == sym) then
                call encode_det(sys%basis, occ_list, f)
                call create_diagonal_density_matrix_particle(f, sys%basis%string_len, &
                                                            sys%basis%tensor_label_len, real_factor, ireplica, spawn)
                ipsip = ipsip + 1
            end if
        end do

    end subroutine init_grand_canonical_ensemble

    subroutine set_level_probabilities(sys, ptrunc_level, max_excit)

        ! Set the probabilities for creating a determinant on a given
        ! excitation level so that get_random_det_truncate_space can be used.
        !
        ! In:
        !    sys: system being studied.
        !    max_excit: maximum excitation of a reference determinant we wish to
        !        consider.
        ! In/Out:
        !    ptrunc_level: array containing probabilities for spawning at a
        !        particular (ialpha, ilevel) excitation level.

        use system, only: sys_t
        use utils, only: binom_r

        type(sys_t), intent(in) :: sys
        integer, intent(in) :: max_excit
        real(dp), intent(inout) :: ptrunc_level(0:sys%nalpha, sys%max_number_excitations)

        integer :: ialpha, ilevel
        real(dp) :: offset

        ! Set up the probabilities for creating a determinant on a given
        ! excitation level.
        ! ptunc_level(ialpha, ilevel) gives the probability of creating a
        ! determinant with excitation level by exciting ialpha
        ! alpha spins (and (ilevel-ialpha) beta spins). We give each possible
        ! realisation of an excitation equal probability and allow each
        ! excitation to occur with equal probability.
        ptrunc_level = 0.0_dp
        do ilevel = 1, max_excit
            do ialpha = max(0,ilevel-sys%nbeta), min(ilevel,sys%nalpha)
                ! Number of ways of exciting ialpha electrons such that
                ! ialpha + nbeta = ilevel.
                ptrunc_level(ialpha, ilevel) = binom_r(ilevel,ialpha)
            end do
        end do

        ! Normalisation for the distribution at this excitation level.
        forall(ilevel=1:max_excit) ptrunc_level(:,ilevel) = ptrunc_level(:,ilevel) / (max_excit*sum(ptrunc_level(:,ilevel)))

        offset = 0.0_dp
        do ilevel = 1, max_excit
            do ialpha = max(0,ilevel-sys%nbeta), min(ilevel,sys%nalpha)
                ptrunc_level(ialpha, ilevel) = offset + ptrunc_level(ialpha, ilevel)
                offset = ptrunc_level(ialpha, ilevel)
            end do
        end do

    end subroutine set_level_probabilities

    subroutine create_diagonal_density_matrix_particle(f, string_len, tensor_label_len, nspawn, particle_type, spawn)

        ! Create a psip on a diagonal element of the density matrix by adding
        ! it to the spawned walkers list. This list can then be sorted correctly
        ! by the direct_annihilation routine.

        ! In:
        !    f: bitstring representation of index of the diagonal element upon
        !        which a new psip shall be placed.
        !    string_len: length of bit array storing a many-particle basis function
        !        (e.g. a determinant or spin product).
        !    tensor_label_len: length of bit array storing the label of a density
        !        matrix element (usually 2xstring_len).
        !    nspawn: the number of particles to be added to this diagonal
        !        element.
        !    particle_type: the label of the replica to which this particle is
        !        to sample.
        ! In/Out:
        !    spawn: spawn_t object to which the spawned particle is added.

        use hashing
        use parallel
        use errors, only: stop_all
        use spawn_data, only: spawn_t
        use spawning, only: add_spawned_particle

        integer, intent(in) :: string_len, tensor_label_len
        integer(i0), intent(in) :: f(string_len)
        integer(int_p), intent(in) :: nspawn
        integer ::particle_type
        integer(i0) :: f_new(tensor_label_len)
        type(spawn_t), intent(inout) :: spawn
#ifndef PARALLEL
        integer, parameter :: iproc_spawn = 0
#else
        integer :: iproc_spawn
#endif

        ! Create the bitstring of the determinant.
        f_new = 0_i0
        f_new(:string_len) = f
        f_new((string_len+1):(tensor_label_len)) = f

#ifdef PARALLEL
        ! Need to determine which processor the spawned psip should be sent to.
        iproc_spawn = modulo(murmurhash_bit_string(f_new, &
                                tensor_label_len, spawn%hash_seed), nprocs)
#endif

        ! spawn%head_start(nthreads-1,i) stores the last entry before the
        ! start of the block of spawned particles to be sent to processor i.
        if (spawn%head(0,iproc_spawn)+1 - spawn%head_start(nthreads-1,iproc_spawn) >= spawn%block_size) &
            call stop_all('create_diagonal_density_matrix_particle', 'There is no space left in the spawning array.')

        call add_spawned_particle(f_new, nspawn, particle_type, iproc_spawn, spawn)

    end subroutine create_diagonal_density_matrix_particle

    subroutine decode_dm_bitstring(basis, f, isym, rdm_info)

        ! This function maps a full DMQMC bitstring to two bitstrings encoding
        ! the subsystem-A RDM bitstrings. These resulting bitstrings are stored
        ! in the end1 and end2 components of rdm_info.

        ! Crucially, the mapping is performed so that, if there are two
        ! subsystems which are equivalent by symmetry, then equivalent sites in
        ! those two subsystems will be mapped to the same RDM bitstrings. This
        ! is clearly a requirement to obtain a correct representation of the RDM
        ! when using translational symmetry.

        ! In:
        !    basis: information about the single-particle basis.
        !    f: bitstring representation of the subsystem-A state.
        !    isym: The label of the symmetry vector being considered.
        ! In/Out:
        !    rdm_info: information about the RDM and subsystem for the RDM being
        !        considered.

        use basis_types, only: basis_t
        use dmqmc_data, only: rdm_t

        type(basis_t), intent(in) :: basis
        integer(i0), intent(in) :: f(:)
        integer, intent(in) :: isym
        type(rdm_t), intent(inout) :: rdm_info

        integer :: i, bit_pos, bit_element

        ! Start from all bits down, so that we can flip bits up one by one.
        rdm_info%end1 = 0_i0
        rdm_info%end2 = 0_i0

        ! Loop over all the sites in the subsystem considered for the reduced
        ! density matrix.
        do i = 1, rdm_info%A_nsites
            ! Find the final bit positions and elements.
            bit_pos = basis%bit_lookup(1,i)
            bit_element = basis%bit_lookup(2,i)

            ! If the spin is up, set the corresponding bit in the first
            ! bitstring.
            if (btest(f(rdm_info%bit_pos(i,isym,2)),rdm_info%bit_pos(i,isym,1))) &
                rdm_info%end1(bit_element) = ibset(rdm_info%end1(bit_element),bit_pos)
            ! Similarly for the second index, by looking at the second end of
            ! the bitstring.
            if (btest(f(rdm_info%bit_pos(i,isym,2)+basis%string_len),rdm_info%bit_pos(i,isym,1))) &
                rdm_info%end2(bit_element) = ibset(rdm_info%end2(bit_element),bit_pos)
        end do

    end subroutine decode_dm_bitstring

    subroutine update_sampling_weights(rng, basis, qmc_in, psip_list, weighted_sampling)

        ! This routine updates the values of the weights used in importance
        ! sampling. It also removes or adds psips from the various excitation
        ! levels accordingly.

        ! In/Out:
        !    rng: random number generator.
        !    psip_list: main particle list.
        ! In:
        !    basis: information about the single-particle basis.
        !    qmc_in: Input options relating to QMC methods.
        ! In/Out:
        !    weighted_sampling: type containing weighted sampling information.

        use annihilation, only: remove_unoccupied_dets
        use basis_types, only: basis_t
        use dSFMT_interface, only: dSFMT_t, get_rand_close_open
        use excitations, only: get_excitation_level
        use fciqmc_data, only: real_factor
        use qmc_data, only: qmc_in_t, particle_t
        use dmqmc_data, only: dmqmc_weighted_sampling_t

        type(dSFMT_t), intent(inout) :: rng
        type(basis_t), intent(in) :: basis
        type(qmc_in_t), intent(in) :: qmc_in
        type(particle_t), intent(inout) :: psip_list
        type(dmqmc_weighted_sampling_t), intent(inout) :: weighted_sampling

        integer :: idet, ireplica, excit_level, nspawn, sign_factor
        real(p) :: new_population_target(psip_list%nspaces)
        integer(int_p) :: old_population(psip_list%nspaces), new_population(psip_list%nspaces)
        real(dp) :: r, pextra

        ! Alter weights for the next iteration.
        weighted_sampling%probs = real(weighted_sampling%probs,p)*weighted_sampling%altering_factors

        ! When the weights for an excitation level are increased by a factor,
        ! the number of psips on that level has to decrease by the same factor,
        ! else the wavefunction which the psips represent will not be the
        ! correct importance sampled wavefunction for the new weights. The code
        ! below loops over every psips and destroys (or creates) it with the
        ! appropriate probability.
        do idet = 1, psip_list%nstates

            excit_level = get_excitation_level(psip_list%states(1:basis%string_len,idet), &
                    psip_list%states(basis%string_len+1:basis%tensor_label_len,idet))

            old_population = abs(psip_list%pops(:,idet))

            ! The new population that we are aiming for. If this is not an
            ! integer then we will have to round up or down to an integer with
            ! an unbiased probability.
            new_population_target = abs(real(psip_list%pops(:,idet),p))/weighted_sampling%altering_factors(excit_level)
            new_population = int(new_population_target, int_p)

            ! If new_population_target is not an integer, round it up or down
            ! with an unbiased probability. Do this for each replica.
            do ireplica = 1, psip_list%nspaces

                pextra = new_population_target(ireplica) - new_population(ireplica)

                if (pextra > depsilon) then
                    r = get_rand_close_open(rng)
                    if (r < pextra) new_population(ireplica) = new_population(ireplica) + 1_int_p
                end if

                ! Finally, update the walker population.
                associate(pops=>psip_list%pops)
                    pops(ireplica,idet) = sign(new_population(ireplica), pops(ireplica,idet))
                end associate

            end do

            ! Update the total number of walkers.
            psip_list%nparticles = psip_list%nparticles + real(new_population - old_population, p)/real_factor

        end do

        ! Call the annihilation routine to update the main walker list, as some
        ! sites will have become unoccupied and so need removing from the
        ! simulation.
        call remove_unoccupied_dets(rng, psip_list, qmc_in%real_amplitudes)

    end subroutine update_sampling_weights

    subroutine output_and_alter_weights(dmqmc_in, max_number_excitations, weighted_sampling)

        ! This routine will alter and output the sampling weights used in
        ! importance sampling. It uses the excitation distribution, calculated
        ! on the beta loop which has just finished, and finds the weights needed
        ! so that each excitation level will have roughly equal numbers of psips
        ! in the next loop. For example, to find the weights of psips on the 1st
        ! excitation level, divide the number of psips on the 1st excitation
        ! level by the number on the 0th level, then multiply the old sampling
        ! weight by this number to give the new weight. This can be used when
        ! the weights are being introduced gradually each beta loop, too. The
        ! weights are output and can then be used in future DMQMC runs.

        ! In:
        !    dmqmc_in: input options relating to DMQMC.
        !    max_number_excitations: maximum number of excitations possible (see
        !       sys_t type in system for details).
        !    vary_weights: vary weights with beta?
        ! In/Out:
        !    weighted_sampling: type containing weighted sampling information.

        use fciqmc_data, only: excit_dist
        use dmqmc_data, only: dmqmc_in_t, dmqmc_weighted_sampling_t
        use parallel

        integer, intent(in) :: max_number_excitations
        type(dmqmc_in_t), intent(inout) :: dmqmc_in
        type(dmqmc_weighted_sampling_t), intent(inout) :: weighted_sampling

        integer :: i, ierr
#ifdef PARALLEL
        real(p) :: merged_excit_dist(0:max_number_excitations)
        call mpi_allreduce(excit_dist, merged_excit_dist, max_number_excitations+1, &
            MPI_PREAL, MPI_SUM, MPI_COMM_WORLD, ierr)

        excit_dist = merged_excit_dist
#endif

        ! It is assumed that there is an even maximum number of excitations.
        do i = 1, (max_number_excitations/2)
            ! Don't include levels where there are very few psips accumulated.
            if (excit_dist(i-1) > 10.0_p .and. excit_dist(i) > 10.0_p) then
                ! Alter the sampling weights using the relevant excitation
                ! distribution.
                dmqmc_in%sampling_probs(i) = dmqmc_in%sampling_probs(i)*&
                    (excit_dist(i)/excit_dist(i-1))
                dmqmc_in%sampling_probs(max_number_excitations+1-i) = dmqmc_in%sampling_probs(i)**(-1)
            end if
        end do

        ! Recalculate weighted_sampling%probs with the new weights.
        do i = 1, max_number_excitations
            weighted_sampling%probs(i) = weighted_sampling%probs(i-1)*dmqmc_in%sampling_probs(i)
        end do

        ! If vary_weights is true then the weights are to be introduced
        ! gradually at the start of each beta loop. This requires redefining
        ! weighted_sampling%altering_factors to coincide with the new sampling weights.
        if (dmqmc_in%vary_weights) then
            weighted_sampling%altering_factors = real(weighted_sampling%probs,dp)**(1/real(dmqmc_in%finish_varying_weights,dp))
            ! Reset the weights for the next loop.
            weighted_sampling%probs = 1.0_p
        end if

        if (parent) then
            ! Print out weights in a form which can be copied into an input
            ! file.
            write(6, '(a31,2X)', advance = 'no') ' # Importance sampling weights:'
            do i = 1, max_number_excitations
                write (6, '(es12.4,2X)', advance = 'no') dmqmc_in%sampling_probs(i)
            end do
            write (6, '()', advance = 'yes')
        end if

    end subroutine output_and_alter_weights

end module dmqmc_procedures

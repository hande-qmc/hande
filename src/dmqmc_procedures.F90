module dmqmc_procedures

use const
implicit none

contains

    subroutine init_dmqmc(sys, qmc_in, dmqmc_in, nreplicas, qs, dmqmc_estimates, weighted_sampling)

        ! In:
        !    nreplicas: number of replicas being used.
        !    sys: system being studied. This should be left in an unmodified
        !       state on output.
        !    qmc_in: Input options relating to QMC methods.
        !    dmqmc_in: Input options relating to DMQMC.
        ! In/Out:
        !    qs: estimators not specific to DMQMC.
        !    dmqmc_estimates: type containing estimates for observables.
        !    weighted_sampling: type containing information for weighted
        !        sampling.

        use calc, only: doing_dmqmc_calc, dmqmc_correlation
        use checking, only: check_allocate
        use system, only: sys_t

        use qmc_data, only: qmc_in_t, qmc_state_t
        use dmqmc_data, only: dmqmc_in_t, dmqmc_estimates_t, dmqmc_weighted_sampling_t

        type(sys_t), intent(in) :: sys
        integer, intent(in) :: nreplicas
        type(qmc_in_t), intent(in) :: qmc_in
        type(dmqmc_in_t), intent(in) :: dmqmc_in
        type(qmc_state_t), intent(inout) :: qs
        type(dmqmc_estimates_t), intent(inout) :: dmqmc_estimates
        type(dmqmc_weighted_sampling_t), intent(inout) :: weighted_sampling

        integer :: ierr, i, bit_position, bit_element

        allocate(dmqmc_estimates%trace(nreplicas), stat=ierr)
        call check_allocate('dmqmc_estimates%trace',size(dmqmc_estimates%trace),ierr)
        dmqmc_estimates%trace = 0.0_p

        dmqmc_estimates%inst_rdm%nrdms = dmqmc_in%rdm%nrdms
        allocate(dmqmc_estimates%inst_rdm%traces(nreplicas, dmqmc_in%rdm%nrdms), stat=ierr)
        call check_allocate('dmqmc_estimates%inst_rdm%traces', size(dmqmc_estimates%inst_rdm%traces),ierr)
        dmqmc_estimates%inst_rdm%traces = 0.0_p

        ! If calculating a correlation function then set up the necessary bit
        ! mask. This has a bit set for each of the two sites/orbitals being
        ! considered in the correlation function.
        if (doing_dmqmc_calc(dmqmc_correlation)) then
            allocate(dmqmc_estimates%correlation_mask(1:sys%basis%string_len), stat=ierr)
            call check_allocate('dmqmc_estimates%correlation_mask',sys%basis%string_len,ierr)
            dmqmc_estimates%correlation_mask = 0_i0
            do i = 1, 2
                bit_position = sys%basis%bit_lookup(1,dmqmc_in%correlation_sites(i))
                bit_element = sys%basis%bit_lookup(2,dmqmc_in%correlation_sites(i))
                dmqmc_estimates%correlation_mask(bit_element) = ibset(dmqmc_estimates%correlation_mask(bit_element), bit_position)
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
        if (dmqmc_in%symmetric) then
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
            qs%init_beta = dmqmc_in%init_beta / sys%ueg%ef
        else
            qs%init_beta = dmqmc_in%init_beta
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
            allocate(weighted_sampling%sampling_probs(sys%max_number_excitations), stat=ierr)
            if (allocated(dmqmc_in%sampling_probs)) then
                weighted_sampling%sampling_probs(:size(dmqmc_in%sampling_probs)) = dmqmc_in%sampling_probs
            end if
            call check_allocate('weighted_sampling%sampling_probs',sys%max_number_excitations,ierr)
            if (dmqmc_in%propagate_to_beta .and. dmqmc_in%symmetric) then
                ! If using symmetric version of ipdmqmc let the last entry contain 0.5*(beta-tau).
                allocate(weighted_sampling%probs(0:sys%max_number_excitations+1), stat=ierr)
                call check_allocate('weighted_sampling%probs',sys%max_number_excitations+2,ierr)
                allocate(weighted_sampling%probs_old(0:sys%max_number_excitations+1), stat=ierr)
                call check_allocate('weighted_sampling%probs_old',sys%max_number_excitations+2,ierr)
            else
                allocate(weighted_sampling%probs(0:sys%max_number_excitations), stat=ierr)
                call check_allocate('weighted_sampling%probs', sys%max_number_excitations+1, ierr)
                allocate(weighted_sampling%probs_old(0:sys%max_number_excitations), stat=ierr)
                call check_allocate('weighted_sampling%probs_old',sys%max_number_excitations+1,ierr)
            end if
            weighted_sampling%probs(0) = 1.0_p
            weighted_sampling%probs_old = 1.0_p
            do i = 1, sys%max_number_excitations
                weighted_sampling%probs(i) = weighted_sampling%probs(i-1)*weighted_sampling%sampling_probs(i)
            end do
            if (dmqmc_in%vary_weights) then
                ! Allocate an array to store the factors by which the weights
                ! will change each iteration.
                allocate(weighted_sampling%altering_factors(0:sys%max_number_excitations), stat=ierr)
                call check_allocate('weighted_sampling%altering_factors',sys%max_number_excitations+1,ierr)
                weighted_sampling%altering_factors = &
                    & real(weighted_sampling%probs(:sys%max_number_excitations),dp)**(1/real(dmqmc_in%finish_varying_weights,dp))
                ! If varying the weights, start the accumulated probabilties
                ! as all 1.0 initially, and then alter them gradually later.
                weighted_sampling%probs = 1.0_p
            end if
        else
            ! If not using the importance sampling procedure, turn it off by
            ! setting all amplitudes to 1.0 in the relevant arrays.
            if (dmqmc_in%propagate_to_beta .and. dmqmc_in%symmetric) then
                ! If using symmetric version of ipdmqmc let the last entry contain 0.5*(beta-tau).
                allocate(weighted_sampling%probs(0:sys%max_number_excitations+1), stat=ierr)
                call check_allocate('weighted_sampling%probs',sys%max_number_excitations+2,ierr)
                allocate(weighted_sampling%probs_old(0:sys%max_number_excitations), stat=ierr)
                call check_allocate('weighted_sampling%probs_old',sys%max_number_excitations+1,ierr)
            else
                allocate(weighted_sampling%probs(0:sys%max_number_excitations), stat=ierr)
                call check_allocate('weighted_sampling%probs',sys%max_number_excitations+1,ierr)
                allocate(weighted_sampling%probs_old(0:sys%max_number_excitations), stat=ierr)
                call check_allocate('weighted_sampling%probs_old',sys%max_number_excitations+1,ierr)
            end if
            weighted_sampling%probs = 1.0_p
            weighted_sampling%probs_old = 1.0_p
        end if

        ! If doing a reduced density matrix calculation, allocate and define
        ! the bit masks that have 1's at the positions referring to either
        ! subsystems A or B.
        if (dmqmc_in%rdm%doing_rdm) call setup_rdm_arrays(sys, .true., dmqmc_estimates%subsys_info, &
                                                          dmqmc_estimates%ground_rdm%rdm, qmc_in, dmqmc_in%rdm, &
                                                          dmqmc_estimates%inst_rdm, nreplicas, qs%psip_list%pop_real_factor)

    end subroutine init_dmqmc

    subroutine setup_rdm_arrays(sys, called_from_dmqmc, subsys_info, ground_rdm, qmc_in, rdm_in, inst_rdms, nreplicas, real_factor)

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
        ! In (optional):
        !    qmc_in: Input options relating to QMC methods.  Only needed for
        !        spawn_cutoff and if calc_inst_rdm is true.
        !    rdm_in: Input options relating to reduced density matrices.
        !    nreplicas: number of replicas being used.  Must be specified if
        !        qmc_in is.
        !    real_factor: The factor by which populations are multiplied to
        !        enable non-integer populations.  Must be specified if
        !        calc_inst_rdm is true.
        ! Out:
        !     ground_rdm: The array used to store the RDM in ground-state
        !         calculations.
        ! In/Out:
        !     subsys_info: information relating to RDM subsystems being studied.
        !     inst_rdms (optional): estimates of instantaneous
        !         (temperature-dependent) reduced density matrices.

        use calc, only: doing_dmqmc_calc, dmqmc_rdm_r2, init_proc_map_t
        use checking, only: check_allocate
        use dmqmc_data, only: subsys_t, dmqmc_inst_rdms_t
        use errors
        use hash_table, only: alloc_hash_table
        use parallel, only: parent
        use spawn_data, only: alloc_spawn_t, proc_map_t
        use system, only: sys_t, heisenberg
        use utils, only: int_fmt

        use dmqmc_data, only: dmqmc_rdm_in_t
        use qmc_data, only: qmc_in_t

        type(sys_t), intent(in) :: sys
        logical, intent(in) :: called_from_dmqmc
        type(subsys_t), intent(inout) :: subsys_info(:)
        real(p), allocatable, intent(out) :: ground_rdm(:,:)
        type(qmc_in_t), intent(in), optional :: qmc_in
        type(dmqmc_rdm_in_t), intent(in), optional :: rdm_in
        type(dmqmc_inst_rdms_t), intent(inout), optional :: inst_rdms
        integer, intent(in), optional :: nreplicas
        integer(int_p), intent(in), optional :: real_factor

        integer :: i, ierr, nrdms
        integer :: size_spawned_rdm, total_size_spawned_rdm
        integer :: nbytes_int, spawn_length_loc
        logical :: calc_ground_rdm
        type(proc_map_t) :: pm_dummy

        ! If this routine was not called from DMQMC then we must be doing a
        ! ground state RDM calculation.
        calc_ground_rdm = .not. called_from_dmqmc
        ! This should only be present if called_from_dmqmc is true.
        if (present(rdm_in)) calc_ground_rdm = rdm_in%calc_ground_rdm

        nrdms = size(subsys_info)

        ! For the Heisenberg model only currently.
        if (sys%system == heisenberg) then
            call find_rdm_masks(sys, subsys_info)
        else
            call stop_all("setup_rdm_arrays","The use of RDMs is currently only implemented for &
                           &the Heisenberg model.")
        end if

        ! Loop over all subsystems for which we are calculating RDMs.
        ! Setup the rdms array.
        do i = 1, nrdms
            ! Initialise the instance of the rdm type for this subsystem.
            subsys_info(i)%string_len = ceiling(real(subsys_info(i)%A_nsites)/i0_length)

            ! With the calc_ground_rdm option, the entire RDM is allocated. If
            ! the following condition is met then the number of rows is greater
            ! than the maximum integer accessible. This would clearly be too
            ! large, so abort in this case.
            if (calc_ground_rdm .and. subsys_info(i)%string_len > 1) call stop_all("setup_rdm_arrays",&
                "A requested RDM is too large for all indices to be addressed by a single integer.")
        end do

        ! For an ms = 0 subspace, assuming less than or exactly half the spins
        ! in the subsystem are in the subsystem, then any combination of spins
        ! can occur in the subsystem, from all spins down to all spins up. Hence
        ! the total size of the reduced density matrix will be 2**(number of
        ! spins in subsystem A).
        if (calc_ground_rdm) then
            if (sys%Ms == 0 .and. subsys_info(1)%A_nsites <= floor(real(sys%lattice%nsites,p)/2.0_p)) then
                allocate(ground_rdm(2**subsys_info(1)%A_nsites,2**subsys_info(1)%A_nsites), stat=ierr)
                call check_allocate('ground_rdm', 2**(2*subsys_info(1)%A_nsites),ierr)
                ground_rdm = 0.0_p
            else
                if (sys%Ms /= 0) then
                    call stop_all("setup_rdm_arrays","Reduced density matrices can only be used for Ms=0 &
                                   &calculations.")
                else if (subsys_info(1)%A_nsites > floor(real(sys%lattice%nsites,p)/2.0_p)) then
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
                allocate(inst_rdms%spawn(nrdms), stat=ierr)
                call check_allocate('inst_rdms%spawn', nrdms, ierr)

                ! If calculating Renyi entropy (S2).
                if (doing_dmqmc_calc(dmqmc_rdm_r2)) then
                    allocate(inst_rdms%renyi_2(nrdms), stat=ierr)
                    call check_allocate('inst_rdms%renyi_2', nrdms, ierr)
                    inst_rdms%renyi_2 = 0.0_p
                end if

                total_size_spawned_rdm = 0
                nbytes_int = bit_size(i)/8

                ! Hard-code default load balancing settings for rdm particles.
                call init_proc_map_t(1, pm_dummy)
                do i = 1, nrdms
                    ! Allocate the spawn_t and hash table instances for this RDM.
                    if (.not.present(qmc_in)) call stop_all('setup_rdm_arrays', 'qmc_in not supplied.')
                    size_spawned_rdm = (subsys_info(i)%string_len*2+nreplicas)*int_s_length/8
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

                    ! Also we use spawn_t only as a lookup/storage/compression device for the RDMs so
                    ! need not worry about hashing identical amounts of data irrespective of DET_SIZE.
                    if (.not.present(real_factor)) call stop_all('setup_rdm_arrays', 'real_factor not supplied.')
                    call alloc_spawn_t(subsys_info(i)%string_len*2, subsys_info(i)%string_len*2*i0_length, nreplicas, .false., &
                                       spawn_length_loc, qmc_in%spawn_cutoff, real_factor, pm_dummy, 27, &
                                       qmc_in%use_mpi_barriers, inst_rdms%spawn(i)%spawn)
                    ! Hard code hash table collision limit for now.  The length of
                    ! the table is three times as large as the spawning arrays and
                    ! each hash value can have 7 clashes. This was found to give
                    ! reasonable performance.
                    call alloc_hash_table(3*spawn_length_loc, 7, subsys_info(i)%string_len*2, &
                                          0, 0, 17, inst_rdms%spawn(i)%ht, inst_rdms%spawn(i)%spawn%sdata)
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

    subroutine find_rdm_masks(sys, subsys_info)

        ! Initialise bit masks for converting a density matrix basis function
        ! into its corresponding reduced density matrix basis function. Bit
        ! masks will be calculated for the subsystems requested by the user, and
        ! also for all subsystems which are equivalent by translational symmetry.

        ! In:
        !    sys: system being studied.

        use checking, only: check_allocate, check_deallocate
        use dmqmc_data, only: subsys_t
        use errors
        use real_lattice, only: find_translational_symmetry_vecs, map_vec_to_cell, enumerate_lattice_vectors
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        type(subsys_t), intent(inout) :: subsys_info(:)

        integer :: i, j, k, l, nrdms, nsym_vecs, ipos, ierr
        integer :: basis_find, bit_position, bit_element
        integer(i0) :: A_mask(sys%basis%string_len)
        real(p), allocatable :: sym_vecs(:,:)
        integer :: r(sys%lattice%ndim)
        integer, allocatable :: lvecs(:,:)

        nrdms = size(subsys_info)

        ! Return all translational symmetry vectors in sym_vecs.
        call find_translational_symmetry_vecs(sys, sym_vecs, nsym_vecs)

        ! Allocate the RDM arrays.
        do i = 1, nrdms
            if (subsys_info(i)%A_nsites == sys%lattice%nsites) then
                call stop_all('find_rdm_masks','You are attempting to use the full density matrix &
                              &as an RDM. This is not supported. You should use the &
                              &dmqmc_full_renyi_2 option to calculate the Renyi entropy of the &
                              &whole lattice.')
            else
                allocate(subsys_info(i)%B_masks(sys%basis%string_len,nsym_vecs), stat=ierr)
                call check_allocate('subsys_info(i)%B_masks', nsym_vecs*sys%basis%string_len,ierr)
                allocate(subsys_info(i)%bit_pos(subsys_info(i)%A_nsites,nsym_vecs,2), stat=ierr)
                call check_allocate('subsys_info(i)%bit_pos', nsym_vecs*subsys_info(i)%A_nsites*2,ierr)
            end if
            subsys_info(i)%B_masks = 0_i0
            subsys_info(i)%bit_pos = 0
        end do

        ! Run through every site on every subsystem and add every translational
        ! symmetry vector.
        allocate(lvecs(sys%lattice%ndim,3**sys%lattice%ndim), stat=ierr)
        call check_allocate('lvecs', size(lvecs), ierr)
        call enumerate_lattice_vectors(sys%lattice, lvecs)
        do i = 1, nrdms ! Over every subsystem.
            do j = 1, nsym_vecs ! Over every symmetry vector.
                A_mask = 0_i0
                do k = 1, subsys_info(i)%A_nsites ! Over every site in the subsystem.
                    r = sys%basis%basis_fns(subsys_info(i)%subsystem_A(k))%l
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
                            subsys_info(i)%bit_pos(k,j,1) = bit_position
                            subsys_info(i)%bit_pos(k,j,2) = bit_element
                        end if
                    end do
                end do

                ! We cannot just flip the mask for system A to get that for
                ! system B, because the trailing bits on the end don't refer to
                ! anything and should be set to 0. So, first set these to 1 and
                ! then flip all the bits.
                subsys_info(i)%B_masks(:,j) = A_mask
                do ipos = 0, i0_end
                    basis_find = sys%basis%basis_lookup(ipos, sys%basis%string_len)
                    if (basis_find == 0) then
                        subsys_info(i)%B_masks(sys%basis%string_len,j) = ibset(subsys_info(i)%B_masks(sys%basis%string_len,j),ipos)
                    end if
                end do
                subsys_info(i)%B_masks(:,j) = not(subsys_info(i)%B_masks(:,j))

            end do
        end do

        deallocate(lvecs, stat=ierr)
        call check_deallocate('lvecs', ierr)
        deallocate(sym_vecs,stat=ierr)
        call check_deallocate('sym_vecs',ierr)

    end subroutine find_rdm_masks

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
        use spawning, only: assign_particle_processor_dmqmc, add_spawned_particle

        integer(i0), intent(in) :: f(:)
        integer, intent(in) :: string_len, tensor_label_len
        integer(int_p), intent(in) :: nspawn
        integer, intent(in) :: particle_type
        type(spawn_t), intent(inout) :: spawn

        integer(i0) :: f_new(tensor_label_len)
#ifndef PARALLEL
        integer, parameter :: iproc_spawn = 0
#else
        integer :: iproc_spawn, slot
#endif

        ! Create the bitstring of the determinant.
        f_new = 0_i0
        f_new(:string_len) = f
        f_new((string_len+1):(tensor_label_len)) = f

#ifdef PARALLEL
        ! Need to determine which processor the spawned psip should be sent to.
        call assign_particle_processor_dmqmc(f_new, spawn%bit_str_nbits, spawn%hash_seed, spawn%hash_shift, spawn%move_freq, &
                                       nprocs, iproc_spawn, slot, spawn%proc_map%map, spawn%proc_map%nslots)
#endif

        call add_spawned_particle(f_new, nspawn, particle_type, iproc_spawn, spawn)

        if (spawn%error) call stop_all('create_diagonal_density_matrix_particle', 'Ran out of space in the spawned list while&
                                  & generating the initial density matrix.')

    end subroutine create_diagonal_density_matrix_particle

    subroutine create_diagonal_density_matrix_particle_initiator(f, string_len, tensor_label_len, nspawn, &
                                                                 & particle_type, initiator_pop, pop_real_factor, spawn)

        ! Create a psip on a diagonal element of the density matrix by adding
        ! it to the spawned walkers list. This list can then be sorted correctly
        ! by the direct_annihilation routine. This routine also checks if the
        ! psip resides on and initiator determinant.

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
        !    pop_real_factor: The factor by which populations are multiplied to
        !        enable non-integer populations.
        ! In/Out:
        !    spawn: spawn_t object to which the spawned particle is added.

        use hashing
        use parallel
        use errors, only: stop_all
        use spawn_data, only: spawn_t
        use spawning, only: assign_particle_processor_dmqmc, add_flagged_spawned_particle
        use idmqmc, only: set_parent_flag_dmqmc

        integer(i0), intent(in) :: f(:)
        integer, intent(in) :: string_len, tensor_label_len
        integer(int_p), intent(in) :: nspawn
        integer, intent(in) :: particle_type
        real(p), intent(in) :: initiator_pop
        integer(int_p), intent(in) :: pop_real_factor
        type(spawn_t), intent(inout) :: spawn

        integer(i0) :: f_new(tensor_label_len)
        integer :: flag
#ifndef PARALLEL
        integer, parameter :: iproc_spawn = 0
#else
        integer :: iproc_spawn, slot
#endif

        ! Create the bitstring of the determinant.
        f_new = 0_i0
        f_new(:string_len) = f
        f_new((string_len+1):(tensor_label_len)) = f

#ifdef PARALLEL
        ! Need to determine which processor the spawned psip should be sent to.
        call assign_particle_processor_dmqmc(f_new, spawn%bit_str_nbits, spawn%hash_seed, spawn%hash_shift, spawn%move_freq, &
                                       nprocs, iproc_spawn, slot, spawn%proc_map%map, spawn%proc_map%nslots)
#endif

        call set_parent_flag_dmqmc(real(nspawn,p)/pop_real_factor, initiator_pop, f_new, f_new, 0, flag)
        call add_flagged_spawned_particle(f_new, nspawn, particle_type, flag, iproc_spawn, spawn)

        if (spawn%error) call stop_all('create_diagonal_density_matrix_particle_initiator', &
                                        'Ran out of space in the spawned list while generating the initial density matrix.')

    end subroutine create_diagonal_density_matrix_particle_initiator

    subroutine decode_dm_bitstring(basis, f, isym, subsys_info, rdm_f1, rdm_f2)

        ! This function maps a full DMQMC bitstring to two bitstrings encoding
        ! the subsystem-A RDM bitstrings.

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
        !    subsys_info: information about the subsystem for the RDM being
        !        considered.
        ! Out:
        !     rdm_f1: the bitstring for the first label of the RDM.
        !     rdm_f2: the bitstring for the second label of the RDM.

        use basis_types, only: basis_t
        use dmqmc_data, only: subsys_t

        type(basis_t), intent(in) :: basis
        integer(i0), intent(in) :: f(:)
        integer, intent(in) :: isym
        type(subsys_t), intent(inout) :: subsys_info
        integer(i0), intent(out) :: rdm_f1(:), rdm_f2(:)

        integer :: i, bit_pos, bit_element

        ! Start from all bits down, so that we can flip bits up one by one.
        rdm_f1 = 0_i0
        rdm_f2 = 0_i0

        ! Loop over all the sites in the subsystem considered for the reduced
        ! density matrix.
        do i = 1, subsys_info%A_nsites
            ! Find the final bit positions and elements.
            bit_pos = basis%bit_lookup(1,i)
            bit_element = basis%bit_lookup(2,i)

            ! If the spin is up, set the corresponding bit in the first
            ! bitstring.
            if (btest(f(subsys_info%bit_pos(i,isym,2)),subsys_info%bit_pos(i,isym,1))) &
                rdm_f1(bit_element) = ibset(rdm_f1(bit_element),bit_pos)
            ! Similarly for the second index, by looking at the second end of
            ! the bitstring.
            if (btest(f(subsys_info%bit_pos(i,isym,2)+basis%string_len),subsys_info%bit_pos(i,isym,1))) &
                rdm_f2(bit_element) = ibset(rdm_f2(bit_element),bit_pos)
        end do

    end subroutine decode_dm_bitstring

    subroutine update_sampling_weights(rng, basis, qmc_in, psip_list, max_number_excitations, weighted_sampling)

        ! This routine updates the values of the weights used in importance
        ! sampling. It also removes or adds psips from the various excitation
        ! levels accordingly.

        ! In/Out:
        !    rng: random number generator.
        !    psip_list: main particle list.
        ! In:
        !    basis: information about the single-particle basis.
        !    qmc_in: Input options relating to QMC methods.
        !    max_number_excitations: maximum number of excitations allowed in the system.
        ! In/Out:
        !    weighted_sampling: type containing weighted sampling information.

        use annihilation, only: remove_unoccupied_dets
        use basis_types, only: basis_t
        use dSFMT_interface, only: dSFMT_t, get_rand_close_open
        use excitations, only: get_excitation_level
        use qmc_data, only: qmc_in_t, particle_t
        use dmqmc_data, only: dmqmc_weighted_sampling_t

        type(dSFMT_t), intent(inout) :: rng
        type(basis_t), intent(in) :: basis
        type(qmc_in_t), intent(in) :: qmc_in
        integer, intent(in) :: max_number_excitations
        type(particle_t), intent(inout) :: psip_list
        type(dmqmc_weighted_sampling_t), intent(inout) :: weighted_sampling

        integer :: idet, ireplica, excit_level
        real(p) :: new_population_target(psip_list%nspaces)
        integer(int_p) :: old_population(psip_list%nspaces), new_population(psip_list%nspaces)
        real(dp) :: r, pextra

        ! Alter weights for the next iteration.
        weighted_sampling%probs(:max_number_excitations) = &
                    & real(weighted_sampling%probs(:max_number_excitations),p)*weighted_sampling%altering_factors

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
            psip_list%nparticles = psip_list%nparticles + real(new_population - old_population, p)/psip_list%pop_real_factor

        end do

        ! Call the annihilation routine to update the main walker list, as some
        ! sites will have become unoccupied and so need removing from the
        ! simulation.
        call remove_unoccupied_dets(rng, psip_list, qmc_in%real_amplitudes)

    end subroutine update_sampling_weights

    subroutine output_and_alter_weights(dmqmc_in, max_number_excitations, excit_dist, weighted_sampling)

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
        ! In/Out:
        !    excit_dist: distribution of particles across excitations levels of
        !        the density matrix.
        !    weighted_sampling: type containing weighted sampling information.

        use dmqmc_data, only: dmqmc_in_t, dmqmc_weighted_sampling_t
        use parallel

        integer, intent(in) :: max_number_excitations
        type(dmqmc_in_t), intent(in) :: dmqmc_in
        real(p), intent(inout) :: excit_dist(0:)
        type(dmqmc_weighted_sampling_t), intent(inout) :: weighted_sampling

        integer :: i
#ifdef PARALLEL
        integer :: ierr
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
                weighted_sampling%sampling_probs(i) = weighted_sampling%sampling_probs(i)*&
                    (excit_dist(i)/excit_dist(i-1))
                weighted_sampling%sampling_probs(max_number_excitations+1-i) = weighted_sampling%sampling_probs(i)**(-1)
            end if
        end do

        ! Recalculate weighted_sampling%probs with the new weights.
        do i = 1, max_number_excitations
            weighted_sampling%probs(i) = weighted_sampling%probs(i-1)*weighted_sampling%sampling_probs(i)
        end do

        ! If vary_weights is true then the weights are to be introduced
        ! gradually at the start of each beta loop. This requires redefining
        ! weighted_sampling%altering_factors to coincide with the new sampling weights.
        if (dmqmc_in%vary_weights) then
            weighted_sampling%altering_factors = &
                & real(weighted_sampling%probs(:max_number_excitations),dp)**(1/real(dmqmc_in%finish_varying_weights,dp))
            ! Reset the weights for the next loop.
            weighted_sampling%probs = 1.0_p
        end if

        if (parent) then
            ! Print out weights in a form which can be copied into an input
            ! file.
            write(6, '(a31,2X)', advance = 'no') ' # Importance sampling weights:'
            do i = 1, max_number_excitations
                write (6, '(es12.4,2X)', advance = 'no') weighted_sampling%sampling_probs(i)
            end do
            write (6, '()', advance = 'yes')
        end if

    end subroutine output_and_alter_weights

end module dmqmc_procedures

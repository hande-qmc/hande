module qmc

! Launcher and initialisation routines for the various QMC algorithms.

use fciqmc_data

implicit none

contains

! --- Initialisation routines ---

    subroutine init_qmc(sys, qmc_in, restart_in, load_bal_in, reference_in, annihilation_flags, qmc_state, dmqmc_in, fciqmc_in)

        ! Initialisation for fciqmc calculations.
        ! Setup the spin polarisation for the system, initialise the RNG,
        ! allocate the required memory for the list of walkers and set the
        ! initial walker.

        ! In:
        !    sys: system being studied.
        !    restart_in: input options for HDF5 restart files.
        !    load_bal_in: input options for load balancing.
        !    fciqmc_in (optional): input options relating to FCIQMC.  Default
        !       fciqmc_in_t settings are used if not present.
        !    referemce_in: reference determinant.  If set (ie components
        !       allocated) then this is copied into qmc_state%ref.
        !       Otherwise a best guess is made based upon symmetry/spin/number
        !       of electrons/etc in set_reference_det.
        ! In/Out:
        !    qmc_in: input options relating to QMC methods.
        !    dmqmc_in (optional): input options relating to DMQMC.
        !    annihilation_flags: calculation specific annihilation flags.
        !    qmc_state: qmc_state_t object.  On output the QMC state is
        !       initialsed (potentially from a restart file) with components
        !       correctly allocated and useful information printed out...

        use checking, only: check_allocate, check_deallocate
        use errors, only: stop_all, warning
        use parallel
        use utils, only: int_fmt

        use basis, only: write_basis_fn
        use calc
        use dmqmc_procedures, only: init_dmqmc
        use determinants, only: decode_det, encode_det, write_det
        use energy_evaluation, only: nparticles_start_ind, calculate_hf_signed_pop
        use qmc_common, only: find_single_double_prob
        use reference_determinant, only: set_reference_det, copy_reference_t
        use particle_t_utils
        use proc_pointers, only: sc0_ptr, op0_ptr, energy_diff_ptr
        use spawn_data, only: alloc_spawn_t
        use spawning, only: assign_particle_processor
        use system
        use symmetry, only: symmetry_orb_list
        use utils, only: factorial_combination_1
        use restart_hdf5, only: read_restart_hdf5, restart_info_t, init_restart_info_t

        use qmc_data, only: qmc_in_t, fciqmc_in_t, restart_in_t, load_bal_in_t, annihilation_flags_t, qmc_state_t, &
                            reference_t, neel_singlet, no_guiding, single_basis
        use dmqmc_data, only: dmqmc_in_t

        type(sys_t), intent(in) :: sys
        type(qmc_in_t), intent(in) :: qmc_in
        type(restart_in_t), intent(in) :: restart_in
        type(load_bal_in_t), intent(in) :: load_bal_in
        type(reference_t), intent(in) :: reference_in
        type(annihilation_flags_t), intent(inout) :: annihilation_flags
        type(qmc_state_t), intent(inout) :: qmc_state
        type(dmqmc_in_t), intent(in), optional :: dmqmc_in
        type(fciqmc_in_t), intent(in), optional :: fciqmc_in

        integer :: ierr
        integer :: i, D0_proc, D0_inv_proc, ipos, occ_list0_inv(sys%nel), slot
        integer :: size_spawned_walker, max_nspawned_states
        integer :: nhash_bits
        integer :: ref_sym ! the symmetry of the reference determinant
        integer(i0) :: f0_inv(sys%basis%string_len)
        real(p) :: spawn_cutoff
        type(fciqmc_in_t) :: fciqmc_in_loc
        type(restart_info_t) :: ri
#ifdef PARALLEL
        integer(int_64) :: tmp_int_64
#endif

        if (present(fciqmc_in)) fciqmc_in_loc = fciqmc_in

        call copy_reference_t(reference_in, qmc_state%ref)

        ! --- Array sizes depending upon QMC algorithms ---

        associate(pl=>qmc_state%psip_list, reference=>qmc_state%ref, spawn=>qmc_state%spawn_store%spawn, &
                  spawn_recv=>qmc_state%spawn_store%spawn_recv)

            if (doing_calc(hfs_fciqmc_calc)) then
                pl%nspaces = pl%nspaces + 1
            else if (present(dmqmc_in)) then
                if (dmqmc_in%replica_tricks) pl%nspaces = pl%nspaces + 1
            end if

            ! Each determinant occupies string_len kind=i0 integers,
            ! pl%nspaces kind=int_p integers, pl%nspaces kind=p reals and one
            ! integer. If the Neel singlet state is used as the reference state for
            ! the projected estimator, then a further 2 reals are used per
            ! determinant.
            if (fciqmc_in_loc%trial_function == neel_singlet) pl%info_size = 2

            annihilation_flags%initiator_approx = qmc_in%initiator_approx
            annihilation_flags%real_amplitudes = qmc_in%real_amplitudes
            if (present(dmqmc_in)) then
                annihilation_flags%propagate_to_beta = dmqmc_in%propagate_to_beta
                annihilation_flags%replica_tricks = dmqmc_in%replica_tricks
            end if
            if (present(fciqmc_in)) then
                annihilation_flags%trial_function = fciqmc_in%trial_function
                qmc_state%trial%wfn = fciqmc_in%trial_function
                qmc_state%trial%guide = fciqmc_in%guiding_function
            end if

            ! Each spawned_walker occupies spawned_size kind=int_s integers.
            if (qmc_in%initiator_approx) then
                size_spawned_walker = (sys%basis%tensor_label_len+pl%nspaces+1)*int_s_length/8
            else
                size_spawned_walker = (sys%basis%tensor_label_len+pl%nspaces)*int_s_length/8
            end if
            max_nspawned_states = qmc_in%spawned_walker_length
            if (max_nspawned_states < 0) then
                ! Given in MB.  Convert.
                ! Note that we store 2 arrays.
                max_nspawned_states = int((-real(max_nspawned_states,p)*10**6)/(2*size_spawned_walker))
            end if

            if (parent) then
                write (6,'(1X,a57,f7.2)') &
                    'Memory allocated per core for spawned walker lists (MB): ', &
                    size_spawned_walker*real(2*max_nspawned_states,p)/10**6
                write (6,'(1X,a51,'//int_fmt(max_nspawned_states,1)//',/)') &
                    'Number of elements per core in spawned walker list:', max_nspawned_states
            end if

            ! --- Memory allocation ---

            call init_parallel_t(pl%nspaces, nparticles_start_ind-1, fciqmc_in_loc%non_blocking_comm, qmc_state%par_info, &
                                 load_bal_in%nslots)

            allocate(reference%f0(sys%basis%string_len), stat=ierr)

            ! Allocate main particle lists.  Include the memory used by semi_stoch_t%determ in the
            ! calculation of memory occupied by the main particle lists.
            call init_particle_t(qmc_in%walker_length, 1, sys%basis%tensor_label_len, qmc_in%real_amplitudes, &
                                 qmc_in%real_amplitude_force_32, pl)

            ! Allocate the shift.
            allocate(qmc_state%shift(pl%nspaces), stat=ierr)
            call check_allocate('qmc_state%shift', size(qmc_state%shift), ierr)
            allocate(qmc_state%vary_shift(pl%nspaces), stat=ierr)
            call check_allocate('qmc_state%vary_shift', size(qmc_state%vary_shift), ierr)
            qmc_state%shift = qmc_in%initial_shift

            ! Allocate spawned particle lists.
            if (mod(max_nspawned_states, nprocs) /= 0) then
                if (parent) write (6,'(1X,a68)') 'spawned_walker_length is not a multiple of the number of processors.'
                max_nspawned_states = ceiling(real(max_nspawned_states)/nprocs)*nprocs
                if (parent) write (6,'(1X,a35,'//int_fmt(max_nspawned_states,1)//',a1,/)') &
                                            'Increasing spawned_walker_length to',max_nspawned_states,'.'
            end if

            ! If not using real amplitudes then we always want spawn_cutoff to be
            ! equal to 1.0, so overwrite the default before creating spawn_t objects.
            spawn_cutoff = qmc_in%spawn_cutoff
            if (.not. qmc_in%real_amplitudes) spawn_cutoff = 0.0_p

            if (doing_calc(dmqmc_calc)) then
                ! Hash the entire first bit array and the minimum number of bits
                ! in the second bit array.
                nhash_bits = sys%basis%nbasis + i0_length*sys%basis%string_len
            else
                nhash_bits = sys%basis%nbasis
            end if
            call alloc_spawn_t(sys%basis%tensor_label_len, nhash_bits, pl%nspaces, qmc_in%initiator_approx, &
                               max_nspawned_states, spawn_cutoff, pl%pop_real_factor, qmc_state%par_info%load%proc_map, 7, &
                               qmc_in%use_mpi_barriers, spawn)
            if (fciqmc_in_loc%non_blocking_comm) then
                call alloc_spawn_t(sys%basis%tensor_label_len, nhash_bits, pl%nspaces, qmc_in%initiator_approx, &
                                   max_nspawned_states, spawn_cutoff, pl%pop_real_factor, qmc_state%par_info%load%proc_map, 7, &
                                   .false., spawn_recv)
            end if

            call check_allocate('reference%f0',sys%basis%string_len,ierr)
            allocate(reference%hs_f0(sys%basis%string_len), stat=ierr)
            call check_allocate('reference%hs_f0', size(reference%hs_f0), ierr)

            ! --- Importance sampling ---

            if (qmc_state%trial%wfn == neel_singlet) then
                ! Calculate all the possible different amplitudes for the Neel singlet state
                ! and store them in an array
                allocate(qmc_state%trial%wfn_dat(-1:(sys%lattice%nsites/2)+1), stat=ierr)
                call check_allocate('qmc_state%trial%wfn_dat',(sys%lattice%nsites/2)+1,ierr)

                qmc_state%trial%wfn_dat(-1) = 0
                qmc_state%trial%wfn_dat((sys%lattice%nsites/2)+1) = 0
                do i=0,(sys%lattice%nsites/2)
                    qmc_state%trial%wfn_dat(i) = factorial_combination_1( (sys%lattice%nsites/2)-i , i )
                    qmc_state%trial%wfn_dat(i) = -(2*mod(i,2)-1) * qmc_state%trial%wfn_dat(i)
                end do
            end if

            ! --- Initial walker distributions ---
            ! Note occ_list could be set and allocated in the input.

            if (restart_in%read_restart) then
                if (.not.allocated(reference%occ_list0)) then
                    allocate(reference%occ_list0(sys%nel), stat=ierr)
                    call check_allocate('reference%occ_list0',sys%nel,ierr)
                end if
                call init_restart_info_t(ri, read_id=restart_in%read_id)
                call read_restart_hdf5(ri, fciqmc_in_loc%non_blocking_comm, qmc_state)
                ! Need to re-calculate the reference determinant data
                call decode_det(sys%basis, reference%f0, reference%occ_list0)
                if (fciqmc_in_loc%trial_function == neel_singlet) then
                    ! Set the Neel state data for the reference state, if it is being used.
                    reference%H00 = 0.0_p
                else
                    reference%H00 = sc0_ptr(sys, reference%f0)
                end if
                if (doing_calc(hfs_fciqmc_calc)) reference%O00 = op0_ptr(sys, reference%f0)
            else

                ! Reference det

                ! Set the reference determinant to be the spin-orbitals with the lowest
                ! single-particle eigenvalues which satisfy the spin polarisation.
                ! Note: this is for testing only!  The symmetry input is currently
                ! ignored.
                if (sys%symmetry < sys%sym_max) then
                    call set_reference_det(sys, reference%occ_list0, .false., sys%symmetry)
                else
                    call set_reference_det(sys, reference%occ_list0, .false.)
                end if

                call encode_det(sys%basis, reference%occ_list0, reference%f0)

                if (allocated(reference%hs_occ_list0)) then
                    call encode_det(sys%basis, reference%hs_occ_list0, reference%hs_f0)
                else
                    allocate(reference%hs_occ_list0(sys%nel), stat=ierr)
                    call check_allocate('reference%hs_occ_list0', size(reference%hs_occ_list0), ierr)
                    reference%hs_occ_list0 = reference%occ_list0
                    reference%hs_f0 = reference%f0
                end if

                ! Energy of reference determinant.
                reference%H00 = sc0_ptr(sys, reference%f0)
                ! Exchange energy of reference determinant.
                select case (sys%system)
                case (ueg, read_in)
                    reference%energy_shift = energy_diff_ptr(sys, reference%occ_list0)
                case default
                    ! [todo] - Implement for all models.
                end select
                if (doing_calc(hfs_fciqmc_calc)) reference%O00 = op0_ptr(sys, reference%f0)

                ! In general FCIQMC, we start with psips only on the
                ! reference determinant, so set pl%nstates = 1 and
                ! initialise pl%pops. For DMQMC, this is
                ! not required, as psips are spawned along the diagonal
                ! initially.
                if (doing_calc(dmqmc_calc)) then
                    ! Initial distribution handled later
                    pl%nstates = 0
                else
                    pl%nstates = 1
                    ! Zero all populations...
                    pl%pops(:,pl%nstates) = 0_int_p
                    ! Set initial population of Hamiltonian walkers.
                    pl%pops(1,pl%nstates) = nint(qmc_in%D0_population)*pl%pop_real_factor
                    ! Set the bitstring of this psip to be that of the
                    ! reference state.
                    pl%states(:,pl%nstates) = reference%f0

                    ! Determine and set properties for the reference state which we start on.
                    ! By definition, when using a single determinant as a reference state:
                    pl%dat(1,pl%nstates) = 0.0_p
                    ! Or if not using a single determinant:
                    if (fciqmc_in_loc%trial_function == neel_singlet) then
                        ! Set the Neel state data for the reference state, if it is being used.
                        pl%dat(1,pl%nstates) = reference%H00
                        reference%H00 = 0.0_p

                        pl%dat(pl%nspaces+1,pl%nstates) = sys%lattice%nsites/2
                        ! For a rectangular bipartite lattice, nbonds = ndim*nsites.
                        ! The Neel state cannot be used for non-bipartite lattices.
                        pl%dat(pl%nspaces+2,pl%nstates) = sys%lattice%ndim*sys%lattice%nsites
                    end if

                    ! Finally, we need to check if the reference determinant actually
                    ! belongs on this processor.
                    ! If it doesn't, set the walkers array to be empty.
                    associate(pm=>spawn%proc_map)
                        call assign_particle_processor(reference%f0, spawn%bit_str_nbits, spawn%hash_seed, spawn%hash_shift, &
                                                       spawn%move_freq, nprocs, D0_proc, slot, pm%map, pm%nslots)
                    end associate
                    if (D0_proc /= iproc) pl%nstates = 0
                end if

                ! For the Heisenberg model and open shell systems, it is often useful to
                ! have psips start on both the reference state and the spin-flipped version.
                if (fciqmc_in_loc%init_spin_inv_D0) then

                    ! Need to handle the Heisenberg model (consisting of spinors on
                    ! lattice sites) and electron systems differently, as the
                    ! Heisenberg model has no concept of unoccupied basis
                    ! functions/holes.
                    select case (sys%system)
                    case (heisenberg)
                        ! Flip all spins in f0 to get f0_inv
                        f0_inv = not(reference%f0)
                        ! In general, the basis bit string has some padding at the
                        ! end which must be unset.  We need to clear this...
                        ! Loop over all bits after the last basis function.
                        i = sys%basis%bit_lookup(2,sys%basis%nbasis)
                        do ipos = sys%basis%bit_lookup(1,sys%basis%nbasis)+1, i0_end
                            f0_inv(i) = ibclr(f0_inv(i), ipos)
                        end do
                    case default
                        ! Swap each basis function for its spin-inverse
                        ! This looks somewhat odd, but relies upon basis
                        ! functions alternating down (Ms=-1) and up (Ms=1).
                        do i = 1, sys%nel
                            if (mod(reference%occ_list0(i),2) == 1) then
                                occ_list0_inv(i) = reference%occ_list0(i) + 1
                            else
                                occ_list0_inv(i) = reference%occ_list0(i) - 1
                            end if
                        end do
                        call encode_det(sys%basis, occ_list0_inv, f0_inv)
                    end select

                    associate(pm=>spawn%proc_map)
                        call assign_particle_processor(f0_inv, spawn%bit_str_nbits, spawn%hash_seed, spawn%hash_shift, &
                                                       spawn%move_freq, nprocs, D0_inv_proc, slot, pm%map, pm%nslots)
                    end associate

                    ! Store if not identical to reference det.
                    if (D0_inv_proc == iproc .and. any(reference%f0 /= f0_inv)) then
                        pl%nstates = pl%nstates + 1
                        ! Zero all populations for this determinant.
                        pl%pops(:,pl%nstates) = 0_int_p
                        ! Set the population for this basis function.
                        pl%pops(1,pl%nstates) = nint(qmc_in%D0_population)*pl%pop_real_factor
                        pl%dat(1,pl%nstates) = sc0_ptr(sys, reference%f0) - reference%H00
                        select case(sys%system)
                        case(heisenberg)
                            if (fciqmc_in_loc%trial_function /= single_basis) then
                                pl%dat(1,pl%nstates) = 0
                            else
                                pl%dat(1,pl%nstates) = sc0_ptr(sys, reference%f0) - reference%H00
                            end if
                        case default
                            pl%dat(1,pl%nstates) = sc0_ptr(sys, reference%f0) - reference%H00
                        end select
                        pl%states(:,pl%nstates) = f0_inv
                        ! If we are using the Neel state as a reference in the
                        ! Heisenberg model, then set the required data.
                        if (fciqmc_in_loc%trial_function == neel_singlet) then
                            pl%dat(pl%nspaces+1,pl%nstates) = 0
                            pl%dat(pl%nspaces+2,pl%nstates) = 0
                        end if
                    end if
                end if

            end if ! End of initialisation of reference state(s)/restarting from previous calculations

            ! Total number of particles on processor.
            ! Probably should be handled more simply by setting it to be either 0 or
            ! D0_population or obtaining it from the restart file, as appropriate.
            forall (i=1:pl%nspaces) pl%nparticles(i) = sum(abs( real(pl%pops(i,:pl%nstates),p)/pl%pop_real_factor))
            ! Should we already be in varyshift mode (e.g. restarting a calculation)?
#ifdef PARALLEL
            do i=1, pl%nspaces
                call mpi_allgather(pl%nparticles(i), 1, MPI_PREAL, pl%nparticles_proc(i,:), 1, MPI_PREAL, MPI_COMM_WORLD, ierr)
            end do
            ! When restarting a non-blocking calculation this sum will not equal
            ! tot_nparticles as some walkers have been communicated around the report
            ! loop. The correct total is in the restart file so get it from there.
            if (.not. (restart_in%read_restart .and. fciqmc_in_loc%non_blocking_comm)) &
                forall(i=1:pl%nspaces) pl%tot_nparticles(i) = sum(pl%nparticles_proc(i,:))
#else
            pl%tot_nparticles = pl%nparticles
            pl%nparticles_proc(:pl%nspaces,1) = pl%nparticles(:pl%nspaces)
#endif

            ! Decide whether the shift should be turned on from the start.
            qmc_state%target_particles = qmc_in%target_particles
            qmc_state%vary_shift = pl%tot_nparticles >= qmc_state%target_particles

            if (doing_calc(hfs_fciqmc_calc)) then
#ifdef PARALLEL
                tmp_int_64 = calculate_hf_signed_pop(pl)
                call mpi_allreduce(tmp_int_64, qmc_state%estimators%hf_signed_pop, pl%nspaces, MPI_INTEGER8, MPI_SUM, &
                                   MPI_COMM_WORLD, ierr)
#else
                qmc_state%estimators%hf_signed_pop = calculate_hf_signed_pop(pl)
#endif
            end if

            ! calculate the reference determinant symmetry
            ref_sym = symmetry_orb_list(sys, reference%occ_list0)

            ! If not set at input, set probability of selecting single or double
            ! excitations based upon the reference determinant and assume other
            ! determinants have a roughly similar ratio of single:double
            ! excitations.
            if (qmc_in%pattempt_single < 0 .or. qmc_in%pattempt_double < 0) then
                call find_single_double_prob(sys, reference%occ_list0, qmc_state%pattempt_single, qmc_state%pattempt_double)
            else
                ! renormalise just in case input wasn't
                qmc_state%pattempt_single = qmc_in%pattempt_single/(qmc_in%pattempt_single+qmc_in%pattempt_double)
                qmc_state%pattempt_double = 1.0_p - qmc_in%pattempt_single
            end if

            ! Set initial values from input
            qmc_state%tau = qmc_in%tau
            qmc_state%estimators%D0_population = qmc_in%D0_population

        end associate

        if (parent) then
            write (6,'(1X,"Reference determinant, |D0> =",1X)',advance='no')
            call write_det(sys%basis, sys%nel, qmc_state%ref%f0, new_line=.true.)
            write (6,'(1X,"E0 = <D0|H|D0> =",f20.12)') qmc_state%ref%H00
            if (doing_calc(hfs_fciqmc_calc)) write (6,'(1X,"O00 = <D0|O|D0> =",f20.12)')  qmc_state%ref%O00
            write(6,'(1X,"Symmetry of reference determinant:")',advance='no')
            select case(sys%system)
            case (hub_k)
                call write_basis_fn(sys, sys%basis%basis_fns(2*ref_sym), new_line=.true., print_full=.false.)
            case default
                write(6,'(1X,i0)') ref_sym
            end select

            if (doing_calc(dmqmc_calc)) then
                write (6,'(1X,"Initial population on the trace of the density matrix:",1X,i0)') int(qmc_in%D0_population,int_64)
            else
                write (6,'(1X,"Initial population on reference determinant:",1X,f11.4)') qmc_in%D0_population
                write (6,'(1X,"Note that the correlation energy is relative to |D0>.")')
            end if
            write (6,'()')

        end if

    end subroutine init_qmc

    subroutine init_proc_pointers(sys, qmc_in, reference, dmqmc_in, fciqmc_in)

        ! Set function pointers for QMC calculations.

        ! In:
        !    sys: system being studied.
        !    qmc_in: input options relating to QMC methods.
        !    reference: reference_t object defining the reference state/determinant.
        !    dmqmc_in (optional): input options relating to DMQMC, only required if doing DMQMC.
        !    fciqmc_in (optional): input options relating to FCIQMC, only required if doing FCIQMC.

        ! System and calculation data
        use calc, only: doing_calc, doing_dmqmc_calc, dmqmc_calc, hfs_fciqmc_calc, &
                        ras, dmqmc_correlation, dmqmc_energy, dmqmc_energy_squared, dmqmc_staggered_magnetisation, &
                        dmqmc_kinetic_energy, dmqmc_H0_energy, dmqmc_potential_energy
        use hfs_data
        use system
        use parallel, only: parent
        use qmc_data, only: qmc_in_t, fciqmc_in_t, reference_t, single_basis, neel_singlet, neel_singlet_guiding, &
                            excit_gen_renorm, excit_gen_no_renorm
        use dmqmc_data, only: dmqmc_in_t, free_electron_dm

        ! Procedures to be pointed to.
        use death, only: stochastic_death
        use determinants
        use dmqmc_estimators
        use dmqmc_procedures
        use energy_evaluation
        use excit_gen_mol
        use excit_gen_op_mol
        use excit_gen_hub_k
        use excit_gen_op_hub_k
        use excit_gen_real_lattice
        use excit_gen_ringium
        use excit_gen_ueg, only: gen_excit_ueg_no_renorm
        use hamiltonian_chung_landau, only: slater_condon0_chung_landau
        use hamiltonian_hub_k, only: slater_condon0_hub_k
        use hamiltonian_hub_real, only: slater_condon0_hub_real
        use hamiltonian_heisenberg, only: diagonal_element_heisenberg, diagonal_element_heisenberg_staggered
        use hamiltonian_molecular, only: slater_condon0_mol, double_counting_correction_mol
        use hamiltonian_ringium, only: slater_condon0_ringium
        use hamiltonian_ueg, only: slater_condon0_ueg, kinetic_energy_ueg, exchange_energy_ueg, potential_energy_ueg
        use heisenberg_estimators
        use importance_sampling
        use operators
        use spawning

        ! Procedure pointers
        use proc_pointers

        ! Utilities
        use errors, only: stop_all

        type(sys_t), intent(in) :: sys
        type(qmc_in_t), intent(in) :: qmc_in
        type(reference_t), intent(in) :: reference
        type(dmqmc_in_t), intent(in), optional :: dmqmc_in
        type(fciqmc_in_t), intent(in), optional :: fciqmc_in

        logical :: truncate_space

        ! 0. In general, use the default spawning routine.
        spawner_ptr => spawn_standard

        ! 1. Set system-specific procedure pointers.
        !     * projected energy estimator
        !     * diagonal hamiltonian matrix element evaluation
        !     * spawning
        !     * excitation generators
        select case(sys%system)
        case(hub_k)

            decoder_ptr => decode_det_spinocc_spinunocc
            update_proj_energy_ptr => update_proj_energy_hub_k
            sc0_ptr => slater_condon0_hub_k
            spawner_ptr => spawn_lattice_split_gen
            select case(qmc_in%excit_gen)
            case(excit_gen_no_renorm)
                gen_excit_ptr%full => gen_excit_hub_k_no_renorm
                gen_excit_ptr%init => gen_excit_init_hub_k_no_renorm
                gen_excit_ptr%finalise => gen_excit_finalise_hub_k_no_renorm
            case(excit_gen_renorm)
                gen_excit_ptr%full => gen_excit_hub_k
                gen_excit_ptr%init => gen_excit_init_hub_k
                gen_excit_ptr%finalise => gen_excit_finalise_hub_k
            case default
                call stop_all('init_proc_pointers', 'Selected excitation generator not implemented.')
            end select

        case(hub_real, chung_landau)

            ! The Hubbard model in a local orbital basis and the Chung--Landau
            ! Hamiltonian have the same off-diagonal operator so we can use the
            ! same excitation generators and just a different function for the
            ! diagonal matrix elements.
            ! Note that the Chung--Landau model (as in Phys Rev B 85 (2012)
            ! 115115) is contains spinless fermions.
            decoder_ptr => decode_det_occ
            update_proj_energy_ptr => update_proj_energy_hub_real
            if (sys%system == hub_real) then
                sc0_ptr => slater_condon0_hub_real
            else
                sc0_ptr => slater_condon0_chung_landau
            end if

            select case(qmc_in%excit_gen)
            case(excit_gen_no_renorm)
                gen_excit_ptr%full => gen_excit_hub_real_no_renorm
            case(excit_gen_renorm)
                gen_excit_ptr%full => gen_excit_hub_real
            case default
                call stop_all('init_proc_pointers', 'Selected excitation generator not implemented.')
            end select

        case(heisenberg)

            ! Only need occupied orbitals list, as for the real Hubbard case.
            decoder_ptr => decode_det_occ
            ! Set which trial wavefunction to use for the energy estimator.
            update_proj_energy_ptr => update_proj_energy_heisenberg_basic
            if (present(fciqmc_in)) then
                select case(fciqmc_in%trial_function)
                case (single_basis)
                    update_proj_energy_ptr => update_proj_energy_heisenberg_basic
                case (neel_singlet)
                    update_proj_energy_ptr => update_proj_energy_heisenberg_neel_singlet
                end select
            end if

            ! Set whether the applied staggered magnetisation is non-zero.
            if (abs(sys%heisenberg%staggered_magnetic_field) > depsilon) then
                sc0_ptr => diagonal_element_heisenberg_staggered
            else
                sc0_ptr => diagonal_element_heisenberg
            end if

            ! Set which guiding wavefunction to use, if requested.
            select case(qmc_in%excit_gen)
            case(excit_gen_no_renorm)
                gen_excit_ptr%full => gen_excit_heisenberg_no_renorm
            case(excit_gen_renorm)
                gen_excit_ptr%full => gen_excit_heisenberg
            case default
                call stop_all('init_proc_pointers', 'Selected excitation generator not implemented.')
            end select
            if (present(fciqmc_in)) then
                select case(fciqmc_in%guiding_function)
                case (neel_singlet_guiding)
                    spawner_ptr => spawn_importance_sampling
                    gen_excit_ptr%trial_fn => neel_trial_state
                end select
            end if

        case(read_in)

            update_proj_energy_ptr => update_proj_energy_mol
            sc0_ptr => slater_condon0_mol
            energy_diff_ptr => double_counting_correction_mol

            select case(qmc_in%excit_gen)
            case(excit_gen_no_renorm)
                gen_excit_ptr%full => gen_excit_mol_no_renorm
                decoder_ptr => decode_det_occ
            case(excit_gen_renorm)
                gen_excit_ptr%full => gen_excit_mol
                decoder_ptr => decode_det_occ_symunocc
            case default
                call stop_all('init_proc_pointers', 'Selected excitation generator not implemented.')
            end select

        case(ueg)

            update_proj_energy_ptr => update_proj_energy_ueg
            sc0_ptr => slater_condon0_ueg
            energy_diff_ptr => exchange_energy_ueg

            gen_excit_ptr%full => gen_excit_ueg_no_renorm
            decoder_ptr => decode_det_occ
            select case(qmc_in%excit_gen)
            case(excit_gen_no_renorm)
            case(excit_gen_renorm)
                if (parent) then
                    write (6,'(1X,"WARNING: renormalised excitation generators not implemented.")')
                    write (6,'(1X,"WARNING: If this upsets you, please send patches.",/)')
                end if
            case default
                call stop_all('init_proc_pointers', 'Selected excitation generator not implemented.')
            end select

        case(ringium)

            update_proj_energy_ptr => update_proj_energy_ringium
            sc0_ptr => slater_condon0_ringium

            gen_excit_ptr%full => gen_excit_ringium_no_renorm
            decoder_ptr => decode_det_occ
            select case(qmc_in%excit_gen)
            case(excit_gen_no_renorm)
            case(excit_gen_renorm)
                if (parent) then
                    write (6,'(1X,"WARNING: renormalised excitation generators not implemented.")')
                    write (6,'(1X,"WARNING: If this upsets you, please send patches.",/)')
                end if
            case default
                call stop_all('init_proc_pointers', 'Selected excitation generator not implemented.')
            end select

        case default
            call stop_all('init_proc_pointers','QMC not implemented for this system yet.')
        end select

        ! 2. Set calculation-specific procedure pointers
        truncate_space = reference%ex_level /= sys%nel

        ! 2: initiator-approximation
        if (qmc_in%initiator_approx) then
            if (all(ras > 0)) then
                create_spawned_particle_ptr => create_spawned_particle_initiator_ras
            else if (truncate_space) then
                create_spawned_particle_ptr => create_spawned_particle_initiator_truncated
            else
                create_spawned_particle_ptr => create_spawned_particle_initiator
            end if
        else
            if (all(ras > 0)) then
                create_spawned_particle_ptr => create_spawned_particle_ras
            else if (truncate_space) then
                create_spawned_particle_ptr => create_spawned_particle_truncated
            else
                create_spawned_particle_ptr => create_spawned_particle
            end if
        end if

        ! 2: density-matrix
        if (doing_calc(dmqmc_calc)) then

            if (.not.present(dmqmc_in)) call stop_all('init_proc_pointers', 'DMQMC options not present.')

            ! Spawned particle creation.
            if (dmqmc_in%half_density_matrix) then
                if (truncate_space) then
                    create_spawned_particle_dm_ptr => create_spawned_particle_truncated_half_density_matrix
                else
                    if (qmc_in%initiator_approx) then
                        create_spawned_particle_dm_ptr => create_spawned_particle_half_density_matrix_initiator
                    else
                        create_spawned_particle_dm_ptr => create_spawned_particle_half_density_matrix
                    end if
                end if
            else
                if (truncate_space) then
                    create_spawned_particle_dm_ptr => create_spawned_particle_truncated_density_matrix
                else
                    if (qmc_in%initiator_approx) then
                        create_spawned_particle_dm_ptr => create_spawned_particle_density_matrix_initiator
                    else
                        create_spawned_particle_dm_ptr => create_spawned_particle_density_matrix
                    end if
                end if
            end if

            ! Weighted importance sampling routines.
            if (dmqmc_in%weighted_sampling) then
                spawner_ptr => spawn_importance_sampling
                gen_excit_ptr%trial_fn => dmqmc_weighting_fn
            end if

            ! Expectation values.
            if (doing_dmqmc_calc(dmqmc_energy)) then
                if (dmqmc_in%propagate_to_beta) then
                    update_dmqmc_energy_and_trace_ptr => dmqmc_energy_and_trace_propagate
                else
                    update_dmqmc_energy_and_trace_ptr => dmqmc_energy_and_trace
                end if
            end if
            select case(sys%system)
            case(heisenberg)
                if (doing_dmqmc_calc(dmqmc_energy_squared)) &
                                         update_dmqmc_energy_squared_ptr => dmqmc_energy_squared_heisenberg
                if (doing_dmqmc_calc(dmqmc_correlation)) update_dmqmc_correlation_ptr => &
                                     dmqmc_correlation_function_heisenberg
                if (doing_dmqmc_calc(dmqmc_staggered_magnetisation)) &
                                         update_dmqmc_stag_mag_ptr => dmqmc_stag_mag_heisenberg
            case(ueg)
                if (dmqmc_in%propagate_to_beta) then
                    if (dmqmc_in%initial_matrix == free_electron_dm) then
                        trial_dm_ptr => kinetic_energy_ueg
                    else
                        trial_dm_ptr => slater_condon0_ueg
                    end if
                end if
                    if (doing_dmqmc_calc(dmqmc_kinetic_energy)) then
                        kinetic_diag_ptr => kinetic_energy_ueg
                        update_dmqmc_kinetic_energy_ptr => dmqmc_kinetic_energy_diag
                    end if
                    if (doing_dmqmc_calc(dmqmc_potential_energy)) then
                        potential_energy_ptr => potential_energy_ueg
                    end if
                    if (doing_dmqmc_calc(dmqmc_H0_energy) .and. .not. dmqmc_in%propagate_to_beta) then
                        ! Assume that we want to evaluate <H_HF> rather than the kinetic energy.
                        trial_dm_ptr => slater_condon0_ueg
                    end if
            case(hub_k)
                if (dmqmc_in%propagate_to_beta) then
                    if (dmqmc_in%initial_matrix == free_electron_dm) then
                        trial_dm_ptr => kinetic0_hub_k
                    else
                        trial_dm_ptr => slater_condon0_hub_k
                    end if
                end if
            case(read_in)
                if (dmqmc_in%propagate_to_beta) then
                    trial_dm_ptr => slater_condon0_mol
                end if
            end select

        end if

        ! 2: Hellmann--Feynman operator sampling
        if (doing_calc(hfs_fciqmc_calc)) then
            select case(hf_operator)
            case(hamiltonian_operator)
                op0_ptr => sc0_ptr
                update_proj_hfs_ptr => update_proj_hfs_hamiltonian
                spawner_hfs_ptr => spawner_ptr
            case(kinetic_operator)
                update_proj_hfs_ptr => update_proj_hfs_diagonal
                spawner_hfs_ptr => spawn_null
                if (sys%system == hub_k) then
                    op0_ptr => kinetic0_hub_k
                else
                    call stop_all('init_proc_pointers','System not yet supported in HFS with operator given.')
                end if
            case(double_occ_operator)
                if (sys%system == hub_k) then
                    ! Shamelessly re-use the Hamiltonian excitation generators.
                    gen_excit_hfs_ptr%full => gen_excit_ptr%full
                    gen_excit_hfs_ptr%init => gen_excit_ptr%init
                    gen_excit_hfs_ptr%finalise => gen_excit_ptr%finalise
                    spawner_hfs_ptr => spawn_lattice_split_gen_importance_sampling
                    ! Scale the Hamiltonian matrix element to obtain the matrix
                    ! element of this operator.
                    gen_excit_hfs_ptr%trial_fn => gen_excit_double_occ_matel_hub_k
                    update_proj_hfs_ptr => update_proj_hfs_double_occ_hub_k
                    op0_ptr => double_occ0_hub_k
                else
                    call stop_all('init_proc_pointers','System not yet supported in HFS with operator given.')
                end if
            case(dipole_operator)
                if (sys%system == read_in) then
                    op0_ptr => one_body0_mol
                    update_proj_hfs_ptr => update_proj_hfs_one_body_mol
                    spawner_hfs_ptr => spawner_ptr

                    select case(qmc_in%excit_gen)
                    case(excit_gen_no_renorm)
                        gen_excit_hfs_ptr%full => gen_excit_one_body_mol_no_renorm
                    case(excit_gen_renorm)
                        gen_excit_hfs_ptr%full => gen_excit_one_body_mol
                    case default
                        call stop_all('init_proc_pointers', 'Selected excitation generator not implemented.')
                    end select
                else
                    call stop_all('init_proc_pointers','System not yet supported in HFS with operator given.')
                end if
            case default
                call stop_all('init_proc_pointers','Operator given is not yet supported')
            end select
        end if

    end subroutine init_proc_pointers

end module qmc

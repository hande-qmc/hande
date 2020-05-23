module qmc

! Launcher and initialisation routines for the various QMC algorithms.

implicit none

contains

! --- Initialisation routines ---

    subroutine init_qmc(sys, qmc_in, restart_in, load_bal_in, reference_in, io_unit, annihilation_flags, qmc_state, uuid_restart, &
                        restart_version_restart, dmqmc_in, fciqmc_in, qmc_state_restart, regenerate_info)

        ! Initialisation for fciqmc calculations.
        ! Setup the spin polarisation for the system, initialise the RNG,
        ! allocate the required memory for the list of walkers and set the
        ! initial walker.

        ! In:
        !    sys: system being studied.
        !    qmc_in: input options relating to QMC methods.
        !    restart_in: input options for HDF5 restart files.
        !    load_bal_in: input options for load balancing.
        !    reference_in: reference determinant.  If set (ie components
        !       allocated) then this is copied into qmc_state%ref.
        !       Otherwise a best guess is made based upon symmetry/spin/number
        !       of electrons/etc in set_reference_det.
        !    io_unit: io unit to write all output to.
        !    dmqmc_in (optional): input options relating to DMQMC.
        !    fciqmc_in (optional): input options relating to FCIQMC.  Default
        !       fciqmc_in_t settings are used if not present.
        ! In/Out:
        !    qmc_state_restart (optional): qmc_state_t object from a calculation
        !       to restart. Deallocated on exit.
        ! Out:
        !    annihilation_flags: calculation specific annihilation flags.
        !    qmc_state: qmc_state_t object.  On output the QMC state is
        !       initialised (potentially from a restart file) with components
        !       correctly allocated and useful information printed out...
        !    uuid_restart: if using a restart file, the UUID of the calculations
        !       that generated it.
        !    restart_version_restart: version of restart that was read in.
        !    regenerate_info (optional): true if additional information within
        !       bit string needs to be regenerated for a read-in restart file.

        use checking, only: check_allocate

        use calc, only: doing_calc, hfs_fciqmc_calc, dmqmc_calc, GLOBAL_META
        use energy_evaluation, only: get_comm_processor_indx, est_buf_data_size, est_buf_n_per_proc
        use load_balancing, only: init_parallel_t
        use particle_t_utils, only: init_particle_t
        use system
        use restart_hdf5, only: read_restart_hdf5, restart_info_t, init_restart_info_t, get_reference_hdf5

        use qmc_data, only: qmc_in_t, fciqmc_in_t, restart_in_t, load_bal_in_t, annihilation_flags_t, qmc_state_t, &
                            neel_singlet
        use reference_determinant, only: reference_t
        use dmqmc_data, only: dmqmc_in_t
        use excit_gens, only: dealloc_excit_gen_data_t
        use const, only: p
        use parallel, only: parent
        use hamiltonian_ueg, only: calc_fock_values_3d_ueg
        use determinants, only: sum_fock_values_occ_list
        use, intrinsic :: iso_fortran_env, only: error_unit
        use utils, only: int_fmt


        type(sys_t), intent(in) :: sys
        type(qmc_in_t), intent(in) :: qmc_in
        type(restart_in_t), intent(in) :: restart_in
        type(load_bal_in_t), intent(in) :: load_bal_in
        type(reference_t), intent(in) :: reference_in
        integer, intent(in) :: io_unit
        type(annihilation_flags_t), intent(out) :: annihilation_flags
        type(qmc_state_t), intent(out) :: qmc_state
        character(36), intent(out) :: uuid_restart
        integer, intent(out) :: restart_version_restart
        type(dmqmc_in_t), intent(in), optional :: dmqmc_in
        type(fciqmc_in_t), intent(in), optional :: fciqmc_in
        type(qmc_state_t), intent(inout), optional :: qmc_state_restart
        logical, intent(out), optional :: regenerate_info

        integer :: ierr, iunit, proc_data_info(2,est_buf_n_per_proc), ntot_proc_data
        type(fciqmc_in_t) :: fciqmc_in_loc
        type(dmqmc_in_t) :: dmqmc_in_loc
        type(restart_info_t) :: ri
        logical :: regenerate_info_loc, qmc_state_restart_loc
        real :: t1, t2, set_up_time

        iunit = 6

        regenerate_info_loc = .false.
        restart_version_restart = 0
        qmc_state_restart_loc = .false.

        if (present(fciqmc_in)) fciqmc_in_loc = fciqmc_in
        if (present(dmqmc_in)) dmqmc_in_loc = dmqmc_in
        if (present(qmc_state_restart)) qmc_state_restart_loc = .true.

        if (restart_in%read_restart) call init_restart_info_t(ri, read_id=restart_in%read_id)

        call init_proc_pointers(sys, qmc_in, reference_in, io_unit, dmqmc_in, fciqmc_in)

        ! Note it is not possible to override a reference if restarting.
        if (restart_in%read_restart) then
            call init_reference_restart(sys, reference_in, ri, qmc_state%ref)
        else if (qmc_state_restart_loc) then
            qmc_state%ref = qmc_state_restart%ref
        else
            call init_reference(sys, reference_in, io_unit, qmc_state%ref)
        end if
        
        ! In read-in systems, sp_fock = sp_eigv. This is different in the case of the UEG, where sp_fock is filled with <i|F|i>.
        ! For the UEG, the Fock values are only implemented for the 3D version!
        ! [todo] - implement 2D, etc.
        allocate(qmc_state%propagator%sp_fock(sys%basis%nbasis), stat=ierr)
        call check_allocate('qmc_state%propagator%sp_fock', sys%basis%nbasis, ierr)
        qmc_state%propagator%sp_fock = sys%basis%basis_fns%sp_eigv
        if ((sys%system == ueg) .and. (sys%lattice%ndim == 3)) then
            call calc_fock_values_3d_ueg(sys, qmc_state%propagator, qmc_state%ref%occ_list0)
        end if
        ! [WARNING - TODO] - ref%fock_sum not initialised in init_reference, etc! 
        qmc_state%ref%fock_sum = sum_fock_values_occ_list(sys, qmc_state%propagator%sp_fock, qmc_state%ref%occ_list0)

        ! --- Allocate psip list ---
        if (doing_calc(hfs_fciqmc_calc)) then
            qmc_state%psip_list%nspaces = qmc_state%psip_list%nspaces + 1
        else if (dmqmc_in_loc%replica_tricks) then
            qmc_state%psip_list%nspaces = qmc_state%psip_list%nspaces + 1
        else if (fciqmc_in_loc%replica_tricks) then
            qmc_state%psip_list%nspaces = qmc_state%psip_list%nspaces * 2
        end if
        if (sys%read_in%comp) then
            qmc_state%psip_list%nspaces = qmc_state%psip_list%nspaces * 2
        end if
        ! Each determinant occupies tot_string_len kind=i0 integers,
        ! qmc_state%psip_list%nspaces kind=int_p integers, qmc_state%psip_list%nspaces kind=p reals and one
        ! integer. If the Neel singlet state is used as the reference state for
        ! the projected estimator, then a further 2 reals are used per
        ! determinant.
        if (fciqmc_in_loc%trial_function == neel_singlet) qmc_state%psip_list%info_size = 2

        ! Allocate main particle lists.  Include the memory used by semi_stoch_t%determ in the
        ! calculation of memory occupied by the main particle lists.
        if (.not. (qmc_state_restart_loc)) then
            call init_particle_t(qmc_in%walker_length, 1, sys%basis%tensor_label_len, qmc_in%real_amplitudes, &
                                 qmc_in%real_amplitude_force_32, qmc_state%psip_list, io_unit=io_unit)
        end if

        call get_comm_processor_indx(qmc_state%psip_list%nspaces, proc_data_info, ntot_proc_data)
        call init_parallel_t(ntot_proc_data, est_buf_data_size, fciqmc_in_loc%non_blocking_comm, &
                             qmc_state%par_info, load_bal_in%nslots)

        uuid_restart = ''
        if (present(qmc_state_restart)) then
            call move_qmc_state_t(qmc_state_restart, qmc_state)
            qmc_state%mc_cycles_done = qmc_state_restart%mc_cycles_done
            uuid_restart = GLOBAL_META%uuid
            ! Prevent pattempt_single to be overwritten by default value in init_excit_gen.
            ! It is only overwritten if user specifies pattempt_single by qmc_in.
            qmc_state%excit_gen_data%p_single_double%pattempt_restart_store = .true.
        else
            call init_spawn_store(qmc_in, qmc_state%psip_list%nspaces, qmc_state%psip_list%pop_real_factor, sys%basis, &
                                  fciqmc_in_loc%non_blocking_comm, qmc_state%par_info%load%proc_map, io_unit, &
                                  qmc_state%spawn_store)
            ! Allocate the shift.
            allocate(qmc_state%shift(qmc_state%psip_list%nspaces), stat=ierr)
            call check_allocate('qmc_state%shift', qmc_state%psip_list%nspaces, ierr)
            allocate(qmc_state%vary_shift(qmc_state%psip_list%nspaces), stat=ierr)
            call check_allocate('qmc_state%vary_shift', qmc_state%psip_list%nspaces, ierr)
            qmc_state%shift = qmc_in%initial_shift
            qmc_state%vary_shift = .false.

            allocate(qmc_state%estimators(qmc_state%psip_list%nspaces))

            ! Initial walker distributions
            if (restart_in%read_restart) then
                call read_restart_hdf5(ri, sys%basis%nbasis, fciqmc_in_loc%non_blocking_comm, sys%basis%info_string_len, &
                                        qmc_state, uuid_restart, regenerate_info_loc, restart_version_restart)
            else if (doing_calc(dmqmc_calc)) then
                ! Initial distribution handled later
                qmc_state%psip_list%nstates = 0
            else
                call initial_distribution(sys, qmc_state%spawn_store%spawn, qmc_in%D0_population, fciqmc_in_loc, &
                                          qmc_state%ref, qmc_state%psip_list)
            end if
        end if
        if (present(regenerate_info)) regenerate_info = regenerate_info_loc

        call init_annihilation_flags(qmc_in, fciqmc_in_loc, dmqmc_in_loc, annihilation_flags)
        call init_trial(sys, fciqmc_in_loc, qmc_state%trial)
        call init_estimators(sys, qmc_in, restart_in%read_restart, qmc_state_restart_loc, fciqmc_in_loc%non_blocking_comm, &
                        qmc_state)
        if (present(qmc_state_restart)) call dealloc_excit_gen_data_t(qmc_state_restart%excit_gen_data)

        if (parent) write(iunit, '(1X, "# Starting the excitation generator initialisation.")')
        call cpu_time(t1)
        call init_excit_gen(sys, qmc_in, qmc_state%ref, qmc_state%vary_shift(1), qmc_state%excit_gen_data)
        call cpu_time(t2)
        set_up_time = t2 - t1
        if (parent) write(iunit, &
            '(1X, "# Finishing the excitation generator initialisation, time taken:",1X,es17.10)') set_up_time

        qmc_state%propagator%quasi_newton = qmc_in%quasi_newton
        if (qmc_in%quasi_newton) then
            if (qmc_in%quasi_newton_threshold < 0.0_p) then ! Not set by user, use auto value.
                ! Assume that fock values are ordered and that the number of basis functions is bigger than the number
                ! of electrons!
                if (parent) then
                    write (error_unit,'(1X,"# Warning in init_qmc: Doing quasi_newton with quasi_newton_threshold not supplied. &
                        &It is now estimated using (a multiple of) the difference in sp_fock energies of the basis functions at &
                        &indices (if CAS specified, then after freezing)",'//int_fmt(sys%nel,1)//', " &
                        &and",'//int_fmt(sys%nel+1,1)//', ". If these are not HOMO and LUMO, specify quasi_newton_threshold &
                        &directly.", /)') sys%nel, sys%nel+1
                end if
                qmc_state%propagator%quasi_newton_threshold = &
                    qmc_state%propagator%sp_fock(sys%nel+1) - qmc_state%propagator%sp_fock(sys%nel)
                if (sys%system == ueg) then
                    ! Know that by symmetry, the sum of Fock values of ref det to next excited det is twice the HOMO LUMO gap.
                    ! Ignore symmetry for the other systems for now...
                    qmc_state%propagator%quasi_newton_threshold = 2.0_p*qmc_state%propagator%quasi_newton_threshold
                end if
            else
                qmc_state%propagator%quasi_newton_threshold = qmc_in%quasi_newton_threshold
            end if
        end if
        if (qmc_in%quasi_newton_value < 0.0_p) then
            ! Default, set equal to quasi_newton_threshold.
            qmc_state%propagator%quasi_newton_value = qmc_state%propagator%quasi_newton_threshold
        else
            qmc_state%propagator%quasi_newton_value = qmc_in%quasi_newton_value
        end if
        if (qmc_state%propagator%quasi_newton) then
            if (qmc_in%quasi_newton_pop_control < 0.0_p) then
                qmc_state%propagator%quasi_newton_pop_control = 1.0_p/qmc_state%propagator%quasi_newton_threshold
            else
                qmc_state%propagator%quasi_newton_pop_control = qmc_in%quasi_newton_pop_control
            end if
        else
            ! Set to 1 if not using QN!
            qmc_state%propagator%quasi_newton_pop_control = 1.0_p
        end if
        ! Need to ensure we end up with a sensible value of shift damping to use.
        ! qmc_state%shift_damping will be set to either its default value or one
        ! read in from a restart file.
        if (qmc_in%shift_damping < huge(1.0_p)) then
            ! If we've passed in a specific value to qmc_in use that.
            qmc_state%shift_damping = qmc_in%shift_damping
        end if

        qmc_state%restart_in = restart_in
    end subroutine init_qmc

    subroutine init_proc_pointers(sys, qmc_in, reference, io_unit, dmqmc_in, fciqmc_in)

        ! Set function pointers for QMC calculations.

        ! In:
        !    sys: system being studied.
        !    qmc_in: input options relating to QMC methods.
        !    reference: reference_t object defining the reference state/determinant.
        !    io_unit: io unit to write all output to.
        !    dmqmc_in (optional): input options relating to DMQMC, only required if doing DMQMC.
        !    fciqmc_in (optional): input options relating to FCIQMC, only required if doing FCIQMC.

        ! System and calculation data
        use calc, only: doing_calc, doing_dmqmc_calc, dmqmc_calc, hfs_fciqmc_calc, &
                        ras, dmqmc_correlation, dmqmc_energy, dmqmc_energy_squared, dmqmc_staggered_magnetisation, &
                        dmqmc_kinetic_energy, dmqmc_H0_energy, dmqmc_potential_energy
        use hfs_data
        use system
        use parallel, only: parent
        use qmc_data, only: qmc_in_t, fciqmc_in_t, single_basis, neel_singlet, neel_singlet_guiding, &
                            excit_gen_renorm, excit_gen_renorm_spin, excit_gen_no_renorm, excit_gen_no_renorm_spin, &
                            excit_gen_power_pitzer_occ, excit_gen_power_pitzer_occ_ij, excit_gen_power_pitzer, &
                            excit_gen_power_pitzer_orderN, excit_gen_cauchy_schwarz_occ, excit_gen_cauchy_schwarz_occ_ij, &
                            excit_gen_heat_bath, excit_gen_heat_bath_uniform, excit_gen_heat_bath_single
        use dmqmc_data, only: dmqmc_in_t, free_electron_dm
        use reference_determinant, only: reference_t

        ! Procedures to be pointed to.
        use death, only: stochastic_death
        use determinants
        use determinant_decoders
        use dmqmc_estimators
        use dmqmc_procedures
        use dmqmc_data, only: hartree_fock_dm
        use energy_evaluation
        use excit_gen_mol
        use excit_gen_power_pitzer_mol, only: gen_excit_mol_power_pitzer_occ
        use excit_gen_power_pitzer_mol, only: gen_excit_mol_power_pitzer_occ_ref, gen_excit_mol_power_pitzer_orderN
        use excit_gen_heat_bath_mol, only: gen_excit_mol_heat_bath, gen_excit_mol_heat_bath_uniform
        use excit_gen_ueg, only:  gen_excit_ueg_power_pitzer
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
        use hamiltonian_molecular, only: slater_condon0_mol, double_counting_correction_mol, hf_hamiltonian_energy_mol, &
                                         slater_condon1_mol_excit, slater_condon2_mol_excit, get_one_e_int_mol, get_two_e_int_mol, & 
                                         create_weighted_excitation_list_mol, abs_hmatel_mol, single_excitation_weight_mol, &
                                         get_two_body_int_cou_mol_real, get_two_body_int_ex_mol_real 
        use hamiltonian_periodic_complex, only: slater_condon0_periodic_complex, slater_condon1_periodic_excit_complex, &
                                                slater_condon2_periodic_excit_complex, &
                                                create_weighted_excitation_list_periodic_complex, abs_hmatel_periodic_complex, &
                                                single_excitation_weight_periodic
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
        integer, intent(in) :: io_unit
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
            ! Quasi-newton requires cdet%fock_sum to be set and this is not provided
            ! by the split excitation generators.
            if (.not.qmc_in%quasi_newton) spawner_ptr => spawn_lattice_split_gen
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
            if (sys%read_in%comp) then
                update_proj_energy_ptr => update_proj_energy_periodic_complex
                sc0_ptr => slater_condon0_periodic_complex
                energy_diff_ptr => null()
                spawner_ptr => spawn_complex
                slater_condon1_excit_ptr => slater_condon1_periodic_excit_complex
                slater_condon2_excit_ptr => slater_condon2_periodic_excit_complex
                create_weighted_excitation_list_ptr => create_weighted_excitation_list_periodic_complex
                abs_hmatel_ptr => abs_hmatel_periodic_complex
                single_excitation_weight_ptr => single_excitation_weight_periodic
            else
                update_proj_energy_ptr => update_proj_energy_mol
                sc0_ptr => slater_condon0_mol
                energy_diff_ptr => double_counting_correction_mol
                slater_condon1_excit_ptr => slater_condon1_mol_excit
                slater_condon2_excit_ptr => slater_condon2_mol_excit
                create_weighted_excitation_list_ptr => create_weighted_excitation_list_mol
                ! Default: we use Power Pitzer type integrals when using on-the-fly weights for doubles.
                get_two_body_int_cou_ex_mol_real_ptr => get_two_body_int_ex_mol_real
                abs_hmatel_ptr => abs_hmatel_mol
                single_excitation_weight_ptr => single_excitation_weight_mol
            end if

            select case(qmc_in%excit_gen)
            case(excit_gen_no_renorm)
                gen_excit_ptr%full => gen_excit_mol_no_renorm
                decoder_ptr => decode_det_occ
            case(excit_gen_no_renorm_spin)
                gen_excit_ptr%full => gen_excit_mol_no_renorm_spin
                decoder_ptr => decode_det_spinocc
            case(excit_gen_renorm)
                gen_excit_ptr%full => gen_excit_mol
                decoder_ptr => decode_det_occ_symunocc
            case(excit_gen_renorm_spin)
                gen_excit_ptr%full => gen_excit_mol_spin
                decoder_ptr => decode_det_spinocc_symunocc
            case(excit_gen_power_pitzer_occ)
                gen_excit_ptr%full => gen_excit_mol_power_pitzer_occ
                decoder_ptr => decode_det_spinocc_spinsymunocc
            case(excit_gen_power_pitzer_occ_ij)
                gen_excit_ptr%full => gen_excit_mol_power_pitzer_occ
                decoder_ptr => decode_det_spinocc_spinsymunocc
            case(excit_gen_cauchy_schwarz_occ)
                gen_excit_ptr%full => gen_excit_mol_power_pitzer_occ
                decoder_ptr => decode_det_spinocc_spinsymunocc
                get_two_body_int_cou_ex_mol_real_ptr => get_two_body_int_cou_mol_real
            case(excit_gen_cauchy_schwarz_occ_ij)
                gen_excit_ptr%full => gen_excit_mol_power_pitzer_occ
                decoder_ptr => decode_det_spinocc_spinsymunocc
                get_two_body_int_cou_ex_mol_real_ptr => get_two_body_int_cou_mol_real
            case(excit_gen_power_pitzer)
                gen_excit_ptr%full => gen_excit_mol_power_pitzer_occ_ref
                decoder_ptr => decode_det_occ
            case(excit_gen_power_pitzer_orderN)
                ! [todo] - check this decoder is correct.
                gen_excit_ptr%full => gen_excit_mol_power_pitzer_orderN
                decoder_ptr => decode_det_occ
            case(excit_gen_heat_bath)
                gen_excit_ptr%full => gen_excit_mol_heat_bath
                decoder_ptr => decode_det_occ
            case(excit_gen_heat_bath_uniform)
                gen_excit_ptr%full => gen_excit_mol_heat_bath_uniform
                decoder_ptr => decode_det_occ_symunocc
            case(excit_gen_heat_bath_single)
                gen_excit_ptr%full => gen_excit_mol_heat_bath_uniform
                ! [todo] - the unocc part is only needed for singles. Too expensive here?
                decoder_ptr => decode_det_occ_unocc
            case default
                call stop_all('init_proc_pointers', 'Selected excitation generator not implemented.')
            end select

            get_one_e_int_ptr => get_one_e_int_mol
            get_two_e_int_ptr => get_two_e_int_mol

        case(ueg)

            update_proj_energy_ptr => update_proj_energy_ueg
            sc0_ptr => slater_condon0_ueg

            gen_excit_ptr%full => gen_excit_ueg_no_renorm
            decoder_ptr => decode_det_occ
            select case(qmc_in%excit_gen)
            case(excit_gen_no_renorm)
            case(excit_gen_power_pitzer)
                gen_excit_ptr%full => gen_excit_ueg_power_pitzer
            case(excit_gen_renorm)
                if (parent) then
                    write (io_unit,'(1X,"WARNING: renormalised excitation generators not implemented.")')
                    write (io_unit,'(1X,"WARNING: If this upsets you, please send patches.",/)')
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
                    write (io_unit,'(1X,"WARNING: renormalised excitation generators not implemented.")')
                    write (io_unit,'(1X,"WARNING: If this upsets you, please send patches.",/)')
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

            ! Reweighting routines for different initial density matrices.
            if (dmqmc_in%ipdmqmc .and. dmqmc_in%grand_canonical_initialisation) then
                select case(sys%system)
                case(ueg)
                    energy_diff_ptr => exchange_energy_ueg
                case(read_in)
                    energy_diff_ptr => double_counting_correction_mol
                case default
                    ! Please implement.
                end select
            end if

            ! Weighted importance sampling routines.
            if (dmqmc_in%weighted_sampling) then
                spawner_ptr => spawn_importance_sampling
                gen_excit_ptr%trial_fn => dmqmc_weighting_fn
            end if

            ! Expectation values.
            if (doing_dmqmc_calc(dmqmc_energy)) then
                if (dmqmc_in%ipdmqmc) then
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
                if (dmqmc_in%ipdmqmc) then
                    if (dmqmc_in%initial_matrix == free_electron_dm) then
                        h0_ptr => kinetic_energy_ueg
                        if (dmqmc_in%symmetric) then
                            if (dmqmc_in%weighted_sampling) then
                                gen_excit_ptr%trial_fn => dmqmc_int_pic_free_importance_sampling
                            else
                                spawner_ptr => spawn_importance_sampling
                                gen_excit_ptr%trial_fn => interaction_picture_reweighting_free
                            endif
                        end if
                    else
                        h0_ptr => slater_condon0_ueg
                        if (dmqmc_in%symmetric) then
                            if (dmqmc_in%weighted_sampling) then
                                gen_excit_ptr%trial_fn => dmqmc_int_pic_hf_importance_sampling
                            else
                                spawner_ptr => spawn_importance_sampling
                                gen_excit_ptr%trial_fn => interaction_picture_reweighting_hartree_fock
                            endif
                        end if
                    end if
                end if
                    if (doing_dmqmc_calc(dmqmc_kinetic_energy)) then
                        kinetic_diag_ptr => kinetic_energy_ueg
                        update_dmqmc_kinetic_energy_ptr => dmqmc_kinetic_energy_diag
                    end if
                    if (doing_dmqmc_calc(dmqmc_potential_energy)) then
                        potential_energy_ptr => potential_energy_ueg
                    end if
                    if (doing_dmqmc_calc(dmqmc_H0_energy) .and. .not. dmqmc_in%ipdmqmc) then
                        ! Assume that we want to evaluate <H_HF> rather than the kinetic energy.
                        h0_ptr => slater_condon0_ueg
                    end if
            case(hub_k)
                if (dmqmc_in%ipdmqmc) then
                    if (dmqmc_in%initial_matrix == free_electron_dm) then
                        h0_ptr => kinetic0_hub_k
                    else
                        h0_ptr => slater_condon0_hub_k
                    end if
                end if
            case(read_in)
                if (dmqmc_in%ipdmqmc) then
                    if (dmqmc_in%initial_matrix == free_electron_dm) then
                        h0_ptr => hf_hamiltonian_energy_mol
                    else
                        h0_ptr => slater_condon0_mol
                    end if
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

    subroutine init_trial(sys, fciqmc_in, trial)

        ! Initialise trial_t object with trial wavefunction details

        ! In:
        !   sys: system being studied
        !   fciqmc_in: input options relating to fciqmc
        ! In/Out:
        !   trial: trial wavefunction for importance sampling

        use checking, only: check_allocate
        use system, only: sys_t
        use qmc_data, only: fciqmc_in_t, trial_t, neel_singlet
        use utils, only: factorial_combination_1

        type(sys_t), intent(in) :: sys
        type(fciqmc_in_t), intent(in) :: fciqmc_in
        type(trial_t), intent(out) :: trial

        integer :: i, ierr

        trial%wfn = fciqmc_in%trial_function
        trial%guide = fciqmc_in%guiding_function

        if (trial%wfn == neel_singlet) then
            ! Calculate all the possible different amplitudes for the Neel singlet state
            ! and store them in an array
            allocate(trial%wfn_dat(-1:(sys%lattice%nsites/2)+1), stat=ierr)
            call check_allocate('qmc_state%trial%wfn_dat',(sys%lattice%nsites/2)+1,ierr)

            trial%wfn_dat(-1) = 0
            trial%wfn_dat((sys%lattice%nsites/2)+1) = 0
            do i=0,(sys%lattice%nsites/2)
                trial%wfn_dat(i) = factorial_combination_1((sys%lattice%nsites/2)-i, i)
                trial%wfn_dat(i) = -(2*mod(i,2)-1) * trial%wfn_dat(i)
            end do
        end if

    end subroutine init_trial

    subroutine init_annihilation_flags(qmc_in, fciqmc_in, dmqmc_in, annihilation_flags)

        ! Initialise annihilation flags

        ! In:
        !   qmc_in: input options relating to qmc
        !   fciqmc_in: input options relating to fciqmc
        !   dmqmc_in: input options relating to dmqmc
        ! Out:
        !   annihilation_flags: calculation specific annihilation flags.

        use qmc_data, only: qmc_in_t, fciqmc_in_t, annihilation_flags_t
        use dmqmc_data, only: dmqmc_in_t

        type(qmc_in_t), intent(in) :: qmc_in
        type(fciqmc_in_t), intent(in) :: fciqmc_in
        type(dmqmc_in_t), intent(in) :: dmqmc_in
        type(annihilation_flags_t), intent(out) :: annihilation_flags

        annihilation_flags%initiator_approx = qmc_in%initiator_approx
        annihilation_flags%real_amplitudes = qmc_in%real_amplitudes
        annihilation_flags%trial_function = fciqmc_in%trial_function
        annihilation_flags%ipdmqmc = dmqmc_in%ipdmqmc
        annihilation_flags%replica_tricks = dmqmc_in%replica_tricks
        annihilation_flags%symmetric = dmqmc_in%symmetric

    end subroutine init_annihilation_flags

    subroutine init_estimators(sys, qmc_in, restart_read_in, qmc_state_restart, fciqmc_non_blocking_comm, qmc_state)

        ! Initialise estimators and related components of qmc_state

        ! In:
        !   sys: system being studied
        !   qmc_in: input options relating to qmc methods
        !   restart_read_in : If true, have restarted and have read in restart file
        !   fciqmc_non_blocking_comm : If true, use non blocking communication and fciqmc
        ! In/Out:
        !   qmc_state: current state of qmc calculation

        use system, only: sys_t
        use qmc_data, only: qmc_in_t, qmc_state_t
        use calc, only: doing_calc, hfs_fciqmc_calc
        use parallel
        use energy_evaluation, only: calculate_hf_signed_pop
        use checking, only: check_allocate

        type(sys_t), intent(in) :: sys
        type(qmc_in_t), intent(in) :: qmc_in
        logical, intent(in) :: restart_read_in, qmc_state_restart, fciqmc_non_blocking_comm
        type(qmc_state_t), intent(inout) :: qmc_state

        logical :: have_tot_nparticles
        integer :: i
#ifdef PARALLEL
        real(dp) :: tmp_dp
        integer :: ierr
#endif
        have_tot_nparticles = .false.
        if ((restart_read_in) .and. (fciqmc_non_blocking_comm)) then
            have_tot_nparticles = .true.
        end if

        associate(pl=>qmc_state%psip_list)
            ! Total number of particles on processor.
            ! Probably should be handled more simply by setting it to be either 0 or
            ! D0_population or obtaining it from the restart file, as appropriate.
            forall (i=1:pl%nspaces) pl%nparticles(i) = sum(abs( real(pl%pops(i,:pl%nstates),p)/pl%pop_real_factor))
            ! Should we already be in varyshift mode (e.g. restarting a calculation)?
#ifdef PARALLEL
            do i=1, pl%nspaces
                call mpi_allgather(pl%nparticles(i), 1, MPI_REAL8, pl%nparticles_proc(i,:), 1, MPI_REAL8, MPI_COMM_WORLD, ierr)
            end do
            ! When restarting a non-blocking calculation this sum will not equal
            ! tot_nparticles as some walkers have been communicated around the report
            ! loop. The correct total is in the restart file so get it from there.
            if (.not. have_tot_nparticles) &
                forall(i=1:pl%nspaces) pl%tot_nparticles(i) = sum(pl%nparticles_proc(i,:))
#else
            pl%tot_nparticles = pl%nparticles
            pl%nparticles_proc(:pl%nspaces,1) = pl%nparticles(:pl%nspaces)
#endif
        end associate

        ! Decide whether the shift should be turned on from the start.
        qmc_state%target_particles = qmc_in%target_particles
        if (qmc_in%vary_shift_present) then
            ! User input overrides other factors determining whether to vary shift
            qmc_state%vary_shift = qmc_in%vary_shift
            if (.not. qmc_in%vary_shift) qmc_state%shift = qmc_in%initial_shift
        end if

        if (doing_calc(hfs_fciqmc_calc)) then
#ifdef PARALLEL
            tmp_dp = calculate_hf_signed_pop(qmc_state%psip_list)
            call mpi_allreduce(tmp_dp, qmc_state%estimators%hf_signed_pop, qmc_state%psip_list%nspaces, mpi_real8, MPI_SUM, &
                               MPI_COMM_WORLD, ierr)
#else
            qmc_state%estimators%hf_signed_pop = calculate_hf_signed_pop(qmc_state%psip_list)
#endif
        end if

        ! Set initial values from input
        qmc_state%tau = qmc_in%tau

    end subroutine init_estimators

    subroutine init_excit_gen(sys, qmc_in, ref, vary_shift_first_el, excit_gen_data)

        ! Initialise pre-computed data for excitation generators

        ! In:
        !   sys: system being studied
        !   qmc_in: input options for qmc
        !   vary_shift_first_el: qmc_state%vary_shift(1), is shift varying in first space?
        !   ref: reference determinant
        ! Out:
        !   excit_gen_data: pre-computed data for fast excitation generation.

        use const, only: p
        use system, only: sys_t, ueg, read_in
        use qmc_data, only: qmc_in_t, excit_gen_renorm_spin, excit_gen_no_renorm_spin, &
                            excit_gen_power_pitzer, excit_gen_power_pitzer_orderN, excit_gen_power_pitzer_occ_ij, &
                            excit_gen_heat_bath, excit_gen_heat_bath_uniform, excit_gen_heat_bath_single, &
                            excit_gen_cauchy_schwarz_occ, excit_gen_cauchy_schwarz_occ_ij
        use excit_gen_power_pitzer_mol, only: init_excit_mol_power_pitzer_occ_ref, init_excit_mol_power_pitzer_orderM_ij, &
                                              init_excit_mol_power_pitzer_orderN
        use excit_gen_heat_bath_mol, only: init_excit_mol_heat_bath
        use excit_gen_ueg, only: init_excit_ueg_power_pitzer
        use excit_gens, only: excit_gen_data_t, zero_p_single_double_coll_t
        use parallel, only: parent
        use ueg_system, only: init_ternary_conserve
        use qmc_common, only: find_single_double_prob, find_parallel_spin_prob_mol
        use reference_determinant, only: reference_t

        type(sys_t), intent(in) :: sys
        type(qmc_in_t), intent(in) :: qmc_in
        type(reference_t), intent(in) :: ref
        logical, intent(in) :: vary_shift_first_el
        type(excit_gen_data_t), intent(inout) :: excit_gen_data
        real :: tw1, tw2, set_up_timew
        integer :: iunit
        
        iunit = 6

        ! Set type of excitation generator to use
        excit_gen_data%excit_gen = qmc_in%excit_gen

        ! If not set at input or available from restart file, set probability of selecting single or double
        ! excitations based upon the reference determinant and assume other
        ! determinants have a roughly similar ratio of single:double
        ! excitations. Input overrides pattempt_single stored in restart file.
        if (qmc_in%pattempt_single < 0 .or. qmc_in%pattempt_double < 0) then
            if (.not. excit_gen_data%p_single_double%pattempt_restart_store) then
                call find_single_double_prob(sys, ref%occ_list0, excit_gen_data%pattempt_single, excit_gen_data%pattempt_double)
            else if (qmc_in%pattempt_zero_accum_data) then
                ! User wants accumulated data to be zeroed (e.g. restarting with a different excitation generator).
                ! pattempt_single and pattempt_double stay what they were.
                call zero_p_single_double_coll_t(excit_gen_data%p_single_double%total)
                excit_gen_data%p_single_double%counter = 1.0_p
            end if
        else
            ! [todo] - is it wise that the user has to specify both *_single and *_double to overwrite?
            ! zero data received from restart file if user overwrittes it by specifying a qmc_in%pattempt_single and *_double.
            call zero_p_single_double_coll_t(excit_gen_data%p_single_double%total)
            excit_gen_data%p_single_double%counter = 1.0_p
            ! renormalise just in case input wasn't
            excit_gen_data%pattempt_single = qmc_in%pattempt_single/(qmc_in%pattempt_single+qmc_in%pattempt_double)
            excit_gen_data%pattempt_double = 1.0_p - qmc_in%pattempt_single
        end if

        ! If not set at input, set probability of selecting parallel ij to the ratio of |Hij->ab| with parallel spins to
        ! total |Hij->ab|.
        ! [todo] - It might be desireable in the future to store pattempt_parallel in the restart file.
        ! [todo] - Currently, we have decided not to store it so a user can easily restart from another system without
        ! [todo] - having to specify that a new pattempt_parallel needs to be calculated. In a sample calculation with 120
        ! [todo] - orbitals it took 14 seconds to initialise (6 mpi procs on three nodes with 12 threads each on ARCHER
        ! [todo] - http://archer.ac.uk/).
        if ((qmc_in%pattempt_parallel < 0.0_p) .and. (sys%system == read_in) .and. &
            ((qmc_in%excit_gen == excit_gen_renorm_spin) .or. (qmc_in%excit_gen == excit_gen_no_renorm_spin))) then
            call find_parallel_spin_prob_mol(sys, excit_gen_data%pattempt_parallel)
        else if (.not.(qmc_in%pattempt_parallel < 0.0_p)) then
            ! pattempt_parallel set by user. check_input.f90 makes sure this is only set if read_in system and excit. gens.
            ! are no_renorm_spin or renorm_spin.
            excit_gen_data%pattempt_parallel = qmc_in%pattempt_parallel
        end if

        ! UEG allowed excitations
        if (sys%system == ueg) call init_ternary_conserve(sys, excit_gen_data%ueg_ternary_conserve)

        ! Init weighted excitation generators if required.
        if (qmc_in%excit_gen==excit_gen_power_pitzer) then
            if (parent) write(iunit, '(1X, "# Starting the Power Pitzer excitation generator initialisation.")')
            call cpu_time(tw1)
            if (sys%system == read_in) then
               excit_gen_data%excit_gen_pp%power_pitzer_min_weight = qmc_in%power_pitzer_min_weight
               call init_excit_mol_power_pitzer_occ_ref(sys, ref, excit_gen_data%excit_gen_pp)
            else if (sys%system == ueg) then 
               call init_excit_ueg_power_pitzer(sys, ref, excit_gen_data%excit_gen_pp)
            end if
            call cpu_time(tw2)
            set_up_timew = tw2 - tw1
            if (parent) write(iunit, &
                '(1X, "# Finishing the Power Pitzer excitation generator initialisation, time taken:",1X,es17.10)') set_up_timew
        end if

        if (qmc_in%excit_gen==excit_gen_power_pitzer_orderN) then
            if (parent) write(iunit, '(1X, "# Starting the P.P. Order N excitation generator initialisation.")')
            call cpu_time(tw1)
            excit_gen_data%excit_gen_pp%power_pitzer_min_weight = qmc_in%power_pitzer_min_weight
            call init_excit_mol_power_pitzer_orderN(sys, ref, excit_gen_data%excit_gen_pp)
            call cpu_time(tw2)
            set_up_timew = tw2 - tw1
            if (parent) write(iunit, &
                '(1X, "# Finishing the P.P. Order N excitation generator initialisation, time taken:",1X,es17.10)') set_up_timew
        end if
        
        if ((qmc_in%excit_gen==excit_gen_power_pitzer_occ_ij) .or. (qmc_in%excit_gen==excit_gen_cauchy_schwarz_occ_ij)) then
            if (parent) write(iunit, '(1X, "# Starting the P.P./C.S. O(M)ij excitation generator initialisation.")')
            call cpu_time(tw1)
            call init_excit_mol_power_pitzer_orderM_ij(sys, ref, excit_gen_data%excit_gen_pp)
            call cpu_time(tw2)
            set_up_timew = tw2 - tw1
            if (parent) write(iunit, &
                '(1X, "# Finishing P.P./C.S. O(M)ij excitation generator initialisation, time taken:",1X,es17.10)') set_up_timew
        end if

        if (qmc_in%excit_gen==excit_gen_heat_bath) then
            if (parent) write(iunit, '(1X, "# Starting the heat bath excitation generator initialisation.")')
            call cpu_time(tw1)
            call init_excit_mol_heat_bath(sys, excit_gen_data%excit_gen_hb, .true.)
            call cpu_time(tw2)
            set_up_timew = tw2 - tw1
            if (parent) write(iunit, &
                '(1X, "# Finishing the heat bath excitation generator initialisation, time taken:",1X,es17.10)') set_up_timew
        end if
        if ((qmc_in%excit_gen==excit_gen_heat_bath_uniform) .or. (qmc_in%excit_gen==excit_gen_heat_bath_single)) then
            if (parent) write(iunit, '(1X, "# Starting the heat bath excitation generator initialisation.")')
            call cpu_time(tw1)
            call init_excit_mol_heat_bath(sys, excit_gen_data%excit_gen_hb, .false.)
            call cpu_time(tw2)
            set_up_timew = tw2 - tw1
            if (parent) write(iunit, &
                '(1X, "# Finishing the heat bath excitation generator initialisation, time taken:",1X,es17.10)') set_up_timew
        end if

        ! Set vary_psingles.
        ! check_input makes sure the pattempt_update cannot be true if we are using the heat bath excitation generator.
        if (.not. (vary_shift_first_el) .and. (qmc_in%pattempt_update)) then
            ! We sample singles with probability pattempt_single. It therefore makes sense to update pattempt_single 
            ! [todo] - DMQMC
            ! for FCIQMC and CCMC on the fly (at least in the beginning of the calculation).
            excit_gen_data%p_single_double%vary_psingles = .true.
        else if ((excit_gen_data%p_single_double%vary_psingles) .and. .not.(qmc_in%pattempt_update)) then
            ! If we are restarting say and vary_psingles is true in the restart file but now pattempt_update is false,
            ! set vary_psingles to false as well.
            excit_gen_data%p_single_double%vary_psingles = .false.
        end if

    end subroutine init_excit_gen

! [review] - AJWT: document occlist.  Does any call to this function actually use it?
    subroutine init_reference(sys, reference_in, io_unit, reference, occlist)

        ! Set the reference determinant from input options

        ! In:
        !   sys: system being studied.
        !   reference_in: reference provided in input, if any.
        !   io_unit: io unit to write any information to.
        ! Out:
        !   reference: reference selected for the qmc calculation.

        use system, only: sys_t, ueg, read_in
        use proc_pointers, only: sc0_ptr, op0_ptr
        use calc, only: doing_calc, hfs_fciqmc_calc
        use reference_determinant, only: reference_t, set_reference_det
        use checking, only: check_allocate
        use determinants, only: encode_det

        type(sys_t), intent(in) :: sys
        type(reference_t), intent(in) :: reference_in
        integer, intent(in) :: io_unit
        type(reference_t), intent(out) :: reference
        integer, allocatable, optional :: occlist(:)

        integer :: ierr
        ! Note occ_list could be set and allocated in the input.
        reference = reference_in

    if (present(occlist)) reference%occ_list0 = occlist

        ! Set the reference determinant to be the spin-orbitals with the lowest
        ! single-particle eigenvalues which satisfy the spin polarisation and, if
        ! specified, the symmetry.

        call set_reference_det(sys, reference%occ_list0, .false., sys%symmetry, io_unit)

        if (.not. allocated(reference%f0)) then
            allocate(reference%f0(sys%basis%tot_string_len), stat=ierr)
            call check_allocate('reference%f0',sys%basis%tot_string_len,ierr)
        end if
        call encode_det(sys%basis, reference%occ_list0, reference%f0)

        if (.not. allocated(reference%hs_f0)) then
            allocate(reference%hs_f0(sys%basis%tot_string_len), stat=ierr)
            call check_allocate('reference%hs_f0', sys%basis%tot_string_len, ierr)
        end if

        ! Set hilbert space reference if not given in input
        if (allocated(reference%hs_occ_list0)) then
            call encode_det(sys%basis, reference%hs_occ_list0, reference%hs_f0)
        else
            allocate(reference%hs_occ_list0(sys%nel), stat=ierr)
            call check_allocate('reference%hs_occ_list0', sys%nel, ierr)
            reference%hs_occ_list0 = reference%occ_list0
            reference%hs_f0 = reference%f0
        end if

        ! Energy of reference determinant.
        reference%H00 = sc0_ptr(sys, reference%f0)
        ! Operators of HFS sampling.
        if (doing_calc(hfs_fciqmc_calc)) reference%O00 = op0_ptr(sys, reference%f0)
        ! [WARNING - TODO] - ref%fock_sum not initialised here! 

    end subroutine init_reference

    subroutine init_reference_restart(sys, reference_in, ri, reference)

        ! Set the reference determinant from a HDF5 restart file.

        ! In:
        !   sys: system being studied.
        !   reference_in: reference provided in input.  Only ex_level is used.
        !   ri: restart information.
        ! Out:
        !   reference: reference selected for the qmc calculation.

        use reference_determinant, only: reference_t
        use system, only: sys_t, ueg
        use restart_hdf5, only: restart_info_t, get_reference_hdf5
        use calc, only: doing_calc, hfs_fciqmc_calc
        use proc_pointers, only: sc0_ptr, op0_ptr
        use checking, only: check_allocate
        use determinants, only: decode_det

        type(sys_t), intent(in) :: sys
        type(reference_t), intent(in) :: reference_in
        type(restart_info_t), intent(in) :: ri
        type(reference_t), intent(out) :: reference

        integer :: ierr

        allocate(reference%f0(sys%basis%tot_string_len), stat=ierr)
        call check_allocate('reference%f0',sys%basis%tot_string_len,ierr)
        allocate(reference%hs_f0(sys%basis%tot_string_len), stat=ierr)
        call check_allocate('reference%hs_f0', sys%basis%tot_string_len, ierr)
        allocate(reference%occ_list0(sys%nel), stat=ierr)
        call check_allocate('reference%occ_list0',sys%nel,ierr)
        allocate(reference%hs_occ_list0(sys%nel), stat=ierr)
        call check_allocate('reference%hs_occ_list0',sys%nel,ierr)
        call get_reference_hdf5(ri, sys%basis%nbasis, sys%basis%info_string_len, reference)

        ! Need to re-calculate the reference determinant data
        call decode_det(sys%basis, reference%f0, reference%occ_list0)
        call decode_det(sys%basis, reference%hs_f0, reference%hs_occ_list0)

        reference%H00 = sc0_ptr(sys, reference%f0)
        if (doing_calc(hfs_fciqmc_calc)) reference%O00 = op0_ptr(sys, reference%f0)

        reference%ex_level = reference_in%ex_level
        ! [WARNING - TODO] - ref%fock_sum not initialised here! 

    end subroutine init_reference_restart

    subroutine init_secondary_reference(sys,reference_in,io_unit,qs)
        ! Set the secondary reference determinant from input options
        ! and use it to set up the maximum considered excitation level
        ! for the calculation.

        ! In:
        !   sys: system being studied.
        !   reference_in: secondary reference provided in input.
        !   io_unit: io unit to write any information to.
        ! In/Out:
        !   qs: qmc_state used in the calculation.
    
        use reference_determinant, only: reference_t
        use system, only: sys_t
        use qmc_data, only: qmc_state_t 
        use excitations, only: get_excitation_level
        type(sys_t), intent(in) :: sys
        type(reference_t), intent(in) :: reference_in
        integer, intent(in) :: io_unit
        type(qmc_state_t), intent(inout) :: qs

        call init_reference(sys, reference_in, io_unit, qs%second_ref)
        qs%ref%max_ex_level = qs%ref%ex_level + get_excitation_level(qs%ref%f0(:sys%basis%bit_string_len), &
                                                                  qs%second_ref%f0(:sys%basis%bit_string_len))
        ! [WARNING - TODO] - ref%fock_sum not initialised here! 
    end subroutine


    subroutine init_spawn_store(qmc_in, nspaces, pop_real_factor, basis, non_blocking_comm, proc_map, io_unit, spawn_store)

        ! Allocate and initialise spawn store

        ! In:
        !   qmc_in: input options relating to qmc methods
        !   nspaces: number of particle types in use.
        !   pop_real_factor: scaling factor to encode real populations.
        !   basis: information on the single particle basis used.
        !   non_blocking_comm: whether non-blocking MPI communication is being used.
        !   proc_map: settings for mapping hash of a determinant to a processor.
        !   io_unit: io unit to write all output to.
        ! Out:
        !   spawn_store: spawned array.  All components allocated/initialised on exit.

        use qmc_data, only: qmc_in_t, spawned_particle_t
        use const
        use parallel, only: parent, nprocs
        use calc, only: doing_calc, dmqmc_calc
        use basis_types, only: basis_t
        use spawn_data, only: proc_map_t, alloc_spawn_t

        type(qmc_in_t), intent(in) :: qmc_in
        integer, intent(in) :: nspaces, io_unit
        integer(int_p), intent(in) :: pop_real_factor
        type(basis_t), intent(in) :: basis
        logical, intent(in) :: non_blocking_comm
        type(proc_map_t), intent(in) :: proc_map
        type(spawned_particle_t), intent(out) :: spawn_store

        integer :: size_spawned_walker, max_nspawned_states, nhash_bits
        real(p) :: spawn_cutoff

        ! Calculate length of spawned particle arrays

        ! Each spawned_walker occupies spawned_size kind=int_s integers.
        if (qmc_in%initiator_approx) then
            size_spawned_walker = (basis%tensor_label_len+nspaces+1)*int_s_length/8
        else
            size_spawned_walker = (basis%tensor_label_len+nspaces)*int_s_length/8
        end if

        max_nspawned_states = qmc_in%spawned_walker_length
        if (max_nspawned_states < 0) then
            ! Given in MB.  Convert.
            ! Note that we store 2 arrays.
            max_nspawned_states = int((-real(max_nspawned_states,p)*10**6)/(2*size_spawned_walker))
        end if
        if (parent) then
            write (io_unit,'(1X,a57,f9.2)') &
                'Memory allocated per core for spawned walker lists (MB): ', &
                size_spawned_walker*real(2*max_nspawned_states,p)/10**6
            write (io_unit,'(1X,a51,1x,i0,/)') &
                'Number of elements per core in spawned walker list:', max_nspawned_states
        end if
        if (mod(max_nspawned_states, nprocs) /= 0) then
            max_nspawned_states = ceiling(real(max_nspawned_states)/nprocs)*nprocs
            if (parent) then
               write (io_unit,'(1X,a68)') 'spawned_walker_length is not a multiple of the number of processors.'
               write (io_unit,'(1X,a35,1x,i0,a1,/)') &
                'Increasing spawned_walker_length to',max_nspawned_states,'.'
            end if
        end if

        ! If not using real amplitudes then we always want spawn_cutoff to be
        ! equal to 1.0, so overwrite the default before creating spawn_t objects.
        spawn_cutoff = qmc_in%spawn_cutoff
        if (.not. qmc_in%real_amplitudes) spawn_cutoff = 0.0_p

        if (doing_calc(dmqmc_calc)) then
            ! Hash the entire first bit array and the minimum number of bits
            ! in the second bit array.
            nhash_bits = basis%nbasis + i0_length*(basis%bit_string_len)
        else
            nhash_bits = basis%nbasis
        end if

        ! Allocate spawned particle lists.
        call alloc_spawn_t(basis%tensor_label_len, nhash_bits, nspaces, qmc_in%initiator_approx, &
                           max_nspawned_states, spawn_cutoff, pop_real_factor, proc_map, 7, &
                           qmc_in%use_mpi_barriers, spawn_store%spawn)
        if (non_blocking_comm) then
            call alloc_spawn_t(basis%tensor_label_len, nhash_bits, nspaces, qmc_in%initiator_approx, &
                               max_nspawned_states, spawn_cutoff, pop_real_factor, proc_map, 7, &
                               .false., spawn_store%spawn_recv)
        end if

    end subroutine init_spawn_store

    subroutine initial_distribution(sys, spawn, D0_pop, fciqmc_in, reference, pl)

        ! Set the initial psip distribution for a CCMC/FCIQMC calculation.

        ! In:
        !   sys: system being studied.
        !   spawn: spawned array.
        !   D0_pop: initial population on the reference determinant(s).
        !   fciqmc_in: input options relating to fciqmc.
        ! In/Out:
        !   reference: reference determinant.  H00 is reset on exit if using neel singlet trial function.
        !   pl: main particle list.  On exit, the initial distribution is set.

        use const, only: i0, i0_end, int_p, p
        use qmc_data, only: particle_t, fciqmc_in_t, neel_singlet, single_basis
        use parallel, only: iproc, nprocs
        use spawn_data, only: spawn_t
        use system, only: sys_t, heisenberg
        use proc_pointers, only: sc0_ptr
        use reference_determinant, only: reference_t
        use spawning, only: assign_particle_processor
        use determinants, only: encode_det

        type(sys_t), intent(in) :: sys
        type(spawn_t), intent(in) :: spawn
        type(fciqmc_in_t), intent(in) :: fciqmc_in
        real(p), intent(in) :: D0_pop
        type(reference_t), intent(inout) :: reference
        type(particle_t), intent(inout) :: pl

        integer :: D0_proc, slot, i, ipos, D0_inv_proc
        integer :: occ_list0_inv(sys%nel)
        integer(i0) :: f0_inv(sys%basis%tot_string_len)

        ! We start with psips only on the reference determinant,
        ! so set pl%nstates = 1 and initialise pl%pops.

        pl%nstates = 1
        ! Set the bitstring of this psip to be that of the
        ! reference state.
        pl%states(:,pl%nstates) = reference%f0
        ! Zero population for all spaces.
        pl%pops(:,pl%nstates) = 0_int_p
        ! Set initial population of Hamiltonian walkers.
        if (sys%read_in%comp) then
            ! Only have population in real spaces
            pl%pops(1::2,pl%nstates) = nint(D0_pop)*pl%pop_real_factor
        else
            pl%pops(:,pl%nstates) = nint(D0_pop)*pl%pop_real_factor
        end if

        ! Determine and set properties for the reference state which we start on.
        ! By definition, when using a single determinant as a reference state:
        pl%dat(1,pl%nstates) = 0.0_p
        ! Or if not using a single determinant:
        if (fciqmc_in%trial_function == neel_singlet) then
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
        call assign_particle_processor(reference%f0, spawn%bit_str_nbits, spawn%hash_seed, spawn%hash_shift, &
                                       spawn%move_freq, nprocs, D0_proc, slot, spawn%proc_map%map, spawn%proc_map%nslots)
        if (D0_proc /= iproc) pl%nstates = 0

        ! For the Heisenberg model and open shell systems, it is often useful to
        ! have psips start on both the reference state and the spin-flipped version.
        if (fciqmc_in%init_spin_inv_D0) then

            ! Need to handle the Heisenberg model (consisting of spinors on
            ! lattice sites) and electron systems differently, as the
            ! Heisenberg model has no concept of unoccupied basis
            ! functions/holes.
            select case (sys%system)
            case (heisenberg)
                ! Flip all spins in f0 to get f0_inv
                f0_inv = not(reference%f0)
                ! Need to account for lengthening of bit string to store ex_lvl for ccmc.
                if (sys%basis%info_string_len /= 0) f0_inv(sys%basis%bit_string_len+1:) = 0_i0
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

            call assign_particle_processor(f0_inv, spawn%bit_str_nbits, spawn%hash_seed, spawn%hash_shift, &
                                           spawn%move_freq, nprocs, D0_inv_proc, slot, spawn%proc_map%map, spawn%proc_map%nslots)

            ! Store if not identical to reference det.
            if (D0_inv_proc == iproc .and. any(reference%f0 /= f0_inv)) then
                pl%nstates = pl%nstates + 1
                ! Zero all populations for this determinant.
                pl%pops(:,pl%nstates) = 0_int_p
                ! Set the population for this basis function.
                pl%pops(1,pl%nstates) = nint(D0_pop)*pl%pop_real_factor
                pl%dat(1,pl%nstates) = sc0_ptr(sys, reference%f0) - reference%H00
                select case(sys%system)
                case(heisenberg)
                    if (fciqmc_in%trial_function /= single_basis) then
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
                if (fciqmc_in%trial_function == neel_singlet) then
                    pl%dat(pl%nspaces+1,pl%nstates) = 0
                    pl%dat(pl%nspaces+2,pl%nstates) = 0
                end if
            end if
        end if

    end subroutine initial_distribution

    subroutine move_qmc_state_t(qmc_state_old, qmc_state_new)

        ! Move components of qmc_state_old to qmc_state_new

        ! In/Out:
        !   qmc_state_old: existing qmc_state.  Deallocated on exit.
        !   qmc_state_new: qmc_state_t.  Components allocated an initialised on exit.

        use dSFMT_interface, only: move_dSFMT_state_t
        use qmc_data, only: qmc_state_t
        use spawn_data, only: move_spawn_t
        use excit_gens, only: move_pattempt_data
        use particle_t_utils, only: move_particle_t
        use errors, only: stop_all

        type(qmc_state_t), intent(inout) :: qmc_state_old, qmc_state_new

        ! Assume if particle_t is correctly present in qmc_state_old, other components are as well. This is a first-order
        ! sanity check on whether qmc_state_old is actually allocated and valid (e.g. has not been already passed in to a previous
        ! calculation and deallocated).
        if (.not.allocated(qmc_state_old%psip_list%states)) call stop_all('move_qmc_state_t', &
                                                                          'Attempting to restart from an invalid qmc_state.')

        call move_spawn_t(qmc_state_old%spawn_store%spawn, qmc_state_new%spawn_store%spawn)
        call move_spawn_t(qmc_state_old%spawn_store%spawn_recv, qmc_state_new%spawn_store%spawn_recv)
        call move_alloc(qmc_state_old%shift, qmc_state_new%shift)
        call move_alloc(qmc_state_old%vary_shift, qmc_state_new%vary_shift)
        call move_alloc(qmc_state_old%estimators, qmc_state_new%estimators)
        call move_dSFMT_state_t(qmc_state_old%rng_state, qmc_state_new%rng_state)
        call move_particle_t(qmc_state_old%psip_list, qmc_state_new%psip_list)
        call move_pattempt_data(qmc_state_old%excit_gen_data, qmc_state_new%excit_gen_data)

    end subroutine move_qmc_state_t

end module qmc

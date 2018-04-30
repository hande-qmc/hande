module check_input
! Parse input options and check input for validity.

implicit none

contains

    subroutine check_sys(sys)

        ! Check the input parameters for the system.

        ! I don't pretend this is the most comprehensive of tests, but at least
        ! make sure a few things are not completely insane.

        ! In:
        !    sys: system object, as set in read_input.

        use parallel, only: default_max_broadcast_chunk
        use system

        use errors, only: stop_all, warning

        type(sys_t), intent(in) :: sys

        integer :: ivec, jvec
        character(*), parameter :: this='check_sys'

        if (sys%system /= heisenberg) then
            if (sys%nel < 0) call stop_all(this,'Number of electrons must be positive.')
            if (sys%nel == 0 .and. sys%Ms /= huge(1)) call stop_all(this, &
                                        'Must specify both or neither of nel and Ms.')
            if (sys%nel /= 0 .and. sys%Ms == huge(1)) call stop_all(this, &
                                        'Must specify both or neither of nel and Ms.')
        end if

        if (sys%system /= ueg .and. sys%system /= read_in .and. sys%system /= ringium) then
            if (.not.(allocated(sys%lattice%lattice))) call stop_all(this, 'Lattice vectors not provided')
            do ivec = 1, sys%lattice%ndim
                do jvec = ivec+1, sys%lattice%ndim
                    if (dot_product(sys%lattice%lattice(:,ivec), sys%lattice%lattice(:,jvec)) /= 0) then
                        call stop_all(this, 'Lattice vectors are not orthogonal.')
                    end if
                end do
            end do
        end if

        if (sys%system == heisenberg) then
            if (abs(sys%heisenberg%staggered_magnetic_field) > depsilon .and. (.not.sys%lattice%bipartite_lattice)) &
                call stop_all(this, 'Cannot set a staggered field for this lattice because it is frustrated.')
            if (abs(sys%heisenberg%staggered_magnetic_field) > depsilon .and. abs(sys%heisenberg%magnetic_field) > depsilon) &
                call stop_all(this, 'Cannot set a uniform and a staggered field at the same time.')
        else if (sys%system == hub_k .or. sys%system == hub_real) then
            if (sys%nel > 2*sys%lattice%nsites) call stop_all(this, 'More than two electrons per site.')
        end if

        if (sys%system /= read_in) then
            if (sys%lattice%ndim > 3) call stop_all(this, 'Limited to 1, 2 or 3 dimensions')
            if (sys%system == ueg .and. sys%lattice%ndim == 1) call stop_all(this, 'UEG only functional in 2D and 3D')
        end if

        if (sys%momentum_space .or. sys%system == read_in) then
            if (sys%real_lattice%finite_cluster) call stop_all(this,'"finite" keyword only valid for calculations in real-space.')
        end if

        if (sys%system == ringium .and. sys%ringium%radius < depsilon) call stop_all(this, 'Ringium must have a positive radius.')

        if (sys%system == read_in .and. sys%read_in%max_broadcast_chunk /= default_max_broadcast_chunk) &
            call warning(this, 'Not using default max_broadcast_chunk. This option is only intended for testing of &
                &broadcasting functionality and gives no benefit to functionality.')

        if (sys%read_in%comp .and. sys%read_in%uhf) call warning(this, &
                            'UHF translational symmetry currently untested. Use with caution!')

    end subroutine check_sys

    subroutine check_fciqmc_opts(sys, fciqmc_in, blocking_in)

        ! Check the FCIQMC specific options.

        ! In:
        !   sys: system being studied.
        !   fciqmc_in: FCIQMC input options.

        use qmc_data, only: fciqmc_in_t, single_basis, no_guiding, neel_singlet_guiding, neel_singlet, blocking_in_t
        use system, only: sys_t, heisenberg
        use errors, only: stop_all

        type(sys_t), intent(in) :: sys
        type(fciqmc_in_t), intent(in) :: fciqmc_in
        type(blocking_in_t), intent(in) :: blocking_in

        character(*), parameter :: this = 'check_fciqmc_opts'

        if (sys%system /= heisenberg) then
            if (fciqmc_in%trial_function /= single_basis) &
                call stop_all(this, 'Only a single determinant can be used as the reference state for this system. &
                                     &Other trial functions are not available.')
            if (fciqmc_in%guiding_function /= no_guiding) &
                call stop_all(this, 'Importance sampling is only avaliable for the Heisenberg model currently.')
         else
            if (fciqmc_in%guiding_function == neel_singlet_guiding .and. fciqmc_in%trial_function /= neel_singlet) &
                call stop_all(this, &
                    'This guiding function is only available when using the Neel singlet state as an energy estimator.')
        end if

        if (fciqmc_in%init_spin_inv_D0 .and. sys%Ms /= 0) then
            call stop_all(this, &
                'Flipping the reference state will give a state which has a different value of Ms and so cannot be used here.')
        end if

        if (fciqmc_in%non_blocking_comm .and. sys%read_in%comp) then
            call stop_all(this, &
                'Complex non-blocking comms not yet implemented.')
        end if

        ! This is mainly a precaution, as the majority of functionality is possible, but totally unused.
        ! Major unresolved issue is likely to be semistochastic restart files.
        if (sys%basis%info_string_len /= 0) call stop_all(this, &
            'Fciqmc is incompatible with additional information being stored in the bit string. Please implement if needed.')

        if (fciqmc_in%replica_tricks .and. blocking_in%blocking_on_the_fly) call stop_all(this, &
            'Blocking on the fly is not currently compatible with replica tricks in fciqmc.')

    end subroutine check_fciqmc_opts

    subroutine check_qmc_opts(qmc_in, sys, need_length, restarting, qmc_state_restart)

        ! Check options common to QMC methods.

        ! In:
        !   qmc_in: QMC input options.
        !   need_length: whether the size of main/spawned particle arrays is needed.
        !       false if using simple fciqmc algorithm or restarting from memory.
        !   restarting: whether the calculation is restarting from a previous one.
        !   qmc_state_restart (optional): qmc_state object being restarted from.

        use qmc_data, only: qmc_in_t, qmc_state_t
        use errors, only: stop_all
        use system, only: sys_t

        type(qmc_in_t), intent(in) :: qmc_in
        type(sys_t), intent(in) :: sys
        logical, intent(in) :: need_length, restarting
        type(qmc_state_t), intent(in), optional :: qmc_state_restart

        character(*), parameter :: this = 'check_qmc_opts'


        if (qmc_in%tau <= 0) call stop_all(this,'Tau must be positive.')
        if (qmc_in%shift_damping <= 0) call stop_all(this,'Shift damping must be positive.')

        if (need_length) then
            if (qmc_in%walker_length == 0) call stop_all(this,'Walker length must not be zero.')
            if (qmc_in%spawned_walker_length == 0) call stop_all(this,'Spawned walker length must not be zero.')
        end if

        if (.not. restarting) then
            if (qmc_in%D0_population <= 0) call stop_all(this, 'Initial population must be positive.')
        end if

        if (sys%read_in%comp .and. qmc_in%quasi_newton) call stop_all(this, 'Quasi-Newton not currently &
            &compatible with complex systems.')

        if (qmc_in%ncycles <= 0) call stop_all(this, 'mc_cycles must be positive.')
        if (qmc_in%nreport <= 0) call stop_all(this, 'nreports must be positive.')

        if (present(qmc_state_restart)) then
            if (sys%basis%tot_string_len /= size(qmc_state_restart%psip_list%states,dim=1)) call stop_all(this, &
                'Attempting to restart calculation within lua using different information string lengths.&
                & This is only possible by producing and reading back in a restart file.')
        end if

    end subroutine check_qmc_opts

    subroutine check_fci_opts(sys, fci_in, lanczos)

        ! Check the input options provided in the fci table.

        ! In:
        !   sys: system being studied.
        !   fci_in: input options for FCI.
        !   lanczos: is this a lanczos or full diagonalisation?

        use fci_utils, only: fci_in_t
        use system, only: sys_t, read_in

        use errors, only: stop_all

        type(sys_t), intent(in) :: sys
        type(fci_in_t), intent(in) :: fci_in
        logical, intent(in) :: lanczos

        character(*), parameter :: this = 'check_fci_opts'

        if (sys%system == read_in) then
            if (fci_in%analyse_fci_wfn /= 0 .and. sys%read_in%dipole_int_file == '') then
                call stop_all(this, 'Cannot analyse FCI wavefunction without a dipole integrals file.')
            end if
        end if

        if (lanczos) then
            if (fci_in%lanczos_string_len <= 0) call stop_all(this,'Lanczos basis not positive.')
            if (fci_in%nlanczos_eigv <= 0) call stop_all(this,'# lanczos eigenvalues not positive.')
        end if
        if (sys%basis%info_string_len /= 0) call stop_all(this, &
            'FCI is incompatible with additional information being stored in the bit string. Please implement if needed.')

    end subroutine check_fci_opts

    subroutine check_load_bal_opts(load_bal_in)

        ! Check load balancing input options.

        ! In:
        !   load_bal_in: load balancing options.

        use qmc_data, only: load_bal_in_t

        use errors, only: stop_all

        type(load_bal_in_t), intent(in) :: load_bal_in

        character(*), parameter :: this = 'check_load_bal_opts'

        if (load_bal_in%nslots < 0) call stop_all(this, 'Number of slots for load balancing is not positive.')
        if (load_bal_in%pop < 0) call stop_all(this, 'Load balancing population must be positive.')
        if (load_bal_in%percent < 0 .or. load_bal_in%percent > 1.0) &
            call stop_all(this, 'Percentage imbalance must be positive and less that 1.')
        if (load_bal_in%max_attempts < 0) call stop_all(this, 'Maximum number of load balancing attempts must be positive')

    end subroutine check_load_bal_opts

    subroutine check_dmqmc_opts(sys, dmqmc_in)

        ! Check validity of dmqmc input options.

        ! In:
        !   sys: system being studied.
        !   dmqmc_in: DMQMC options.

        use system, only: sys_t, heisenberg, ueg, read_in
        use dmqmc_data, only: dmqmc_in_t, hartree_fock_dm
        use calc, only: dmqmc_rdm_r2, doing_dmqmc_calc, dmqmc_full_r2
        use calc, only: dmqmc_staggered_magnetisation, dmqmc_energy_squared, dmqmc_correlation
        use calc, only: dmqmc_potential_energy, dmqmc_kinetic_energy, dmqmc_HI_energy

        use errors, only: stop_all
        use const, only: depsilon, p

        type(sys_t), intent(in) :: sys
        type(dmqmc_in_t), intent(in) :: dmqmc_in

        character(*), parameter :: this = 'check_dmqmc_opts'

        if (dmqmc_in%rdm%calc_inst_rdm .and. dmqmc_in%rdm%spawned_length == 0) call stop_all(this,'Spawned RDM length zero.')

        if (doing_dmqmc_calc(dmqmc_staggered_magnetisation) .and. (.not.sys%lattice%bipartite_lattice)) &
            call stop_all(this,'Staggered magnetisation calculation is only supported on a bipartite lattice.')

        if (.not. sys%system == heisenberg) then
            if (doing_dmqmc_calc(dmqmc_staggered_magnetisation)) &
                call stop_all(this,'The staggered magnetisation operator is not supported for this system.')
            if (doing_dmqmc_calc(dmqmc_energy_squared)) &
                call stop_all(this,'The energy squared operator is not supported for this system.')
            if (doing_dmqmc_calc(dmqmc_correlation)) &
                call stop_all(this,'The correlation function operator is not supported for this system.')
            if (dmqmc_in%rdm%calc_ground_rdm .or. dmqmc_in%rdm%calc_inst_rdm) &
                call stop_all(this,'The calculation of reduced density matrices is not supported for this system.')
        end if
        if (.not. sys%system == ueg) then
            if (doing_dmqmc_calc(dmqmc_kinetic_energy)) &
                call stop_all(this,'The kinetic energy operator is not supported for this system.')
            if (doing_dmqmc_calc(dmqmc_potential_energy)) &
                call stop_all(this,'The potential energy operator is not supported for this system.')
        end if

        if (allocated(dmqmc_in%correlation_sites)) then
            if (size(dmqmc_in%correlation_sites) /= 2) &
                call stop_all(this, 'You must enter exactly two sites for the correlation function option.')
        end if

        if (dmqmc_in%find_weights .and. dmqmc_in%calc_excit_dist) &
            call stop_all(this, 'find_weights and excit_dist options cannot be used together.')

        if (doing_dmqmc_calc(dmqmc_rdm_r2) .and. (.not. dmqmc_in%replica_tricks)) &
            call stop_all(this, 'The replica_tricks option must be used in order to calculate the Renyi-2 entropy.')
        if (doing_dmqmc_calc(dmqmc_rdm_r2) .and. (.not. dmqmc_in%rdm%calc_inst_rdm)) &
            call stop_all(this, 'The instantaneous_rdm option must be used in order to calculate the Renyi-2 entropy.')
        if (doing_dmqmc_calc(dmqmc_full_r2) .and. (.not. dmqmc_in%replica_tricks)) &
            call stop_all(this, 'The replica_tricks option must be used in order to calculate the Renyi-2 entropy.')
        if (doing_dmqmc_calc(dmqmc_HI_energy) .and. (.not. dmqmc_in%symmetric)) &
            call stop_all(this, 'Evaluation of interaction picture Hamiltonian only possible when using symmetric algorithm.')

        if (dmqmc_in%all_spin_sectors) then
            if (abs(sys%heisenberg%magnetic_field) > depsilon .or. &
                abs(sys%heisenberg%staggered_magnetic_field) > depsilon) &
                call stop_all(this, 'The use of all spin sectors simultaneously is not supported with magnetic fields.')
            if (dmqmc_in%rdm%calc_ground_rdm) &
                call stop_all(this, 'The use of all spin sectors simultaneously is not supported with ground-state RDMs.')
        end if

        if (dmqmc_in%vary_weights .and. (.not. dmqmc_in%weighted_sampling)) then
            call stop_all(this, 'The vary_weights option can only be used together with the weighted_sampling option.')
        end if

        if (dmqmc_in%ipdmqmc .and. .not. dmqmc_in%grand_canonical_initialisation &
            & .and.  dmqmc_in%metropolis_attempts == 0) then
            call stop_all(this, 'metropolis_attempts must be non-zero to sample the correct initial density matrix&
                                 & if not using grand_canonical_initialisation.')
        end if
        if (dmqmc_in%grand_canonical_initialisation .and. dmqmc_in%initial_matrix == hartree_fock_dm &
            .and. sys%system /= ueg .and. sys%system /= read_in) then
            call stop_all(this, 'Reweighting of initial matrix not supported for this system. Please implement.')
        end if
        if (dmqmc_in%symmetric .and. dmqmc_in%ipdmqmc .and. sys%system /= ueg) then
            call stop_all(this, 'Symmetric propagation is only implemented for the UEG. Please implement.')
        end if

        if (dmqmc_in%target_beta < depsilon .and. dmqmc_in%grand_canonical_initialisation) then
            call stop_all(this, 'target_beta must be greater than zero if grand_canonical_initialisation is to be used.')
        else if (dmqmc_in%target_beta < 0) then
            call stop_all(this, 'target_beta must be non-negative.')
        end if

        if (dmqmc_in%metropolis_attempts < 0) then
            call stop_all(this, 'metropolis_attempts must be greater than zero.')
        end if

        if (sys%read_in%comp) call stop_all(this, 'Complex DMQMC not yet implemented')

        if (sys%basis%info_string_len /= 0) call stop_all(this, &
            'DMQMC is incompatible with additional information being stored in the bit string. Please implement if needed.')

    end subroutine check_dmqmc_opts

    subroutine check_ccmc_opts(sys, ccmc_in)

        ! Check the CCMC input options

        ! In:
        !   sys: system being studied.
        !   ccmc_in: CCMC options


        use qmc_data, only: ccmc_in_t
        use system, only: sys_t, read_in
        use errors, only: stop_all

        type(sys_t), intent(in) :: sys
        type(ccmc_in_t), intent(in) :: ccmc_in

        character(*), parameter :: this = 'check_ccmc_opts'

        if (ccmc_in%move_freq >= 32) then
            call stop_all(this, "move_frequency must be less than 32")
        else if (ccmc_in%move_freq < 0) then
            call stop_all(this, "move_frequency must be non-negative")
        end if

        if (ccmc_in%cluster_multispawn_threshold <= 0) then
            call stop_all(this, "cluster_multispawn_threshold must be positive")
        end if

        if (ccmc_in%density_matrices .and. sys%system /= read_in) then
            call stop_all(this, "CCMC density matrices not implemented for this system type.")
        end if

        if (sys%read_in%comp .and. ccmc_in%linked) call stop_all(this, 'Linked complex CCMC not yet implemented')
        if (sys%read_in%comp .and. ccmc_in%density_matrices) &
                call stop_all(this, 'RDMS not yet implemented for complex CCMC.')

        if (ccmc_in%even_selection) then
            if (.not. ccmc_in%full_nc) call &
                stop_all(this, "Even selection routines currently only work with full noncomposite enabled.")
            if (ccmc_in%linked) call stop_all(this, &
                "Even selection routines are not currently confirmed to work with linked ccmc, so have been disabled")
        end if

        if (sys%basis%info_string_len /= 0 .and. .not. ccmc_in%even_selection) call stop_all(this, &
            'Additional space allocated in bit strings for no reason. Something has gone wrong.')

    end subroutine check_ccmc_opts

    subroutine check_blocking_opts(sys, blocking_in, restart_in)

        ! Check block on the fly input options.
        ! The majority of these can't be sanity-checked as a nonsensical (ie. negative)
        ! value is used to indicate that we should default to suitable values.

        ! In:
        !   sys: system being studied.
        !   blocking_in: Reblocking on the fly input options.
        !   restart_in: Calculation restarting input options.

        use qmc_data, only: blocking_in_t, restart_in_t
        use errors, only: warning, stop_all
        use system, only: sys_t
        use const, only: p

        type(sys_t), intent(in) :: sys
        type(blocking_in_t), intent(in) :: blocking_in
        type(restart_in_t), intent(in) :: restart_in

        character(*), parameter :: this = 'check_blocking_opts'

        if (blocking_in%blocking_on_the_fly) then
            if (restart_in%read_restart) then
                call warning(this, 'Restarting from restart file while using reblocking. &
                        &No reblocking data from before restart has been saved. If this annoys you send patches.')
            end if

            if (restart_in%write_restart) then
                call warning(this, 'Writing a restart file from a calculation using on-the-fly reblocking. &
                        &No reblocking data has been saved. If this annoys you send patches.')
            end if
        end if

        if (blocking_in%auto_shift_damping .and. .not. blocking_in%blocking_on_the_fly) then
            call stop_all(this, 'Automatic shift damping optimisation requires blocking on the fly to be activated. &
                        & Please restart calculation with blocking enabled.')
        end if

        if (blocking_in%auto_shift_damping .and. blocking_in%shift_damping_precision < 1.5_p) then
            call stop_all(this, 'Automatic shift damping optimisation precision should be a positive float greater &
                       &than or equal to 1.5 for reasonable results.')
        end if

        if (sys%read_in%comp .and. blocking_in%blocking_on_the_fly) then
            call stop_all(this, "Blocking on the fly is not currently compatible with complex calculations.")
        end if

    end subroutine check_blocking_opts

end module check_input

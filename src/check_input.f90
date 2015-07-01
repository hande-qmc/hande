module check_input
! Parse input options and check input for validity.

implicit none

contains

    subroutine check_sys(sys)

        ! Check the input parameters for the system.

        ! I don't pretend this is the most comprehensive of tests, but at least
        ! make sure a few things are not completely insane.

        ! In:
        !    sys: system object, as set in read_input

        use system

        use errors, only: stop_all, warning

        type(sys_t), intent(in) :: sys

        integer :: ivec, jvec
        character(*), parameter :: this='check_sys'

        if (sys%system /= heisenberg) then
            if (sys%nel <= 0) call stop_all(this,'Number of electrons must be positive.')
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
            if (sys%heisenberg%staggered_magnetic_field /= 0.0_p .and. (.not.sys%lattice%bipartite_lattice)) &
                call stop_all(this, 'Cannot set a staggered field for this lattice because it is frustrated.')
            if (sys%heisenberg%staggered_magnetic_field /= 0.0_p .and. sys%heisenberg%magnetic_field /= 0.0_p) &
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


    end subroutine check_sys

    subroutine check_fciqmc_opts(sys, fciqmc_in)

        ! Check the FCIQMC specific options.

        ! In:
        !   sys: system being studied
        !   fciqmc_in: FCIQMC input options

        use qmc_data, only: fciqmc_in_t, single_basis, no_guiding, neel_singlet_guiding, neel_singlet
        use system, only: sys_t, heisenberg
        use errors, only: stop_all

        type(sys_t), intent(in) :: sys
        type(fciqmc_in_t), intent(in) :: fciqmc_in

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

    end subroutine check_fciqmc_opts

    subroutine check_qmc_opts(qmc_in, simple_fciqmc)

        ! Check options common to QMC methods.

        ! In:
        !   qmc_in: QMC input options
        !   simple_fciqmc: true if using the simple FCIQMC algorithm.

        use qmc_data, only: qmc_in_t
        use errors, only: stop_all

        type(qmc_in_t), intent(in) :: qmc_in
        logical, intent(in) :: simple_fciqmc

        character(*), parameter :: this = 'check_qmc_opts'

        if (qmc_in%tau <= 0) call stop_all(this,'Tau not positive.')
        if (qmc_in%shift_damping <= 0) call stop_all(this,'Shift damping not positive.')

        if (.not. simple_fciqmc) then
            if (qmc_in%walker_length == 0) call stop_all(this,'Walker length zero.')
            if (qmc_in%spawned_walker_length == 0) call stop_all(this,'Spawned walker length zero.')
        end if

    end subroutine check_qmc_opts

    subroutine check_fci_opts(sys, fci_in, lanczos)

        ! Check the input options provided in the fci table

        ! In:
        !   sys: system being studied
        !   fci_in: Input options for FCI
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

    end subroutine check_fci_opts

    subroutine check_load_bal_opts(load_bal_in)

        ! Check load balancing input options

        ! In:
        !   load_bal_in: load balancing options

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

        ! Check validity of dmqmc input options

        ! In:
        !   sys: system being studied
        !   dmqmc_in: DMQMC options

        use system, only: sys_t, heisenberg
        use dmqmc_data, only: dmqmc_in_t
        use calc, only: dmqmc_calc_type, dmqmc_energy, dmqmc_rdm_r2, dmqmc_staggered_magnetisation, doing_dmqmc_calc

        use errors, only: stop_all
        use const, only: depsilon

        type(sys_t), intent(in) :: sys
        type(dmqmc_in_t), intent(in) :: dmqmc_in

        character(*), parameter :: this = 'check_dmqmc_opts'

        if (dmqmc_in%rdm%calc_inst_rdm .and. dmqmc_in%rdm%spawned_length == 0) call stop_all(this,'Spawned RDM length zero.')

        if (sys%system == heisenberg) then
            if (doing_dmqmc_calc(dmqmc_staggered_magnetisation) .and. (.not.sys%lattice%bipartite_lattice)) then
                call stop_all(this,'Staggered magnetisation can only be calculated on a bipartite lattice.')
            end if
        else if (dmqmc_calc_type /= dmqmc_energy) then
            if (dmqmc_in%all_spin_sectors) call stop_all(this, &
                'The option to use all symmetry sectors at the same time is only available for the Heisenberg model.')
            call stop_all(this, 'The observable requested is not currently implemented for this Hamiltonian.')
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

        if (dmqmc_in%all_spin_sectors) then
            if (abs(sys%heisenberg%magnetic_field) > depsilon .or. &
                abs(sys%heisenberg%staggered_magnetic_field) > depsilon) &
                call stop_all(this, 'The use_all_spin_sectors option cannot be used with magnetic fields.')
            if (dmqmc_in%rdm%calc_ground_rdm) &
                call stop_all(this, 'The use_all_spin_sectors and ground_state_rdm options cannot be used together.')
        end if

        if (dmqmc_in%vary_weights .and. (.not. dmqmc_in%weighted_sampling)) then
            call stop_all(this, 'The vary_weights option can only be used together with the weighted_sampling option.')
        end if
    end subroutine check_dmqmc_opts

end module check_input

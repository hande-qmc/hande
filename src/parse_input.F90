module parse_input
! Parse input options and check input for validity.

implicit none

contains

    subroutine check_input(sys, qmc_in, fciqmc_in, ccmc_in, semi_stoch_in, restart_in, reference, load_bal_in, dmqmc_in, fci_in)

        ! I don't pretend this is the most comprehensive of tests, but at least
        ! make sure a few things are not completely insane.

        ! In/Out:
        !    sys: system object, as set in read_input (invalid settings are overridden).
        !    qmc_in: input options relating to QMC methods.
        !    fciqmc_in: input options relating to FCIQMC.
        !    ccmc_in: input options relating to CCMC.
        !    load_bal_in: input options for load balancing.
        !    fci_in: input options relating to FCI.
        ! In:
        !    semi_stoch_in: Input options for the semi-stochastic adaptation.
        !    restart_in: input options for HDF5 restart files.
        !    reference: reference determinant.
        !    dmqmc_in: input options relating to DMQMC.

        use calc
        use const
        use qmc_data, only: qmc_in_t, fciqmc_in_t, ccmc_in_t, semi_stoch_in_t
        use qmc_data, only: restart_in_t, reference_t, load_bal_in_t, empty_determ_space
        use qmc_data, only: neel_singlet, neel_singlet_guiding, no_guiding, single_basis
        use dmqmc_data, only: dmqmc_in_t
        use fci_utils, only: fci_in_t
        use system
        use parallel, only: parent

        use errors, only: stop_all, warning

        type(sys_t), intent(inout) :: sys
        type(qmc_in_t), intent(inout) :: qmc_in
        type(fciqmc_in_t), intent(inout) :: fciqmc_in
        type(ccmc_in_t), intent(inout) :: ccmc_in
        type(semi_stoch_in_t), intent(in) :: semi_stoch_in
        type(restart_in_t), intent(in) :: restart_in
        type(reference_t), intent(in) :: reference
        type(load_bal_in_t), intent(in) :: load_bal_in
        type(dmqmc_in_t), intent(in) :: dmqmc_in
        type(fci_in_t), intent(inout) :: fci_in

        integer :: ivec, jvec
        character(*), parameter :: this='check_input'

        if (sys%system /= heisenberg) then
            if (sys%nel <= 0) call stop_all(this,'Number of electrons must be positive.')
            if (fciqmc_in%trial_function /= single_basis) call stop_all(this, 'Only a single determinant can be used as the&
                                                     & reference state for this system. Other trial functions are not avaliable.')
            if (fciqmc_in%guiding_function /= no_guiding) &
                call stop_all(this, 'Importance sampling is only avaliable for the Heisenberg model&
                                         & currently.')
            if (dmqmc_in%all_spin_sectors) call stop_all(this,'The option to use all symmetry sectors at the same time is only&
                                         & available for the Heisenberg model.')
        end if

        if (sys%system == read_in) then

            if (fci_in%analyse_fci_wfn /= 0 .and. sys%read_in%dipole_int_file == '') then
                if (parent) call warning(this, 'Cannot analyse FCI wavefunction without a dipole &
                             &integrals file.  Turning fci_in%analyse_fci_wfn option off...')
                fci_in%analyse_fci_wfn = 0
            end if

        else

            if (sys%system /= ueg .and. sys%system /= ringium) then
                if (.not.(allocated(sys%lattice%lattice))) call stop_all(this, 'Lattice vectors not provided')
                do ivec = 1, sys%lattice%ndim
                    do jvec = ivec+1, sys%lattice%ndim
                        if (dot_product(sys%lattice%lattice(:,ivec), sys%lattice%lattice(:,jvec)) /= 0) then
                            call stop_all(this, 'Lattice vectors are not orthogonal.')
                        end if
                    end do
                end do
                if (doing_calc(mc_canonical_kinetic_energy)) call stop_all(this, 'estimate_canonical_kinetic_energy&
                                                                           & only implemented for the UEG.')
            end if

            if (sys%system == heisenberg) then
                if (.not. dmqmc_in%all_spin_sectors) then
                    if (sys%Ms > sys%lattice%nsites) call stop_all(this,'Value of Ms given is too large for this lattice.')
                    if ((-sys%Ms) > sys%lattice%nsites) call stop_all(this,'Value of Ms given is too small for this lattice.')
                    if (mod(abs(sys%Ms),2) /=  mod(sys%lattice%nsites,2)) call stop_all(this, 'Ms value specified is not&
                                                                                          & possible for this lattice.')
                end if
                if (sys%heisenberg%staggered_magnetic_field /= 0.0_p .and. (.not.sys%lattice%bipartite_lattice)) &
                    call stop_all(this, 'Cannot set a staggered field&
                                        & for this lattice because it is frustrated.')
                if (sys%heisenberg%staggered_magnetic_field /= 0.0_p .and. sys%heisenberg%magnetic_field /= 0.0_p) &
                    call stop_all(this, 'Cannot set a uniform and a staggered field at the same time.')
                if ((fciqmc_in%guiding_function==neel_singlet_guiding) .and. fciqmc_in%trial_function /= neel_singlet) &
                    call stop_all(this, 'This guiding function is only avaliable when using the Neel singlet state &
                                        &as an energy estimator.')
                if (doing_dmqmc_calc(dmqmc_staggered_magnetisation) .and. (.not.sys%lattice%bipartite_lattice)) then
                    if (parent) call warning(this,'Staggered magnetisation can only be calculated on a bipartite lattice.&
                                          & This is not a bipartite lattice. Changing options so that it will not be calculated.')
                    dmqmc_calc_type = dmqmc_calc_type - dmqmc_staggered_magnetisation
                end if
            else if (sys%system == hub_k .or. sys%system == hub_real) then
                if (sys%nel > 2*sys%lattice%nsites) call stop_all(this, 'More than two electrons per site.')
            end if

            if (sys%lattice%ndim > 3) call stop_all(this, 'Limited to 1,  2 or 3 dimensions')
            if (sys%system == ueg .and. sys%lattice%ndim == 1) call stop_all(this, 'UEG only functional in 2D and 3D')
            if (sys%system == ringium .and. sys%lattice%ndim /= 1) then
                if (parent) call warning(this, 'Ringium must be 1D')
                sys%lattice%ndim = 1
            end if
            if (sys%system == ringium .and. mod(sys%nel + sys%ringium%maxlz, 2) == 0) call stop_all(this, 'Max lz must have &
                                        &opposite parity to the number of electrons in ringium')
            ! [review] - JSS: danger in equality comparison with floats.  Use abs(radius)<depsilon instead.
            if (sys%system == ringium .and. sys%ringium%radius == 0.0) call stop_all(this, 'Ringium must have a finite radius')
            if (sys%system == ringium .and. sys%nel /= sys%ms) then
                if (parent) call warning(this, 'Ringium must have spin and number of electrons equal')
                sys%ms = sys%nel
            end if

        end if

        ! Real amplitude checks.
        if (qmc_in%real_amplitudes) then
            if (doing_calc(ct_fciqmc_calc) .or. doing_calc(hfs_fciqmc_calc)) then
                call stop_all(this, 'The real_amplitudes option is not implemented with the method you have requested.')
            end if
        end if

        ! Semi-stochastic checks.
        if (semi_stoch_in%start_iter /= 0 .and. semi_stoch_in%space_type == empty_determ_space .and. parent) &
            call warning(this,'You have specified an iteration to turn semi-stochastic on but have not &
                         &specified a deterministic space to use.')
        if (semi_stoch_in%space_type /= empty_determ_space .and. (doing_calc(dmqmc_calc) .or. &
                                   doing_calc(ct_fciqmc_calc) .or. doing_calc(hfs_fciqmc_calc))) &
              call stop_all(this, 'Semi-stochastic is only implemented with the FCIQMC method.')

        if (fciqmc_in%init_spin_inv_D0 .and. sys%Ms /= 0) then
            if (parent) call warning(this, 'Flipping the reference state will give &
                                            &a state which has a different value of Ms and so cannot be used here.')
            fciqmc_in%init_spin_inv_D0 = .false.
        end if

        if (allocated(dmqmc_in%correlation_sites) .and. size(dmqmc_in%correlation_sites) /= 2) &
                            call stop_all(this, 'You must enter exactly two sites for the correlation function option.')

          if (dmqmc_in%find_weights .and. dmqmc_in%calc_excit_dist) call stop_all(this, 'DMQMC_FIND_WEIGHTS and OUTPUT_EXCITATION&
              &_DISTRIBUTION options cannot be used together.')

        ! Calculation specific checking.
        if (doing_calc(lanczos_diag)) then
            if (fci_in%lanczos_string_len <= 0) call stop_all(this,'Lanczos basis not positive.')
            if (fci_in%nlanczos_eigv <= 0) call stop_all(this,'# lanczos eigenvalues not positive.')
        end if

        if (.not.doing_calc(dmqmc_calc) .and. dmqmc_calc_type /= 0 .and. parent) call warning(this,&
               'You are not performing a DMQMC calculation but have requested DMQMC options to be calculated.')
        if (doing_calc(fciqmc_calc)) then
            if (.not.doing_calc(simple_fciqmc_calc)) then
                if (qmc_in%walker_length == 0) call stop_all(this,'Walker length zero.')
                if (qmc_in%spawned_walker_length == 0) call stop_all(this,'Spawned walker length zero.')
            end if
            if (dmqmc_in%rdm%calc_inst_rdm .and. dmqmc_in%rdm%spawned_length == 0) call stop_all(this,'Spawned RDM length zero.')
            if (qmc_in%tau <= 0) call stop_all(this,'Tau not positive.')
            if (qmc_in%shift_damping <= 0) call stop_all(this,'Shift damping not positive.')
            if (allocated(reference%occ_list0)) then
                if (size(reference%occ_list0) /= sys%nel) then
                    if (sys%system /= heisenberg) then
                        call stop_all(this,'Number of electrons specified is different from &
                        &number of electrons used in the reference determinant.')
                    end if
                end if
            end if
            if (load_bal_in%nslots < 0) call stop_all(this, 'Number of slots for load balancing is not positive.')
            if (load_bal_in%pop < 0) call stop_all(this, 'Load balancing population must be positive.')
            if (load_bal_in%percent < 0 .or. load_bal_in%percent > 1.0) &
                call stop_all(this, 'Percentage imbalance must be positive and less that 1.')
            if (load_bal_in%max_attempts < 0) call stop_all(this, 'Maximum number of load balancing attempts must be positive')
        end if
        if (doing_calc(ct_fciqmc_calc)) qmc_in%ncycles = 1

        if (doing_dmqmc_calc(dmqmc_rdm_r2) .and. (.not. dmqmc_in%replica_tricks)) call stop_all(this,&
                    'The replica_tricks option must be used in order to calculate the Renyi-2 entropy.')
        if (doing_dmqmc_calc(dmqmc_rdm_r2) .and. (.not. dmqmc_in%rdm%calc_inst_rdm)) call stop_all(this,&
                    'The instantaneous_rdm option must be used in order to calculate the Renyi-2 entropy.')

        ! If the FINITE_CLUSTER keyword was detected then make sure that
        ! we are doing a calculation in real-space. If we're not then
        ! unset finite cluster,tell the user and carry on
        if(sys%momentum_space) then
            if (sys%real_lattice%finite_cluster .and. parent) call warning(this,'FINITE_CLUSTER keyword only valid for hubbard&
                                      & calculations in real-space: ignoring keyword')
            sys%real_lattice%finite_cluster = .false.
        end if

        if (dmqmc_in%all_spin_sectors) then
            if (.not. doing_calc(dmqmc_calc)) call stop_all(this, 'The use_all_spin_sectors option can only be used in&
                                                                   & DMQMC calculations.')
            if (abs(sys%heisenberg%magnetic_field) > depsilon .or. &
                abs(sys%heisenberg%staggered_magnetic_field) > depsilon) &
                    call stop_all(this, 'The use_all_spin_sectors option cannot be used with magnetic fields.')
            if (dmqmc_in%rdm%calc_ground_rdm) call stop_all(this, 'The use_all_spin_sectors and ground_state_rdm options&
                                                      & cannot be used together.')
        end if

        if (dmqmc_in%vary_weights .and. (.not. dmqmc_in%weighted_sampling)) then
            call stop_all(this, 'The vary_weights option can only be used together with the weighted_sampling option.')
        end if
        if (sys%system /= heisenberg .and. dmqmc_calc_type > dmqmc_energy) then
            call stop_all(this, 'The observable requested is not currently implemented for this Hamiltonian.')
        end if

        if (parent) write (6,'(/,1X,13("-"),/)')

    end subroutine check_input

end module parse_input

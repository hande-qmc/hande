module parse_input
! Parse input options and check input for validity.

use const

use parallel, only: parent, block_size
use errors
use hilbert_space
use calc
use lanczos
use determinants
use determinant_enumeration, only: write_determinants, determinant_file
use fciqmc_data
use restart_hdf5, only: restart_info_global, restart_info_global_shift
use hfs_data
use semi_stoch

implicit none

contains

    subroutine read_input(sys)

        ! Read input options from a file (if specified on the command line) or via
        ! STDIN.

        ! In/Out:
        !   sys: system being studied.  Parameters specified in the input file
        !        are set directly in the system object, components which are not
        !        mentioned in the input file are not altered.

! nag doesn't automatically bring in command-line option handling.
#ifdef NAGF95
        use f90_unix_env
#endif

        use system

        use input
        use utils, only: get_free_unit
        use checking, only: check_allocate

#ifdef NAGF95
        use f90_unix_env, ONLY: getarg,iargc
#else
        integer :: iargc ! External function.
#endif

        type (sys_t), intent(inout) :: sys

        character(255) :: cInp
        character(100) :: w
        integer :: ios
        logical :: eof, t_exists
        integer :: ivec, i, j, ierr, nweights

        if (iargc() > 0) then
            ! Input file specified on the command line.
            ir = get_free_unit()
            call GetArg(1, cInp)
            inquire(file=cInp, exist=t_exists)
            if (.not.t_exists) then
                call stop_all('read_input','File does not exist:'//trim(cInp))
            end if
            open(ir, file=cInp, status='old', form='formatted', iostat=ios)
        else
            if (parent) write (6,'(a19)') 'Reading from STDIN'
            ir = 5
            ios = 0
        end if

        if (parent) write (6,'(a14,/,1X,13("-"),/)') 'Input options'
        call input_options(echo_lines=parent, skip_blank_lines=.true.)

        do ! loop over lines in input file.
            call read_line(eof)
            if (eof) exit
            call readu(w)
            select case(w)

            ! System type
            case('HUBBARD_REAL')
                sys%system = hub_real
            case('HUBBARD_K','HUBBARD_MOMENTUM')
                sys%system = hub_k
            case('HEISENBERG')
                sys%system = heisenberg
            case('READ')
                sys%system = read_in
                if (item /= nitems) call reada(sys%read_in%fcidump)
            case('UEG')
                sys%system = ueg
            case('CHUNG-LANDAU')
                sys%system = chung_landau

                ! System information.
            case('LATTICE')
                ! Lattice block
                call read_line(eof)
                if (eof) call stop_all('read_input','Unexpected end of file reading lattice vectors.')
                ! nitems gives the number of items in the line, and thus the number
                ! of dimensions...
                sys%lattice%ndim = nitems
                allocate(sys%lattice%lattice(sys%lattice%ndim,sys%lattice%ndim), stat=ierr)
                call check_allocate('sys%lattice%lattice',sys%lattice%ndim*sys%lattice%ndim,ierr)
                do ivec = 1, sys%lattice%ndim
                    if (nitems /= sys%lattice%ndim) &
                        call stop_all('read_input', 'Do not understand lattice vector.')
                    do i = 1, sys%lattice%ndim
                        call readi(sys%lattice%lattice(i, ivec))
                    end do
                    if (ivec /= sys%lattice%ndim) then
                        call read_line(eof)
                        if (eof) call stop_all('read_input', 'Unexpected end of file reading lattice vectors.')
                    end if
                end do
            case('2D')
                if (sys%lattice%ndim > 0) then
                    if (parent) call warning('read_input','Dimension already set; ignoring keyword '//w)
                else
                    sys%lattice%ndim = 2
                end if
            case('3D')
                if (sys%lattice%ndim > 0) then
                    if (parent) call warning('read_input','Dimension already set; ignoring keyword '//w)
                else
                    sys%lattice%ndim = 3
                end if
            case('NEL', 'ELECTRONS')
                if (sys%system == heisenberg) &
                     call stop_all('read_input', 'Cannot set electron number for Heisenberg. &
                     &Please enter a Ms value instead.')
                call readi(sys%nel)

            ! Hubbard-specific system info
            case('T')
                call readf(sys%hubbard%t)
            case('U')
                call readf(sys%hubbard%u)

            ! Heisenberg-specific system info.
            case('J')
                call readf(sys%heisenberg%J)
            case('MAGNETIC_FIELD')
                call readf(sys%heisenberg%magnetic_field)
            case('STAGGERED_MAGNETIC_FIELD')
                call readf(sys%heisenberg%staggered_magnetic_field)

            ! UEG-specific system info.
            case('RS','DENSITY')
                call readf(sys%ueg%r_s)
            case('ECUTOFF')
                call readf(sys%ueg%ecutoff)

            ! Additional information for systems read in (i.e. molecular)
            case('DIPOLE_INTEGRALS')
                call reada(sys%read_in%dipole_int_file)

            ! Select symmetry of wavefunction.
            case('MS')
                call readi(ms_in)
            case('SYM','SYMMETRY')
                call readi(sym_in)
            case("LZ")
                sys%read_in%useLz=.true.
            case('CAS')
                do i = 1,2
                    call readi(sys%CAS(i))
                end do
            case('RAS')
                do i = 1,2
                    call readi(RAS(i))
                end do
                call readi(ras3_max)
            case('TWIST')
                allocate(sys%k_lattice%ktwist(nitems-item), stat=ierr)
                call check_allocate('sys%k_lattice%ktwist',nitems-item,ierr)
                do i = 1, nitems-item
                    call readf(sys%k_lattice%ktwist(i))
                end do

            ! Calculation type.
            case('EXACT','FCI')
                calc_type = calc_type + exact_diag
            case('LANCZOS_DIRECT')
                calc_type = calc_type + lanczos_diag
                direct_lanczos = .true.
            case('LANCZOS')
                calc_type = calc_type + lanczos_diag
            case('SIMPLE_FCIQMC')
                calc_type = calc_type + fciqmc_calc + simple_fciqmc_calc
            case('FCIQMC')
                calc_type = calc_type + fciqmc_calc
            case('IFCIQMC')
                calc_type = calc_type + fciqmc_calc
                initiator_approximation = .true.
            case('CT_FCIQMC')
                calc_type = calc_type + ct_fciqmc_calc
            case('DMQMC')
                calc_type = calc_type + dmqmc_calc
            case('CCMC')
                calc_type = calc_type + ccmc_calc
            case('ICCMC')
                calc_type = calc_type + ccmc_calc
                initiator_approximation = .true.
            case('HELLMANN-FEYNMAN')
                calc_type = calc_type + hfs_fciqmc_calc
            case('ESTIMATE_HILBERT_SPACE')
                calc_type = calc_type + mc_hilbert_space
                call readi(nhilbert_cycles)

            case('REPLICA_TRICKS')
                replica_tricks = .true.
            case('PROPAGATE_TO_BETA')
                propagate_to_beta = .true.
            case('INIT_BETA')
                call readf(init_beta)
            case('METROPOLIS_ATTEMPTS')
                call readi(metropolis_attempts)
            case('FREE_ELECTRON_TRIAL')
                free_electron_trial = .true.
            case('CHEM_POT')
                call readf(chem_pot)
            case('GRAND_CANONICAL_ENSEMBLE')
                grand_canonical_ensemble = .true.

            case('CCMC_FULL_NC')
                ccmc_full_nc = .true.
            case('CCMC_LINKED')
                linked_ccmc = .true.

            case('REAL_AMPLITUDES')
                real_amplitudes = .true.
            case('SPAWN_CUTOFF')
                call readf(spawn_cutoff)

            ! Semi-stochastic options.
            case('SEMI_STOCH_ITERATION')
                call readi(semi_stoch_start_iter)
                real_amplitudes = .true.
            case('SEMI_STOCH_SHIFT_START')
                call readi(semi_stoch_shift_iter)
                semi_stoch_start_iter = -1
                real_amplitudes = .true.
            ! Deterministic spaces.
            case('SEMI_STOCH_HIGH_POP')
                determ_space_type = high_pop_determ_space
                call readi(determ_target_size)
            case('SEMI_STOCH_READ')
                determ_space_type = read_determ_space
                ! Not needed.
                determ_target_size = -1
            case('WRITE_DETERM_SPACE')
                write_determ_space = .true.
            case('SEMI_STOCH_COMBINE_ANNIHIL')
                separate_determ_annihil = .false.

            ! DMQMC expectation values to be calculated.
            case('DMQMC_FULL_RENYI_2')
                dmqmc_calc_type = dmqmc_calc_type + dmqmc_full_r2
            case('DMQMC_ENERGY')
                dmqmc_calc_type = dmqmc_calc_type + dmqmc_energy
            case('DMQMC_ENERGY_SQUARED')
                dmqmc_calc_type = dmqmc_calc_type + dmqmc_energy_squared
            case('DMQMC_CORRELATION_FUNCTION')
                dmqmc_calc_type = dmqmc_calc_type + dmqmc_correlation
                allocate(correlation_sites(nitems-1), stat=ierr)
                call check_allocate('correlation_sites',nitems-1,ierr)
                do i = 1, nitems-1
                    call readi(correlation_sites(i))
                end do
            case('DMQMC_STAGGERED_MAGNETISATION')
                dmqmc_calc_type = dmqmc_calc_type + dmqmc_staggered_magnetisation
            case('DMQMC_WEIGHTED_SAMPLING')
                call readi(nweights)
                allocate(dmqmc_sampling_probs(nweights), stat=ierr)
                call check_allocate('dmqmc_sampling_probs', nweights, ierr)
                dmqmc_weighted_sampling = .true.
                call read_line(eof)
                if (eof) call stop_all('read_input', 'Unexpected end of file reading DMQMC weights.')
                do i = 1, nweights
                    call readf(dmqmc_sampling_probs(i))
                end do
            case('DMQMC_VARY_WEIGHTS')
                call readi(finish_varying_weights)
                dmqmc_vary_weights = .true.
            case('DMQMC_FIND_WEIGHTS')
                dmqmc_find_weights = .true.
                dmqmc_weighted_sampling = .true.
            case('OUTPUT_EXCITATION_DISTRIBUTION')
                calculate_excit_distribution = .true.
            case('USE_ALL_SYM_SECTORS')
                all_sym_sectors = .true.
            case('USE_ALL_MOM_SECTORS')
                all_mom_sectors = .true.
            case('REDUCED_DENSITY_MATRIX')
                call readi(nrdms)
                allocate(rdms(nrdms), stat=ierr)
                call check_allocate('rdms', nrdms, ierr)
                doing_reduced_dm = .true.
                do i = 1, nrdms
                    call read_line(eof)
                    if (eof) call stop_all('read_input','Unexpected end of file reading reduced density matrices.')
                    rdms(i)%A_nsites = nitems
                    allocate(rdms(i)%subsystem_A(nitems))
                    call check_allocate('rdms(i)%subsystem_A',nitems,ierr)
                    do j = 1, nitems
                        call readi(rdms(i)%subsystem_A(j))
                    end do
                end do
            case('GROUND_STATE_RDM')
                calc_ground_rdm = .true.
            case('INSTANTANEOUS_RDM')
                calc_inst_rdm = .true.
            case('OUTPUT_RDM')
                output_rdm = .true.
            case('EXACT_RDM_EIGENVALUES')
                doing_exact_rdm_eigv = .true.
                calc_ground_rdm = .true.
            case('CONCURRENCE')
                doing_concurrence = .true.
            case('VON_NEUMANN_ENTROPY')
                doing_von_neumann_entropy = .true.
            case('RENYI_ENTROPY_2')
                dmqmc_calc_type = dmqmc_calc_type + dmqmc_rdm_r2
            case('START_AVERAGING')
                call readi(start_averaging)
            ! calculation options: DMQMC
            case('TRUNCATION_LEVEL')
                truncate_space = .true.
                call readi(truncation_level)
            case('HALF_DENSITY_MATRIX')
                half_density_matrix = .true.

            ! Calculation options: lanczos.
            case('LANCZOS_BASIS')
                call readi(lanczos_string_len)
            case('LANCZOS_SOLUTIONS','LANCZOS_SOLNS')
                call readi(nlanczos_eigv)
            case('SPARSE_HAMILTONIAN')
                use_sparse_hamil = .true.

            ! Calculation options: lanczos/exact diagonalisation.
            case('PRINT_FCI_WFN')
                print_fci_wfn = -1
                if (item /= nitems) call readi(print_fci_wfn)
                if (item /= nitems) call reada(print_fci_wfn_file)
            case('ANALYSE_FCI_WFN')
                analyse_fci_wfn = -1
                if (item /= nitems) call readi(analyse_fci_wfn)

            ! Calculation options: fciqmc.
            case('MC_CYCLES')
                call readi(ncycles)
            case('NREPORTS')
                call readi(nreport)
                if (nreport < 0) nreport = huge(nreport)
            case('BETA_LOOPS')
                call readi(beta_loops)
            case('WALKER_LENGTH')
                call readi(walker_length)
                if (item /= nitems) then
                    call readu(w)
                    if (w == 'MB') then
                        walker_length = -walker_length
                    else
                        call report('Keyword '//trim(w)//' not recognized.', .true.)
                    end if
                end if
            case('SPAWNED_WALKER_LENGTH')
                call readi(spawned_walker_length)
                if (item /= nitems) then
                    call readu(w)
                    if (w == 'MB') then
                        spawned_walker_length = -spawned_walker_length
                    else
                        call report('Keyword '//trim(w)//' not recognized.', .true.)
                    end if
                end if
            case('SPAWNED_RDM_LENGTH')
                call readi(spawned_rdm_length)
                if (item /= nitems) then
                    call readu(w)
                    if (w == 'MB') then
                        spawned_rdm_length = -spawned_rdm_length
                    else
                        call report('Keyword '//trim(w)//' not recognized.', .true.)
                    end if
                end if
            case('TAU')
                call readf(tau)
            case('TAU_SEARCH')
                tau_search = .true.
            case('INITIAL_SHIFT')
                call readf(initial_shift)
                ! We assume the user is sensible/knows what he/she is doing if
                ! initial_shift and vary_shift_from are set.
                vary_shift_from = initial_shift
            case('VARY_SHIFT_FROM')
                call readu(w)
                if (w == 'PROJE') then
                    vary_shift_from_proje = .true.
                else
                    call reread(0)
                    call readf(vary_shift_from)
                end if
            case('DMQMC_AVERAGE_SHIFT')
                call readi(average_shift_until)
            case('VARYSHIFT_TARGET')
                call readf(target_particles)
            case('INIT_POP')
                call readf(D0_population)
            case('REFERENCE_DET')
                allocate(occ_list0(nitems-1), stat=ierr)
                call check_allocate('occ_list0',nitems-1,ierr)
                do i = 1, nitems-1
                    call readi(occ_list0(i))
                end do
            case('HS_REFERENCE_DET')
                allocate(hs_occ_list0(nitems-1), stat=ierr)
                call check_allocate('hs_occ_list0',nitems-1,ierr)
                do i = 1, nitems-1
                    call readi(hs_occ_list0(i))
                end do
            case('NO_RENORM')
                no_renorm = .true.
            case('SELECT_REFERENCE_DET')
                select_ref_det_every_nreports = 20
                if (item /= nitems) call readi(select_ref_det_every_nreports)
                if (item /= nitems) call readf(ref_det_factor)
            case('ATTEMPT_SPAWN_PROB')
                call readf(pattempt_single)
                call readf(pattempt_double)

            ! Calculation options: CCMC.
            case('move_freq')
                call readi(ccmc_move_freq)

            ! use a negative number to indicate that the restart numbers have
            ! been fixed.
            case('RESTART')
                restart = .true.
                if (item /= nitems) then
                    call readi(restart_info_global%read_id)
                    restart_info_global%read_id = -restart_info_global%read_id -1
                end if
            case('DUMP_RESTART')
                if (item /= nitems) then
                    call readu(w)
                    ! Do we want to dump a restart file, when the shift turns on.
                    if(w == 'SHIFT') then
                        dump_restart_file_shift = .true.
                        ! Do we have a restart number for when the shift turns on.
                        if (item /= nitems) then
                            call readi(restart_info_global_shift%write_id)
                            restart_info_global_shift%write_id = -restart_info_global_shift%write_id-1
                        end if
                    ! Otherwise we have read a restart number.
                    else
                        call reread(0)
                        call readi(restart_info_global%write_id)
                        restart_info_global%write_id = -restart_info_global%write_id-1
                        dump_restart_file = .true.
                    end if
                else
                dump_restart_file = .true.
                end if
                
                ! If semi-stochastic is being used then a semi-stoch file will
                ! automatically be dumped when using this option.
                write_determ_space = .true.
            case('DUMP_RESTART_FREQUENCY')
                call readi(restart_info_global%write_restart_freq)
            case('SEED')
                call readi(seed)
            case('SHIFT_DAMPING')
                call readf(shift_damping)
            case('CLUSTER_MULTISPAWN_THRESHOLD')
                call readf(cluster_multispawn_threshold)
            case('INIT_SPIN_INVERSE_REFERENCE_DET')
                init_spin_inv_D0 = .true.

            ! Calculation options: initiator-fciqmc.
            case('INITIATOR_POPULATION')
                call readf(initiator_population)

            ! Calculation options: operators sampled using Hellmann--Feynman.
            case('OPERATOR')
                call readu(w)
                select case(w)
                case('HAMILTONIAN')
                    hf_operator = hamiltonian_operator
                case('KINETIC')
                    hf_operator = kinetic_operator
                case('DOUBLE_OCCUPANCY')
                    hf_operator = double_occ_operator
                case('DIPOLE')
                    hf_operator = dipole_operator
                end select
            ! Integral file for dipole moment.
            case('ALPHA0')
                call readi(alpha0)

            ! Output information.
            case('HAMIL','HAMILTONIAN')
                write_hamiltonian = .true.
                if (item /= nitems) call reada(hamiltonian_file)
            case('DET','DETERMINANTS')
                write_determinants = .true.
                if (item /= nitems) call reada(determinant_file)

            ! Parameters for parallel calculations.
            case('BLOCK_SIZE')
                call readi(block_size)
            case('NON_BLOCKING_COMM')
                non_blocking_comm = .true.
            case('LOAD_BALANCING')
                doing_load_balancing = .true.
            case('LOAD_BALANCING_SLOTS')
                call readi(par_info%load%nslots)
            case('LOAD_BALANCING_POP')
                call readli(par_info%load%pop)
            case('PERCENT_IMBAL')
                call readf(par_info%load%percent)
            case('MAX_LOAD_ATTEMPTS')
                call readi(par_info%load%max_attempts)
            case('WRITE_LOAD_INFO')
                par_info%load%write_info = .true.
            case('USE_MPI_BARRIERS')
                use_mpi_barriers = .true.

            case('FINITE_CLUSTER')
                ! this will be checked in check_input to ensure that it
                ! is only used when we are formulating the calculation
                ! in real-space
                sys%real_lattice%finite_cluster = .true.
            case('TRIANGULAR_LATTICE')
                sys%lattice%triangular_lattice = .true.

            case('NEEL_SINGLET_ESTIMATOR')
                trial_function = neel_singlet
            case('NEEL_SINGLET_GUIDING')
                guiding_function = neel_singlet_guiding
                trial_function = neel_singlet

            case('END')
                exit
            case default
                call report('Keyword '//trim(w)//' not recognized.', .true.)
            end select
        end do ! end reading of input.

        close(ir, status='keep')
        if (ios.gt.0) call stop_all('read_input','Problem reading input.')

    end subroutine read_input

    subroutine check_input(sys)

        ! I don't pretend this is the most comprehensive of tests, but at least
        ! make sure a few things are not completely insane.

        ! In/Out:
        !    sys: system object, as set in read_input (invalid settings are overridden).

        use const
        use system

        type(sys_t), intent(inout) :: sys

        integer :: ivec, jvec
        character(*), parameter :: this='check_input'

        if (sys%system /= heisenberg) then
            if (sys%nel <= 0) call stop_all(this,'Number of electrons must be positive.')
            if (trial_function /= single_basis) call stop_all(this, 'Only a single determinant can be used as the reference&
                                                     & state for this system. Other trial functions are not avaliable.')
            if (guiding_function /= no_guiding) &
                call stop_all(this, 'Importance sampling is only avaliable for the Heisenberg model&
                                         & currently.')
            if (all_sym_sectors) call stop_all(this,'The option to use all symmetry sectors at the same time is only&
                                         & available for the Heisenberg model.')
        end if

        if (sys%system == read_in) then

            if (analyse_fci_wfn /= 0 .and. sys%read_in%dipole_int_file == '') then
                if (parent) call warning(this, 'Cannot analyse FCI wavefunction without a dipole &
                             &integrals file.  Turning analyse_fci_wfn option off...')
                analyse_fci_wfn = 0
            end if

        else

            if (sys%system /= ueg) then
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
                if (ms_in > sys%lattice%nsites .and. (.not. all_sym_sectors)) call stop_all(this,'Value of Ms given is&
                                                                             & too large for this lattice.')
                if ((-ms_in) > sys%lattice%nsites) call stop_all(this,'Value of Ms given is too small for this lattice.')
                if (mod(abs(ms_in),2) /=  mod(sys%lattice%nsites,2)) call stop_all(this, 'Ms value specified is not&
                                                                              & possible for this lattice.')
                if (sys%heisenberg%staggered_magnetic_field /= 0.0_p .and. (.not.sys%lattice%bipartite_lattice)) &
                    call stop_all(this, 'Cannot set a staggered field&
                                        & for this lattice because it is frustrated.')
                if (sys%heisenberg%staggered_magnetic_field /= 0.0_p .and. sys%heisenberg%magnetic_field /= 0.0_p) &
                    call stop_all(this, 'Cannot set a uniform and a staggered field at the same time.')
                if ((guiding_function==neel_singlet_guiding) .and. trial_function /= neel_singlet) call stop_all(this, 'This &
                                                         &guiding function is only avaliable when using the Neel singlet state &
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

        end if

        ! Real amplitude checks.
        if (real_amplitudes) then
            if (doing_calc(ct_fciqmc_calc) .or. doing_calc(hfs_fciqmc_calc)) then
                call stop_all(this, 'The real_amplitudes option is not implemented with the method you have requested.')
            end if
            if (bit_size(0_int_p) == 32 .and. parent) then
                call warning(this,'You are using 32-bit walker populations with real amplitudes. The maximum &
                             &population size on a given determinant is 2^20=1048576. Errors will occur if this is &
                             &exceeded. Compile HANDE with the CPPFLAG -DPOP_SIZE=64 to use 64-bit populations.')
            end if
        end if

        ! Semi-stochastic checks.
        if (semi_stoch_start_iter /= 0 .and. determ_space_type == empty_determ_space .and. parent) &
            call warning(this,'You have specified an iteration to turn semi-stochastic on but have not &
                         &specified a deterministic space to use.')
        if (determ_space_type /= empty_determ_space .and. (doing_calc(dmqmc_calc) .or. doing_calc(ct_fciqmc_calc) .or. &
              doing_calc(hfs_fciqmc_calc))) &
              call stop_all(this, 'Semi-stochastic is only implemented with the FCIQMC method.')

        if (init_spin_inv_D0 .and. ms_in /= 0) then
            if (parent) call warning(this, 'Flipping the reference state will give &
                                            &a state which has a different value of Ms and so cannot be used here.')
            init_spin_inv_D0 = .false.
        end if

        if (allocated(correlation_sites) .and. size(correlation_sites) /= 2) call stop_all(this, 'You must enter exactly two &
               &sites for the correlation function option.')

          if (dmqmc_find_weights .and. calculate_excit_distribution) call stop_all(this, 'DMQMC_FIND_WEIGHTS and OUTPUT_EXCITATION&
              &_DISTRIBUTION options cannot be used together.')

        ! Calculation specific checking.
        if (doing_calc(lanczos_diag)) then
            if (lanczos_string_len <= 0) call stop_all(this,'Lanczos basis not positive.')
            if (nlanczos_eigv <= 0) call stop_all(this,'# lanczos eigenvalues not positive.')
        end if

        if (.not.doing_calc(dmqmc_calc) .and. dmqmc_calc_type /= 0 .and. parent) call warning(this,&
               'You are not performing a DMQMC calculation but have requested DMQMC options to be calculated.')
        if (doing_calc(fciqmc_calc)) then
            if (.not.doing_calc(simple_fciqmc_calc)) then
                if (walker_length == 0) call stop_all(this,'Walker length zero.')
                if (spawned_walker_length == 0) call stop_all(this,'Spawned walker length zero.')
            end if
            if (calc_inst_rdm .and. spawned_rdm_length == 0) call stop_all(this,'Spawned RDM length zero.')
            if (tau <= 0) call stop_all(this,'Tau not positive.')
            if (shift_damping <= 0) call stop_all(this,'Shift damping not positive.')
            if (allocated(occ_list0)) then
                if (size(occ_list0) /= sys%nel) then
                    if (sys%system /= heisenberg) then
                        call stop_all(this,'Number of electrons specified is different from &
                        &number of electrons used in the reference determinant.')
                    end if
                end if
            end if
            if (par_info%load%nslots < 0) call stop_all(this, 'Number of slots for load balancing is not positive.')
            if (par_info%load%pop < 0) call stop_all(this, 'Load balancing population must be positive.')
            if (par_info%load%percent < 0 .or. par_info%load%percent > 1.0) &
                call stop_all(this, 'Percentage imbalance must be positive and less that 1.')
            if (par_info%load%max_attempts < 0) call stop_all(this, 'Maximum number of load balancing attempts must be positive')
        end if
        if (doing_calc(ct_fciqmc_calc)) ncycles = 1

        if (doing_dmqmc_calc(dmqmc_rdm_r2) .and. (.not. replica_tricks)) call stop_all(this,&
                    'The replica_tricks option must be used in order to calculate the Renyi-2 entropy.')
        if (doing_dmqmc_calc(dmqmc_rdm_r2) .and. (.not. calc_inst_rdm)) call stop_all(this,&
                    'The instantaneous_rdm option must be used in order to calculate the Renyi-2 entropy.')

        ! If the FINITE_CLUSTER keyword was detected then make sure that
        ! we are doing a calculation in real-space. If we're not then
        ! unset finite cluster,tell the user and carry on
        if(sys%momentum_space) then
            if (sys%real_lattice%finite_cluster .and. parent) call warning(this,'FINITE_CLUSTER keyword only valid for hubbard&
                                      & calculations in real-space: ignoring keyword')
            sys%real_lattice%finite_cluster = .false.
        end if

        if (all_sym_sectors) then
            if (.not. doing_calc(dmqmc_calc)) call stop_all(this, 'The use_all_sym_sectors option can only be used in&
                                                                   & DMQMC calculations.')
            if (abs(sys%heisenberg%magnetic_field) > depsilon .or. &
                abs(sys%heisenberg%staggered_magnetic_field) > depsilon) &
                    call stop_all(this, 'The use_all_sym_sectors option cannot be used with magnetic fields.')
            if (calc_ground_rdm) call stop_all(this, 'The use_all_sym_sectors and ground_state_rdm options cannot be&
                                                      & used together.')
        end if

        if (dump_restart_file_shift) then
             if (restart_info_global_shift%write_id<0 .and. restart_info_global%write_id<0 &
                 .and. restart_info_global%write_id == restart_info_global_shift%write_id) &
                 call stop_all(this, 'The ids of the restart files are the same.')
             if (restart_info_global_shift%write_id<0 .and. restart_info_global%write_restart_freq /= huge(0) )&
                 call stop_all(this, 'The ids of the restart files could be the same')
        end if   
        if (dmqmc_vary_weights .and. (.not. dmqmc_weighted_sampling)) then
            call stop_all(this, 'The dmqmc_vary_weights option can only be used together with the dmqmc_weighted_sampling option.')
        end if
        if (sys%system /= heisenberg .and. dmqmc_calc_type > dmqmc_energy) then
            call stop_all(this, 'The observable requested is not currently implemented for this Hamiltonian.')
        end if

        if (parent) write (6,'(/,1X,13("-"),/)')

    end subroutine check_input

    subroutine distribute_input(sys)

        ! Distribute the data read in by the parent processor to all other
        ! processors.

        ! Completely empty (courtesy of C pre-processing) when compiled in
        ! serial.

        ! In/Out:
        !    sys: object describing the system.  All parameters which can be set
        !       in the input file are distributed to other processors.

#ifndef PARALLEL

        use system, only: sys_t

        type(sys_t), intent(inout) :: sys

#else

        use mpi
        use parallel
        use checking, only: check_allocate

        use system

        type(sys_t), intent(inout) :: sys

        integer :: i, ierr, occ_list_size
        logical :: option_set

        call mpi_bcast(sys%system, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(sys%read_in%fcidump, len(sys%read_in%fcidump), mpi_character, 0, mpi_comm_world, ierr)
        call mpi_bcast(sym_in, 1, mpi_integer, 0, mpi_comm_world, ierr)

        call mpi_bcast(sys%lattice%ndim, 1, mpi_integer, 0, mpi_comm_world, ierr)
        if (parent) option_set = allocated(sys%lattice%lattice)
        call mpi_bcast(option_set, 1, mpi_logical, 0, mpi_comm_world, ierr)
        if (option_set) then
            if (.not.parent) then
                allocate(sys%lattice%lattice(sys%lattice%ndim,sys%lattice%ndim), stat=ierr)
                call check_allocate('sys%lattice%lattice',sys%lattice%ndim*sys%lattice%ndim,ierr)
            end if
            call mpi_bcast(sys%lattice%lattice, sys%lattice%ndim*sys%lattice%ndim, mpi_integer, 0, mpi_comm_world, ierr)
        end if
        call mpi_bcast(real_amplitudes, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(spawn_cutoff, 1, mpi_preal, 0, mpi_comm_world, ierr)
        call mpi_bcast(semi_stoch_start_iter, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(semi_stoch_shift_iter, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(determ_space_type, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(determ_target_size, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(write_determ_space, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(separate_determ_annihil, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(ccmc_full_nc, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(linked_ccmc, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(replica_tricks, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(sys%real_lattice%finite_cluster, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(sys%lattice%triangular_lattice, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(trial_function, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(guiding_function, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(sys%nel, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(sys%hubbard%t, 1, mpi_preal, 0, mpi_comm_world, ierr)
        call mpi_bcast(sys%hubbard%u, 1, mpi_preal, 0, mpi_comm_world, ierr)
        call mpi_bcast(sys%heisenberg%J, 1, mpi_preal, 0, mpi_comm_world, ierr)
        call mpi_bcast(sys%heisenberg%magnetic_field, 1, mpi_preal, 0, mpi_comm_world, ierr)
        call mpi_bcast(sys%heisenberg%staggered_magnetic_field, 1, mpi_preal, 0, mpi_comm_world, ierr)
        if (parent) option_set = allocated(sys%k_lattice%ktwist)
        call mpi_bcast(option_set, 1, mpi_logical, 0, mpi_comm_world, ierr)
        if (option_set) then
            if (.not.parent) then
                allocate(sys%k_lattice%ktwist(sys%lattice%ndim), stat=ierr)
                call check_allocate('sys%k_lattice%ktwist',sys%lattice%ndim,ierr)
            end if
            call mpi_bcast(sys%k_lattice%ktwist, sys%lattice%ndim, mpi_preal, 0, mpi_comm_world, ierr)
        end if
        call mpi_bcast(sys%ueg%r_s, 1, mpi_preal, 0, mpi_comm_world, ierr)
        call mpi_bcast(sys%ueg%ecutoff, 1, mpi_preal, 0, mpi_comm_world, ierr)
        call mpi_bcast(sys%read_in%dipole_int_file, len(sys%read_in%dipole_int_file), mpi_character, 0, mpi_comm_world, ierr)
        call mpi_bcast(select_ref_det_every_nreports, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(ref_det_factor, 1, mpi_preal, 0, mpi_comm_world, ierr)
        call mpi_bcast(sys%cas, 2, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(ras, 2, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(ras3_max, 1, mpi_integer, 0, mpi_comm_world, ierr)

        call mpi_bcast(ms_in, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(sym_in, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(sys%read_in%useLz, 1, mpi_logical, 0, mpi_comm_world, ierr)

        call mpi_bcast(calc_type, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(dmqmc_calc_type, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(direct_lanczos, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(nhilbert_cycles, 1, mpi_integer, 0, mpi_comm_world, ierr)

        call mpi_bcast(lanczos_string_len, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(nlanczos_eigv, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(use_sparse_hamil, 1, mpi_logical, 0, mpi_comm_world, ierr)

        call mpi_bcast(print_fci_wfn, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(analyse_fci_wfn, 1, mpi_integer, 0, mpi_comm_world, ierr)

        call mpi_bcast(ncycles, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(nreport, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(beta_loops, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(walker_length, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(spawned_walker_length, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(spawned_rdm_length, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(tau, 1, mpi_preal, 0, mpi_comm_world, ierr)
        call mpi_bcast(tau_search, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(initial_shift, 1, mpi_preal, 0, mpi_comm_world, ierr)
        call mpi_bcast(vary_shift_from, 1, mpi_preal, 0, mpi_comm_world, ierr)
        call mpi_bcast(vary_shift_from_proje, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(average_shift_until, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(target_particles, 1, mpi_integer8, 0, mpi_comm_world, ierr)
        call mpi_bcast(doing_reduced_dm, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(calc_ground_rdm, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(calc_inst_rdm, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(output_rdm, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(nrdms, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(doing_exact_rdm_eigv, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(doing_von_neumann_entropy, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(doing_concurrence, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(start_averaging, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(dmqmc_weighted_sampling, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(dmqmc_vary_weights, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(dmqmc_find_weights, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(all_sym_sectors, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(all_mom_sectors, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(finish_varying_weights, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(propagate_to_beta, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(free_electron_trial, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(init_beta, 1, mpi_preal, 0, mpi_comm_world, ierr)
        call mpi_bcast(half_density_matrix, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(calculate_excit_distribution, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(metropolis_attempts, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(chem_pot, 1, mpi_preal, 0, mpi_comm_world, ierr)
        call mpi_bcast(grand_canonical_ensemble, 1, mpi_logical, 0, mpi_comm_world, ierr)
        option_set = .false.
        if (parent) option_set = allocated(dmqmc_sampling_probs)
        call mpi_bcast(option_set, 1, mpi_logical, 0, mpi_comm_world, ierr)
        if (option_set) then
            occ_list_size = size(dmqmc_sampling_probs)
            call mpi_bcast(occ_list_size, 1, mpi_integer, 0, mpi_comm_world, ierr)
            if (.not.parent) then
                allocate(dmqmc_sampling_probs(occ_list_size), stat=ierr)
                call check_allocate('dmqmc_sampling_probs',occ_list_size,ierr)
            end if
            call mpi_bcast(dmqmc_sampling_probs, occ_list_size, mpi_preal, 0, mpi_comm_world, ierr)
        end if
        option_set = .false.
        if (parent) option_set = allocated(rdms)
        call mpi_bcast(option_set, 1, mpi_logical, 0, mpi_comm_world, ierr)
        if (option_set) then
            occ_list_size = size(rdms)
            call mpi_bcast(occ_list_size, 1, mpi_integer, 0, mpi_comm_world, ierr)
            if (.not.parent) then
                allocate(rdms(occ_list_size), stat=ierr)
                call check_allocate('rdms',occ_list_size,ierr)
            end if
            do i = 1, occ_list_size
                call mpi_bcast(rdms(i)%A_nsites, 1, mpi_integer, 0, mpi_comm_world, ierr)
                if (.not.parent) then
                    allocate(rdms(i)%subsystem_A(rdms(i)%A_nsites), stat=ierr)
                    call check_allocate('rdms(i)%subsystem_A',rdms(i)%A_nsites,ierr)
                end if
                call mpi_bcast(rdms(i)%subsystem_A, rdms(i)%A_nsites, mpi_integer, 0, mpi_comm_world, ierr)
            end do
        end if
        option_set = .false.
        if (parent) option_set = allocated(correlation_sites)
        call mpi_bcast(option_set, 1, mpi_logical, 0, mpi_comm_world, ierr)
        if (option_set) then
            occ_list_size = size(correlation_sites)
            call mpi_bcast(occ_list_size, 1, mpi_integer, 0, mpi_comm_world, ierr)
            if (.not.parent) then
                allocate(correlation_sites(occ_list_size), stat=ierr)
                call check_allocate('correlation_sites',occ_list_size,ierr)
            end if
            call mpi_bcast(correlation_sites, occ_list_size, mpi_integer, 0, mpi_comm_world, ierr)
        end if
        option_set = .false.
        call mpi_bcast(truncate_space, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(truncation_level, 1, mpi_integer, 0, mpi_comm_world, ierr)
        if (parent) option_set = allocated(occ_list0)
        call mpi_bcast(option_set, 1, mpi_logical, 0, mpi_comm_world, ierr)
        if (option_set) then
            ! Have not yet set sys%nel in the Heisenberg model.
            occ_list_size = size(occ_list0)
            call mpi_bcast(occ_list_size, 1, mpi_integer, 0, mpi_comm_world, ierr)
            if (.not.parent) then
                allocate(occ_list0(occ_list_size), stat=ierr)
                call check_allocate('occ_list0', occ_list_size, ierr)
            end if
            call mpi_bcast(occ_list0, occ_list_size, mpi_integer, 0, mpi_comm_world, ierr)
        end if
        if (parent) option_set = allocated(hs_occ_list0)
        call mpi_bcast(option_set, 1, mpi_logical, 0, mpi_comm_world, ierr)
        if (option_set) then
            ! Have not yet set sys%nel in the Heisenberg model.
            occ_list_size = size(hs_occ_list0)
            call mpi_bcast(occ_list_size, 1, mpi_integer, 0, mpi_comm_world, ierr)
            if (.not.parent) then
                allocate(hs_occ_list0(occ_list_size), stat=ierr)
                call check_allocate('hs_occ_list0', occ_list_size, ierr)
            end if
            call mpi_bcast(hs_occ_list0, occ_list_size, mpi_integer, 0, mpi_comm_world, ierr)
        end if
        call mpi_bcast(restart, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(dump_restart_file, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(dump_restart_file_shift, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(restart_info_global%read_id, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(restart_info_global%write_id, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(restart_info_global%write_restart_freq, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(restart_info_global_shift%write_id, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(seed, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(shift_damping, 1, mpi_preal, 0, mpi_comm_world, ierr)
        call mpi_bcast(cluster_multispawn_threshold, 1, mpi_preal, 0, mpi_comm_world, ierr)
        call mpi_bcast(D0_population, 1, mpi_preal, 0, mpi_comm_world, ierr)
        call mpi_bcast(no_renorm, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(pattempt_single, 1, mpi_preal, 0, mpi_comm_world, ierr)
        call mpi_bcast(pattempt_double, 1, mpi_preal, 0, mpi_comm_world, ierr)

        call mpi_bcast(ccmc_move_freq, 1, mpi_integer, 0, mpi_comm_world, ierr)

        call mpi_bcast(init_spin_inv_D0, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(initiator_population, 1, mpi_preal, 0, mpi_comm_world, ierr)
        call mpi_bcast(initiator_approximation, 1, mpi_logical, 0, mpi_comm_world, ierr)

        call mpi_bcast(hf_operator, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(alpha0, 1, mpi_integer, 0, mpi_comm_world, ierr)

        call mpi_bcast(write_hamiltonian, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(write_determinants, 1, mpi_logical, 0, mpi_comm_world, ierr)

        call mpi_bcast(block_size, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(non_blocking_comm, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(doing_load_balancing, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(par_info%load%nslots, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(par_info%load%pop, 1, mpi_integer8, 0, mpi_comm_world, ierr)
        call mpi_bcast(par_info%load%percent, 1, mpi_preal, 0, mpi_comm_world, ierr)
        call mpi_bcast(par_info%load%max_attempts, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(par_info%load%write_info, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(use_mpi_barriers, 1, mpi_logical, 0, mpi_comm_world, ierr)

#endif

    end subroutine distribute_input

end module parse_input

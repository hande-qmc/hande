module parse_input
! Parse input options and check input for validity.

use const

use parallel, only: parent, block_size
use errors
use hilbert_space
use canonical_kinetic_energy, only: nkinetic_cycles
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

    subroutine read_input(sys, qmc_in, fciqmc_in, ccmc_in, semi_stoch_in, restart_in, reference, load_bal_in, dmqmc_in)

        ! Read input options from a file (if specified on the command line) or via
        ! STDIN.

        ! In/Out:
        !    sys: system being studied.  Parameters specified in the input file
        !         are set directly in the system object, components which are not
        !         mentioned in the input file are not altered.
        !    qmc_in: input options relating to QMC methods.
        !    fciqmc_in: input options relating to FCIQMC.
        !    ccmc_in: input options relating to CCMC.
        !    semi_stoch_in: Input options for the semi-stochastic adaptation.
        !    restart_in: input options for HDF5 restart files.
        !    load_bal_in: input options for load balancing.
        !    dmqmc_in: input options relating to DMQMC.

! nag doesn't automatically bring in command-line option handling.
#ifdef NAGF95
        use f90_unix_env
#endif

        use qmc_data, only: qmc_in_t, fciqmc_in_t, ccmc_in_t, semi_stoch_in_t
        use qmc_data, only: restart_in_t, reference_t, load_bal_in_t
        use qmc_data, only: read_determ_space, high_pop_determ_space
        use dmqmc_data, only: dmqmc_in_t
        use system

        use input
        use utils, only: get_free_unit
        use checking, only: check_allocate

#ifdef NAGF95
        use f90_unix_env, ONLY: getarg,iargc
#else
        integer :: iargc ! External function.
#endif

        type(sys_t), intent(inout) :: sys
        type(qmc_in_t), intent(inout) :: qmc_in
        type(fciqmc_in_t), intent(inout) :: fciqmc_in
        type(ccmc_in_t), intent(inout) :: ccmc_in
        type(semi_stoch_in_t), intent(inout) :: semi_stoch_in
        type(restart_in_t), intent(inout) :: restart_in
        type(reference_t), intent(inout) :: reference
        type(load_bal_in_t), intent(inout) :: load_bal_in
        type(dmqmc_in_t), intent(inout) :: dmqmc_in

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
                qmc_in%initiator_approx = .true.
            case('CT_FCIQMC')
                calc_type = calc_type + ct_fciqmc_calc
            case('DMQMC')
                calc_type = calc_type + dmqmc_calc
            case('CCMC')
                calc_type = calc_type + ccmc_calc
            case('ICCMC')
                calc_type = calc_type + ccmc_calc
                qmc_in%initiator_approx = .true.
            case('HELLMANN-FEYNMAN')
                calc_type = calc_type + hfs_fciqmc_calc
            case('ESTIMATE_HILBERT_SPACE')
                calc_type = calc_type + mc_hilbert_space
                call readi(nhilbert_cycles)
            case('ESTIMATE_CANONICAL_KINETIC_ENERGY')
                calc_type = calc_type + mc_canonical_kinetic_energy
            case('NUM_KINETIC_CYCLES')
                call readi(nkinetic_cycles)

            case('REPLICA_TRICKS')
                dmqmc_in%replica_tricks = .true.
            case('PROPAGATE_TO_BETA')
                dmqmc_in%propagate_to_beta = .true.
            case('INIT_BETA')
                call readf(dmqmc_in%init_beta)
            case('METROPOLIS_ATTEMPTS')
                call readi(dmqmc_in%metropolis_attempts)
            case('MAX_METROPOLIS_MOVE')
                call readi(dmqmc_in%max_metropolis_move)
            case('FREE_ELECTRON_TRIAL')
                dmqmc_in%free_electron_trial = .true.
            case('CHEM_POT')
                call readf(sys%ueg%chem_pot)
            case('GRAND_CANONICAL_INITIALISATION')
                dmqmc_in%grand_canonical_initialisation = .true.
            case('FERMI_TEMPERATURE')
                dmqmc_in%fermi_temperature = .true.

            case('CCMC_FULL_NC')
                ccmc_in%full_nc = .true.
            case('CCMC_LINKED')
                ccmc_in%linked = .true.

            case('REAL_AMPLITUDES')
                qmc_in%real_amplitudes = .true.
            case('SPAWN_CUTOFF')
                call readf(qmc_in%spawn_cutoff)

            ! Semi-stochastic options.
            case('SEMI_STOCH_ITERATION')
                call readi(semi_stoch_in%start_iter)
                qmc_in%real_amplitudes = .true.
            case('SEMI_STOCH_SHIFT_START')
                call readi(semi_stoch_in%shift_iter)
                semi_stoch_in%start_iter = -1
                qmc_in%real_amplitudes = .true.
            ! Deterministic spaces.
            case('SEMI_STOCH_HIGH_POP')
                semi_stoch_in%determ_space_type = high_pop_determ_space
                call readi(semi_stoch_in%target_size)
            case('SEMI_STOCH_READ')
                semi_stoch_in%determ_space_type = read_determ_space
                ! Not needed.
                semi_stoch_in%target_size = -1
            case('WRITE_DETERM_SPACE')
                semi_stoch_in%write_determ_space = .true.
            case('SEMI_STOCH_COMBINE_ANNIHIL')
                semi_stoch_in%separate_annihil = .false.

            ! DMQMC expectation values to be calculated.
            case('DMQMC_FULL_RENYI_2')
                dmqmc_calc_type = dmqmc_calc_type + dmqmc_full_r2
            case('DMQMC_ENERGY')
                dmqmc_calc_type = dmqmc_calc_type + dmqmc_energy
            case('DMQMC_ENERGY_SQUARED')
                dmqmc_calc_type = dmqmc_calc_type + dmqmc_energy_squared
            case('DMQMC_CORRELATION_FUNCTION')
                dmqmc_calc_type = dmqmc_calc_type + dmqmc_correlation
                allocate(dmqmc_in%correlation_sites(nitems-1), stat=ierr)
                call check_allocate('dmqmc_in%correlation_sites',nitems-1,ierr)
                do i = 1, nitems-1
                    call readi(dmqmc_in%correlation_sites(i))
                end do
            case('DMQMC_STAGGERED_MAGNETISATION')
                dmqmc_calc_type = dmqmc_calc_type + dmqmc_staggered_magnetisation
            case('DMQMC_WEIGHTED_SAMPLING')
                call readi(nweights)
                allocate(dmqmc_in%sampling_probs(nweights), stat=ierr)
                call check_allocate('dmqmc_in%sampling_probs', nweights, ierr)
                dmqmc_in%weighted_sampling = .true.
                call read_line(eof)
                if (eof) call stop_all('read_input', 'Unexpected end of file reading DMQMC weights.')
                do i = 1, nweights
                    call readf(dmqmc_in%sampling_probs(i))
                end do
            case('DMQMC_VARY_WEIGHTS')
                call readi(dmqmc_in%finish_varying_weights)
                dmqmc_in%vary_weights = .true.
            case('DMQMC_FIND_WEIGHTS')
                dmqmc_in%find_weights = .true.
                dmqmc_in%weighted_sampling = .true.
            case('OUTPUT_EXCITATION_DISTRIBUTION')
                dmqmc_in%calc_excit_dist = .true.
            case('USE_ALL_SYM_SECTORS')
                dmqmc_in%all_sym_sectors = .true.
            case('USE_ALL_SPIN_SECTORS')
                dmqmc_in%all_spin_sectors = .true.
            case('REDUCED_DENSITY_MATRIX')
                call readi(dmqmc_in%rdm%nrdms)
                allocate(rdms(dmqmc_in%rdm%nrdms), stat=ierr)
                call check_allocate('rdms', dmqmc_in%rdm%nrdms, ierr)
                dmqmc_in%rdm%doing_rdm = .true.
                do i = 1, dmqmc_in%rdm%nrdms
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
                dmqmc_in%rdm%calc_ground_rdm = .true.
            case('INSTANTANEOUS_RDM')
                calc_inst_rdm = .true.
            case('OUTPUT_RDM')
                output_rdm = .true.
            case('EXACT_RDM_EIGENVALUES')
                doing_exact_rdm_eigv = .true.
                dmqmc_in%rdm%calc_ground_rdm = .true.
            case('CONCURRENCE')
                doing_concurrence = .true.
            case('VON_NEUMANN_ENTROPY')
                doing_vn_entropy = .true.
            case('RENYI_ENTROPY_2')
                dmqmc_calc_type = dmqmc_calc_type + dmqmc_rdm_r2
            case('START_AVERAGING_EXCITATION_DIST')
                call readi(dmqmc_in%start_av_excit_dist)
            case('START_AVERAGING_RDM')
                call readi(dmqmc_in%start_av_rdm)
            ! calculation options: DMQMC
            case('TRUNCATION_LEVEL')
                truncate_space = .true.
                call readi(truncation_level)
            case('HALF_DENSITY_MATRIX')
                dmqmc_in%half_density_matrix = .true.

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
                call readi(qmc_in%ncycles)
            case('NREPORTS')
                call readi(qmc_in%nreport)
                if (qmc_in%nreport < 0) qmc_in%nreport = huge(qmc_in%nreport)
            case('BETA_LOOPS')
                call readi(dmqmc_in%beta_loops)
            case('WALKER_LENGTH')
                call readi(qmc_in%walker_length)
                if (item /= nitems) then
                    call readu(w)
                    if (w == 'MB') then
                        qmc_in%walker_length = -qmc_in%walker_length
                    else
                        call report('Keyword '//trim(w)//' not recognized.', .true.)
                    end if
                end if
            case('SPAWNED_WALKER_LENGTH')
                call readi(qmc_in%spawned_walker_length)
                if (item /= nitems) then
                    call readu(w)
                    if (w == 'MB') then
                        qmc_in%spawned_walker_length = -qmc_in%spawned_walker_length
                    else
                        call report('Keyword '//trim(w)//' not recognized.', .true.)
                    end if
                end if
            case('SPAWNED_RDM_LENGTH')
                call readi(dmqmc_in%rdm%spawned_length)
                if (item /= nitems) then
                    call readu(w)
                    if (w == 'MB') then
                        dmqmc_in%rdm%spawned_length = -dmqmc_in%rdm%spawned_length
                    else
                        call report('Keyword '//trim(w)//' not recognized.', .true.)
                    end if
                end if
            case('TAU')
                call readf(qmc_in%tau)
            case('TAU_SEARCH')
                qmc_in%tau_search = .true.
            case('INITIAL_SHIFT')
                call readf(qmc_in%initial_shift)
                ! We assume the user is sensible/knows what he/she is doing if
                ! initial_shift and vary_shift_from are set.
                qmc_in%vary_shift_from = qmc_in%initial_shift
            case('VARY_SHIFT_FROM')
                call readu(w)
                if (w == 'PROJE') then
                    qmc_in%vary_shift_from_proje = .true.
                else
                    call reread(0)
                    call readf(qmc_in%vary_shift_from)
                end if
            case('VARYSHIFT_TARGET')
                call readf(qmc_in%target_particles)
            case('INIT_POP')
                call readf(qmc_in%D0_population)
            case('REFERENCE_DET')
                allocate(reference%occ_list0(nitems-1), stat=ierr)
                call check_allocate('reference%occ_list0',nitems-1,ierr)
                do i = 1, nitems-1
                    call readi(reference%occ_list0(i))
                end do
            case('HS_REFERENCE_DET')
                allocate(reference%hs_occ_list0(nitems-1), stat=ierr)
                call check_allocate('reference%hs_occ_list0',nitems-1,ierr)
                do i = 1, nitems-1
                    call readi(reference%hs_occ_list0(i))
                end do
            case('NO_RENORM')
                qmc_in%no_renorm = .true.
            case('SELECT_REFERENCE_DET')
                fciqmc_in%select_ref_det_every_nreports = 20
                if (item /= nitems) call readi(fciqmc_in%select_ref_det_every_nreports)
                if (item /= nitems) call readf(fciqmc_in%ref_det_factor)
            case('ATTEMPT_SPAWN_PROB')
                call readf(qmc_in%pattempt_single)
                call readf(qmc_in%pattempt_double)

            ! Calculation options: CCMC.
            case('move_freq')
                call readi(ccmc_in%move_freq)

            ! use a negative number to indicate that the restart numbers have
            ! been fixed.
            case('RESTART')
                restart_in%read_restart = .true.
                if (item /= nitems) then
                    call readi(restart_info_global%read_id)
                    restart_info_global%read_id = -restart_info_global%read_id -1
                end if
            case('DUMP_RESTART')
                if (item /= nitems) then
                    call readu(w)
                    ! Do we want to dump a restart file, when the shift turns on.
                    if(w == 'SHIFT') then
                        restart_in%dump_restart_file_shift = .true.
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
                        restart_in%dump_restart = .true.
                    end if
                else
                restart_in%dump_restart = .true.
                end if
                
                ! If semi-stochastic is being used then a semi-stoch file will
                ! automatically be dumped when using this option.
                semi_stoch_in%write_determ_space = .true.
            case('DUMP_RESTART_FREQUENCY')
                call readi(restart_info_global%write_restart_freq)
            case('SEED')
                call readi(qmc_in%seed)
            case('SHIFT_DAMPING')
                call readf(qmc_in%shift_damping)
            case('CLUSTER_MULTISPAWN_THRESHOLD')
                call readf(ccmc_in%cluster_multispawn_threshold)
            case('INIT_SPIN_INVERSE_REFERENCE_DET')
                fciqmc_in%init_spin_inv_D0 = .true.

            ! Calculation options: initiator-fciqmc.
            case('INITIATOR_POPULATION')
                call readf(qmc_in%initiator_pop)

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
                fciqmc_in%non_blocking_comm = .true.
            case('LOAD_BALANCING')
                fciqmc_in%doing_load_balancing = .true.
            case('LOAD_BALANCING_SLOTS')
                call readi(load_bal_in%nslots)
            case('LOAD_BALANCING_POP')
                call readli(load_bal_in%pop)
            case('PERCENT_IMBAL')
                call readf(load_bal_in%percent)
            case('MAX_LOAD_ATTEMPTS')
                call readi(load_bal_in%max_attempts)
            case('WRITE_LOAD_INFO')
                load_bal_in%write_info = .true.
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

    subroutine check_input(sys, qmc_in, fciqmc_in, ccmc_in, semi_stoch_in, restart_in, reference, load_bal_in, dmqmc_in)

        ! I don't pretend this is the most comprehensive of tests, but at least
        ! make sure a few things are not completely insane.

        ! In/Out:
        !    sys: system object, as set in read_input (invalid settings are overridden).
        !    qmc_in: input options relating to QMC methods.
        !    fciqmc_in: input options relating to FCIQMC.
        !    ccmc_in: input options relating to CCMC.
        !    load_bal_in: input options for load balancing.
        ! In:
        !    semi_stoch_in: Input options for the semi-stochastic adaptation.
        !    restart_in: input options for HDF5 restart files.
        !    reference: reference determinant.
        !    dmqmc_in: input options relating to DMQMC.

        use const
        use qmc_data, only: qmc_in_t, fciqmc_in_t, ccmc_in_t, semi_stoch_in_t
        use qmc_data, only: restart_in_t, reference_t, load_bal_in_t, empty_determ_space
        use dmqmc_data, only: dmqmc_in_t
        use system

        type(sys_t), intent(inout) :: sys
        type(qmc_in_t), intent(inout) :: qmc_in
        type(fciqmc_in_t), intent(inout) :: fciqmc_in
        type(ccmc_in_t), intent(inout) :: ccmc_in
        type(semi_stoch_in_t), intent(in) :: semi_stoch_in
        type(restart_in_t), intent(in) :: restart_in
        type(reference_t), intent(in) :: reference
        type(load_bal_in_t), intent(in) :: load_bal_in
        type(dmqmc_in_t), intent(in) :: dmqmc_in

        integer :: ivec, jvec
        character(*), parameter :: this='check_input'

        if (sys%system /= heisenberg) then
            if (sys%nel <= 0) call stop_all(this,'Number of electrons must be positive.')
            if (trial_function /= single_basis) call stop_all(this, 'Only a single determinant can be used as the reference&
                                                     & state for this system. Other trial functions are not avaliable.')
            if (guiding_function /= no_guiding) &
                call stop_all(this, 'Importance sampling is only avaliable for the Heisenberg model&
                                         & currently.')
            if (dmqmc_in%all_spin_sectors) call stop_all(this,'The option to use all symmetry sectors at the same time is only&
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
                if (doing_calc(mc_canonical_kinetic_energy)) call stop_all(this, 'estimate_canonical_kinetic_energy&
                                                                           & only implemented for the UEG.')
            end if

            if (sys%system == heisenberg) then
                if (ms_in > sys%lattice%nsites .and. (.not. dmqmc_in%all_spin_sectors)) call stop_all(this,'Value of Ms given is&
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
        if (qmc_in%real_amplitudes) then
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
        if (semi_stoch_in%start_iter /= 0 .and. semi_stoch_in%determ_space_type == empty_determ_space .and. parent) &
            call warning(this,'You have specified an iteration to turn semi-stochastic on but have not &
                         &specified a deterministic space to use.')
        if (semi_stoch_in%determ_space_type /= empty_determ_space .and. (doing_calc(dmqmc_calc) .or. &
                                   doing_calc(ct_fciqmc_calc) .or. doing_calc(hfs_fciqmc_calc))) &
              call stop_all(this, 'Semi-stochastic is only implemented with the FCIQMC method.')

        if (fciqmc_in%init_spin_inv_D0 .and. ms_in /= 0) then
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
            if (lanczos_string_len <= 0) call stop_all(this,'Lanczos basis not positive.')
            if (nlanczos_eigv <= 0) call stop_all(this,'# lanczos eigenvalues not positive.')
        end if

        if (.not.doing_calc(dmqmc_calc) .and. dmqmc_calc_type /= 0 .and. parent) call warning(this,&
               'You are not performing a DMQMC calculation but have requested DMQMC options to be calculated.')
        if (doing_calc(fciqmc_calc)) then
            if (.not.doing_calc(simple_fciqmc_calc)) then
                if (qmc_in%walker_length == 0) call stop_all(this,'Walker length zero.')
                if (qmc_in%spawned_walker_length == 0) call stop_all(this,'Spawned walker length zero.')
            end if
            if (calc_inst_rdm .and. dmqmc_in%rdm%spawned_length == 0) call stop_all(this,'Spawned RDM length zero.')
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
            if (load_bal_in%percent < 0 .or. par_info%load%percent > 1.0) &
                call stop_all(this, 'Percentage imbalance must be positive and less that 1.')
            if (load_bal_in%max_attempts < 0) call stop_all(this, 'Maximum number of load balancing attempts must be positive')
        end if
        if (doing_calc(ct_fciqmc_calc)) qmc_in%ncycles = 1

        if (doing_dmqmc_calc(dmqmc_rdm_r2) .and. (.not. dmqmc_in%replica_tricks)) call stop_all(this,&
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

        if (dmqmc_in%all_spin_sectors) then
            if (.not. doing_calc(dmqmc_calc)) call stop_all(this, 'The use_all_spin_sectors option can only be used in&
                                                                   & DMQMC calculations.')
            if (abs(sys%heisenberg%magnetic_field) > depsilon .or. &
                abs(sys%heisenberg%staggered_magnetic_field) > depsilon) &
                    call stop_all(this, 'The use_all_spin_sectors option cannot be used with magnetic fields.')
            if (dmqmc_in%rdm%calc_ground_rdm) call stop_all(this, 'The use_all_spin_sectors and ground_state_rdm options&
                                                      & cannot be used together.')
        end if

        if (restart_in%dump_restart_file_shift) then
             if (restart_info_global_shift%write_id<0 .and. restart_info_global%write_id<0 &
                 .and. restart_info_global%write_id == restart_info_global_shift%write_id) &
                 call stop_all(this, 'The ids of the restart files are the same.')
             if (restart_info_global_shift%write_id<0 .and. restart_info_global%write_restart_freq /= huge(0) )&
                 call stop_all(this, 'The ids of the restart files could be the same')
        end if   
        if (dmqmc_in%vary_weights .and. (.not. dmqmc_in%weighted_sampling)) then
            call stop_all(this, 'The vary_weights option can only be used together with the weighted_sampling option.')
        end if
        if (sys%system /= heisenberg .and. dmqmc_calc_type > dmqmc_energy) then
            call stop_all(this, 'The observable requested is not currently implemented for this Hamiltonian.')
        end if

        if (parent) write (6,'(/,1X,13("-"),/)')

    end subroutine check_input

    subroutine distribute_input(sys, qmc_in, fciqmc_in, ccmc_in, semi_stoch_in, restart_in, load_bal_in, reference, dmqmc_in)

        ! Distribute the data read in by the parent processor to all other
        ! processors.

        ! Completely empty (courtesy of C pre-processing) when compiled in
        ! serial.

        ! In/Out:
        !    sys: object describing the system.  All parameters which can be set
        !       in the input file are distributed to other processors.
        !    fciqmc_in: input options relating to FCIQMC.
        !    ccmc_in: input options relating to CCMC.
        !    qmc_in: input options relating to QMC methods.
        !    semi_stoch_in: Input options for the semi-stochastic adaptation.
        !    restart_in: input options for HDF5 restart files.
        !    reference: current reference determinant.
        !    dmqmc_in: input options relating to DMQMC.

        use qmc_data, only: qmc_in_t, fciqmc_in_t, ccmc_in_t, semi_stoch_in_t
        use qmc_data, only: restart_in_t, load_bal_in_t, reference_t
        use dmqmc_data, only: dmqmc_in_t

#ifndef PARALLEL

        use system, only: sys_t

        type(sys_t), intent(inout) :: sys
        type(qmc_in_t), intent(inout) :: qmc_in
        type(fciqmc_in_t), intent(inout) :: fciqmc_in
        type(ccmc_in_t), intent(inout) :: ccmc_in
        type(semi_stoch_in_t), intent(inout) :: semi_stoch_in
        type(restart_in_t), intent(inout) :: restart_in
        type(load_bal_in_t), intent(inout) :: load_bal_in
        type(reference_t), intent(inout) :: reference
        type(dmqmc_in_t), intent(inout) :: dmqmc_in

#else

        use mpi
        use parallel
        use checking, only: check_allocate

        use system

        type(sys_t), intent(inout) :: sys
        type(qmc_in_t), intent(inout) :: qmc_in
        type(fciqmc_in_t), intent(inout) :: fciqmc_in
        type(ccmc_in_t), intent(inout) :: ccmc_in
        type(semi_stoch_in_t), intent(inout) :: semi_stoch_in
        type(restart_in_t), intent(inout) :: restart_in
        type(load_bal_in_t), intent(inout) :: load_bal_in
        type(reference_t), intent(inout) :: reference
        type(dmqmc_in_t), intent(inout) :: dmqmc_in

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
        call mpi_bcast(qmc_in%real_amplitudes, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(qmc_in%spawn_cutoff, 1, mpi_preal, 0, mpi_comm_world, ierr)
        call mpi_bcast(semi_stoch_in%start_iter, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(semi_stoch_in%shift_iter, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(semi_stoch_in%determ_space_type, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(semi_stoch_in%target_size, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(semi_stoch_in%write_determ_space, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(semi_stoch_in%separate_annihil, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(ccmc_in%full_nc, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(ccmc_in%linked, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(dmqmc_in%replica_tricks, 1, mpi_logical, 0, mpi_comm_world, ierr)
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
        call mpi_bcast(fciqmc_in%select_ref_det_every_nreports, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(fciqmc_in%ref_det_factor, 1, mpi_preal, 0, mpi_comm_world, ierr)
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
        call mpi_bcast(nkinetic_cycles, 1, mpi_integer, 0, mpi_comm_world, ierr)

        call mpi_bcast(lanczos_string_len, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(nlanczos_eigv, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(use_sparse_hamil, 1, mpi_logical, 0, mpi_comm_world, ierr)

        call mpi_bcast(print_fci_wfn, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(analyse_fci_wfn, 1, mpi_integer, 0, mpi_comm_world, ierr)

        call mpi_bcast(qmc_in%ncycles, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(qmc_in%nreport, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(dmqmc_in%beta_loops, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(qmc_in%walker_length, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(qmc_in%spawned_walker_length, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(dmqmc_in%rdm%spawned_length, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(qmc_in%tau, 1, mpi_preal, 0, mpi_comm_world, ierr)
        call mpi_bcast(qmc_in%tau_search, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(qmc_in%initial_shift, 1, mpi_preal, 0, mpi_comm_world, ierr)
        call mpi_bcast(qmc_in%vary_shift_from, 1, mpi_preal, 0, mpi_comm_world, ierr)
        call mpi_bcast(qmc_in%vary_shift_from_proje, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(qmc_in%target_particles, 1, mpi_preal, 0, mpi_comm_world, ierr)
        call mpi_bcast(dmqmc_in%rdm%doing_rdm, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(dmqmc_in%rdm%calc_ground_rdm, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(calc_inst_rdm, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(output_rdm, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(dmqmc_in%rdm%nrdms, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(doing_exact_rdm_eigv, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(doing_vn_entropy, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(doing_concurrence, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(dmqmc_in%start_av_excit_dist, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(dmqmc_in%start_av_rdm, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(dmqmc_in%weighted_sampling, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(dmqmc_in%vary_weights, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(dmqmc_in%find_weights, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(dmqmc_in%all_sym_sectors, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(dmqmc_in%all_spin_sectors, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(dmqmc_in%finish_varying_weights, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(dmqmc_in%propagate_to_beta, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(dmqmc_in%free_electron_trial, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(dmqmc_in%init_beta, 1, mpi_preal, 0, mpi_comm_world, ierr)
        call mpi_bcast(dmqmc_in%half_density_matrix, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(dmqmc_in%calc_excit_dist, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(dmqmc_in%metropolis_attempts, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(dmqmc_in%max_metropolis_move, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(sys%ueg%chem_pot, 1, mpi_preal, 0, mpi_comm_world, ierr)
        call mpi_bcast(dmqmc_in%grand_canonical_initialisation, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(dmqmc_in%fermi_temperature, 1, mpi_logical, 0, mpi_comm_world, ierr)
        option_set = .false.
        if (parent) option_set = allocated(dmqmc_in%sampling_probs)
        call mpi_bcast(option_set, 1, mpi_logical, 0, mpi_comm_world, ierr)
        if (option_set) then
            occ_list_size = size(dmqmc_in%sampling_probs)
            call mpi_bcast(occ_list_size, 1, mpi_integer, 0, mpi_comm_world, ierr)
            if (.not.parent) then
                allocate(dmqmc_in%sampling_probs(occ_list_size), stat=ierr)
                call check_allocate('dmqmc_in%sampling_probs',occ_list_size,ierr)
            end if
            call mpi_bcast(dmqmc_in%sampling_probs, occ_list_size, mpi_preal, 0, mpi_comm_world, ierr)
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
        if (parent) option_set = allocated(dmqmc_in%correlation_sites)
        call mpi_bcast(option_set, 1, mpi_logical, 0, mpi_comm_world, ierr)
        if (option_set) then
            occ_list_size = size(dmqmc_in%correlation_sites)
            call mpi_bcast(occ_list_size, 1, mpi_integer, 0, mpi_comm_world, ierr)
            if (.not.parent) then
                allocate(dmqmc_in%correlation_sites(occ_list_size), stat=ierr)
                call check_allocate('dmqmc_in%correlation_sites',occ_list_size,ierr)
            end if
            call mpi_bcast(dmqmc_in%correlation_sites, occ_list_size, mpi_integer, 0, mpi_comm_world, ierr)
        end if
        option_set = .false.
        call mpi_bcast(truncate_space, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(truncation_level, 1, mpi_integer, 0, mpi_comm_world, ierr)
        if (parent) option_set = allocated(reference%occ_list0)
        call mpi_bcast(option_set, 1, mpi_logical, 0, mpi_comm_world, ierr)
        if (option_set) then
            ! Have not yet set sys%nel in the Heisenberg model.
            occ_list_size = size(reference%occ_list0)
            call mpi_bcast(occ_list_size, 1, mpi_integer, 0, mpi_comm_world, ierr)
            if (.not.parent) then
                allocate(reference%occ_list0(occ_list_size), stat=ierr)
                call check_allocate('reference%occ_list0', occ_list_size, ierr)
            end if
            call mpi_bcast(reference%occ_list0, occ_list_size, mpi_integer, 0, mpi_comm_world, ierr)
        end if
        if (parent) option_set = allocated(reference%hs_occ_list0)
        call mpi_bcast(option_set, 1, mpi_logical, 0, mpi_comm_world, ierr)
        if (option_set) then
            ! Have not yet set sys%nel in the Heisenberg model.
            occ_list_size = size(reference%hs_occ_list0)
            call mpi_bcast(occ_list_size, 1, mpi_integer, 0, mpi_comm_world, ierr)
            if (.not.parent) then
                allocate(reference%hs_occ_list0(occ_list_size), stat=ierr)
                call check_allocate('reference%hs_occ_list0', occ_list_size, ierr)
            end if
            call mpi_bcast(reference%hs_occ_list0, occ_list_size, mpi_integer, 0, mpi_comm_world, ierr)
        end if
        call mpi_bcast(restart_in%read_restart, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(restart_in%dump_restart, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(restart_in%dump_restart_file_shift, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(restart_info_global%read_id, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(restart_info_global%write_id, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(restart_info_global%write_restart_freq, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(restart_info_global_shift%write_id, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(qmc_in%seed, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(qmc_in%shift_damping, 1, mpi_preal, 0, mpi_comm_world, ierr)
        call mpi_bcast(ccmc_in%cluster_multispawn_threshold, 1, mpi_preal, 0, mpi_comm_world, ierr)
        call mpi_bcast(qmc_in%D0_population, 1, mpi_preal, 0, mpi_comm_world, ierr)
        call mpi_bcast(qmc_in%no_renorm, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(qmc_in%pattempt_single, 1, mpi_preal, 0, mpi_comm_world, ierr)
        call mpi_bcast(qmc_in%pattempt_double, 1, mpi_preal, 0, mpi_comm_world, ierr)

        call mpi_bcast(ccmc_in%move_freq, 1, mpi_integer, 0, mpi_comm_world, ierr)

        call mpi_bcast(fciqmc_in%init_spin_inv_D0, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(qmc_in%initiator_pop, 1, mpi_preal, 0, mpi_comm_world, ierr)
        call mpi_bcast(qmc_in%initiator_approx, 1, mpi_logical, 0, mpi_comm_world, ierr)

        call mpi_bcast(hf_operator, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(alpha0, 1, mpi_integer, 0, mpi_comm_world, ierr)

        call mpi_bcast(write_hamiltonian, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(write_determinants, 1, mpi_logical, 0, mpi_comm_world, ierr)

        call mpi_bcast(block_size, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(fciqmc_in%non_blocking_comm, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(fciqmc_in%doing_load_balancing, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(load_bal_in%nslots, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(load_bal_in%pop, 1, mpi_integer8, 0, mpi_comm_world, ierr)
        call mpi_bcast(load_bal_in%percent, 1, mpi_preal, 0, mpi_comm_world, ierr)
        call mpi_bcast(load_bal_in%max_attempts, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(load_bal_in%write_info, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(use_mpi_barriers, 1, mpi_logical, 0, mpi_comm_world, ierr)

#endif

    end subroutine distribute_input

end module parse_input

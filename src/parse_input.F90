module parse_input
! Parse input options and check input for validity.

use const

use parallel, only: parent
use errors
use hilbert_space
use system
use calc
use lanczos
use determinants
use fciqmc_data
use fciqmc_restart, only: read_restart_number, write_restart_number,&
                          binary_fmt_in, binary_fmt_out 
use hubbard_real, only: finite_cluster
use hfs_data, only: lmag2

implicit none

contains

    subroutine read_input()

        ! Read input options from a file (if specified on the command line) or via
        ! STDIN.

! nag doesn't automatically bring in command-line option handling.
#ifdef NAGF95
        use f90_unix_env
#endif
    
        use input
        use utils, only: get_free_unit
        use checking, only: check_allocate

#ifdef NAGF95
        use f90_unix_env, ONLY: getarg,iargc
#else
        integer :: iargc ! External function.
#endif

        character(255) :: cInp
        character(100) :: w
        integer :: ios
        logical :: eof, t_exists
        integer :: ivec, i, ierr

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
                system_type = hub_real
            case('HUBBARD_K','HUBBARD_MOMENTUM')
                system_type = hub_k
            case('HEISENBERG')
                system_type = heisenberg
            case('READ')
                system_type = read_in
                if (item /= nitems) call reada(fcidump)
            case('UEG')
                system_type = ueg

            ! System information.
            case('LATTICE')
                ! Lattice block
                call read_line(eof)
                if (eof) call stop_all('read_input','Unexpected end of file reading lattice vectors.')
                ! nitems gives the number of items in the line, and thus the number
                ! of dimensions...
                ndim = nitems
                allocate(lattice(ndim,ndim), stat=ierr)
                call check_allocate('lattice',ndim*ndim,ierr)
                do ivec = 1, ndim
                    if (nitems /= ndim) call stop_all('read_input', 'Do not understand lattice vector.')
                    do i = 1, ndim
                        call readi(lattice(i, ivec))
                    end do
                    if (ivec /= ndim) then
                        call read_line(eof)
                        if (eof) call stop_all('read_input', 'Unexpected end of file reading lattice vectors.')
                    end if
                end do
            case('2D')
                if (ndim > 0) then
                    call warning('read_input','Dimension already set; ignoring keyword '//w)
                else
                    ndim = 2
                end if
            case('3D')
                if (ndim > 0) then
                    call warning('read_input','Dimension already set; ignoring keyword '//w)
                else
                    ndim = 3
                end if
            case('NEL', 'ELECTRONS')
                if (system_type == heisenberg) &
                     call stop_all('read_input', 'Cannot set electron number for Heisenberg. &
                     &Please enter a Ms value instead.')
                call readi(nel)

            ! Hubbard-specific system info
            case('T')
                call readf(hubt)
            case('U')
                call readf(hubu)

            ! Heisenberg-specific system info.
            case('J')
                call readf(J_coupling)
            case('MAGNETIC_FIELD')
                call readf(magnetic_field)
            case('STAGGERED_MAGNETIC_FIELD')
                call readf(staggered_magnetic_field)

            ! UEG-specific system info.
            case('RS','DENSITY')
                call readf(r_s)
            case('ECUTOFF')
                call readf(ueg_ecutoff)

            ! Select symmetry of wavefunction.
            case('MS')
                call readi(ms_in)
            case('SYM','SYMMETRY')
                call readi(sym_in)

            case('SEPARATE_STRINGS')
                separate_strings = .true.
            case('CAS')
                do i = 1,2
                    call readi(CAS(i))
                end do
            case('TWIST')
                allocate(ktwist(nitems-item), stat=ierr)
                call check_allocate('ktwist',nitems-item,ierr)
                do i = 1, nitems-item
                    call readf(ktwist(i))
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
                calc_type = calc_type + initiator_fciqmc
            case('CT_FCIQMC')
                calc_type = calc_type + ct_fciqmc_calc
            case('DMQMC')
                calc_type = calc_type + dmqmc_calc
            case('HELLMANN-FEYNMAN')
                calc_type = calc_type + hfs_fciqmc_calc
            case('ESTIMATE_HILBERT_SPACE')
                calc_type = calc_type + mc_hilbert_space
                call readi(nhilbert_cycles)
            case('FOLDED_SPECTRUM')
                calc_type = calc_type + folded_spectrum

            ! DMQMC expectation values to be calculated
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
            print *, correlation_sites
            case('DMQMC_STAGGERED_MAGNETISATION')
                dmqmc_calc_type = dmqmc_calc_type + dmqmc_staggered_magnetisation
            ! Calculate a reduced density matrix
            case('REDUCED_DENSITY_MATRIX')
                doing_reduced_dm = .true.
            case('SUBSYSTEM_SITES')
                allocate(subsystem_A_list(nitems-1), stat=ierr)
                call check_allocate('subsystem_A_list',nitems-1,ierr)
                do i = 1, nitems-1
                    call readi(subsystem_A_list(i))
                end do
            ! calculation options: DMQMC
            case('TRUNCATION_LEVEL')
                truncate_space = .true.
                call readi(truncation_level)

            ! Calculation options: lanczos.
            case('LANCZOS_BASIS')
                call readi(lanczos_basis_length)
            case('LANCZOS_SOLUTIONS','LANCZOS_SOLNS')
                call readi(nlanczos_eigv)

            ! Calculation options: lanczos/exact diagonalisation.
            case('PRINT_GROUND_STATE')
                print_ground_state = .true.
            case('ANALYSE_GROUND_STATE')
                analyse_ground_state = .true.

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
            case('TAU')
                call readf(tau)
            case('INITIAL_SHIFT')
                call readf(shift)
                ! We assume the user is sensible/knows what he/she is doing if
                ! initial_shift and vary_shift_from are set.
                vary_shift_from = shift 
            case('VARY_SHIFT_FROM')
                call readu(w)
                if (w == 'PROJE') then
                    vary_shift_from_proje = .true.
                else
                    call reread(0)
                    call readf(vary_shift_from)
                end if
            case('VARYSHIFT_TARGET')
                call readli(target_particles)
            case('INIT_POP')
                call readf(D0_population)
            case('REFERENCE_DET')
                allocate(occ_list0(nitems-1), stat=ierr)
                call check_allocate('occ_list0',nitems-1,ierr)
                do i = 1, nitems-1
                    call readi(occ_list0(i))
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

            ! Calculation options: Folded spectrum.
            case('FOLD_LINE')
                call readf(fold_line)
            case('P__')
                call readf(P__)
            case('Po_')
                call readf(Po_)
            case('P_o')
                call readf(P_o)

            ! use a negative number to indicate that the restart numbers have
            ! been fixed.
            case('RESTART')
                restart = .true.
                if (item /= nitems) then
                    call readi(read_restart_number)
                    read_restart_number = -read_restart_number-1
                end if
            case('DUMP_RESTART')
                dump_restart_file = .true.
                if (item /= nitems) then
                    call readi(write_restart_number)
                    write_restart_number = -write_restart_number-1
                end if
            case('ASCII_FORMAT_IN')
                binary_fmt_in = .false.
            case('ASCII_FORMAT_OUT')
                binary_fmt_out = .false.
            case('ASCII_FORMAT')
                binary_fmt_in = .false.
                binary_fmt_out = .false.
            case('SEED')
                call readi(seed)
            case('SHIFT_DAMPING')
                call readf(shift_damping)
            case('INIT_SPIN_INVERSE_REFERENCE_DET')
                init_spin_inv_D0 = .true.

            ! Calculation options: initiator-fciqmc.
            case('INITIATOR_CAS')
                call readi(initiator_cas(1))
                call readi(initiator_cas(2))
            case('INITIATOR_POPULATION')
                call readi(initiator_population)

            ! Calculation options: operators sampled using Hellmann--Feynman.
            case('L2')
                ! Set value of |l|^2 which is used
                call readi(lmag2)

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
             
            case('FINITE_CLUSTER')
                ! this will be checked in check_input to ensure that it 
                ! is only used when we are formulating the calculation
                ! in real-space
                finite_cluster = .true.   
            case('TRIANGULAR_LATTICE')
                triangular_lattice = .true.
            
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

    subroutine check_input()

        ! I don't pretend this is the most comprehensive of tests, but at least
        ! make sure a few things are not completely insane.

        use const
        
        integer :: ivec, jvec
        character(*), parameter :: this='check_input'

        if (system_type /= heisenberg) then
            if (nel <= 0) call stop_all(this,'Number of electrons must be positive.')
                if (trial_function /= single_basis) call stop_all(this, 'Only a single determinant can be used as the reference&
                                                     & state for this system. Other trial functions are not avaliable.')
                if (guiding_function /= no_guiding) &
                    call stop_all(this, 'Importance sampling is only avaliable for the Heisenberg model&
                                         & currently.')
        end if

        if (system_type /= read_in) then

            if (.not.(allocated(lattice))) call stop_all(this, 'Lattice vectors not provided')

            if (system_type == heisenberg) then
                if (ms_in > nsites) call stop_all(this,'Value of Ms given is too large for this lattice.')
                if ((-ms_in) > nsites) call stop_all(this,'Value of Ms given is too small for this lattice.')
                if (mod(abs(ms_in),2) /=  mod(nsites,2)) call stop_all(this, 'Ms value specified is not&
                                                                              & possible for this lattice.')
                if (staggered_magnetic_field /= 0.0_p .and. (.not.bipartite_lattice)) &
                    call stop_all(this, 'Cannot set a staggered field&
                                        & for this lattice because it is frustrated.')
                if (staggered_magnetic_field /= 0.0_p .and. magnetic_field /= 0.0_p) &
                    call stop_all(this, 'Cannot set a uniform and a staggered field at the same time.')
                if ((guiding_function==neel_singlet_guiding) .and. trial_function /= neel_singlet) call stop_all(this, 'This &
                                                         &guiding function is only avaliable when using the Neel singlet state &
                                                         &as an energy estimator.') 
                if (doing_dmqmc_calc(dmqmc_staggered_magnetisation) .and. (.not.bipartite_lattice)) then
                    call warning('check_input','Staggered magnetisation can only be calculated on a bipartite lattice.&
                                          & This is not a bipartite lattice. Changing options so that it will not be calculated.')
                    dmqmc_calc_type = dmqmc_calc_type - dmqmc_staggered_magnetisation
                end if
            else if (system_type == hub_k .or. system_type == hub_real) then
                if (nel > 2*nsites) call stop_all(this, 'More than two electrons per site.')
            end if

            if (ndim > 3) call stop_all(this, 'Limited to 1,  2 or 3 dimensions')
            if (system_type == ueg .and. ndim == 1) call stop_all(this, 'UEG only functional in 2D and 3D')
        
            if (triangular_lattice .and. (.not.bipartite_lattice) .and. (.not.finite_cluster)) then
                call warning('check_input','Periodic boundary conditions may not be applied for these particular &
                               &triangular lattice. Periodic boundary conditions are being turned off.')
                finite_cluster = .true.
            end if

            do ivec = 1, ndim
                do jvec = ivec+1, ndim
                    if (dot_product(lattice(:,ivec), lattice(:,jvec)) /= 0) then
                        call stop_all(this, 'Lattice vectors are not orthogonal.')
                    end if
                end do
            end do

        end if

        if (init_spin_inv_D0 .and. ms_in /= 0) then
            call warning(this, 'Flipping the reference state will give &
                                            &a state which has a different value of Ms and so cannot be used here.')
            init_spin_inv_D0 = .false.
        end if

        if (allocated(correlation_sites) .and. size(correlation_sites) /= 2) call stop_all(this, 'You must enter exactly two &
               &sites for the correlation function option.')

        ! Calculation specific checking.
        if (doing_calc(lanczos_diag)) then
            if (lanczos_basis_length <= 0) call stop_all(this,'Lanczos basis not positive.')
            if (nlanczos_eigv <= 0) call stop_all(this,'# lanczos eigenvalues not positive.')
        end if

        if ((.not.doing_calc(dmqmc_calc)) .and. dmqmc_calc_type /= 0) call warning('check_input',&
               'You are not performing a DMQMC calculation but have requested DMQMC options to be calculated.')

        if (doing_calc(fciqmc_calc)) then
            if (.not.doing_calc(simple_fciqmc_calc)) then
                if (walker_length == 0) call stop_all(this,'Walker length zero.')
                if (spawned_walker_length == 0) call stop_all(this,'Spawned walker length zero.')
            end if
            if (tau <= 0) call stop_all(this,'Tau not positive.')
            if (shift_damping <= 0) call stop_all(this,'Shift damping not positive.')
            if (allocated(occ_list0)) then
                if (size(occ_list0) /= nel) then
                    if (system_type /= heisenberg) then
                        call stop_all(this,'Number of electrons specified is different from &
                        &number of electrons used in the reference determinant.')
                    end if
                end if
            end if
            if (any(initiator_CAS < 0)) call stop_all(this,'Initiator CAS space must be non-negative.')
        end if
        if (doing_calc(ct_fciqmc_calc)) ncycles = 1
         
        ! If the FINITE_CLUSTER keyword was detected then make sure that 
        ! we are doing a calculation in real-space. If we're not then
        ! unset finite cluster,tell the user and carry on
        if(momentum_space) then
            if (finite_cluster .and. parent) call warning('check_input','FINITE_CLUSTER keyword only valid for hubbard&
                                      & calculations in real-space: ignoring keyword')
            if (separate_strings .and. parent) call warning('check_input','SEPARATE_STRINGS keyword only valid for hubbard&
                                      & calculations in real-space: ignoring keyword')
            finite_cluster = .false.
            separate_strings = .false.
        end if

        if (separate_strings) then
            if (system_type.ne.hub_real) then
                separate_strings = .false.
                if (parent) call warning('check_input','SEPARATE_STRINGS keyword only valid for hubbard&
                                      & calculations in real-space: ignoring keyword')
            else if (ndim /= 1) then
                separate_strings = .false.
                if (parent) call warning('check_input','SEPARATE_STRINGS keyword only valid for 1D&
                                      & calculations in real-space: ignoring keyword')
            end if
        end if
        
        if (parent) write (6,'(/,1X,13("-"),/)') 

    end subroutine check_input

    subroutine distribute_input()

        ! Distribute the data read in by the parent processor to all other
        ! processors.

        ! Completely empty (courtesy of C pre-processing) when compiled in
        ! serial.

#ifdef PARALLEL

        use mpi
        use parallel
        use checking, only: check_allocate

        integer :: ierr, occ_list_size, subsystem_size
        logical :: option_set

        call mpi_bcast(system_type, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(fcidump, len(fcidump), mpi_character, 0, mpi_comm_world, ierr)
        call mpi_bcast(sym_in, 1, mpi_integer, 0, mpi_comm_world, ierr)

        call mpi_bcast(ndim, 1, mpi_integer, 0, mpi_comm_world, ierr)
        if (.not.parent) then
            allocate(lattice(ndim,ndim), stat=ierr)
            call check_allocate('lattice',ndim*ndim,ierr)
        end if
        call mpi_bcast(lattice, ndim*ndim, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(finite_cluster, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(triangular_lattice, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(trial_function, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(guiding_function, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(nel, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(hubt, 1, mpi_preal, 0, mpi_comm_world, ierr)
        call mpi_bcast(hubu, 1, mpi_preal, 0, mpi_comm_world, ierr)
        call mpi_bcast(J_coupling, 1, mpi_preal, 0, mpi_comm_world, ierr)
        call mpi_bcast(magnetic_field, 1, mpi_preal, 0, mpi_comm_world, ierr)
        call mpi_bcast(staggered_magnetic_field, 1, mpi_preal, 0, mpi_comm_world, ierr)
        if (parent) option_set = allocated(ktwist)
        call mpi_bcast(option_set, 1, mpi_logical, 0, mpi_comm_world, ierr)
        if (option_set) then
            if (.not.parent) then
                allocate(ktwist(ndim), stat=ierr)
                call check_allocate('ktwist',ndim,ierr)
            end if
            call mpi_bcast(ktwist, ndim, mpi_preal, 0, mpi_comm_world, ierr)
        end if
        call mpi_bcast(r_s, 1, mpi_preal, 0, mpi_comm_world, ierr)
        call mpi_bcast(ueg_ecutoff, 1, mpi_preal, 0, mpi_comm_world, ierr)
        call mpi_bcast(separate_strings, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(select_ref_det_every_nreports, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(ref_det_factor, 1, mpi_preal, 0, mpi_comm_world, ierr)
        call mpi_bcast(cas, 2, mpi_integer, 0, mpi_comm_world, ierr)

        call mpi_bcast(ms_in, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(sym_in, 1, mpi_integer, 0, mpi_comm_world, ierr)

        call mpi_bcast(calc_type, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(dmqmc_calc_type, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(direct_lanczos, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(nhilbert_cycles, 1, mpi_integer, 0, mpi_comm_world, ierr)

        call mpi_bcast(lanczos_basis_length, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(nlanczos_eigv, 1, mpi_integer, 0, mpi_comm_world, ierr)

        call mpi_bcast(print_ground_state, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(analyse_ground_state, 1, mpi_logical, 0, mpi_comm_world, ierr)

        call mpi_bcast(ncycles, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(nreport, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(beta_loops, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(walker_length, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(spawned_walker_length, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(tau, 1, mpi_preal, 0, mpi_comm_world, ierr)
        call mpi_bcast(shift, 1, mpi_preal, 0, mpi_comm_world, ierr)
        call mpi_bcast(vary_shift_from, 1, mpi_preal, 0, mpi_comm_world, ierr)
        call mpi_bcast(vary_shift_from_proje, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(target_particles, 1, mpi_integer8, 0, mpi_comm_world, ierr)
        call mpi_bcast(doing_reduced_dm, 1, mpi_logical, 0, mpi_comm_world, ierr)
        option_set = .false.
        if (parent) option_set = allocated(subsystem_A_list)
        call mpi_bcast(option_set, 1, mpi_logical, 0, mpi_comm_world, ierr)
        if (option_set) then
            occ_list_size = size(subsystem_A_list)
            call mpi_bcast(subsystem_A_list, 1, mpi_integer, 0, mpi_comm_world, ierr)
            if (.not.parent) then
                allocate(subsystem_A_list(occ_list_size), stat=ierr)
                call check_allocate('subsystem_A_list',occ_list_size,ierr)
            end if
            call mpi_bcast(subsystem_A_list, occ_list_size, mpi_integer, 0, mpi_comm_world, ierr)
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
            ! Have not yet set nel in the Heisenberg model.
            occ_list_size = size(occ_list0)
            call mpi_bcast(occ_list_size, 1, mpi_integer, 0, mpi_comm_world, ierr)
            if (.not.parent) then
                allocate(occ_list0(occ_list_size), stat=ierr)
                call check_allocate('occ_list0', occ_list_size, ierr)
            end if
            call mpi_bcast(occ_list0, occ_list_size, mpi_integer, 0, mpi_comm_world, ierr)
        end if
        call mpi_bcast(restart, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(dump_restart_file, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(read_restart_number, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(write_restart_number, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(seed, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(shift_damping, 1, mpi_preal, 0, mpi_comm_world, ierr)
        call mpi_bcast(D0_population, 1, mpi_preal, 0, mpi_comm_world, ierr)
        call mpi_bcast(no_renorm, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(pattempt_single, 1, mpi_preal, 0, mpi_comm_world, ierr)
        call mpi_bcast(pattempt_double, 1, mpi_preal, 0, mpi_comm_world, ierr)

        call mpi_bcast(init_spin_inv_D0, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(initiator_CAS, 2, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(initiator_population, 1, mpi_integer, 0, mpi_comm_world, ierr)

        call mpi_bcast(lmag2, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(write_hamiltonian, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(write_determinants, 1, mpi_logical, 0, mpi_comm_world, ierr)

        call mpi_bcast(block_size, 1, mpi_integer, 0, mpi_comm_world, ierr)

        call mpi_bcast(fold_line, 1, mpi_preal, 0, mpi_comm_world, ierr)
        call mpi_bcast(P__, 1, mpi_preal, 0, mpi_comm_world, ierr)
        call mpi_bcast(Po_, 1, mpi_preal, 0, mpi_comm_world, ierr)
        call mpi_bcast(P_o, 1, mpi_preal, 0, mpi_comm_world, ierr)

#endif

    end subroutine distribute_input

end module parse_input

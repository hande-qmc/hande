module parse_input
! Parse input options and check input for validity.

use parallel, only: parent
use errors
use hilbert_space
use system
use calc
use lanczos
use determinants
use fciqmc_data
use fciqmc_restart, only: read_restart_number, write_restart_number 
use hubbard_real, only: finite_cluster

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
            case('REAL_SPACE')
                system_type = hub_real
            case('K_SPACE','MOMENTUM_SPACE')
                system_type = hub_k

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
            case('NEL', 'ELECTRONS')
                call readi(nel)
            case('T')
                call readf(hubt)
            case('U')
                call readf(hubu)
            case('TWIST')
                allocate(ktwist(nitems-item), stat=ierr)
                call check_allocate('ktwist',nitems-item,ierr)
                do i = 1, nitems-item
                    call readf(ktwist(i))
                end do

            ! Select symmetry of wavefunction.
            case('MS')
                call readi(ms_in)
            case('SYM','SYMMETRY')
                call readi(sym_in)

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
            case('HELLMANN-FEYNMAN')
                calc_type = calc_type + hfs_fciqmc_calc
            case('ESTIMATE_HILBERT_SPACE')
                calc_type = calc_type + mc_hilbert_space
                call readi(nhilbert_cycles)

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
            case('WALKER_LENGTH')
                call readi(walker_length)
            case('SPAWNED_WALKER_LENGTH')
                call readi(spawned_walker_length)
            case('TAU')
                call readf(tau)
            case('INITIAL_SHIFT')
                call readf(shift)
            case('VARYSHIFT_TARGET')
                call readi(target_particles)
            case('REFERENCE_DET')
                allocate(occ_list0(nitems-1), stat=ierr)
                call check_allocate('occ_list0',nitems-1,ierr)
                do i = 1, nitems-1
                    call readi(occ_list0(i))
                end do
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
            case('SEED')
                call readi(seed)
            case('SHIFT_DAMPING')
                call readf(shift_damping)
            case('REFERENCE_DET_POPULATION')
                call readi(D0_population)

            ! Calculation options: initiator-fciqmc.
            case('CAS')
                call readi(CAS(1))
                call readi(CAS(2))
            case('INITIATOR_POPULATION')
                call readi(initiator_population)

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

        if (.not.(allocated(lattice))) call stop_all(this, 'Lattice vectors not provided')

        if (ndim > 3) call stop_all(this, 'Limited to 1,  2 or 3 dimensions')

        if (nel <= 0) call stop_all(this,'Number of electrons must be positive.')
        if (nel > 2*nsites) call stop_all(this, 'More than two electrons per site.')

        do ivec = 1, ndim
            do jvec = ivec+1, ndim
                if (dot_product(lattice(:,ivec), lattice(:,jvec)) /= 0) then
                    call stop_all(this, 'Lattice vectors are not orthogonal.')
                end if
            end do
        end do

        if (doing_calc(lanczos_diag)) then
            if (lanczos_basis_length <= 0) call stop_all(this,'Lanczos basis not positive.')
            if (nlanczos_eigv <= 0) call stop_all(this,'# lanczos eigenvalues not positive.')
        end if

        if (doing_calc(fciqmc_calc)) then
            if (.not.doing_calc(simple_fciqmc_calc)) then
                if (walker_length <= 0) call stop_all(this,'Walker length not positive.')
                if (spawned_walker_length <= 0) call stop_all(this,'Spawned walker length not positive.')
            end if
            if (tau <= 0) call stop_all(this,'Tau not positive.')
            if (shift_damping <= 0) call stop_all(this,'Shift damping not positive.')
            if (allocated(occ_list0)) then
                if (size(occ_list0) /= nel) call stop_all(this,'Number of electrons specified is different from &
                                                           &number of electrons used in the reference determinant.')
            end if
            if (any(CAS < 0)) call stop_all(this,'CAS space must be non-negative.')
        end if
         
        ! If the FINITE_CLUSTER keyword was detected then make sure that 
        ! we are doing a calculation in real-space. If we're not then
        ! unset finite cluster,tell the user and carry on
        if(finite_cluster .and. (system_type .ne. hub_real)) then
            finite_cluster = .false.    
            if (parent) call warning('check_input','FINITE_CLUSTER keyword only valid for hubbard&
                                      & calculations in real-space: ignoring keyword')
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

        integer :: ierr
        logical :: option_set

        call mpi_bcast(system_type, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(sym_in, 1, mpi_integer, 0, mpi_comm_world, ierr)

        call mpi_bcast(ndim, 1, mpi_integer, 0, mpi_comm_world, ierr)
        if (.not.parent) then
            allocate(lattice(ndim,ndim), stat=ierr)
            call check_allocate('lattice',ndim*ndim,ierr)
        end if
        call mpi_bcast(lattice, ndim*ndim, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(finite_cluster, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(nel, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(hubt, 1, mpi_preal, 0, mpi_comm_world, ierr)
        call mpi_bcast(hubu, 1, mpi_preal, 0, mpi_comm_world, ierr)
        if (parent) option_set = allocated(ktwist)
        call mpi_bcast(option_set, 1, mpi_logical, 0, mpi_comm_world, ierr)
        if (option_set) then
            if (.not.parent) then
                allocate(ktwist(ndim), stat=ierr)
                call check_allocate('ktwist',ndim,ierr)
            end if
            call mpi_bcast(ktwist, ndim, mpi_preal, 0, mpi_comm_world, ierr)
        end if

        call mpi_bcast(ms_in, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(sym_in, 1, mpi_integer, 0, mpi_comm_world, ierr)

        call mpi_bcast(calc_type, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(direct_lanczos, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(nhilbert_cycles, 1, mpi_integer, 0, mpi_comm_world, ierr)

        call mpi_bcast(lanczos_basis_length, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(nlanczos_eigv, 1, mpi_integer, 0, mpi_comm_world, ierr)

        call mpi_bcast(print_ground_state, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(analyse_ground_state, 1, mpi_logical, 0, mpi_comm_world, ierr)

        call mpi_bcast(ncycles, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(nreport, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(walker_length, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(spawned_walker_length, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(tau, 1, mpi_preal, 0, mpi_comm_world, ierr)
        call mpi_bcast(shift, 1, mpi_preal, 0, mpi_comm_world, ierr)
        call mpi_bcast(target_particles, 1, mpi_integer, 0, mpi_comm_world, ierr)
        if (parent) option_set = allocated(occ_list0)
        call mpi_bcast(option_set, 1, mpi_logical, 0, mpi_comm_world, ierr)
        if (option_set) then
            if (.not.parent) then
                allocate(occ_list0(nel), stat=ierr)
                call check_allocate('occ_list0',nel,ierr)
            end if
            call mpi_bcast(occ_list0, nel, mpi_integer, 0, mpi_comm_world, ierr)
        end if
        call mpi_bcast(restart, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(dump_restart_file, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(read_restart_number, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(write_restart_number, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(seed, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(shift_damping, 1, mpi_preal, 0, mpi_comm_world, ierr)
        call mpi_bcast(D0_population, 1, mpi_integer, 0, mpi_comm_world, ierr)

        call mpi_bcast(CAS, 2, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(initiator_population, 1, mpi_integer, 0, mpi_comm_world, ierr)

        call mpi_bcast(write_hamiltonian, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(write_determinants, 1, mpi_logical, 0, mpi_comm_world, ierr)

        call mpi_bcast(block_size, 1, mpi_integer, 0, mpi_comm_world, ierr)

#endif

    end subroutine distribute_input

end module parse_input

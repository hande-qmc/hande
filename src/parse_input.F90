module parse_input
! Parse input options and check input for validity.

use parallel, only: parent
use errors
use system
use calc
use lanczos
use determinants
use fciqmc_data
use fciqmc_restart, only: read_restart_number, write_restart_number 

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
            ir = 1
            call GetArg(1, cInp)
            inquire(file=cInp, exist=t_exists)
            if (.not.t_exists) then
                write (6,'(a21,1X,a)') 'File does not exist:',trim(cInp)
                stop
            end if
            open(1, file=cInp, status='old', form='formatted', iostat=ios)
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
                t_exact = .true.
            case('LANCZOS_DIRECT')
                t_lanczos = .true.
                direct_lanczos = .true.
            case('LANCZOS')
                t_lanczos = .true.
            case('SIMPLE_FCIQMC')
                tsimple = .true.
                t_fciqmc = .true.
            case('FCIQMC')
                t_fciqmc = .true.
            case('IFCIQMC')
                t_fciqmc = .true.
                initiator = .true.

            ! Calculation options: lanczos.
            case('LANCZOS_BASIS')
                call readi(lanczos_basis_length)
            case('LANCZOS_SOLUTIONS','LANCZOS_SOLNS')
                call readi(nlanczos_eigv)

            ! Calculation options: lanczos/exact diagonalisation.
            case('EIGENVALUES')
                find_eigenvectors = .false.
            case('EIGENVECTORS')
                find_eigenvectors = .true.

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
                do i = 1, nitems-1
                    call readi(occ_list0(i))
                end do
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

            case('END')
                exit
            case default
                call report('Keyword '//trim(w)//' not recognized.', .true.)
            end select
        end do ! end reading of input.

        if (ios.gt.0) then
            if (parent) write (6,*) 'Problem reading input.'
            stop
        end if

    end subroutine read_input

    subroutine check_input()

        ! I don't pretend this is the most comprehensive of tests, but at least
        ! make sure a few things are not completely insane.

        use const

        integer :: ivec, jvec

        if (.not.(allocated(lattice))) call stop_all('check_input', 'Lattice vectors not provided')

        if (ndim > 3) call stop_all('check_input', 'Limited to 1,  2 or 3 dimensions')

        if (nel <= 0) call stop_all('check_input','Number of electrons must be positive.')

        do ivec = 1, ndim
            do jvec = ivec+1, ndim
                if (dot_product(lattice(:,ivec), lattice(:,jvec)) /= 0) then
                    call stop_all('check_input', 'Lattice vectors are not orthogonal.')
                end if
            end do
        end do

        if (nel > 2*nsites) call stop_all('check_input', 'More than two electrons per site.')

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

        integer :: ierr
        logical :: option_set

        call mpi_bcast(system_type, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(sym_in, 1, mpi_integer, 0, mpi_comm_world, ierr)

        call mpi_bcast(ndim, 1, mpi_integer, 0, mpi_comm_world, ierr)
        if (.not.parent) allocate(lattice(ndim,ndim), stat=ierr)
        call mpi_bcast(lattice, ndim*ndim, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(nel, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(hubt, 1, mpi_preal, 0, mpi_comm_world, ierr)
        call mpi_bcast(hubu, 1, mpi_preal, 0, mpi_comm_world, ierr)
        if (.not.parent) allocate(lattice(ndim,ndim), stat=ierr)
        if (parent) option_set = allocated(ktwist)
        call mpi_bcast(option_set, 1, mpi_logical, 0, mpi_comm_world, ierr)
        if (option_set) then
            if (.not.parent) allocate(ktwist(ndim), stat=ierr)
            call mpi_bcast(ktwist, ndim, mpi_preal, 0, mpi_comm_world, ierr)
        end if

        call mpi_bcast(ms_in, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(sym_in, 1, mpi_integer, 0, mpi_comm_world, ierr)

        call mpi_bcast(t_exact, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(t_lanczos, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(direct_lanczos, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(tsimple, 1, mpi_logical, 0, mpi_comm_world, ierr)
        call mpi_bcast(t_fciqmc, 1, mpi_logical, 0, mpi_comm_world, ierr)

        call mpi_bcast(lanczos_basis_length, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_bcast(nlanczos_eigv, 1, mpi_integer, 0, mpi_comm_world, ierr)

        call mpi_bcast(find_eigenvectors, 1, mpi_logical, 0, mpi_comm_world, ierr)

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
            if (.not.parent) allocate(occ_list0(nel), stat=ierr)
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

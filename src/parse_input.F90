module parse_input
! Parse input options and check input for validity.

use parallel, only: iproc, parent
use errors
use system
use hamiltonian

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
            if (iproc == parent) write (6,'(a19)') 'Reading from STDIN'
            ir = 5
            ios = 0
        end if

        if (iproc == parent) write (6,'(a14,/,1X,13("-"),/)') 'Input options'
        call input_options(echo_lines=iproc==parent, skip_blank_lines=.true.)

        do
            call read_line(eof)
            if (eof) exit
            call readu(w)
            select case(w)
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
            case('EXACT','FCI')
                t_exact = .true.
            case('LANCZOS')
                t_lanczos = .true.
            case('EIGENVALUES')
                find_eigenvectors = .false.
            case('EIGENVECTORS')
                find_eigenvectors = .true.
            case('HAMIL','HAMILTONIAN')
                write_hamiltonian = .true.
                if (item /= nitems) call reada(hamiltonian_file)
            case('DET','DETERMINANTS')
                write_determinants = .true.
                if (item /= nitems) call reada(determinant_file)
            case('END')
                exit
            case default
                call report('Keyword '//trim(w)//' not recognized.', .true.)
            end select
        end do

        if (ios.gt.0) then
            if (iproc == parent) write (6,*) 'Problem reading input.'
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

        if (mod(nel, 2) /= 0) call stop_all('check_input', 'Odd number of electrons => open shell.')

        do ivec = 1, ndim
            do jvec = ivec+1, ndim
                if (dot_product(lattice(:,ivec), lattice(:,jvec)) /= 0) then
                    call stop_all('check_input', 'Lattice vectors are not orthogonal.')
                end if
            end do
        end do

        if (nel > 2*nsites) call stop_all('check_input', 'More than two electrons per site.')

        if (iproc == parent) write (6,'(/,1X,13("-"),/)') 

    end subroutine check_input

end module parse_input

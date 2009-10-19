module parse_input
!= Parse input options, set defaults and check input for validity.

use parallel

implicit none

contains

    subroutine read_input()

#ifdef NAGF95
        use f90_unix_env
#endif
    
        use input

! nag doesn't automatically bring in command-line option handling.
#ifdef NAGF95
        use f90_unix_env, ONLY: getarg,iargc
#else
        integer :: iargc ! External function.
#endif

        character(255) :: cInp
        character(100) :: w
        integer :: ios
        logical :: eof, t_exists


        !call set_defaults()

        if (iargc().gt.0) then
            ! Input file specified on the command line.
            ir=1
            call GetArg(1,cInp)
            inquire(file=cInp,exist=t_exists)
            if (.not.t_exists) then
                write (6,'(a21,1X,a)') 'File does not exist:',trim(cInp)
                stop
            end if
            open(1,file=cInp,status='old',form='formatted',iostat=ios)
        else
            if (i_proc.eq.0) write (6,'(a19)') 'Reading from STDIN'
            ir=5
            ios=0
        end if

        if (i_proc.eq.0) write (6,'(a14,/,1X,13("-"),/)') 'Input options'
        call input_options(echo_lines=i_proc.eq.0,skip_blank_lines=.true.)

        do
            call read_line(eof)
            if (eof) exit
            call readu(w)
            select case(w)

            case('END')
                exit
            case default
                call report('Keyword '//trim(w)//' not recognized.',.true.)
            end select
        end do

        if (ios.gt.0) then
            if (i_proc.eq.0) write (6,*) 'Problem reading input.'
            stop
        end if

        if (i_proc.eq.0) write (6,*) ! Formatting.

    end subroutine read_input

end module parse_input

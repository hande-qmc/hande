module errors
! Module for printing out warnings/errors.

implicit none

contains

    subroutine stop_all(sub_name,error_msg)
        ! Stop calculation due to an error.
        ! Exit with code 999.
        !
        ! In:
        !    sub_name:  calling subroutine name.
        !    error_msg: error message.

        use, intrinsic :: iso_fortran_env, only: error_unit

#ifdef PARALLEL
        use mpi
#endif

        character(*), intent(in) :: sub_name,error_msg

        ! It seems that giving STOP a string is far more portable.
        ! mpi_abort requires an integer though.
        character(3), parameter :: error_str='999'
#ifdef PARALLEL
        integer, parameter :: error_code=999
        integer :: ierr
#endif

        write (error_unit,'(/a7)') 'ERROR.'
        write (error_unit,'(1X,a)') 'HANDE stops in subroutine: '//adjustl(sub_name)//'.'
        write (error_unit,'(a9,a)') 'Reason: ',adjustl(error_msg)
        write (error_unit,'(1X,a10)') 'EXITING...'

        flush(error_unit)

        ! Abort all processors.
        ! error code is given to mpi_abort which (apparently) returns it to the invoking environment.
#if PARALLEL
        call mpi_abort(mpi_comm_world, error_code, ierr)
#endif

        stop error_str

        return

    end subroutine stop_all

    subroutine warning(sub_name,error_msg,blank_lines)
        ! Print a warning message in a (helpfully consistent) format.
        ! I was bored of typing the same formatting in different places. ;-)
        !
        ! In:
        !    sub_name:  calling subroutine name.
        !    error_msg: error message.
        !    blank_lines (optional): if 0, print a blank line either side of the
        !        warning message.  This is the default behaviour. If 1, a blank
        !        line is only printed before the warning message.  If 2, a blank
        !        line is only printed after the warning message.  No blank lines
        !        are printed for any other value.

        use, intrinsic :: iso_fortran_env, only: error_unit

        character(*), intent(in) :: sub_name,error_msg
        integer, optional :: blank_lines

        call write_blank(1)
        write (error_unit,'(1X,a)') 'WARNING: error in '//adjustl(sub_name)//'.'
        write (error_unit,'(1X,a)') adjustl(error_msg)
        call write_blank(2)

        return

        contains

            subroutine write_blank(point)

                integer :: point

                if (present(blank_lines)) then
                    if (blank_lines == 0) then
                        write (error_unit,'()')
                    else if (blank_lines == point) then
                        write (error_unit,'()')
                    end if
                else
                    write (error_unit,'()')
                end if

            end subroutine write_blank

    end subroutine warning

    subroutine quiet_stop(msg)
        ! Exit without making any noise.  Useful for when there's no error, but you
        ! still want to exit midway through a calculation (e.g. for testing purposes,
        ! or for use with the SOFTEXIT functionality).
        ! In:
        !    msg (optional) : Print msg before exiting if msg is present.

#ifdef PARALLEL
        use mpi
#endif

        character(*), intent(in), optional :: msg
#ifdef PARALLEL
        integer :: ierr
#endif

        if (present(msg)) then
            write (6,'(1X,a)') adjustl(msg)
            flush(6)
        end if

        ! Abort all processors.
        ! error code is given to mpi_abort which (apparently) returns it to the invoking environment.
#if PARALLEL
        call mpi_abort(mpi_comm_world, 0, ierr)
#endif
        stop

    end subroutine quiet_stop

end module errors

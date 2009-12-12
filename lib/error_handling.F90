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

#ifdef _PARALLEL
        use mpi
#endif

        character(*), intent(in) :: sub_name,error_msg

        ! It seems that giving STOP a string is far more portable.
        ! mpi_abort requires an integer though.
        integer, parameter :: error_code=999
        character(3), parameter :: error_str='999'
        integer :: ierr

        write (6,'(/a7)') 'ERROR.'
        write (6,'(a30,a)') 'hubbard stops in subroutine: ',adjustl(sub_name)
        write (6,'(a9,18X,a)') 'Reason: ',adjustl(error_msg)
        write (6,'(a11)') 'EXITING...'

        call flush(6)

        ! Abort all processors.
        ! error code is given to mpi_abort which (apparently) returns it to the invoking environment.
#if _PARALLEL
        call mpi_abort(mpi_comm_world, error_code, ierr)
#endif

        stop error_str

        return

    end subroutine stop_all

    subroutine warning(sub_name,error_msg)
        ! Print a warning message in a (helpfully consistent) format.
        ! I was bored of typing the same formatting in different places. ;-)
        !
        ! In:
        !    sub_name:  calling subroutine name.
        !    error_msg: error message.

        character(*), intent(in) :: sub_name,error_msg

        write (6,'(/,1X,a)') 'WARNING: error in '//adjustl(sub_name)//'.'
        write (6,'(1X,a/)') adjustl(error_msg)

        return

    end subroutine warning

    subroutine quiet_stop(msg)
        ! Exit without making any noise.  Useful for when there's no error, but you
        ! still want to exit midway through a calculation (e.g. for testing purposes,
        ! or for use with the SOFTEXIT functionality).
        ! In:
        !    msg (optional) : Print msg before exiting if msg is present.

#ifdef _PARALLEL
        use mpi
#endif

        character(*), intent(in), optional :: msg
        integer :: ierr

        if (present(msg)) then
            write (6,'(1X,a)') adjustl(msg)
            call flush(6)
        end if

        ! Abort all processors.
        ! error code is given to mpi_abort which (apparently) returns it to the invoking environment.
#if _PARALLEL
        call mpi_abort(mpi_comm_world, 0, ierr)
#endif
        stop

    end subroutine quiet_stop

end module errors

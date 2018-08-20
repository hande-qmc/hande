module report

! James Spencer, Imperial College London.
!
! Copyright (c) 2012 James Spencer.
!
! Permission is hereby granted, free of charge, to any person
! obtaining a copy of this software and associated documentation
! files (the "Software"), to deal in the Software without
! restriction, including without limitation the rights to use,
! copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the
! Software is furnished to do so, subject to the following
! conditions:
!
! The above copyright notice and this permission notice shall be
! included in all copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
! EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
! OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
! NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
! HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
! WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
! FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
! OTHER DEALINGS IN THE SOFTWARE.

use git_info, only: HANDE_VCS_VERSION => GIT_COMMIT_HASH

implicit none

! Global uuid
character(36) :: GLOBAL_UUID

contains

    subroutine host_writer(message, io_unit)

      use, intrinsic :: iso_c_binding, only: c_char, c_int, c_null_char

      character(kind=c_char), intent(in) :: message(*)
      integer(c_int), intent(in) :: io_unit
      integer(c_int) :: length, i

      length = 0
      do
         if (message(length + 1) == c_null_char) exit
         length = length + 1
      end do

      ! Let whoever sends a message decide its newlines
      write (unit=io_unit, fmt='(1000A)', advance="no") (message(i), i=1, length)

    end subroutine host_writer

    subroutine environment_report(io)

        ! In:
        !    io (optional): unit to which the environment information is written.
        !                   Default: 6.
        !
        ! Print out a summary of the environment:
        !   * when the code was compiled;
        !   * the VCS BASE repository version (there's no guarantee that the code wasn't
        !     changed!);
        !   * whether the working directory contains local changes;
        !   * the working directory;
        !   * the host computer.

        use, intrinsic :: iso_c_binding, only: c_char, c_ptr, c_size_t, c_int, c_associated, c_funloc, c_funptr
        use utils, only: carray_to_fstring
        use const, only: i0, int_p

        ! Accessing the hostname and working directory directly from Fortran are
        ! (admittedly commonly implemented) extensions so we cannot rely on them
        ! existing *cough*ibmandnag*cough*.  It is safer to use POSIX-standard
        ! functions.
        interface
            subroutine print_info(writer, io_unit) bind(C, name="print_info")
              import
              type(c_funptr), intent(in), value :: writer
              integer(c_int), intent(in) :: io_unit
            end subroutine print_info
        end interface

        integer, intent(in), optional :: io

        integer :: io_unit
        integer :: date_values(8)

        if (present(io)) then
            io_unit = io
        else
            io_unit = 6
        end if

        write (io_unit,'(1X,64("="))')
        flush(io_unit)

        call print_info(c_funloc(host_writer), io_unit)

        call date_and_time(VALUES=date_values)

        write (io_unit,'(a18,1X,i2.2,"/",i2.2,"/",i4.4,1X,a2,1X,i2.2,2(":",i2.2))') &
                   "Started running on", date_values(3:1:-1), "at", date_values(5:7)
        write (io_unit,'(1X,"Calculation UUID:",1X,a36,".")') GLOBAL_UUID

        write (io_unit,'(1X,64("="),/)')

    end subroutine environment_report

    subroutine comm_global_uuid()

        ! Send UUID from root to all other processors.

#ifdef PARALLEL
        use parallel
        integer :: ierr
        call mpi_bcast(GLOBAL_UUID, len(GLOBAL_UUID), mpi_character, 0, mpi_comm_world, ierr)
#endif

    end subroutine comm_global_uuid

    subroutine get_uuid(uuid)

        ! Out:
        !     UUID: a (reasonably!) unique identification string.

        use, intrinsic :: iso_c_binding
        use utils, only: carray_to_fstring

        implicit none

        character(36), intent(out) :: uuid

#ifdef DISABLE_UUID
        uuid = 'UNKNOWN: UUID GENERATION DISABLED.  '
#else
        character(c_char) :: uuid_bin(16)
        character(c_char), target :: uuid_str(37)
        type(c_ptr) :: ptr

        interface
            subroutine uuid_generate(uu) bind(c)
                use, intrinsic :: iso_c_binding
                implicit none
                character(c_char), intent(out) :: uu(16)
            end subroutine
            subroutine uuid_unparse(uu, uu_str) bind(c)
                use, intrinsic :: iso_c_binding
                implicit none
                character(c_char), intent(in) :: uu(16)
                character(c_char), intent(out) :: uu_str(37)
            end subroutine uuid_unparse
        end interface

        call uuid_generate(uuid_bin)
        ptr = c_loc(uuid_str)
        call uuid_unparse(uuid_bin, uuid_str)
        uuid = carray_to_fstring(uuid_str)
#endif

    end subroutine get_uuid

    subroutine end_report(wall_time, cpu_time_used, io)

        ! Print out date at end of calculation and how long it took.

        ! In:
        !    wall_time: number of seconds between the start and end of the
        !        calculation.
        !    cpu_time_used: number of seconds took by process and any and all
        !       child processes.
        !    io (optional): unit to which the environment information is written.
        !        Default: 6.

        real, intent(in) :: wall_time, cpu_time_used
        integer, intent(in), optional :: io
        integer :: date_values(8), io_unit

        if (present(io)) then
            io_unit = io
        else
            io_unit = 6
        end if

        write (io_unit,'(1X,64("="))')

        call date_and_time(VALUES=date_values)

        write (io_unit,'(1X,a19,1X,i2.2,"/",i2.2,"/",i4.4,1X,a2,1X,i2.2,2(":",i2.2))') &
                   "Finished running on", date_values(3:1:-1), "at", date_values(5:7)
        write (io_unit,'(1X,a20,17X,f14.2)') "Wall time (seconds):", wall_time
        write (io_unit,'(1X,a34,3X,f14.2)') "CPU time (per processor, seconds):", cpu_time_used

        write (io_unit,'(1X,64("="),/)')

    end subroutine end_report

    subroutine write_date_time_close(iunit, date_values)

        ! Write final message to a given io unit, if open, and close

        integer, intent(in) :: iunit, date_values(8)
        logical :: io_open

        inquire(unit=iunit, opened = io_open)

        if (io_open) then
            write (iunit,'(1X,64("="))')
            write (iunit,'(1X,a19,1X,i2.2,"/",i2.2,"/",i4.4,1X,a2,1X,i2.2,2(":",i2.2))') &
                      "Finished running on", date_values(3:1:-1), "at", date_values(5:7)
            write (iunit,'(1X,64("="))')

            close(iunit, status='keep')
        end if

    end subroutine write_date_time_close

end module report

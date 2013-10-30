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

implicit none

! Global uuid
character(36) :: GLOBAL_UUID

contains

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

        ! The VCS information is passed to environment_report via c-preprocessing.  The
        ! following definitions are used:
        !
        ! * __DATE__ and __TIME__: date and time of compilation.
        ! * _VCS_LOCAL_CHANGES: define if the local source code directory contains local
        !   (uncommitted) changes.
        ! * _VCS_VERSION: set to a (quoted!) string containing the VCS revision id of the current
        !   commit.
        !
        ! The VCS information can be simply obtained via shell commands in the makefile:
        ! see the example provided for the commands needed for git and subversion
        ! repositories.  For instance, if a subversion repository is being used, then the
        ! version id can be obtained using::
        !
        !     VCS_VERSION:=$(shell echo -n '"'$(svn info | awk '/Revision/{print $NF}')'"')
        !
        ! and if the working directory contains local changes, then the command::
        !
        !     svn st -q | xargs -i test -z {}
        !
        ! gives a non-zero return code.  Hence the variable VCS_LOCAL_CHANGES is
        ! equal to -D_VCS_LOCAL_CHANGES or is a null string depending on whether there
        ! are local changes or not::
        !
        !    VCS_LOCAL_CHANGES:=$(shell svn st -q | xargs -i test -z {} || echo -n "-D_VCS_LOCAL_CHANGES")
        !
        ! environment_report.F90 can thus be appropriately compiled with the command::
        !
        !    $(FC) $(VCS_LOCAL_CHANGES) -D_VCS_VERSION='$(VCS_VERSION)' -c environment_report.F90 -o environment_report.o
        !
        ! where $(FC) is defined to be the desired fortran compiler.
        !
        ! Similarly for git::
        !
        !     VCS_VERSION:=$(shell echo -n '"' && git log --max-count=1 --pretty=format:%H && echo -n '"')
        !     VCS_LOCAL_CHANGES := $(shell git diff-index --quiet --cached HEAD
        !                             --ignore-submodules -- && git diff-files --quiet
        !                             --ignore-submodules || echo -n "-D_VCS_LOCAL_CHANGES") # on one line.
        !
        ! environment_report has been tested with gfortran, g95, ifort, pgf90, nag and
        ! pathscale.  Note that hostnm and getcwd are intrinsic functions (at least
        ! according to gnu documentation), unless you have to use nag.  There should also
        ! exist the intrinsic subroutines, but pgf90 and ifort gave segmentation faults
        ! when they were used, so we use intrinsic functions unless using nag.

        use, intrinsic :: iso_c_binding, only: c_char, c_ptr, c_size_t, c_int, c_associated
        use utils, only: carray_to_fstring

        ! Accessing the hostname and working directory directly from Fortran are
        ! (admittedly commonly implemented) extensions so we cannot rely on them
        ! existing *cough*ibmandnag*cough*.  It is safer to use POSIX-standard
        ! functions.
        interface
            function gethostname(hostname, len) result(stat) bind(c)
                use, intrinsic :: iso_c_binding, only: c_int, c_char, c_size_t
                integer(c_int) :: stat
                integer(c_size_t), intent(in), value :: len
                character(kind=c_char), intent(inout) :: hostname(len)
            end function gethostname
            function getcwd(cwd, len) result(path) bind(c)
                use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t, c_char
                type(c_ptr) :: path
                integer(c_size_t), intent(in), value :: len
                character(kind=c_char), intent(inout) :: cwd(len)
            end function getcwd
        end interface

        integer, intent(in), optional :: io
        integer :: io_unit
        integer :: date_values(8)
        integer(c_size_t), parameter :: str_len = 255
        character(kind=c_char) :: str(str_len)
        type(c_ptr) :: path
        integer(c_int) :: stat

! Set defaults in case not provided by the preprocessor.
#ifndef __DATE__
#define __DATE__ 'unknown'
#endif
#ifndef __TIME__
#define __TIME__ 'unknown'
#endif
#ifndef _CONFIG
#define _CONFIG 'unknown'
#endif
#ifndef _VCS_VERSION
#define _VCS_VERSION 'unknown'
#define _VCS_LOCAL_CHANGES
#endif

        if (present(io)) then
            io_unit = io
        else
            io_unit = 6
        end if

        write (io_unit,'(1X,64("="))')

        write (io_unit,'(a13,a,a4,a)') 'Compiled on ',__DATE__,'at ',__TIME__
        write (io_unit,'(a16,a)') 'Compiled using ', _CONFIG

        write (io_unit,'(a29,/,5X,a)') 'VCS BASE repository version:',_VCS_VERSION
#ifdef _VCS_LOCAL_CHANGES
        write (io_unit,'(a46)') 'Source code directory contains local changes.'
#endif

        stat = gethostname(str, str_len)
        if (stat == 0) then
            write (io_unit, '(a10)') 'Hostname:'
        else
            write (io_unit, '(a53)') 'Hostname exceeds 255 characters; truncated hostname:'
        end if
        write (io_unit,'(5X,a)') trim(carray_to_fstring(str))

        path = getcwd(str, str_len)
        if (c_associated(path)) then
            write (io_unit,'(a20)') 'Working directory: '
        else
            write (io_unit,'(a63)') 'Working directory exceeds 255 characters; truncated directory:'
        end if
        write (io_unit,'(5X,a)') trim(carray_to_fstring(str))

        call date_and_time(VALUES=date_values)

        write (6,'(1X,a18,1X,i2.2,"/",i2.2,"/",i4.4,1X,a2,1X,i2.2,2(":",i2.2))') &
                   "Started running on", date_values(3:1:-1), "at", date_values(5:7)
        call get_uuid(GLOBAL_UUID)
        write (io_unit,'(1X,"Calculation UUID:",1X,a36,".")') GLOBAL_UUID

        write (io_unit,'(1X,64("="),/)')

        return

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
        write (io_unit,'(1X,a10,17X,f14.2,a1)') "Wall time:", wall_time, "s"
        write (io_unit,'(1X,a25,2X,f14.2,a1)') "CPU time (per processor):", cpu_time_used, "s"

        write (io_unit,'(1X,64("="),/)')

    end subroutine end_report

end module report

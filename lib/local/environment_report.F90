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

#ifndef _VCS_VERSION
#define _VCS_VERSION 'unknown'
#endif
character(*), parameter :: VCS_VERSION = _VCS_VERSION

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

        ! The VCS information is passed to environment_report via C-preprocessing.  The
        ! following definitions are used:
        !
        ! * __DATE__ and __TIME__: date and time of compilation.
        ! * _VCS_VERSION: set to a (quoted!) string containing the VCS revision id of the current
        !   commit.
        ! * _CONFIG: set to quoted string containing the configuration (e.g compilers) used.
        !
        ! The VCS information can be simply obtained via shell commands in the makefile, e.g. for git:
        !
        !    GIT_SHA1 := $(shell git rev-parse HEAD 2> /dev/null || echo "unknown")
        !    GIT_SHA1 := $(GIT_SHA1)$(shell test -z "$$(git status --porcelain -- $(SRCDIRS))" || echo -dirty)
        !
        ! environment_report.F90 can thus be appropriately compiled with the command::
        !
        !    $(FC) -D_VCS_VERSION="$(GIT_SHA1)" -D_CONFIG="<config>" environment_report.F90 -o environment_report.o
        !
        ! where $(FC) is defined to be the desired fortran compiler.

        use, intrinsic :: iso_c_binding, only: c_char, c_ptr, c_size_t, c_int, c_associated
        use utils, only: carray_to_fstring
        use const, only: i0, int_p

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
            function getcwd_c(cwd, len) result(path) bind(c, name='getcwd')
                use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t, c_char
                type(c_ptr) :: path
                integer(c_size_t), intent(in), value :: len
                character(kind=c_char), intent(inout) :: cwd(len)
            end function getcwd_c
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

        if (present(io)) then
            io_unit = io
        else
            io_unit = 6
        end if

        write (io_unit,'(1X,64("="))')

        write (io_unit,'(a13,a,a4,a)') 'Compiled on ',__DATE__,'at ',__TIME__
        write (io_unit,'(a16,a)') 'Compiled using ', _CONFIG

        write (io_unit,'(a29,/,5X,a)') 'VCS BASE repository version:',VCS_VERSION

        stat = gethostname(str, str_len)
        if (stat == 0) then
            write (io_unit, '(a10)') 'Hostname:'
        else
            write (io_unit, '(a53)') 'Hostname exceeds 255 characters; truncated hostname:'
        end if
        write (io_unit,'(5X,a)') trim(carray_to_fstring(str))

        path = getcwd_c(str, str_len)
        if (c_associated(path)) then
            write (io_unit,'(a20)') 'Working directory: '
        else
            write (io_unit,'(a63)') 'Working directory exceeds 255 characters; truncated directory:'
        end if
        write (io_unit,'(5X,a)') trim(carray_to_fstring(str))

        call date_and_time(VALUES=date_values)

        write (io_unit,'(1X,a18,1X,i2.2,"/",i2.2,"/",i4.4,1X,a2,1X,i2.2,2(":",i2.2))') &
                   "Started running on", date_values(3:1:-1), "at", date_values(5:7)
        call get_uuid(GLOBAL_UUID)
        write (io_unit,'(1X,"Calculation UUID:",1X,a36,".")') GLOBAL_UUID

        write (io_unit,'(1X,"Preprocessor settings:")')
        ! Sadly preprocessor does not retokenize the output, so have to do this by hand...
#ifdef DISABLE_HDF5
        write (io_unit,'(5X,"DISABLE_HDF5 defined.  HDF5 disabled.")')
#else
        write (io_unit,'(5X,"DISABLE_HDF5 not defined.  HDF5 enabled.")')
#endif
#ifdef DISABLE_LANCZOS
        write (io_unit,'(5X,"DISABLE_LANCZOS defined.  Lanczos disabled.")')
#else
        write (io_unit,'(5X,"DISABLE_LANCZOS not defined.  Lanczos enabled.")')
#endif
#ifdef DISABLE_UUID
        write (io_unit,'(5X,"DISABLE_UUID defined.  UUID disabled.")')
#else
        write (io_unit,'(5X,"DISABLE_UUID not defined.  UUID enabled.")')
#endif
#ifdef PARALLEL
        write (io_unit,'(5X,"PARALLEL defined.  MPI parallelization enabled.")')
#else
        write (io_unit,'(5X,"PARALLEL not defined.  MPI parallelization disabled.")')
#endif
#ifdef DISABLE_SCALAPACK
        write (io_unit,'(5X,"DISABLE_SCALAPACK defined.  ScaLAPACK disabled.")')
#else
        write (io_unit,'(5X,"DISABLE_SCALAPACK not defined.  ScaLAPACK enabled.")')
#endif
#ifdef SINGLE_PRECISION
        write (io_unit,'(5X,"SINGLE_PRECISION defined.  Single precision used where relevant.")')
#else
        write (io_unit,'(5X,"SINGLE_PRECISION not defined.  Double precision used throughout.")')
#endif
#ifdef USE_POPCNT
        write (io_unit,'(5X,"USE_POPCNT defined.  Fortran 2003 POPCNT procedure used.")')
#else
        write (io_unit,'(5X,"USE_POPCNT not defined.  Internal POPCNT procedure used.")')
#endif
        write (io_unit,'(5X,"DET_SIZE = ",i2,".")') bit_size(0_i0)
        write (io_unit,'(5X,"POP_SIZE = ",i2,".")') bit_size(0_int_p)

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

end module report

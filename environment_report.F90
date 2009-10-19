module report

! James Spencer, CUC3, University of Cambridge.
!
! Copyright (c) 2009 James Spencer.
!
! environment_report prints out a summary of the environment:
!   * when the code was compiled;
!   * the VCS BASE repository version (there's no guarantee that the code wasn't
!     changed!);
!   * whether the working directory contains local changes;
!   * the working directory;
!   * the host computer.
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
    ! * NAGF95: define if using the nag compiler so that certain nag-only modules are
    !   used.
    ! * _WORKING_DIR_CHANGES: define if the working directory contains local
    !   (uncommitted) changes.
    ! * _VCS_VER: set to a (quoted!) string containing the VCS revision id of the current
    !   commit.
    ! 
    ! The VCS information can be simply obtained via shell commands in the makefile:
    ! see the example provided for the commands needed for git and subversion
    ! repositories.  For instance, if a subversion repository is being used, then the
    ! version id can be obtained using::
    !
    !     VCS_VER:=$(shell echo -n \"`svn info | grep 'Revision'| sed -e 's/Revision: //'`\")
    ! 
    ! and if the working directory contains local changes, then the command::
    ! 
    !     svn st -q | xargs -i test -z {}
    ! 
    ! gives a non-zero return code.  Hence the variable WORKING_DIR_CHANGES is 
    ! equal to -D_WORKING_DIR_CHANGES or is a null string depending on whether there
    ! are local changes or not::
    ! 
    !    WORKING_DIR_CHANGES:=$(shell svn st -q | xargs -i test -z {} || echo -n "-D_WORKING_DIR_CHANGES")
    ! 
    ! environment_report.F90 can thus be appropriately compiled with the command::
    ! 
    !    $(FC) $(WORKING_DIR_CHANGES) -D_VCS_VER='$(VCS_VER)' -c environment_report.F90 -o environment_report.o
    !   
    ! where $(FC) is defined to be the desired fortran compiler.
    ! 
    ! Similarly for git::
    ! 
    !     VCS_VER:=$(shell echo -n \" && git log --max-count=1 --pretty=format:%H && echo -n \")
    !     WORKING_DIR_CHANGES := $(shell git diff-index --quiet --cached HEAD
    !                             --ignore-submodules -- && git diff-files --quiet
    !                             --ignore-submodules || echo -n "-D_WORKING_DIR_CHANGES") # on one line. 
    ! 
    ! environment_report has been tested with gfortran, g95, ifort, pgf90, nag and
    ! pathscale.  Note that hostnm and getcwd are intrinsic functions (at least
    ! according to gnu documentation), unless you have to use nag.  There should also
    ! exist the intrinsic subroutines, but pgf90 and ifort gave segmentation faults
    ! when they were used, so we use intrinsic functions unless using nag.

#ifdef NAGF95
    use f90_unix_dir, only: getcwd
    use f90_unix_env, only: gethostname
#endif

    implicit none
    integer, intent(in), optional :: io
    integer :: stat,io_unit
#ifndef NAGF95
    integer :: getcwd,hostnm
#endif
    character(255) :: dirname,host

    if (present(io)) then
        io_unit = io
    else
        io_unit = 6
    end if

    write (io_unit,'(1X,64("="))')

#ifndef NAGF95
    ! Stupid nag: it uses fpp, which doesn't define __DATE__ and __TIME__ like cpp
    ! does.
    write (io_unit,'(a13,a,a4,a)') 'Compiled on ',__DATE__,'at ',__TIME__
#endif

    write (io_unit,'(a29,/,5X,a)') 'VCS BASE repository version:',_VCS_VER
#ifdef _WORKING_DIR_CHANGES
    write (io_unit,'(a42)') 'Working directory contains local changes.'
#endif

#ifdef NAGF95
    call getcwd(dirname,errno=stat)
#else
    stat=getcwd(dirname)
#endif

    if (stat.eq.0) then
        write (io_unit,'(a20)') 'Working directory: '
        write (io_unit,'(5X,a)') trim(dirname)
    end if

#ifdef NAGF95
    call gethostname(host)
    stat=0
#else
    stat=hostnm(host)
#endif

    if (stat.eq.0) then
        write (io_unit,'(a13,a)') 'Running on: ',trim(host)
    end if
    write (io_unit,'(1X,64("="),/)')

    return 

    end subroutine environment_report

end module report

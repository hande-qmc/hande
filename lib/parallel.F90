module parallel

! Wrappers for parallel tasks.
! This allows (in principle) much of the code to be identical for both the 
! serial and parallel cases and avoids the rest of the code being littered
! with preprocessing statements.

#ifdef _PARALLEL
use mpi
#endif

implicit none

! Processor id/rank.
integer :: iproc

! Number of processors.
integer :: nprocs

! True if the processor is the root (i.e. iproc==0) processor.
! In particular, output should only be performed on the parent processor.
logical :: parent

contains

    subroutine init_parallel()

        ! Initialise the parallel environment.
        ! If in serial mode, then the appropriate dummy module variables are
        ! set.

        use errors
        integer :: ierr

#ifdef _PARALLEL
        call mpi_init(ierr)
        if (ierr /= mpi_success) call stop_all('init_parallel','Error initialising MPI.')
        call mpi_comm_size(mpi_comm_world, nprocs, ierr)
        call mpi_comm_rank(mpi_comm_world, iproc, ierr)
#else
        iproc = 0
        nprocs = 1
#endif

        parent = iproc == 0 

    end subroutine init_parallel

    subroutine end_parallel()

        ! Terminate the parallel environment.
        ! This is just a empty procedure in serial mode.
        use errors
        integer :: ierr

#if _PARALLEL
        call mpi_finalize(ierr)
        if (ierr /= mpi_success) call stop_all('end_parallel','Error terminating MPI.')
#endif

    end subroutine end_parallel

    subroutine stop_all_processors(error_code, error_string)

        ! Abort all processors.
        ! This is just a empty procedure in serial mode.
        ! In:
        !    error_code: error code given to mpi_abort which returns it to the
        !        invoking environment.
        !    error_string: string (length 3) containing error code.  Used as an
        !        argument to stop.

        character(3), intent(in) :: error_string
        integer :: error_code, ierr

#if _PARALLEL
        call mpi_abort(mpi_comm_world, error_code, ierr)
        stop error_string
#endif

    end subroutine stop_all_processors

end module parallel

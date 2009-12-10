module parallel
! Just a dummy module for now.

#ifdef _PARALLEL
use mpi
#endif

implicit none

integer :: iproc, nprocs
logical :: parent

contains

    subroutine init_parallel()

        use errors

#ifdef _PARALLEL
        call mpi_init(ierr)
        if (ierr /= MPI_SUCCESS) call stop_all('init_parallel','Error initialising MPI.')
        call mpi_comm_size(mpi_comm_world, nprocs, ierr)
        call mpi_comm_rank(mpi_comm_world, iproc, ierr)
#else
        iproc = 0
        nprocs = 1
#endif

        parent = iproc == 0 

    end subroutine init_parallel

end module parallel

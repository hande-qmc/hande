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

! Processor grid dimensions (for use with scalapack and blacs)
integer :: nproc_cols, nproc_rows

! True if the processor is the root (i.e. iproc==0) processor.
! In particular, output should only be performed on the parent processor.
logical :: parent

! blacs and scalapack split a matrix up into n x n blocks which are then
! distributed around the processors in a cyclic fashion.
! The block size is critical to performance.  64 seems to be a good value (see
! scalapack documentation).
integer, parameter :: block_size = 64

! Type for storing information about a processor as used in blacs and scalapack.
type blacs_info
    ! Location of the processor within the grid.
    integer :: procx, procy
    ! Number of rows and columns of the global matrix stored on the processor.
    integer :: nrows, ncols
    ! blacs and scalapack use a 9 element integer array as a description of how
    ! a matrix is distributed throughout the processor grid.
    ! Store two descriptors: desca refers to the matrix and descz to the output
    ! eigenvector matrix.
    integer :: my_desca(9), my_descz(9)
end type blacs_info

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

    function get_blacs_info(matrix_size) result(my_blacs_info)

        ! In:
        !    matrix_size: leading dimension of the square array to be
        !                 distributed amongst the processor mesh.
        ! Returns:
        !    my_blacs_info: derived type containing information about the
        !                   processor within the blacs mesh.
        !                   In serial mode a dummy result is returned containing
        !                   the appropriate information apart from the
        !                   descriptor arrays, which are set to zero.

        type(blacs_info) :: my_blacs_info
        integer, intent(in) :: matrix_size
        integer :: context
        integer :: i, j, k
        integer :: numroc ! scalapack function 
        integer :: procy, procx, nrows, ncols
        integer :: desca(9), descz(9)
        integer :: ierr

        ! Set processor grid dimensions.
        ! Square is good.  Make it as "square" as possible.
        i = int(sqrt(real(nprocs,8)))
        do j = i, 1, -1
            k = nprocs/j
            if (j*k == nprocs) then
                nproc_cols = j
                nproc_rows = k
                exit
            end if
        end do

#ifdef _PARALLEL
        ! Initialise a single BLACS context.
        call blacs_get(0, 0, context)
        ! Initialise processor grid, nproc_rows x nproc_cols
        call blacs_gridinit(context, 'Row-major', nproc_rows, nproc_cols)

        ! Get the information about this processor: what is its coordinate
        ! (procx, procy) in the processor grid?
        call blacs_gridinfo(context, nproc_rows, nproc_cols, procx, procy)

        ! Find how many rows and columns are owned by the processor.
        nrows = numroc(matrix_size, block_size, procx, 0, nproc_rows)
        ncols = numroc(matrix_size, block_size, procy, 0, nproc_cols)

        ! Initialise the descriptor vectors needed for scalapack procedures.
        call descinit(desca, matrix_size, matrix_size, block_size, block_size, 0, 0, context, nrows, ierr)
        call descinit(descz, matrix_size, matrix_size, block_size, block_size, 0, 0, context, nrows, ierr)

        my_blacs_info = blacs_info(procx, procy, nrows, ncols, desca, descz)
#else
        desca = 0
        descz = 0
        my_blacs_info = blacs_info(0, 0, 1, 1, desca, descz)
#endif

    end function get_blacs_info

end module parallel

module parallel

! Wrappers for parallel tasks.
! This allows (in principle) much of the code to be identical for both the
! serial and parallel cases and avoids the rest of the code being littered
! with preprocessing statements.

#ifdef PARALLEL
use mpi
#endif

use const

implicit none

!=== Processor ID ===

!--- MPI ---

! Processor id/rank.
integer :: iproc

! Number of processors.
integer :: nprocs

! Choose root processor to have rank 0.
integer, parameter :: root = 0

! True if the processor is the root (i.e. iproc==root) processor.
! In particular, output should only be performed on the parent processor.
logical :: parent

!--- OpenMP ---

! WARNING: we assume that we *never* do any I/O that is not either in a (e.g.)
! master block or outside of the OpenMP region.
integer :: nthreads

!=== Comms info ===

!--- BLACS/ScaLAPACK ---

! Type for storing information about a processor as used in blacs and scalapack.
! Conveniently filled by the function get_blacs_info.
type blacs_info
    ! Location of the processor within the grid.  Negative if the processor
    ! isn't involved in the grid.
    integer :: procx, procy
    ! Processor grid dimensions
    integer :: nproc_cols, nproc_rows
    ! Number of rows and columns of the global matrix stored on the processor.
    ! Zero if no rows/columns stored on the processor.
    integer :: nrows, ncols
    ! size of n x n blocks of a array which are cyclicly distributed around processors.
    integer :: block_size
    ! blacs and scalapack use a 9 element integer array as a description of how
    ! a matrix is distributed throughout the processor grid.
    ! Descriptor for distributing an NxN matrix:
    integer :: desc_m(9)
    ! Descriptor for distributing vector of length N:
    integer :: desc_v(9)
end type blacs_info

!--- Information about MPI timing ---

! Type for holding MPI barrier and communication times.
type parallel_timing_t
    ! If true then check the MPI barrier time before communication.
    logical :: check_barrier_time = .false.
    ! The amount of time spent waiting for other processes before an MPI call
    ! can be made.
    real(p) :: barrier_time = 0.0_p
    ! The total time taken for MPI communications.
    real(p) :: comm_time = 0.0_p
end type parallel_timing_t

!--- MPI Communication data types ---

! MPI data type for 32-bit or 64-bit integer used in determinant bit arrays.
#ifdef PARALLEL
#if DET_SIZE == 64
integer, parameter :: mpi_det_integer = MPI_INTEGER8
#else
! Use 32-bits by default.  (Note that const.F90 throws a compile-time error if DET_SIZE is not 32 or 64 bits...)
integer, parameter :: mpi_det_integer = MPI_INTEGER
#endif

! MPI data type for 32-bit or 64-bit integers used to store walker populations.
#if POP_SIZE == 32
integer, parameter :: mpi_pop_integer = MPI_INTEGER
#elif POP_SIZE == 64
integer, parameter :: mpi_pop_integer = MPI_INTEGER8
#else
! Use 32-bit integers by default.
integer, parameter :: mpi_pop_integer = MPI_INTEGER
#endif

! MPI data type for 32-bit or 64-bit integers used in the sdata component of spawn_t.
#if POP_SIZE == 64 || DET_SIZE == 64
integer, parameter :: mpi_sdata_integer = MPI_INTEGER8
#else
integer, parameter :: mpi_sdata_integer = MPI_INTEGER
#endif

! MPI data type for reals of single/double precision (as chosen at
! compile-time).
#ifdef SINGLE_PRECISION
integer, parameter :: mpi_preal = MPI_REAL
#else
integer, parameter :: mpi_preal = MPI_REAL8
#endif
#endif

contains

    subroutine init_parallel()

        ! Initialise the parallel environment.
        ! If in serial mode, then the appropriate dummy module variables are
        ! set.

        use omp_lib
        use errors

#ifdef PARALLEL
        integer :: ierr

        call mpi_init(ierr)
        if (ierr /= mpi_success) call stop_all('init_parallel','Error initialising MPI.')
        call mpi_comm_size(mpi_comm_world, nprocs, ierr)
        call mpi_comm_rank(mpi_comm_world, iproc, ierr)
#else
        iproc = 0
        nprocs = 1
#endif
#ifdef _OPENMP
        nthreads = omp_get_max_threads()
#else
        nthreads = 1
#endif

        parent = iproc == root

    end subroutine init_parallel

    subroutine parallel_report()

        ! Print out information about the parallel environment.

        use utils, only: int_fmt

        character(4) :: fmt1, fmt2

        fmt1 = int_fmt(nprocs,1)
        fmt2 = int_fmt(nthreads,1)
        if (nprocs == 1) then
            write (6,'(1X,"Number of MPI processes running on:"'//fmt1//')') nprocs
            if (nthreads == 1) then
                write (6,'(1X,"Running with"'//fmt2//'" thread.",/)') nthreads
            else
                write (6,'(1X,"Running with"'//fmt2//'" threads.",/)') nthreads
            end if
        else
            write (6,'(1X,"Number of MPI processes running on:"'//fmt1//')') nprocs
            if (nthreads == 1) then
                write (6,'(1X,"Running with"'//fmt2//'" thread per MPI process.",/)') nthreads
            else
                write (6,'(1X,"Running with"'//fmt2//'" threads per MPI process.",/)') nthreads
            end if
        end if

    end subroutine parallel_report

    subroutine end_parallel()

        ! Terminate the parallel environment.
        ! This is just a empty procedure in serial mode.
        use errors

#ifdef PARALLEL
        integer :: ierr

        call mpi_finalize(ierr)
        if (ierr /= mpi_success) call stop_all('end_parallel','Error terminating MPI.')
#endif

    end subroutine end_parallel

    function get_blacs_info(matrix_size, block_size, proc_grid) result(my_blacs_info)

        ! In:
        !    matrix_size: leading dimension of the square array to be
        !                 distributed amongst the processor mesh.
        !   block_size: size of n x n blocks of a array which are cyclicly
        !                 distributed around processors.
        !   proc_grid(2) (optional): set the processor mesh to be
        !                 (/ nproc_cols, nproc_rows /).  If not present then
        !                 a mesh as close to square as possible is used.
        ! Returns:
        !    my_blacs_info: derived type containing information about the
        !                   processor within the blacs mesh.
        !                   In serial mode a dummy result is returned containing
        !                   the appropriate information apart from the
        !                   descriptor arrays, which are set to zero.

        type(blacs_info) :: my_blacs_info
        integer, intent(in) :: matrix_size, block_size
        integer, intent(in), optional :: proc_grid(2)
        integer :: i, j, k, nproc_rows, nproc_cols
        integer :: procy, procx, nrows, ncols
        integer :: desc_m(9), desc_v(9)
#ifdef PARALLEL
        integer :: numroc ! scalapack function
        integer :: ierr
        integer :: context
#endif

        ! Set processor grid dimensions.
        if (present(proc_grid)) then
            nproc_cols = proc_grid(1)
            nproc_rows = proc_grid(2)
        else
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
        end if

#ifdef PARALLEL
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
        call descinit(desc_m, matrix_size, matrix_size, block_size, block_size, 0, 0, context, max(1,nrows), ierr)
        call descinit(desc_v, matrix_size, 1, block_size, block_size, 0, 0, context, max(1,nrows), ierr)
#else
        procx = 0
        procy = 0
        nproc_cols = 1
        nproc_rows = 1
        nrows = matrix_size
        ncols = matrix_size
        desc_m = 0
        desc_v = 0
#endif

        my_blacs_info = blacs_info(procx, procy, nproc_cols, nproc_rows, nrows, ncols, block_size, desc_m, desc_v)

    end function get_blacs_info

    function get_thread_id() result(id)

        ! Returns:
        !    Thread index.  This is just 0 outside an OpenMP region,
        !    omp_get_thread_num() is returned inside an OpenMP region.

        use omp_lib

        integer :: id

#ifdef _OPENMP
        id = omp_get_thread_num()
#else
        id = 0
#endif

    end function get_thread_id

end module parallel

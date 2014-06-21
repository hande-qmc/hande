module csr

! Handling of sparse matrices in CSR format.
! This is not intended to be a complete implementation, but rather only include
! procedures as required in HANDE.

! Parallelisation (beyond perhaps OpenMP) is beyond the current scope and not
! all procedures have threading so they can be called from regions already
! parallelised using MPI and/or OpenMP.

! We use the abbreviations nrow and nnz to stand for the number of rows of and
! the number of non-zero elements of a matrix, respectively.

use const, only: p

implicit none

! real(p) matrix
type csrp_t
    ! For matrices which will use the symmetric routines, this stores the
    ! non-zero elements of upper *or* lower triangle of matrix.
    ! For general matrices, all of the non-zero elements are stored.
    real(p), allocatable :: mat(:) ! (nnz)
    ! Column index of values stored in mat.
    ! If M_{ij} is stored in mat(k), then col_ind(k) = j.
    integer, allocatable :: col_ind(:) ! (nnz)
    ! row_ptr(i) gives the location in mat of the first non-zero element in row
    ! i, i.e. if mat(k) = M_{ij}, then row_ptr(i) <= k < row_ptr(i+1).  By
    ! convention, row_ptr(nrow+1) = nnz+1.
    integer, allocatable :: row_ptr(:) ! (nrow+1)
    ! WARNING: if the matrix is symmetric, it is assumed that the programmer
    ! only stores the uppper *or* lower triangle of the matrix.
    logical :: symmetric = .false.
end type csrp_t

contains

    subroutine init_csrp(spm, nrow, nnz, symmetric)

        ! Initialise a real(p) sparse symmetrix matrix in csr format.

        ! In:
        !    nrow: number of rows in matrix.
        !    nnz: number of non-zero elements in upper/lower triangle of matrix.
        !    symmetric: (optional, default: false).
        ! Out:
        !    spm: csrp_t object with components correctly allocated and
        !         row_ptr(nrow+1) set to nnz+1.

        use checking, only: check_allocate

        type(csrp_t), intent(out) :: spm
        integer, intent(in) :: nrow, nnz
        logical, intent(in), optional :: symmetric

        integer :: ierr

        allocate(spm%mat(nnz), stat=ierr)
        call check_allocate('spm%mat', nnz, ierr)
        allocate(spm%col_ind(nnz), stat=ierr)
        call check_allocate('spm%col_ind', nnz, ierr)
        allocate(spm%row_ptr(nrow+1), stat=ierr)
        call check_allocate('spm%row_ptr', nrow+1, ierr)
        spm%row_ptr(nrow+1) = nnz+1

        if (present(symmetric)) spm%symmetric = symmetric

    end subroutine init_csrp

    subroutine end_csrp(spm)

        ! Deallocate components of a real(p) sparse matrix in csr format.

        ! In/Out:
        !    spm: csrp_t object with all components deallocated on exit.

        use checking, only: check_deallocate

        type(csrp_t), intent(inout) :: spm

        integer :: ierr

        deallocate(spm%mat, stat=ierr)
        call check_deallocate('spm%mat', ierr)
        deallocate(spm%col_ind, stat=ierr)
        call check_deallocate('spm%col_ind', ierr)
        deallocate(spm%row_ptr, stat=ierr)
        call check_deallocate('spm%row_ptr', ierr)

    end subroutine end_csrp

    subroutine csrpsymv_symmetric(spm, x, y)

        ! Calculate y = m*x, where m is a sparse symmetric matrix and x and
        ! y are dense vectors.

        ! In:
        !   spm: sparse symmetric matrix (real(p) in csr format.  See
        !        module-level notes about storage format.
        !   x: dense vector.  Number of elements must be at least the number of
        !      rows in spm.
        ! Out:
        !   y: dense vector.  Holds m*x on exit..  Number of elements must be at
        !      least the number of rows in spm, with all additional elements set to
        !      0.

        ! WARNING: if the matrix is symmetric, it is assumed that the programmer
        ! only stores the uppper *or* lower triangle of the matrix.

        use errors, only: stop_all

        type(csrp_t), intent(in) :: spm
        real(p), intent(in) :: x(:)
        real(p), intent(out) :: y(:)

        integer :: irow, icol, iz
        real(p) :: rowx

        if (.not.spm%symmetric) call stop_all('csrpsymv', 'Sparse matrix not symmetric.')
        
        y = 0.0_p
        ! Avoid overhead of creating thread pool.
        ! However, by not just doing "!$omp parallel do", we must be *very*
        ! careful when accessing or updating any shared value.
        !$omp parallel
        do irow = 1, size(spm%row_ptr)-1
            !$omp master
            rowx = 0.0_p
            !$omp end master
            !$omp barrier
            ! OpenMP chunk size determined completely empirically from a single
            ! test.  Please feel free to improve...
            !$omp do private(icol) reduction(+:rowx) schedule(dynamic, 200)
            do iz = spm%row_ptr(irow), spm%row_ptr(irow+1)-1
                icol = spm%col_ind(iz)
                y(icol) = y(icol) + spm%mat(iz)*x(irow)
                if (icol /= irow) rowx = rowx + spm%mat(iz)*x(icol)
            end do
            !$omp end do
            !$omp master
            y(irow) = y(irow) + rowx
            !$omp end master
        end do
        !$omp end parallel

    end subroutine csrpsymv_symmetric

    subroutine csrpsymv_general(spm, x, y)

        ! Calculate y = m*x, where m is a sparse matrix and x and y are dense
        ! vectors.

        ! In:
        !   spm: sparse matrix (real(p) in csr format. See module-level notes
        !        about storage format.
        !   x: dense vector.  Number of elements must be at least the number of
        !      columns in spm.
        ! Out:
        !   y: dense vector.  Holds m*x on exit..  Number of elements must be at
        !      least the number of rows in spm, with all additional elements set to
        !      0.

        use errors, only: stop_all

        type(csrp_t), intent(in) :: spm
        real(p), intent(in) :: x(:)
        real(p), intent(out) :: y(:)

        integer :: irow, icol, iz

        ! This routine should not be used for symmetric matrices where only the
        ! upper or lower halves of the matrix are stored.
        if (spm%symmetric) call stop_all('csrpsymv', 'Sparse matrix is symmetric.')
        
        y = 0.0_p
        do irow = 1, size(spm%row_ptr)-1
            do iz = spm%row_ptr(irow), spm%row_ptr(irow+1)-1
                icol = spm%col_ind(iz)
                y(irow) = y(irow) + spm%mat(iz)*x(icol)
            end do
        end do

    end subroutine csrpsymv_general

    subroutine csrpsymv_row(spm, x, y_irow, irow)

        ! Calculate a single value in the vector y = m*x, where m is a sparse
        ! matrix and x and y are dense vectors.

        ! In:
        !   spm: sparse matrix (real(p) in csr format. See module-level notes
        !        about storage format.
        !   x: dense vector.  Number of elements must be at least the number of
        !      columns in spm.
        !   irow: The index of the row of the Hamiltonian to multiply with.
        ! Out:
        !   y_irow: Holds \sum_j m_{irow,j}*x_j on exit.

        use errors, only: stop_all

        type(csrp_t), intent(in) :: spm
        real(p), intent(in) :: x(:)
        real(p), intent(out) :: y_irow
        integer, intent(in) :: irow

        integer :: icol, iz

        y_irow = 0.0_p
        do iz = spm%row_ptr(irow), spm%row_ptr(irow+1)-1
            icol = spm%col_ind(iz)
            y_irow = y_irow + spm%mat(iz)*x(icol)
        end do

    end subroutine csrpsymv_row

end module csr

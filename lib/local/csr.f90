module csr

! Handling of sparse matrices in CSR format.
! This is not intended to be a complete implementation, but rather only include
! procedures as required in HANDE.

! Parallelisation (beyond perhaps OpenMP) is beyond the current scope.

! We use the abbreviations nrow and nnz to stand for the number of rows of and
! the number of non-zero elements of a matrix, respectively.

use const, only: p

implicit none

! Symmetric real(p) matrix
type csrpsy_t
    ! non-zero elements of upper *or* lower triangle of matrix.
    ! WARNING: it is assumed that the programmer only stores one-half of the
    ! symmetric matrix.
    real(p), allocatable :: mat(:) ! (nnz)
    ! Column index of values stored in mat.
    ! If M_{ij} is stored in mat(k), then col_ind(k) = j.
    integer, allocatable :: col_ind(:) ! (nnz)
    ! row_ptr(i) gives the location in mat of the first non-zero element in row
    ! i, i.e. if val(k) = M_{ij}, then row_ptr(i) <= k < row_ptr(i).  By
    ! convention, row_ptr(nrow+1) = nnz+1.
    integer, allocatable :: row_ptr(:) ! (nrow+1)
end type csrpsy_t

contains

    subroutine init_csrpsy(spm, nrow, nnz)

        ! Initialise a real(p) sparse symmetrix matrix in csr format.

        ! In:
        !    nrow: number of rows/columns in matrix.
        !    nnz: number of non-zero elements in upper/lower triangle of matrix.
        ! Out:
        !    spm: csrpsy_t object with components correctly allocated and
        !         row_ptr(nrow+1) set to nnz+1.

        use checking, only: check_allocate

        type(csrpsy_t), intent(out) :: spm
        integer, intent(in) :: nrow, nnz

        integer :: ierr

        allocate(spm%mat(nnz), stat=ierr)
        call check_allocate('spm%mat', nnz, ierr)
        allocate(spm%col_ind(nnz), stat=ierr)
        call check_allocate('spm%col_ind', nnz, ierr)
        allocate(spm%row_ptr(nrow+1), stat=ierr)
        call check_allocate('spm%row_ptr', nrow+1, ierr)
        spm%row_ptr(nrow+1) = nnz+1

    end subroutine init_csrpsy

    subroutine end_csrpsy(spm)

        ! Deallocate components of a real(p) sparse symmetric matrix in csr
        ! format.

        ! In/Out:
        !    spm: csrpsy_t object with all components deallocated on exit.

        use checking, only: check_deallocate

        type(csrpsy_t), intent(inout) :: spm

        integer :: ierr

        deallocate(spm%mat, stat=ierr)
        call check_deallocate('spm%mat', ierr)
        deallocate(spm%col_ind, stat=ierr)
        call check_deallocate('spm%col_ind', ierr)
        deallocate(spm%row_ptr, stat=ierr)
        call check_deallocate('spm%row_ptr', ierr)

    end subroutine end_csrpsy

    pure subroutine csrpsymv(spm, x, y)

        ! Calculate y = m*x, where m is a sparse matrix and x and y are dense
        ! vectors.

        ! In:
        !   spm: sparse symmetric matrix (real(p) in csr format.  See
        !        module-level notes about storage format.
        !   x: dense vector.  Number of elements must be at least the number of
        !      rows in spm.
        ! Out:
        !   y: dense vector.  Holds m*x on exit..  Number of elements must be at
        !      least the number of rows in spm, with all additional elements set to
        !      0.

        type(csrpsy_t), intent(in) :: spm
        real(p), intent(in) :: x(:)
        real(p), intent(out) :: y(:)

        integer :: irow, icol, iz
        
        y = 0.0_p
        do irow = 1, size(spm%row_ptr)-1
            ! TODO: OpenMP.
            do iz = spm%row_ptr(irow), spm%row_ptr(irow+1)-1
                icol = spm%col_ind(iz)
                y(icol) = y(icol) + spm%mat(iz)*x(irow)
                if (icol /= irow) y(irow) = y(irow) + spm%mat(iz)*x(icol)
            end do
        end do

    end subroutine csrpsymv

end module csr

module fci_utils

implicit none

contains

    subroutine generate_hamil(sys, ndets, dets, hamil, hamil_csr, proc_blacs_info, full_mat)

        ! Generate a symmetry block of the Hamiltonian matrix, H = < D_i | H | D_j >.
        ! The list of determinants, {D_i}, is grouped by symmetry and contains
        ! only determinants of a specified spin.
        ! Only generate the upper diagonal for use with (sca)lapack and Lanczos routines.
        ! In:
        !    sys: system to be studied.
        !    ndets: number of determinants in the Hilbert space.
        !    dets: list of determinants in the Hilbert space (bit-string representation).
        !    sparse_mode: if true, generate the Hamiltonian in sparse format (in hamil_csr)
        !        rather than in non-sparse format (in hamil).
        !    distribute_mode (optional): flag for determining how the
        !        Hamiltonian matrix is distributed among processors.  It is
        !        irrelevant if only one processor is used: the distribution schemes
        !        all reduce to the storing the entire matrix on the single
        !        processor.  Can take the values given by the distribute_off,
        !        distribute_blocks and distribute_cols parameters.  See above
        !        for descriptions of the different behaviours.
        !    full_mat (optional): if present and true generate the full matrix rather than
        !        just storing one triangle.

        use checking, only: check_allocate, check_deallocate
        use csr, only: init_csrp, end_csrp, csrp_t
        use utils, only: get_free_unit
        use errors, only: stop_all
        use parallel

        use hamiltonian, only: get_hmatel
        use real_lattice
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        integer, intent(in) :: ndets
        integer(i0), intent(in) :: dets(:,:)
        real(p), intent(out), allocatable, optional :: hamil(:,:)
        type(csrp_t), intent(out), optional :: hamil_csr
        type(blacs_info), intent(in), optional :: proc_blacs_info
        logical, intent(in), optional :: full_mat
        integer :: ierr, iunit, n1, n2
        integer :: i, j, ii, jj, ilocal, iglobal, jlocal, jglobal, nnz, imode
        logical :: sparse_mode
        real(p) :: hmatel

        sparse_mode = present(hamil_csr) .and. .not.present(hamil)
        if (.not.present(hamil) .and. .not.present(hamil_csr)) &
            call stop_all('generate_hamil', 'Must supply either hamil or hamil_csr in argument list.')

        if (sparse_mode .and. present(proc_blacs_info)) then
            call stop_all('generate_hamil', &
                'Sparse distributed matrices are not currently implemented.  &
                &If this is disagreeable to you, please contribute patches resolving this situation.')
        end if

        if (present(hamil)) then
            if (present(proc_blacs_info)) then
                allocate(hamil(proc_blacs_info%nrows,proc_blacs_info%ncols), stat=ierr)
            else
                allocate(hamil(ndets, ndets), stat=ierr)
            end if
            call check_allocate('hamil', size(hamil), ierr)
        end if

        ! Form the Hamiltonian matrix < D_i | H | D_j >.
        if (present(proc_blacs_info)) then
            ! blacs divides the matrix up into sub-matrices of size block_size x block_size.
            ! The blocks are distributed in a cyclic fashion amongst the
            ! processors.
            ! Each processor stores a total of nrows of the matrix.
            ! i gives the index of the first row in the current block.
            ! ii gives the index of the current row within the current block.
            ! The local i index is thus i-1+ii.  This is used to refer to the
            ! matrix element as stored (continuously) on the processor.
            ! The global i index is given by the sum of the rows held on
            ! preceeding proecssors for the previous loops over processors (as
            ! done in the loop over i), the rows held on other processors during
            ! the corrent loop over processors and the position within the
            ! current block.
            ! Similarly for the other indices.
            do i = 1, proc_blacs_info%nrows, block_size
                do ii = 1, min(block_size, proc_blacs_info%nrows - i + 1)
                    ilocal = i - 1 + ii
                    iglobal = (i-1)*nproc_rows + proc_blacs_info%procx*block_size + ii
                    do j = 1, proc_blacs_info%ncols, block_size
                        do jj = 1, min(block_size, proc_blacs_info%ncols - j + 1)
                            jlocal = j - 1 + jj
                            jglobal = (j-1)*nproc_cols + proc_blacs_info%procy*block_size + jj
                            hamil(ilocal, jlocal) = get_hmatel(sys, dets(:,iglobal), dets(:,jglobal))
                        end do
                    end do
                end do
            end do
        else if (sparse_mode) then
                ! First, find out how many non-zero elements there are.
                ! We'll be naive and not just test for non-zero by symmetry, but
                ! actually calculate all matrix elements for now.
                ! Then, store the Hamiltonian.
                do imode = 1, 2
                    nnz = 0
                    !$omp parallel
                    do i = 1, ndets
                        ii = i
                        if (present(full_mat)) then
                            if (full_mat) ii = 1
                        end if
                        ! OpenMP chunk size determined completely empirically
                        ! from a single test.  Please feel free to improve...
                        !$omp do private(j, hmatel) schedule(dynamic, 200)
                        do j = ii, ndets
                            hmatel = get_hmatel(sys, dets(:,i), dets(:,j))
                            if (abs(hmatel) > depsilon) then
                                !$omp critical
                                nnz = nnz + 1
                                if (imode == 2) then
                                    hamil_csr%mat(nnz) = hmatel
                                    hamil_csr%col_ind(nnz) = j
                                    if (hamil_csr%row_ptr(i) == 0) hamil_csr%row_ptr(i) = nnz
                                end if
                                !$omp end critical
                            end if
                        end do
                        !$omp end do
                    end do
                    !$omp end parallel
                    if (imode == 1) then
                        write (6,'(1X,a50,i8/)') 'Number of non-zero elements in Hamiltonian matrix:', nnz
                        call init_csrp(hamil_csr, ndets, nnz, .true.)
                        hamil_csr%row_ptr(1:ndets) = 0
                    else
                        ! Any element not set in row_ptr means that the
                        ! corresponding row has no non-zero elements.
                        ! Set it to be identical to the next row (this avoids
                        ! looping over the zero-row).
                        do i = ndets, 1, -1
                            if (hamil_csr%row_ptr(i) == 0) hamil_csr%row_ptr(i) = hamil_csr%row_ptr(i+1)
                        end do
                    end if
                end do
        else
            !$omp parallel
            do i = 1, ndets
                ii = i
                if (present(full_mat)) then
                    if (full_mat) ii = 1
                end if
                !$omp do private(j) schedule(dynamic, 200)
                do j = ii, ndets
                    hamil(i,j) = get_hmatel(sys,dets(:,i),dets(:,j))
                end do
                !$omp end do
            end do
            !$omp end parallel
        end if

    end subroutine generate_hamil

end module fci_utils

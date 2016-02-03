module fci_lanczos

use const

implicit none

contains

    subroutine do_fci_lanczos(sys, fci_in, ref_in, sparse_hamil)

        ! Perform an FCI calculation via the Lanczos algorithm.

        ! In:
        !    sys: system of interest.
        !    fci_in: fci input options.
        !    ref_in: reference determinant defining (if relevant) a
        !        truncated Hilbert space.

        use fci_utils, only: fci_in_t, init_fci, generate_hamil, write_hamil, hamil_t
        use hamiltonian, only: get_hmatel
        use qmc_data, only: reference_t
        use reference_determinant, only: copy_reference_t
        use system, only: sys_t, copy_sys_spin_info

        use checking, only: check_allocate
        use csr, only: csrp_t
        use errors, only: warning, stop_all
        use parallel, only: parent, nprocs, blacs_info, get_blacs_info
        use utils, only: int_fmt
        use check_input, only: check_fci_opts

        type(sys_t), intent(inout) :: sys
        type(fci_in_t), intent(in) :: fci_in
        type(reference_t), intent(in) :: ref_in
        logical, intent(in) :: sparse_hamil

        type(sys_t) :: sys_bak
        type(reference_t) :: ref
        integer(i0), allocatable :: dets(:,:)
        real(p), allocatable :: eigv(:)
        integer :: ndets, ierr, i, nfound, block_size
        type(blacs_info) :: proc_blacs_info
        type(hamil_t), allocatable :: hamil

        if (parent) call check_fci_opts(sys, fci_in, .true.)

        call copy_sys_spin_info(sys, sys_bak)
        call copy_reference_t(ref_in, ref)

        if (nprocs > 1) then
            if (fci_in%direct_lanczos) call stop_all('do_fci_lanczos', &
                                                            'Direct Lanczos not implemented with MPI parallelism.')
            if (sparse_hamil) call stop_all('do_fci_lanczos', &
                                            'Lanczos with sparse matrices not implemented with MPI parallelism.')
        end if

        call init_fci(sys, fci_in, ref, ndets, dets)

        allocate(eigv(fci_in%nlanczos_eigv), stat=ierr)
        call check_allocate('eigv', fci_in%nlanczos_eigv, ierr)

        ! TRLan assumes that the only the rows of the matrix are distributed.  Furthermore,
        ! it seems TRLan assumes all processors store at least some of the matrix.
        ! Useful to have this for consistency, even as a dummy object in serial.
        block_size = fci_in%block_size
        if (nprocs*block_size > ndets) then
            if (parent) then
                call warning('do_fci_lanczos','Reducing block size so that all processors contain at least a single row.',3)
                write (6,'(1X,"Consider running on fewer processors or reducing block size in input.")')
                write (6,'(1X,"Old block size was:"'//int_fmt(block_size,1)//')') block_size
            end if
            block_size = ndets/nprocs
            if (parent) write (6,'(1X,"New block size is:"'//int_fmt(block_size,1)//')') block_size
        end if
        proc_blacs_info = get_blacs_info(ndets, block_size, [1, nprocs])

        if (ndets == 1) then
            ! The trivial case seems to trip things up...
            eigv(1) = get_hmatel(sys, dets(:,1), dets(:,1))
            nfound = 1
        else
            if (nprocs == 1) then
                if (sparse_hamil) then
                    call generate_hamil(sys, ndets, dets, hamil, use_sparse_hamil = sparse_hamil)
                else
                    call generate_hamil(sys, ndets, dets, hamil)
                end if
            else
                call generate_hamil(sys, ndets, dets, hamil, proc_blacs_info=proc_blacs_info)
            end if
            if (fci_in%write_hamiltonian) call write_hamil(fci_in%hamiltonian_file, ndets, proc_blacs_info, hamil)
            call lanczos_diagonalisation(sys, fci_in, dets, proc_blacs_info, nfound, eigv, hamil)
        end if

        if (parent) then
            write (6,'(1X,"Lanczos diagonalisation results")')
            write (6,'(1X,"^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^",/)')
            write (6,'(1X," State",1X,4X,"Energy")')
            do i = 1, nfound
                write (6,'(1X,i6,1X,f18.12)') i, eigv(i)
            end do
            write (6,'()')
        end if

        ! Return sys in an unaltered state.
        call copy_sys_spin_info(sys_bak, sys)

    end subroutine do_fci_lanczos

    subroutine lanczos_diagonalisation(sys, fci_in, dets, proc_blacs_info, nfound, eigv, hamil)

        ! Perform a Lanczos diagonalisation of the current (spin) block of the
        ! Hamiltonian matrix.
        ! In:
        !    sys: system to be studied.
        !    fci_in: fci input options.
        !    dets: list of determinants in the Hilbert space in the bit
        !        string representation.
        !    proc_blacs_info: BLACS description of distribution of the Hamiltonian.
        !    hamil (optional): Hamiltonian matrix.
        !    hamil_csr (optional): Hamiltonian matrix in sparse format.
        ! Out:
        !    nfound: number of solutions found from this block.  This is
        !        min(number of determinants with current spin, fci_in%nlanczos_eigv).
        !    eigv(:nfound): Lanczos eigenvalues of the current block of the
        !        Hamiltonian matrix.

#ifndef DISABLE_LANCZOS
        use trl_info
        use trl_interface
        use checking, only: check_allocate, check_deallocate
        use parallel, only: parent

        use operators
#endif

        use csr, only: csrp_t
        use errors, only: stop_all
        use system, only: sys_t
        use fci_utils, only: fci_in_t, hamil_t
        use parallel, only: blacs_info

        type(sys_t), intent(in) :: sys
        type(fci_in_t), intent(in) :: fci_in
        integer(i0), intent(in) :: dets(:,:)
        type(blacs_info), intent(in) :: proc_blacs_info
        integer, intent(out) :: nfound
        real(p), intent(out) :: eigv(:)
        type(hamil_t), intent(in) :: hamil
        logical :: sparse_hamil

        integer :: ndets

#ifdef DISABLE_LANCZOS
        call stop_all('lanczos_diagonalisation','Lanczos diagonalisation disabled at compile-time.')
        ! Avoid compile-time warnings about unset intent(out) variables...
        nfound = 0
        eigv = huge(0.0_p)
#else

        integer, parameter :: lohi = -1
        integer :: mev
        real(dp), allocatable :: eval(:) ! (mev)
        real(dp), allocatable :: evec(:,:) ! (ndets, mev)
        integer :: ierr, nrows, i, nwfn
        type(trl_info_t) :: info

#ifdef SINGLE_PRECISION
        ! TRLan requires a double precision interface; elsewhere we need single precision
        real(p), allocatable :: evec_copy(:,:)
#endif

        if (.not.present(hamil%rmat) .and. .not.present(hamil%mat_sparse)) then
            call stop_all('lanczos_diagonalisation', 'No Hamiltonian supplied.')
        end if
        sparse_hamil = present(hamil%mat_sparse) .and. .not.present(hamil%rmat)

        ndets = ubound(dets, dim=2)

        ! mev: number of eigenpairs that can be stored in eval and evec.
        ! twice the number of eigenvalues to be found is a reasonable default.
        mev = min(2*fci_in%nlanczos_eigv, ndets)

        if (parent) then
            if (fci_in%direct_lanczos) then
                write (6,'(/,1X,a44,/)') 'Performing direct Lanczos diagonalisation...'
            else
                write (6,'(/,1X,a37,/)') 'Performing Lanczos diagonalisation...'
            end if
        end if

        nrows = proc_blacs_info%nrows

        ! Initialise trlan.
        ! info: type(trl_info_t).  Used by trl to store calculation info.
        ! nrows: number of rows of matrix on processor.
        ! fci_in%lanczos_string_len: maximum Lanczos basis size.
        ! lohi: -1 means calculate the smallest eigenvalues first (1 to calculate
        !       the largest).
        ! fci_in%nlanczos_eigv: number of eigenvalues to compute.
        call trl_init_info(info, nrows, fci_in%lanczos_string_len, lohi, min(fci_in%nlanczos_eigv,ndets))

        allocate(eval(mev), stat=ierr)
        call check_allocate('eval',mev,ierr)

        allocate(evec(nrows,mev), stat=ierr)
        call check_allocate('evec',nrows*mev,ierr)

        ! Call Lanczos diagonalizer.
        ! hamil_vector: matrix-vector multiplication routine.
        ! info: created in trl_init_info.
        ! nrows: number of rows of matrix on processor.
        ! mev: number of eigenpairs that can be stored in eval and evec.
        ! eval: array to store eigenvalue
        ! evec: array to store the eigenvectors
        ! lde: the leading dimension of evec (in serial case: ndets).
        if (fci_in%direct_lanczos) then
            call trlan(HPsi_direct, info, nrows, mev, eval, evec, nrows)
        else
            call trlan(HPsi, info, nrows, mev, eval, evec, nrows)
        end if

#ifdef SINGLE_PRECISION
        allocate(evec_copy(nrows,mev), stat=ierr)
        call check_allocate('evec_copy',nrows*mev,ierr)
        evec_copy=evec
#endif

        ! Get info...
        ! dsymv uses 2(N^2+2N) floating point operations.
        ! Ref: Benchmark of the Extended Basic Linear Algebra Subprograms on
        ! the NEC SX-2 Supercomputer, by R. M. Dubash, J. L. Fredin and O. G.
        ! Johnson, Supercomputing, 1st International Conference, Athens,
        ! Greece, June 8-12, 1987, Proceedings, Lecture Notes in Computer
        ! Science, 297 (1987) 894-913.
        ! trl_print_info gathers information from the processors so must be
        ! a global call but only prints information from one processor.
        if (parent) write (6,'(1X,a28,/)') 'TRLan (Lanczos) information:'
        call trl_print_info(info, 2*(ndets**2+2*ndets))
        if (parent) write (6,'()')

        nfound = min(fci_in%nlanczos_eigv,ndets)
        eigv(1:nfound) = real(eval(1:nfound),p)

        nwfn = fci_in%analyse_fci_wfn
        if (nwfn < 0 .or. nwfn > nfound) nwfn = nfound
        do i = 1, nwfn
#ifdef SINGLE_PRECISION
            call analyse_wavefunction(sys, evec_copy(:,i), dets, proc_blacs_info)
#else
            call analyse_wavefunction(sys, evec(:,i), dets, proc_blacs_info)
#endif
        end do
        nwfn = fci_in%print_fci_wfn
        if (nwfn < 0 .or. nwfn > nfound) nwfn = nfound
        do i = 1, nwfn
#ifdef SINGLE_PRECISION
            call print_wavefunction(fci_in%print_fci_wfn_file, dets, proc_blacs_info, evec_copy(:,i))
#else
            call print_wavefunction(fci_in%print_fci_wfn_file, dets, proc_blacs_info, evec(:,i))
#endif
        end do

        deallocate(eval, stat=ierr)
        call check_deallocate('eval',ierr)
        deallocate(evec, stat=ierr)
        call check_deallocate('evec',ierr)
#ifdef SINGLE_PRECISION
        deallocate(evec_copy, stat=ierr)
        call check_deallocate('evec_copy',ierr)
#endif
! end ifdef DISABLE_LANCZOS
#endif

    contains

        subroutine HPsi(nrow, ncol, xin, ldx, yout, ldy)

            ! Matrix-vector multiplication procedure for use with trlan.
            ! In:
            !    nrow: the number of rows on this processor if the problem is distributed
            !        using MPI, otherwise the number of total rows in a Lanczos vector.
            !    ncol: the number of vectors (columns in xin and yout) to be multiplied.
            !    xin: the array to store the input vectors to be multiplied.
            !    ldx: the leading dimension of the array xin when it is declared as
            !       two-dimensional array.
            !    ldy: the leading dimension of the array yout when it is declared as
            !       two-dimensional array.
            ! Out:
            !    yout: the array to store results of the multiplication.

            use csr, only: csrpsymv
            use parallel, only: nprocs

            integer, intent(in) :: nrow, ncol, ldx, ldy
            real(dp), intent(in) :: xin(ldx,ncol)
            real(dp), intent(out) :: yout(ldy,ncol)
#ifdef SINGLE_PRECISION
            real(p) :: xin_p(ldx,ncol)
            real(p) :: yout_p(ldy,ncol)
#endif
            ! local variables
            integer :: i

            ! TRLan requires an interface with xin and yout being double precision.
            ! This is not conducive with running in single precision, where we only
            ! store the Hamiltonian in single precision.  Calling BLAS with
            ! arguments of the wrong kind really screws things up (as you can
            ! imagine...) so we work round it by doing copies to/from the
            ! single-precision equivalents of the input/output variables.
#ifdef SINGLE_PRECISION
            xin_p = xin
#endif

            if (nprocs == 1) then
                ! Use blas to perform matrix-vector multiplication.
                do i = 1, ncol
                    ! y = H x,
                    ! where H is the Hamiltonian matrix, x is the input Lanczos vector
                    ! and y the output Lanczos vector.
#ifdef SINGLE_PRECISION
                    if (sparse_hamil) then
                        ! This could be improved by multiplying the sparse
                        ! hamiltonian matrix by the dense xin matrix, rather than
                        ! doing one vector at a time.
                        call csrpsymv(hamil%mat_sparse, xin_p(:,i), yout_p(:,i))
                    else
                        call ssymv('U', nrow, 1.0_p, hamil%rmat, nrow, xin_p(:,i), 1, 0.0_p, yout_p(:,i), 1)
                    end if
#else
                    if (sparse_hamil) then
                        call csrpsymv(hamil%mat_sparse, xin(:,i), yout(:,i))
                    else
                        call dsymv('U', nrow, 1.0_dp, hamil%rmat, nrow, xin(:,i), 1, 0.0_dp, yout(:,i), 1)
                    end if
#endif
                end do
            else
#ifdef PARALLEL
                if (sparse_hamil) write (6,'(1X, a81)') 'WARNING.  Sparsity not implemented in parallel.  &
                                                  &This is not going to end well...'
                ! Use pblas to perform matrix-vector multiplication.
                if (proc_blacs_info%nrows > 0 .and. proc_blacs_info%ncols > 0) then
                    ! Some of the matrix on this processor.
                    do i = 1, ncol
                        ! y = H x,
                        ! where H is the Hamiltonian matrix, x is the input Lanczos vector
                        ! and y the output Lanczos vector.
#ifdef SINGLE_PRECISION
                        call pssymv('U', ndets, 1.0_p, hamil%rmat, 1, 1,                   &
                                    proc_blacs_info%desc_m, xin_p(:,i), 1, 1,         &
                                    proc_blacs_info%desc_v, 1, 0.0_p, yout_p(:,i), 1, &
                                    1, proc_blacs_info%desc_v, 1)
#else
                        call pdsymv('U', ndets, 1.0_dp, hamil%rmat, 1, 1,                 &
                                    proc_blacs_info%desc_m, xin(:,i), 1, 1,          &
                                    proc_blacs_info%desc_v, 1, 0.0_dp, yout(:,i), 1, &
                                    1, proc_blacs_info%desc_v, 1)
#endif
                    end do
                else
#ifdef SINGLE_PRECISION
                    yout_p = 0.0_p
#else
                    yout = 0.0_dp
#endif
                end if
#endif
            end if

#ifdef SINGLE_PRECISION
            yout = yout_p
#endif

        end subroutine HPsi

        subroutine HPsi_direct(nrow, ncol, xin, ldx, yout, ldy)

            ! Matrix-vector multiplication procedure for use with trlan.
            ! The Hamiltonian matrix elements are calculated directly ('on-the-fly')
            ! rather than looked up in the precomputed array hamil.
            ! In:
            !    nrow: the number of rows on this processor if the problem is distributed
            !        using MPI, otherwise the number of total rows in a Lanczos vector.
            !    ncol: the number of vectors (columns in xin and yout) to be multiplied.
            !    xin: the array to store the input vectors to be multiplied.
            !    ldx: the leading dimension of the array xin when it is declared as
            !       two-dimensional array.
            !    ldy: the leading dimension of the array yout when it is declared as
            !       two-dimensional array.
            ! Out:
            !    yout: the array to store results of the multiplication.
            
            ! NOTE: sys is used from the scope of lanczos_diagonalisation.

            use hamiltonian, only: get_hmatel

            integer, intent(in) :: nrow, ncol, ldx, ldy
            real(dp), intent(in) :: xin(ldx,ncol)
            real(dp), intent(out) :: yout(ldy,ncol)
            integer :: i, j, k
            real(p) :: tmp, hmatel

            yout = 0.0_dp
            do k = 1, ncol
                ! y = H x,
                ! where H is the Hamiltonian matrix, x is the input Lanczos vector
                ! and y the output Lanczos vector.
                ! Borrowing from ideas in dsymv, we can perform this only using one
                ! triangle of the Hamiltonian matrix.
                do j = 1, nrow ! Identical to ndets in serial.
                    tmp = 0.0_p
                    do i = 1, j-1
                        hmatel = get_hmatel(sys, dets(:,i), dets(:,j))
                        yout(i,k) = yout(i,k) + hmatel*xin(j, k)
                        tmp = tmp + hmatel*xin(i,k)
                    end do
                    yout(j,k) = get_hmatel(sys, dets(:,j), dets(:,j))*xin(j,k) + tmp
                end do
            end do

        end subroutine HPsi_direct

    end subroutine lanczos_diagonalisation

end module fci_lanczos

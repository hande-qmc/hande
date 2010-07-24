module lanczos

use const

implicit none

! Number of Lanczos eigenpairs to find.
integer :: nlanczos_eigv = 5

! Size of Lanczos basis.
integer :: lanczos_basis_length = 40

logical :: direct_lanczos = .false.

! Procedure to multiply the hamiltonian matrix by a Lanczos vector.
private :: HPsi, HPsi_direct

contains

    subroutine lanczos_diagonalisation(nfound, eigv)

        ! Perform a Lanczos diagonalisation of the current (spin) block of the
        ! Hamiltonian matrix.
        ! Out:
        !    nfound: number of solutions found from this block.  This is
        !        min(number of determinants with current spin, nlanczos_eigv).
        !    eigv(:nfound): Lanczos eigenvalues of the current block of the
        !        Hamiltonian matrix.

        use trl_info
        use trl_interface
        use checking, only: check_allocate, check_deallocate
        use errors, only: stop_all
        use parallel, only: parent, nprocs, get_blacs_info

        use calc
        use determinants, only: ndets
        use operators

        integer, intent(out) :: nfound
        real(dp), intent(out) :: eigv(ndets)
        
        integer, parameter :: lohi = -1
        integer :: mev
        real(dp), allocatable :: eval(:) ! (mev)
        real(dp), allocatable :: evec(:,:) ! (ndets, mev)
        type(trl_info_t) :: info
        integer :: ierr, nrows
        integer :: i

        ! mev: number of eigenpairs that can be stored in eval and evec.
        ! twice the number of eigenvalues to be found is a reasonable default.
        mev = min(2*nlanczos_eigv, ndets)
       
        if (parent) then
            if (direct_lanczos) then
                write (6,'(/,1X,a44,/)') 'Performing direct Lanczos diagonalisation...'
            else
                write (6,'(/,1X,a37,/)') 'Performing Lanczos diagonalisation...'
            end if
        end if
       
        if (distribute /= distribute_off .and. distribute /= distribute_cols) then
            call stop_all('exact_diagonalisation','Incorrect distribution mode used.')
        end if

        if (.not.doing_calc(exact_diag) .and. direct_lanczos) then
            ! Didn't execute generate_hamil so need to set up
            ! the processor grid.
            proc_blacs_info = get_blacs_info(ndets, (/1, nprocs/))
        end if

        if (direct_lanczos .and. nprocs > 1) then
            call stop_all('lanczos_diagonalisation','Direct Lanczos not implemented in parallel.')
        end if

        nrows = proc_blacs_info%nrows

        ! Initialise trlan.
        ! info: type(trl_info_t).  Used by trl to store calculation info.
        ! nrows: number of rows of matrix on processor.
        ! lanczos_basis_length: maximum Lanczos basis size.
        ! lohi: -1 means calculate the smallest eigenvalues first (1 to calculate
        !       the largest).
        ! nlanczos_eigv: number of eigenvalues to compute.
        call trl_init_info(info, nrows, lanczos_basis_length, lohi, min(nlanczos_eigv,ndets))
       
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
        if (direct_lanczos) then
            call trlan(HPsi_direct, info, nrows, mev, eval, evec, nrows)
        else
            call trlan(HPsi, info, nrows, mev, eval, evec, nrows)
        end if

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

        nfound = min(nlanczos_eigv,ndets)
        eigv(1:nfound) = eval(1:nfound)

        if (find_eigenvectors) then
            do i = 1, 1
                call analyse_wavefunction(evec(:,i))
            end do
        end if

        deallocate(eval, stat=ierr)
        call check_deallocate('eval',ierr)
        deallocate(evec, stat=ierr)
        call check_deallocate('evec',ierr)

    end subroutine lanczos_diagonalisation
       
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

        use parallel, only: nprocs

        use calc, only: hamil, proc_blacs_info
        use determinants, only: ndets
 
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
                call ssymv('U', nrow, 1.0_p, hamil, nrow, xin_p(:,i), 1, 0.0_p, yout_p(:,i), 1)
#else
                call dsymv('U', nrow, 1.0_dp, hamil, nrow, xin(:,i), 1, 0.0_dp, yout(:,i), 1)
#endif
            end do
        else
#ifdef PARALLEL
            ! Use pblas to perform matrix-vector multiplication.
            if (proc_blacs_info%nrows > 0 .and. proc_blacs_info%ncols > 0) then
                ! Some of the matrix on this processor.
                do i = 1, ncol
                    ! y = H x,
                    ! where H is the Hamiltonian matrix, x is the input Lanczos vector
                    ! and y the output Lanczos vector.
#ifdef SINGLE_PRECISION
                    call pssymv('U', ndets, 1.0_p, hamil, 1, 1,                   &
                                proc_blacs_info%desc_m, xin_p(:,i), 1, 1,         &
                                proc_blacs_info%desc_v, 1, 0.0_p, yout_p(:,i), 1, &
                                1, proc_blacs_info%desc_v, 1)
#else
                    call pdsymv('U', ndets, 1.0_dp, hamil, 1, 1,                 &
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

        use determinants, only: ndets
        use hamiltonian, only: get_hmatel
 
        integer, intent(in) :: nrow, ncol, ldx, ldy
        real(dp), intent(in) :: xin(ldx,ncol)
        real(dp), intent(out) :: yout(ldy,ncol)
        integer :: i, j, k
        real(dp) :: tmp, hmatel
 
        yout = 0.0_dp
        do k = 1, ncol
            ! y = H x,
            ! where H is the Hamiltonian matrix, x is the input Lanczos vector
            ! and y the output Lanczos vector.
            ! Borrowing from ideas in dsymv, we can perform this only using one
            ! triangle of the Hamiltonian matrix.
            do j = 1, nrow ! Identical to ndets in serial.
                tmp = 0.0_dp
                do i = 1, j-1
                    hmatel = get_hmatel(i,j) 
                    yout(i,k) = yout(i,k) + hmatel*xin(j, k)
                    tmp = tmp + hmatel*xin(i,k)
                end do
                yout(j,k) = get_hmatel(j,j)*xin(j,k) + tmp
            end do
        end do

    end subroutine HPsi_direct
    
end module lanczos

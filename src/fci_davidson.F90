module fci_davidson

! Complete diagonalisation of Hamiltonian matrix using lapack/scalapack
! libraries.

use const

implicit none

contains

   subroutine do_fci_davidson(sys, fci_in, ref_in)

        ! Perform an FCI calculation via the approximate iterative Davidson's method.

        ! In:
        !    sys: system of interest.
        !    fci_in: fci input options.
        !    ref_in: reference determinant defining (if relevant) a
        !        truncated Hilbert space.

        use dmqmc_procedures, only: setup_rdm_arrays
        use fci_utils, only: fci_in_t, init_fci, generate_hamil, write_hamil, hamil_t
        use hamiltonian, only: get_hmatel
        use reference_determinant, only: reference_t
        use system, only: sys_t, copy_sys_spin_info, heisenberg

        use checking, only: check_allocate
        use errors, only: warning
        use parallel, only: parent, nprocs, blacs_info, get_blacs_info
        use check_input, only: check_fci_opts
        use hamiltonian_data

        type(sys_t), intent(inout) :: sys
        type(fci_in_t), intent(inout) :: fci_in
        type(reference_t), intent(in) :: ref_in

        type(sys_t) :: sys_bak
        type(reference_t) :: ref
        integer(i0), allocatable :: dets(:,:)
        real(p), allocatable :: eigv(:), rdm_eigv(:), rdm(:,:)
        integer :: ndets, ierr, i, rdm_size
        type(blacs_info) :: proc_blacs_info
        type(hamil_t) :: hamil
        type(hmatel_t) :: hmatel
        integer :: iunit

        iunit = 6

        if (parent) call check_fci_opts(sys, fci_in)

        call copy_sys_spin_info(sys, sys_bak)
        ref = ref_in

        call init_fci(sys, fci_in, ref, ndets, dets)

        allocate(eigv(fci_in%ndavidson_eigv), stat=ierr)
        call check_allocate('eigv',ndets,ierr)

        if (ndets == 1) then
            ! The trivial case seems to trip things up...
            hmatel = get_hmatel(sys, dets(:,1), dets(:,1))
            eigv(1) = hmatel%r
        else
            if (nprocs == 1) then
                call generate_hamil(sys, ndets, dets, hamil, full_mat=.true.)
            else
                ! Use as square a processor grid as possible.
                proc_blacs_info = get_blacs_info(ndets, fci_in%block_size)
                call generate_hamil(sys, ndets, dets, hamil, proc_blacs_info=proc_blacs_info)
            end if
            if (fci_in%write_hamiltonian) then
                call write_hamil(fci_in%hamiltonian_file, ndets, proc_blacs_info, hamil)
            end if
            call davidson_diagonalisation(sys, fci_in, dets, proc_blacs_info, hamil, eigv)
        end if

        if (parent) then
            write (iunit,'(1X,"Davidson diagonalisation results")')
            write (iunit,'(1X,"^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^",/)')
            write (iunit,'(1X," State",5X,"Energy")')
            do i = 1, fci_in%ndavidson_eigv
                write (iunit,'(1X,i6,1X,f18.12)') i, eigv(i)
            end do
            write (iunit,'()')
        end if

        ! Return sys in an unaltered state.
        call copy_sys_spin_info(sys_bak, sys)

    end subroutine do_fci_davidson

    subroutine davidson_diagonalisation(sys, fci_in, dets, proc_blacs_info, hamil, eigv)

        ! Perform a Davidson diagonalisation of the current (spin) block of the
        ! Hamiltonian matrix.
        ! Note that this destroys the Hamiltonian matrix stored in hamil.
        ! Brief notes on Davidson diagonalisation:
        !   To diagonalise a sparse, (ndets x ndets) matrix A and get nEig lowest eigenpairs, we form ntrial guesses, 
        !   where ntrial is conventionally about double nEigs, and guess vectors are usually just ntrial lowest unit vectors. 
        !   We collect these vectors into a rectangular (ndets x ntrial) matrix V.
        !   Further parameters include:
        !       maxguess: max number of guess vectors V can hold, after it's exceeded we need to collapse the subspace (see below)
        !       maxiter: max number of iterations of the process
        !       nactive: number of guess vectors currently being held in the matrix V
        !   Algorithm:
        !       initialise V
        !       nactive = ntrial
        !       if iter < maxiter do       
        !       1. Orthonormalise V(:,1:nactive) (we use the QR decomposition)
        !       2. Compute T = V^T A V, the Rayleigh matrix / subspace Hamiltonian
        !       3. Diagonalise the subspace Hamiltonian: TC = Ct
        !       4. if nactive < maxguess - ntrial, then
        !           for i in {1..ntrial}
        !           4.1 r^{(i)} = V C^{(i)} where r is the Ritz vector
        !           4.2 w^{(i)} = A r^{(i)}, the residual vector
        !           4.3 Convergence test: if norm(w) < tol then this residual has converged
        !           4.4 q^{(i)} = w^{(i)} / (t_j - A(j,j)), with the denominator the preconditioner and q the new guess vector
        !           4.5 V(:,nactive+i) = q^{(i)}
        !           4.6 nactive += ntrial
        !       5. Else, collapse subspace
        !           5.1 V(:,:ntrial) = V(:,:maxguess) T(:maxguess,:ntrial)
        !           5.2 nactive = ntrial
        !       
        ! In:
        !    sys: system being studied.  Only used if the wavefunction is
        !        analysed.
        !    fci_in: fci input options.
        !    dets: list of determinants in the Hilbert space in the bit
        !        string representation.
        !    proc_blacs_info: BLACS description of distribution of the Hamiltonian.
        !        Used only if running on multiple processors.
        ! In/Out:
        !    hamil: hamil_t derived type containing Hamiltonian matrix of the system
        !        in the Hilbert space used in appropriate format.
        !        On output contains the eigenvectors, if requested, or is
        !        otherwise destroyed.
        ! Out:
        !    eigv: fci_in%ndavidson_eigv lowest eigenvalues of the current block of the
        !        Hamiltonian matrix.

        use checking, only: check_allocate, check_deallocate
        use linalg, only: syev_wrapper, heev_wrapper, gemm, gemv, qr_wrapper, psyev_wrapper, pheev_wrapper, pgemm, pqr_wrapper
        use linalg, only: qr_wrapper, pqr_wrapper
        use parallel, only: parent, nprocs, blacs_info
        use system, only: sys_t
        use errors, only: stop_all, warning

        use fci_utils, only: fci_in_t, hamil_t
        use operators

        type(sys_t), intent(in) :: sys
        type(fci_in_t), intent(in) :: fci_in
        integer(i0), intent(in) :: dets(:,:)
        type(blacs_info), intent(in) :: proc_blacs_info
        type(hamil_t), intent(inout) :: hamil
        real(p), intent(out) :: eigv(:)
        real(p), allocatable :: rwork(:), eigvec(:,:)
        integer :: info, ierr, i, j, nwfn, ndets
        integer :: iunit, maxguess, nactive
        real(p), allocatable :: V(:,:), theta(:), theta_old(:), tmp(:,:), tmpV(:), T(:,:), w(:), V_coll(:,:)
        logical, allocatable :: normconv(:)
        logical :: conv
        real(p) :: norm

        iunit = 6

        ndets = ubound(dets, dim=2)

        if (parent) then
            write (iunit,'(1X,A,/)') 'Performing Davidson diagonalisation...'
        end if

        if (nprocs > 1) then
            call stop_all('davidson_diagonalisation', 'ScaLAPACK not yet supported for Davidson diagonalisation!')
        end if

        if (hamil%comp) then
            call stop_all('davidson_diagonalisation', 'Complex Hamiltonian not yet supported for Davidson diagonalisation!')
        end if

        associate(nEig=>fci_in%ndavidson_eigv, ntrial=>fci_in%ndavidson_trialvec, maxsize=>fci_in%davidson_maxsize, A=>hamil%rmat, &
                  maxiter=>fci_in%davidson_maxiter, tol=>fci_in%davidson_tol)

        if (nprocs == 1) then

            allocate(theta_old(nEig), source=0.0_p, stat=ierr)
            call check_allocate('theta_old', nEig, ierr)
            allocate(normconv(ntrial))
            call check_allocate('normconv',ntrial, ierr)

            ! Integer division always rounds towards zero, i.e 8*50/8 = 48
            maxguess = ntrial*(maxsize/ntrial)
            allocate(theta(maxguess))
            allocate(V(ndets,maxguess+ntrial),source=0.0_p,stat=ierr)
            call check_allocate('V',ndets*(maxguess+ntrial),ierr)

            allocate(V_coll(ndets,ntrial),source=0.0_p,stat=ierr)
            call check_allocate('V',ndets*ntrial,ierr)

            allocate(T(maxguess,maxguess),source=0.0_p,stat=ierr)
            call check_allocate('T',maxguess**2,ierr)

            allocate(tmp,source=V,stat=ierr)
            call check_allocate('tmp',ndets*maxguess,ierr)

            allocate(tmpV(ndets),source=0.0_p,stat=ierr)
            call check_allocate('tmpV',ndets,ierr)

            allocate(w(ndets),source=0.0_p,stat=ierr)
            call check_allocate('w',ndets,ierr)

            ! Initial guesses are ntrial lowest unit vectors
            do i = 1, ntrial
                V(i,i) = 1.0_p
            end do

            nactive = ntrial
            conv = .false.
            
            do i = 1, maxiter
                ! Orthonormalise the current set of guess vectors
                call qr_wrapper(ndets, nactive, V, ndets, info)
                
                ! Form the subspace Hamiltonian / Rayleigh matrix
                ! T = V^T A V
                call gemm('N', 'N', ndets, nactive, ndets, 1.0_p, A, ndets, V, ndets, 0.0_p, tmp, ndets)
                call gemm('T', 'N', nactive, nactive, ndets, 1.0_p, V, ndets, tmp, ndets, 0.0_p, T, nactive)
                ! Diagonalise the subspace Hamiltonian
                ! T C = C t (T gets overwritten by C in syev, theta is the diagonal of t)
                call syev_wrapper('V', 'U', nactive, T, nactive, theta, info)

                ! We don't use this norm for convergence checking as after each subspace collapse the change in norm is 
                ! essentially zero, but we report it nonetheless as during each restart-block it is still 
                ! a useful measure of convergence.
                norm = sqrt(sum((theta(1:nEig)-theta_old)**2))
                write(iunit,'(1X, A, I3, A, I3, A, ES15.6)') 'Iteration ', i, ', basis size ', nactive, ', rmsE ', norm
                theta_old = theta(1:nEig)

                if (all(normconv)) then
                    write(iunit,'(1X, A, ES10.4, A)') 'Residue tolerance of ', tol,' reached, printing results...'
                    conv = .true.
                    exit
                end if

                if (nactive <= (maxguess-ntrial)) then
                    ! If the number of guess vectors can be grown by at least another lot of ntrial
                    do j = 1, ntrial
                        ! Residue vector w = (A-theta(j)*I) V T(:,j)
                        ! Storing a diagonal matrix as large as A is obviously a bad idea, so we use a tmp vector
                        ! Technically speaking, tmpV is the 'Ritz vector' and w is the residue vector
                        call gemv('N', ndets, nactive, 1.0_p, V, ndets, T(:,j), 1, 0.0_p, tmpV, 1)
                        call gemv('N', ndets, ndets, 1.0_p, A, ndets, tmpV, 1, 0.0_p, w, 1)
                        w = w - theta(j)*tmpV
                        if (sqrt(sum(w**2)) < tol) normconv(j) = .true.
                        ! Precondition the residue vector to form the correction vector,
                        ! if preconditioner = 1, we recover the Lanczos algorithm.
                        if (abs((theta(j)-A(j,j))) < 1e-6) then
                            w = w/(theta(j) - A(j,j) + 0.01)
                        else
                            w = w/(theta(j)-A(j,j))
                        end if
                        V(:,(nactive+j)) = w
                    end do
                    nactive = nactive + ntrial
                else
                    ! We need to collapse the subspace into ntrial best guesses and restart the iterations
                    ! V holds the approximate eigenvectors and T holds the CI coefficients,
                    ! so one call to gemm gives us the actual guess vectors.
                    write(iunit, '(1X, A)') 'Collapsing subspace...'
                    ! Matrix multiplication requires the matrices throughout the calculation, so a temporary array is needed
                    ! to avoid compilers deciding to overwrite the original V array during the multiplication.
                    call gemm('N', 'N', ndets, ntrial, maxguess, 1.0_p, V,&
                              ndets, T, maxguess, 0.0_p, V_coll, ndets)
                    V(:,:ntrial) = V_coll
                    nactive = ntrial
                end if

            end do

            if (.not. conv) then
                call warning('davidson_diagonalisation',&
                    'Davidson diagonalisation did not converge, choose a larger maxiter or tolerance.')
            end if

            eigv = theta(1:nEig)
            deallocate(V, V_coll, T, tmpV, w, tmp, theta)
        end if
        end associate

    end subroutine davidson_diagonalisation

end module fci_davidson

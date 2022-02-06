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

        if (parent) call check_fci_opts(sys, fci_in, .true.)

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
                call generate_hamil(sys, ndets, dets, hamil)
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
        integer :: iunit, maxEig, nactive
        real(p), allocatable :: V(:,:), theta(:), theta_old(:), tmp(:,:), tmpV(:), T(:,:), w(:)
        logical, allocatable :: normconv(:)
        logical :: conv
        real(p) :: norm
        integer(8) :: t0, t1, count_rate, count_max, t_qr, t_dgemm1, t_dgemm2, t_diag,t_dgemv1,t_dgemv2, t_collapse, t_start, t_end


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

        associate(nEig=>fci_in%ndavidson_eigv, nTrial=>fci_in%ndavidson_trialvec, maxSize=>fci_in%davidson_maxsize, A=>hamil%rmat, &
                  maxiter=>fci_in%davidson_maxiter, tol=>fci_in%davidson_tol)

        if (nprocs == 1) then

            allocate(theta_old(nEig), source=0.0_p)
            allocate(normconv(nTrial))

            ! Integer division always rounds towards zero, i.e 8*50/8 = 48
            maxEig = nTrial*(maxSize/nTrial)
            allocate(theta(maxEig))
            allocate(V(ndets,maxEig+nTrial),source=0.0_p,stat=ierr)
            call check_allocate('V',ndets*(maxEig+nTrial),ierr)

            allocate(T(maxEig,maxEig),source=0.0_p,stat=ierr)
            call check_allocate('T',maxEig**2,ierr)

            allocate(tmp,source=V,stat=ierr)
            call check_allocate('tmp',ndets*maxEig,ierr)

            allocate(tmpV(ndets),source=0.0_p,stat=ierr)
            call check_allocate('tmpV',ndets,ierr)

            allocate(w(ndets),source=0.0_p,stat=ierr)
            call check_allocate('w',ndets,ierr)

            ! For now Davidson requires the full matrix, [todo] - make it compatible with triangulars
            ! inflate A from upper triangular to full matrix
            do i = 1, ndets
                do j = i, ndets
                    A(j,i) = A(i,j)
                end do
            end do
            
            ! Initial guesses are nTrial lowest unit vectors
            do i = 1, nTrial
                V(i,i) = 1.0_p
            end do
            nactive = nTrial
            conv = .false.
            
            t_qr=0; t_dgemm1=0; t_dgemm2=0; t_diag=0; t_dgemv1=0; t_dgemv2=0; t_collapse=0
            call system_clock(t_start, count_rate, count_max)

            do i = 1, maxiter
                call system_clock(t0)
                ! Orthonormalise the current set of guess vectors
                call qr_wrapper(ndets, nactive, V, ndets, info)
                
                call system_clock(t1)
                t_qr = t_qr + t1 - t0
                t0 = t1

                
                ! Form the subspace Hamiltonian / Rayleigh matrix
                ! T = V^T A V
                call gemm('N', 'N', ndets, nactive, ndets, 1.0_p, A, ndets, V, ndets, 0.0_p, tmp, ndets)
                call system_clock(t1)
                t_dgemm1 = t_dgemm1 + t1 - t0
                t0 = t1
                call gemm('T', 'N', nactive, nactive, ndets, 1.0_p, V, ndets, tmp, ndets, 0.0_p, T, nactive)
                call system_clock(t1)
                t_dgemm2 = t_dgemm2 + t1 - t0
                t0 = t1

                ! Diagonalise the subspace Hamiltonian
                ! T C = C t (T gets overwritten by C in syev, theta is the diagonal of t)
                call syev_wrapper('V', 'U', nactive, T, nactive, theta, info)
                call system_clock(t1)
                t_diag = t_diag + t1 - t0
                t0 = t1

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

                if (nactive <= (maxEig-nTrial)) then
                    ! If the number of guess vectors can be grown by at least another lot of nTrial
                    do j = 1, nTrial
                        ! Residue vector w = (A-theta(j)*I) V T(:,j)
                        ! Storing a diagonal matrix as large as A is obviously a bad idea, so we use a tmp vector
                        ! Technically speaking, tmpV is the 'Ritz vector' and w is the residue vector
                        call gemv('N', ndets, nactive, 1.0_p, V, ndets, T(:,j), 1, 0.0_p, tmpV, 1)
                        call system_clock(t1)
                        t_dgemv1 = t_dgemv1 + t1 - t0
                        t0 = t1
                        call gemv('N', ndets, ndets, 1.0_p, A, ndets, tmpV, 1, 0.0_p, w, 1)
                        call system_clock(t1)
                        t_dgemv2 = t_dgemv2 + t1 - t0
                        t0 = t1
                        w = w - theta(j)*tmpV
                        if (sqrt(sum(w**2)) < tol) normconv(j) = .true.
                        ! Precondition the residue vector to form the correction vector,
                        ! if preconditioner = 1, we recover the Lanczos algorithm.
                        w = w/(theta(j)-A(j,j))
                        V(:,(nactive+j)) = w
                    end do
                    nactive = nactive + nTrial
                else
                    ! We need to collapse the subspace into nTrial best guesses and restart the iterations
                    ! V holds the approximate eigenvectors and T holds the CI coefficients,
                    ! so one call to gemm gives us the actual guess vectors.
                    write(iunit, '(1X, A)') 'Collapsing subspace...'
                    call gemm('N', 'N', ndets, nTrial, maxEig, 1.0_p, V,&
                              ndets, T, maxEig, 0.0_p, V, ndets)
                    call system_clock(t1)
                    t_collapse = t_collapse + t1 - t0
                    t0 = t1
                    nactive = nTrial
                end if

            end do

            if (conv .eqv. .false.) then
                call warning('davidson_diagonalisation',&
                    'Davidson diagonalisation did not converge, choose a larger maxiter or tolerance.')
            end if

            eigv = theta(1:nEig)
        else
            continue
        end if
        call system_clock(t_end)
        write(iunit, '(1X, A, F16.8)') 'Total runtime: ', real((t_end-t_start),kind=p)/count_rate
        write(iunit, '(1X, A, F16.8)') 'QR: ', real(t_qr,kind=p)/count_rate
        write(iunit, '(1X, A, F16.8)') 'dgemm1: ', real(t_dgemm1,kind=p)/count_rate
        write(iunit, '(1X, A, F16.8)') 'dgemm2: ', real(t_dgemm2,kind=p)/count_rate
        write(iunit, '(1X, A, F16.8)') 'diag: ', real(t_diag,kind=p)/count_rate
        write(iunit, '(1X, A, F16.8)') 'dgemv1: ', real(t_dgemv1,kind=p)/count_rate
        write(iunit, '(1X, A, F16.8)') 'dgemv2: ', real(t_dgemv2,kind=p)/count_rate
        write(iunit, '(1X, A, F16.8)') 'collapse: ', real(t_collapse,kind=p)/count_rate
        end associate

    end subroutine davidson_diagonalisation

end module fci_davidson

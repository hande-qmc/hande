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
        use linalg, only: syev_wrapper, heev_wrapper, gemm, qr_wrapper, psyev_wrapper, pheev_wrapper, pgemm, pqr_wrapper
        use linalg, only: qr_wrapper, pqr_wrapper
        use parallel, only: parent, nprocs, blacs_info
        use system, only: sys_t
        use errors, only: stop_all

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
        character(1) :: job
        integer :: iunit, maxEig, nactive
        real(p), allocatable :: V(:,:), theta(:), theta_old(:), tmp(:,:), tmpV(:,:), eye(:,:), T(:,:), w(:,:), tmpdiag(:,:)
        logical, allocatable :: normconv(:)
        logical :: conv
        real(p) :: norm

        iunit = 6

        ndets = ubound(dets, dim=2)

        if (parent) then
            write (iunit,'(1X,A,/)') 'Performing Davidson diagonalisation...'
        end if

        ! Davidson needs the eigenvectors to build the residuals
        job = 'V'

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

            allocate(tmpV(ndets,1),source=0.0_p,stat=ierr)
            call check_allocate('tmpV',ndets,ierr)

            allocate(w(ndets,1),source=0.0_p,stat=ierr)
            call check_allocate('w',ndets,ierr)
            
            ! Initial guesses are nTrial lowest unit vectors
            do i = 1, nTrial
                V(i,i) = 1.0_p
            end do
            nactive = nTrial
            conv = .false.

            conv_while: do while (conv .eqv. .false.)
            maxiter_do: do i = 1, maxiter
                call qr_wrapper(ndets, nactive, V, ndets, info)

                ! Likewise, tmp is designed to hold the largest tmp matrix but only the relevant slice will be written to
                call gemm('N', 'N', ndets, nactive, ndets, 1.0_p, A, ndets, V, ndets, 0.0_p, tmp, ndets)
                call gemm('T', 'N', nactive, nactive, ndets, 1.0_p, V, ndets, tmp, ndets, 0.0_p, T, nactive)

                call syev_wrapper(job, 'U', nactive, T, nactive, theta, info)

                norm = sqrt(sum((theta(1:nEig)-theta_old)**2))
                write(iunit,'(1X, A, I0, A, I0, A, E15.6)') 'Iteration ', i, ', basis size ', nactive, ', rmsE ', norm

                if (all(normconv)) then
                    write(iunit,'(1X, A, E10.4, A)') 'Residue tolerance of ', tol,' reached, printing results...'
                    conv = .true.
                    exit
                end if

                if (nactive <= (maxEig-nTrial)) then
                    do j = 1, nTrial
                        call gemm('N', 'N', ndets, 1, nactive, 1.0_p, V, ndets, T, nactive, 0.0_p, tmpV, ndets)
                        call gemm('N', 'N', ndets, 1, ndets, 1.0_p, A, ndets, tmpV, ndets, 0.0_p, w, ndets)
                        w(:,1) = (w(:,1)-theta(j)*tmpV(:,1))/(theta(j)-A(j,j))
                        if (sqrt(sum(w(:,1)**2)) < tol) normconv(j) = .true.
                        V(:,(nactive+j)) = w(:,1)
                    end do
                    nactive = nactive + nTrial
                else
                    write(iunit, '(1X, A)') 'Collapsing subspace...'
                    call gemm('N', 'N', ndets, nTrial, maxEig, 1.0_p, V, ndets, T, maxEig, 0.0_p, V, ndets)
                    nactive = nTrial
                end if

                theta_old = theta(1:nEig)
            end do maxiter_do
            end do conv_while

            eigv = theta(1:nEig)
        else
            continue
        end if
        end associate

    end subroutine davidson_diagonalisation

end module fci_davidson

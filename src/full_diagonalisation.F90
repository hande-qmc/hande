module full_diagonalisation

! Complete diagonalisation of Hamiltonian matrix using lapack/scalapack
! libraries.

use const

implicit none

contains

    subroutine exact_diagonalisation(sys, dets, eigv)

        ! Perform an exact diagonalisation of the current (spin) block of the
        ! Hamiltonian matrix.
        ! Note that this destroys the Hamiltonian matrix stored in hamil.
        ! In:
        !    sys: system being studied.  Only used if the wavefunction is
        !        analysed.
        !    dets: list of determinants in the Hilbert space in the bit
        !        string representation.
        ! Out:
        !    eigv(ndets): Lanczos eigenvalues of the current block of the
        !        Hamiltonian matrix.

        use checking, only: check_allocate, check_deallocate
        use errors, only: stop_all
        use parallel, only: parent, nprocs
        use system, only: sys_t

        use calc

        use operators

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: dets(:,:)
        real(p), intent(out) :: eigv(:)
        real(p), allocatable :: work(:), eigvec(:,:)
        integer :: info, ierr, lwork, i, nwfn, ndets
        character(1) :: job

        ndets = ubound(dets, dim=2)

        if (parent) then
            write (6,'(/,1X,a35,/)') 'Performing exact diagonalisation...'
        end if

        if (analyse_fci_wfn /= 0 .or. print_fci_wfn /= 0 .or. doing_exact_rdm_eigv) then
            job = 'V'
        else
            job = 'N'
        end if

        if (distribute /= distribute_off .and. distribute /= distribute_blocks) then
            call stop_all('exact_diagonalisation','Incorrect distribution mode used.')
        end if

        if (nprocs > 1) then
            allocate(eigvec(proc_blacs_info%nrows, proc_blacs_info%ncols), stat=ierr)
            call check_allocate('eigvec',proc_blacs_info%nrows*proc_blacs_info%ncols,ierr)
        end if

        ! Find the optimal size of the workspace.
        allocate(work(1), stat=ierr)
        call check_allocate('work',1,ierr)
        if (nprocs == 1) then
#ifdef SINGLE_PRECISION
            call ssyev(job, 'U', ndets, hamil, ndets, eigv, work, -1, info)
#else
            call dsyev(job, 'U', ndets, hamil, ndets, eigv, work, -1, info)
#endif
        else
#ifdef PARALLEL
#ifdef SINGLE_PRECISION
            call pssyev(job, 'U', ndets, hamil, 1, 1,               &
                        proc_blacs_info%desc_m, eigv, eigvec, 1, 1, &
                        proc_blacs_info%desc_m, work, -1, info)
#else
            call pdsyev(job, 'U', ndets, hamil, 1, 1,               &
                        proc_blacs_info%desc_m, eigv, eigvec, 1, 1, &
                        proc_blacs_info%desc_m, work, -1, info)
#endif
#endif
        end if

        lwork = nint(work(1))
        deallocate(work)
        call check_deallocate('work',ierr)

        ! Now perform the diagonalisation.
        allocate(work(lwork), stat=ierr)
        call check_allocate('work',lwork,ierr)

        if (nprocs == 1) then
            ! Use lapack.
#ifdef SINGLE_PRECISION
            call ssyev(job, 'U', ndets, hamil, ndets, eigv, work, lwork, info)
#else
            call dsyev(job, 'U', ndets, hamil, ndets, eigv, work, lwork, info)
#endif
        else
#ifdef PARALLEL
            ! Use scalapack to do the diagonalisation in parallel.
#ifdef SINGLE_PRECISION
            call pssyev(job, 'U', ndets, hamil, 1, 1,               &
                        proc_blacs_info%desc_m, eigv, eigvec, 1, 1, &
                        proc_blacs_info%desc_m, work, lwork, info)
#else
            call pdsyev(job, 'U', ndets, hamil, 1, 1,               &
                        proc_blacs_info%desc_m, eigv, eigvec, 1, 1, &
                        proc_blacs_info%desc_m, work, lwork, info)
#endif
#endif
        end if

        deallocate(work, stat=ierr)
        call check_deallocate('work',ierr)

        nwfn = analyse_fci_wfn
        if (nwfn < 0) nwfn = ndets
        do i = 1, nwfn
            if (nprocs == 1) then
                call analyse_wavefunction(sys, hamil(:,i), dets)
            else
                call analyse_wavefunction(sys, eigvec(:,i), dets)
            end if
        end do
        nwfn = print_fci_wfn
        if (nwfn < 0) nwfn = ndets
        do i = 1, nwfn
            if (nprocs == 1) then
                call print_wavefunction(print_fci_wfn_file, hamil(:,i), dets)
            else
                call print_wavefunction(print_fci_wfn_file, eigvec(:,i), dets)
            end if
        end do

        if (nprocs > 1) then
            deallocate(eigvec, stat=ierr)
            call check_deallocate('eigvec',ierr)
        end if

    end subroutine exact_diagonalisation

end module full_diagonalisation

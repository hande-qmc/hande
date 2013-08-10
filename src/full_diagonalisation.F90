module full_diagonalisation

! Complete diagonalisation of Hamiltonian matrix using lapack/scalapack
! libraries.

use const

implicit none

contains

    subroutine exact_diagonalisation(eigv)

        ! Perform an exact diagonalisation of the current (spin) block of the
        ! Hamiltonian matrix.
        ! Note that this destroys the Hamiltonian matrix stored in hamil.
        ! Out:
        !    eigv(ndets): Lanczos eigenvalues of the current block of the
        !        Hamiltonian matrix.

        use checking, only: check_allocate, check_deallocate
        use errors, only: stop_all
        use fciqmc_data, only: doing_exact_rdm_eigv
        use parallel, only: parent, nprocs

        use calc
        use determinant_enumeration, only: ndets

        use operators

        real(p), intent(out) :: eigv(ndets)
        real(p), allocatable :: work(:), eigvec(:,:)
        integer :: info, ierr, lwork
        character(1) :: job

        if (parent) then
            write (6,'(/,1X,a35,/)') 'Performing exact diagonalisation...'
        end if

        if (analyse_ground_state .or. print_ground_state .or. doing_exact_rdm_eigv) then
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
            if (proc_blacs_info%nrows > 0 .and. proc_blacs_info%ncols > 0) then
                ! Part of matrix on this processor.
#ifdef SINGLE_PRECISION
                call pssyev(job, 'U', ndets, hamil, 1, 1,               &
                            proc_blacs_info%desc_m, eigv, eigvec, 1, 1, &
                            proc_blacs_info%desc_m, work, -1, info)
#else
                call pdsyev(job, 'U', ndets, hamil, 1, 1,               &
                            proc_blacs_info%desc_m, eigv, eigvec, 1, 1, &
                            proc_blacs_info%desc_m, work, -1, info)
#endif
            end if
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
            if (proc_blacs_info%nrows > 0 .and. proc_blacs_info%ncols > 0) then
                ! Part of matrix on this processor.
#ifdef SINGLE_PRECISION
                call pssyev(job, 'U', ndets, hamil, 1, 1,               &
                            proc_blacs_info%desc_m, eigv, eigvec, 1, 1, &
                            proc_blacs_info%desc_m, work, lwork, info)
#else
                call pdsyev(job, 'U', ndets, hamil, 1, 1,               &
                            proc_blacs_info%desc_m, eigv, eigvec, 1, 1, &
                            proc_blacs_info%desc_m, work, lwork, info)
#endif
            end if
#endif
        end if

        deallocate(work, stat=ierr)
        call check_deallocate('work',ierr)

        if (analyse_ground_state) then
            if (nprocs == 1) then
                call analyse_wavefunction(hamil(:,1))
            else
                call analyse_wavefunction(eigvec(:,1))
            end if
        else if (print_ground_state) then
            if (nprocs == 1) then
                call print_wavefunction('GROUND_STATE_WFN', hamil(:,1))
            else
                call print_wavefunction('GROUND_STATE_WFN', eigvec(:,1))
            end if
        end if

        if (nprocs > 1) then
            deallocate(eigvec, stat=ierr)
            call check_deallocate('eigvec',ierr)
        end if

    end subroutine exact_diagonalisation

end module full_diagonalisation

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

        use errors, only: stop_all
        use parallel, only: parent, nprocs

        use calc
        use determinants, only: ndets

        real(dp), intent(out) :: eigv(ndets)
        real(dp), allocatable :: work(:), eigvec(:,:)
        integer :: info, ierr, lwork
        integer :: i
        character(1) :: job

        if (parent) then
            write (6,'(/,1X,a35,/)') 'Performing exact diagonalisation...'
        end if

        if (find_eigenvectors) then
            job = 'V'
        else
            job = 'N'
        end if

        if (distribute /= distribute_off .and. distribute /= distribute_blocks) then
            call stop_all('exact_diagonalisation','Incorrect distribution mode used.')
        end if
        
        allocate(eigvec(proc_blacs_info%nrows, proc_blacs_info%ncols), stat=ierr)

        ! Find the optimal size of the workspace.
        allocate(work(1), stat=ierr)
        if (nprocs == 1) then
            call dsyev(job, 'U', ndets, hamil, ndets, eigv, work, -1, info)
        else
#ifdef _PARALLEL
            if (proc_blacs_info%nrows > 0 .and. proc_blacs_info%ncols > 0) then
                ! Part of matrix on this processor.
                call pdsyev(job, 'U', ndets, hamil, 1, 1,               &
                            proc_blacs_info%desc_m, eigv, eigvec, 1, 1, &
                            proc_blacs_info%desc_m, work, -1, info)
            end if
#endif
        end if

        lwork = work(1)
        deallocate(work)

        ! Now perform the diagonalisation.
        allocate(work(lwork), stat=ierr)

        if (nprocs == 1) then
            ! Use lapack.
            call dsyev(job, 'U', ndets, hamil, ndets, eigv, work, lwork, info)
        else
#ifdef _PARALLEL
            ! Use scalapack to do the diagonalisation in parallel.
            if (proc_blacs_info%nrows > 0 .and. proc_blacs_info%ncols > 0) then
                ! Part of matrix on this processor.
                call pdsyev(job, 'U', ndets, hamil, 1, 1,               &
                            proc_blacs_info%desc_m, eigv, eigvec, 1, 1, &
                            proc_blacs_info%desc_m, work, lwork, info)
            end if
#endif
        end if

        deallocate(work, stat=ierr)
        deallocate(eigvec, stat=ierr)

    end subroutine exact_diagonalisation

end module full_diagonalisation

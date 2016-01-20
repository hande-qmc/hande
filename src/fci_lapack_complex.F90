module fci_lapack_complex

! Complete diagonalisation of complex Hamiltonian matrix using lapack/scalapack
! libraries.
! Use separate algorithm so 

use const

implicit none

contains

    subroutine do_fci_lapack_complex(sys, fci_in, ref_in)

        ! Perform an FCI calculation using LAPACK/ScaLAPACK diagonalisation.

        ! In:
        !    sys: system of interest.
        !    fci_in: fci input options.
        !    ref_in: reference determinant defining (if relevant) a
        !        truncated Hilbert space.

        use dmqmc_procedures, only: setup_rdm_arrays
        use fci_utils, only: fci_in_t, init_fci, generate_hamil, write_hamil
        use hamiltonian, only: get_hmatel
        use qmc_data, only: reference_t
        use reference_determinant, only: copy_reference_t
        use system, only: sys_t, copy_sys_spin_info, heisenberg
        use checking, only: check_allocate
        use errors, only: warning
        use parallel, only: parent, nprocs, blacs_info, get_blacs_info
        use check_input, only: check_fci_opts

        type(sys_t), intent(inout) :: sys
        type(fci_in_t), intent(inout) :: fci_in
        type(reference_t), intent(in) :: ref_in

        type(sys_t) :: sys_bak
        type(reference_t) :: ref
        integer(i0), allocatable :: dets(:,:)
        real(p), allocatable :: eigv(:)
        complex(p), allocatable :: hamil(:,:)
        integer :: ndets, ierr, i
        type(blacs_info) :: proc_blacs_info

        if (parent) call check_fci_opts(sys, fci_in, .false.)

        call copy_sys_spin_info(sys, sys_bak)
        call copy_reference_t(ref_in, ref)

        call init_fci(sys, fci_in, ref, ndets, dets)

        allocate(eigv(ndets), stat=ierr)
        call check_allocate('eigv',ndets,ierr)

        if (ndets == 1) then
            ! The trivial case seems to trip things up...
            eigv(1) = get_hmatel(sys, dets(:,1), dets(:,1))
        else
            if (nprocs == 1) then
                call generate_hamil(sys, ndets, dets, hamil_comp=hamil)
            else
                ! Use as square a processor grid as possible.
                proc_blacs_info = get_blacs_info(ndets, fci_in%block_size)
                call generate_hamil(sys, ndets, dets, hamil_comp=hamil, proc_blacs_info=proc_blacs_info)
            end if
            if (fci_in%write_hamiltonian) then
                call write_hamil(fci_in%hamiltonian_file, ndets, proc_blacs_info, hamil_comp=hamil)
            end if
            call lapack_diagonalisation(sys, fci_in, dets, proc_blacs_info, hamil, eigv)
        end if

        if (parent) then
            write (6,'(1X,"LAPACK diagonalisation results")')
            write (6,'(1X,"^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^",/)')
            write (6,'(1X," State",5X,"Energy")')
            do i = 1, ndets
                write (6,'(1X,i6,1X,f18.12)') i, eigv(i)
            end do
        end if

        ! If requested, calculate and print eigenvalues for an RDM.
        if (allocated(fci_in%subsys_info)) then
            if (parent) call warning('diagonalise','RDM calculations aren"t implemented in complex. Skipping.',3)
        end if

        ! Return sys in an unaltered state.
        call copy_sys_spin_info(sys_bak, sys)

    end subroutine do_fci_lapack_complex

    subroutine lapack_diagonalisation(sys, fci_in, dets, proc_blacs_info, hamil, eigv)

        ! Perform an exact diagonalisation of the current (spin) block of the
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
        !    hamil: Hamiltonian matrix of the system in the Hilbert space used.
        !        On output contains the eigenvectors, if requested, or is
        !        otherwise destroyed.
        ! Out:
        !    eigv(ndets): Lanczos eigenvalues of the current block of the
        !        Hamiltonian matrix.

        use checking, only: check_allocate, check_deallocate
        use parallel, only: parent, nprocs, blacs_info
        use system, only: sys_t

        use fci_utils, only: fci_in_t
        use operators

        type(sys_t), intent(in) :: sys
        type(fci_in_t), intent(in) :: fci_in
        integer(i0), intent(in) :: dets(:,:)
        type(blacs_info), intent(in) :: proc_blacs_info
        complex(p), intent(inout) :: hamil(:,:)
        real(p), intent(out) :: eigv(:)
        complex(p), allocatable :: work(:)
        real(p), allocatable :: rwork(:)
        complex(p), allocatable :: eigvec(:,:)
        integer :: info, ierr, lwork, i, nwfn, ndets
        character(1) :: job

        ndets = ubound(dets, dim=2)

        if (parent) then
            write (6,'(/,1X,a35,/)') 'Performing exact diagonalisation...'

        end if

        if (fci_in%analyse_fci_wfn /= 0 .or. fci_in%print_fci_wfn /= 0 &
                    .or. allocated(fci_in%subsys_info)) then
            job = 'V'
        else
            job = 'N'
        end if

        if (nprocs > 1) then
            allocate(eigvec(proc_blacs_info%nrows, proc_blacs_info%ncols), stat=ierr)
            call check_allocate('eigvec',proc_blacs_info%nrows*proc_blacs_info%ncols,ierr)
        end if

        
        ! Find the optimal size of the workspace.
        allocate(work(1), stat=ierr)
        call check_allocate('work',1,ierr)
        allocate(rwork(max(1,3*ndets-2)), stat = ierr)
        call check_allocate('rwork',max(1, 3*ndets - 2), ierr)
        if (nprocs == 1) then
#ifdef SINGLE_PRECISION
            call cheev(job, 'U', ndets, hamil, ndets, eigv, work, -1, rwork, info)
#else
            call zheev(job, 'U', ndets, hamil, ndets, eigv, work, -1, rwork, info)
#endif
        else
#ifdef PARALLEL
#ifdef SINGLE_PRECISION
            call pcheev(job, 'U', ndets, hamil, 1, 1,               &
                        proc_blacs_info%desc_m, eigv, eigvec, 1, 1, &
                        proc_blacs_info%desc_m, work, -1, rwork, info)
#else
            call pzheev(job, 'U', ndets, hamil, 1, 1,               &
                        proc_blacs_info%desc_m, eigv, eigvec, 1, 1, &
                        proc_blacs_info%desc_m, work, -1, rwork, info)
#endif
#endif
        end if

        lwork = nint(real(work(1)))
        deallocate(work)
        call check_deallocate('work',ierr)
        ! Now perform the diagonalisation.
        allocate(work(lwork), stat=ierr)
        call check_allocate('work',lwork,ierr)

        if (nprocs == 1) then
            ! Use lapack.
#ifdef SINGLE_PRECISION
            call cheev(job, 'U', ndets, hamil, ndets, eigv, work, lwork, rwork, info)
#else
            call zheev(job, 'U', ndets, hamil, ndets, eigv, work, lwork, rwork, info)
#endif
        else
#ifdef PARALLEL
            ! Use scalapack to do the diagonalisation in parallel.
#ifdef SINGLE_PRECISION
            call pcheev(job, 'U', ndets, hamil, 1, 1,               &
                        proc_blacs_info%desc_m, eigv, eigvec, 1, 1, &
                        proc_blacs_info%desc_m, work, lwork, rwork, info)
#else
            call pzheev(job, 'U', ndets, hamil, 1, 1,               &
                        proc_blacs_info%desc_m, eigv, eigvec, 1, 1, &
                        proc_blacs_info%desc_m, work, lwork, rwork, info)
#endif
#endif
        end if
        deallocate(work, stat=ierr)
        call check_deallocate('work',ierr)
        deallocate(rwork, stat=ierr)
        call check_deallocate('rwork',ierr)
        nwfn = fci_in%print_fci_wfn
        if (nwfn < 0) nwfn = ndets
        do i = 1, nwfn
            if (nprocs == 1) then
                call print_wavefunction(fci_in%print_fci_wfn_file, dets, proc_blacs_info, wfn_complex=hamil(:,i))
            else
                call print_wavefunction(fci_in%print_fci_wfn_file, dets, proc_blacs_info, wfn_complex=eigvec(:,i))
            end if
        end do

        if (nprocs > 1) then
            deallocate(eigvec, stat=ierr)
            call check_deallocate('eigvec',ierr)
        end if

    end subroutine lapack_diagonalisation

end module fci_lapack_complex

module fci_lapack

! Complete diagonalisation of Hamiltonian matrix using lapack/scalapack
! libraries.

use const

implicit none

contains

    subroutine do_fci_lapack(sys, ref_in)

        use dmqmc_procedures, only: setup_rdm_arrays
        use fci_utils, only: init_fci, generate_hamil
        use hamiltonian, only: get_hmatel
        use qmc_data, only: reference_t
        use reference_determinant, only: copy_reference_t
        use system, only: sys_t, copy_sys_spin_info

        use checking, only: check_allocate
        use errors, only: warning
        use parallel, only: parent, nprocs, get_blacs_info

        use calc, only: hamil, proc_blacs_info, fci_rdm, doing_exact_rdm_eigv, fci_rdm_info

        type(sys_t), intent(inout) :: sys
        type(reference_t), intent(in) :: ref_in

        type(sys_t) :: sys_bak
        type(reference_t) :: ref
        integer(i0), allocatable :: dets(:,:)
        real(p), allocatable :: eigv(:), rdm_eigv(:)
        integer :: ndets, ierr, i, rdm_size

        call copy_sys_spin_info(sys, sys_bak)
        call copy_reference_t(ref_in, ref)

        call init_fci(sys, ref, ndets, dets)

        allocate(eigv(ndets), stat=ierr)
        call check_allocate('eigv',ndets,ierr)

        if (ndets == 1) then
            ! The trivial case seems to trip things up...
            eigv(1) = get_hmatel(sys, dets(:,1), dets(:,1))
        else
            if (nprocs == 1) then
                call generate_hamil(sys, ndets, dets, hamil)
            else
                ! Use as square a processor grid as possible.
                proc_blacs_info = get_blacs_info(ndets)
                call generate_hamil(sys, ndets, dets, hamil, proc_blacs_info=proc_blacs_info)
            end if
            call lapack_diagonalisation(sys, dets, eigv)
        end if

        ! [todo] - call from lapack_diagonalisation
        if (doing_exact_rdm_eigv) then
            if (nprocs > 1) then
                if (parent) call warning('diagonalise','RDM eigenvalue calculation is only implemented in serial. Skipping.', 3)
            else
                write(6,'(1x,a46)') "Performing reduced density matrix calculation."
                call setup_rdm_arrays(sys, .false., fci_rdm_info, fci_rdm)
                rdm_size = size(fci_rdm, 1)
                allocate(rdm_eigv(rdm_size), stat=ierr)
                call check_allocate('rdm_eigv',rdm_size,ierr)
                call get_rdm_eigenvalues(sys%basis, fci_rdm_info, ndets, dets, rdm_eigv)
            end if
        end if

        if (parent) then
            write (6,'(1X,"LAPACK diagonalisation results")')
            write (6,'(1X,"^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^",/)')
            write (6,'(1X," State",1X,4X,"Energy")')
            do i = 1, ndets
                write (6,'(1X,i6,1X,f18.12)') i, eigv(i)
            end do
            write (6,'()')
        end if

        ! Return sys in an unaltered state.
        call copy_sys_spin_info(sys_bak, sys)

    end subroutine do_fci_lapack

    subroutine lapack_diagonalisation(sys, dets, eigv)

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

    end subroutine lapack_diagonalisation

    subroutine get_rdm_eigenvalues(basis, rdm_info, ndets, dets, rdm_eigenvalues)

        use basis_types, only: basis_t
        use calc, only: fci_rdm, hamil
        use checking, only: check_allocate, check_deallocate
        use dmqmc_data, only: rdm_t
        use dmqmc_procedures, only: decode_dm_bitstring

        type(basis_t), intent(in) :: basis
        type(rdm_t), intent(inout) :: rdm_info(:)
        integer, intent(in) :: ndets
        integer(i0), intent(in) :: dets(:,:)
        real(p), intent(out) :: rdm_eigenvalues(size(fci_rdm,1))
        integer(i0) :: f1(basis%string_len), f2(basis%string_len)
        integer(i0) :: f3(2*basis%string_len)
        integer :: i, j, rdm_size, info, ierr, lwork
        integer(i0) :: end1, end2
        real(p), allocatable :: work(:)
        real(p) :: rdm_element

        write(6,'(1x,a36)') "Setting up reduced density matrix..."

        ! Loop over all elements of the density matrix and add all contributing elements to the RDM.
        do i = 1, ndets
            do j = 1, ndets
                f1 = iand(rdm_info(1)%B_masks(:,1),dets(:,i))
                f2 = iand(rdm_info(1)%B_masks(:,1),dets(:,j))
                ! If the two bitstrings are the same after bits corresponding to subsystem B have
                ! been unset, then these two bitstrings contribute to the RDM.
                if (sum(abs(f1-f2)) == 0) then
                    ! In f3, concatenate the two bitstrings.
                    f3(1:basis%string_len) = dets(:,i)
                    f3(basis%string_len+1:basis%string_len*2) = dets(:,j)

                    ! Get the position in the RDM of this density matrix element.
                    call decode_dm_bitstring(basis, f3, 1, rdm_info(1))
                    rdm_info(1)%end1 = rdm_info(1)%end1 + 1
                    rdm_info(1)%end2 = rdm_info(1)%end2 + 1

                    ! The ground state wave function is stored in hamil(:,1).
                    rdm_element = hamil(i,1)*hamil(j,1)
                    ! Finally add in the contribution from this density matrix element.
                    fci_rdm(rdm_info(1)%end1(1),rdm_info(1)%end2(1)) = &
                        fci_rdm(rdm_info(1)%end1(1),rdm_info(1)%end2(1)) + rdm_element
                end if
            end do
        end do

        ! Now the RDM is completley calculated.

        ! Calculate the eigenvalues:
        write(6,'(1x,a39,/)') "Diagonalising reduced density matrix..."

        rdm_size = size(fci_rdm, 1)
        ! Find the optimal size of the workspace.
        allocate(work(1), stat=ierr)
        call check_allocate('work',1,ierr)
#ifdef SINGLE_PRECISION
        call ssyev('N', 'U', rdm_size, fci_rdm, rdm_size, rdm_eigenvalues, work, -1, info)
#else
        call dsyev('N', 'U', rdm_size, fci_rdm, rdm_size, rdm_eigenvalues, work, -1, info)
#endif
        lwork = nint(work(1))
        deallocate(work)
        call check_deallocate('work',ierr)

        ! Perform the diagonalisation.
        allocate(work(lwork), stat=ierr)
        call check_allocate('work',lwork,ierr)

#ifdef SINGLE_PRECISION
        call ssyev('N', 'U', rdm_size, fci_rdm, rdm_size, rdm_eigenvalues, work, lwork, info)
#else
        call dsyev('N', 'U', rdm_size, fci_rdm, rdm_size, rdm_eigenvalues, work, lwork, info)
#endif

    end subroutine get_rdm_eigenvalues

end module fci_lapack

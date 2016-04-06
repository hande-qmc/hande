module fci_lapack

! Complete diagonalisation of Hamiltonian matrix using lapack/scalapack
! libraries.

use const

implicit none

contains

    subroutine do_fci_lapack(sys, fci_in, ref_in)

        ! Perform an FCI calculation using LAPACK/ScaLAPACK diagonalisation.

        ! In:
        !    sys: system of interest.
        !    fci_in: fci input options.
        !    ref_in: reference determinant defining (if relevant) a
        !        truncated Hilbert space.

        use dmqmc_procedures, only: setup_rdm_arrays
        use fci_utils, only: fci_in_t, init_fci, generate_hamil, write_hamil, hamil_t
        use hamiltonian, only: get_hmatel
        use reference_determinant, only: reference_t, copy_reference_t
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
        real(p), allocatable :: eigv(:), rdm_eigv(:), rdm(:,:)
        integer :: ndets, ierr, i, rdm_size
        type(blacs_info) :: proc_blacs_info
        type(hamil_t) :: hamil

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
                call generate_hamil(sys, ndets, dets, hamil)
            else
                ! Use as square a processor grid as possible.
                proc_blacs_info = get_blacs_info(ndets, fci_in%block_size)
                call generate_hamil(sys, ndets, dets, hamil, proc_blacs_info=proc_blacs_info)
            end if
            if (fci_in%write_hamiltonian) then
                call write_hamil(fci_in%hamiltonian_file, ndets, proc_blacs_info, hamil)
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
            write (6,'()')
        end if

        ! If requested, calculate and print eigenvalues for an RDM.
        if (allocated(fci_in%subsys_info)) then
            if (nprocs > 1) then
                if (parent) call warning('diagonalise','RDM eigenvalue calculation is only implemented in serial. Skipping.',3)
            else if (sys%system /= heisenberg) then
                if (parent) call warning('diagonalise','RDM calculations are only implemented for Heisenberg systems. Skipping.',3)
            else
                write(6,'(1x,a46)') "Performing reduced density matrix calculation."
                call setup_rdm_arrays(sys, .false., fci_in%subsys_info, rdm)
                rdm_size = size(rdm, 1)
                allocate(rdm_eigv(rdm_size), stat=ierr)
                call check_allocate('rdm_eigv',rdm_size,ierr)
                call get_rdm_eigenvalues(sys%basis, fci_in%subsys_info, ndets, dets, hamil%rmat, rdm, rdm_eigv)

                write (6,'(1X,"RDM eigenvalues")')
                write (6,'(1X,"^^^^^^^^^^^^^^^",/)')
                write (6,'(1X," State",5X,"RDM eigenvalue")')
                do i = 1, rdm_size
                    write (6,'(1X,i6,1X,f18.12)') i, rdm_eigv(i)
                end do
                write (6,'()')
            end if
        end if

        ! Return sys in an unaltered state.
        call copy_sys_spin_info(sys_bak, sys)

    end subroutine do_fci_lapack

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
        !    hamil: hamil_t derived type containing Hamiltonian matrix of the system
        !        in the Hilbert space used in appropriate format.
        !        On output contains the eigenvectors, if requested, or is
        !        otherwise destroyed.
        ! Out:
        !    eigv(ndets): Lanczos eigenvalues of the current block of the
        !        Hamiltonian matrix.

        use checking, only: check_allocate, check_deallocate
        use linalg, only: syev_wrapper, heev_wrapper, psyev_wrapper, pheev_wrapper
        use parallel, only: parent, nprocs, blacs_info
        use system, only: sys_t

        use fci_utils, only: fci_in_t, hamil_t
        use operators

        type(sys_t), intent(in) :: sys
        type(fci_in_t), intent(in) :: fci_in
        integer(i0), intent(in) :: dets(:,:)
        type(blacs_info), intent(in) :: proc_blacs_info
        type(hamil_t), intent(inout) :: hamil
        real(p), intent(out) :: eigv(:)
        real(p), allocatable :: rwork(:), eigvec(:,:)
        complex(p), allocatable :: ceigvec(:,:)
        integer :: info, ierr, i, nwfn, ndets
        character(1) :: job

        ndets = ubound(dets, dim=2)

        if (parent) then
            write (6,'(1X,a35,/)') 'Performing exact diagonalisation...'
        end if

        if (fci_in%analyse_fci_wfn /= 0 .or. fci_in%print_fci_wfn /= 0 &
                    .or. allocated(fci_in%subsys_info)) then
            job = 'V'
        else
            job = 'N'
        end if

        if (nprocs > 1) then
            if (hamil%comp) then
                allocate(ceigvec(proc_blacs_info%nrows, proc_blacs_info%ncols), stat=ierr)
                call check_allocate('ceigvec',proc_blacs_info%nrows*proc_blacs_info%ncols,ierr)
            else
                allocate(eigvec(proc_blacs_info%nrows, proc_blacs_info%ncols), stat=ierr)
                call check_allocate('eigvec',proc_blacs_info%nrows*proc_blacs_info%ncols,ierr)
            end if
        end if

        ! Find the optimal size of the workspace and then perform the diagonalisation.
        if (hamil%comp) then
            allocate(rwork(max(1,3*ndets-2)), stat = ierr)
            call check_allocate('rwork',max(1, 3*ndets - 2), ierr)
            if (nprocs == 1) then
                call heev_wrapper(job, 'U', ndets, hamil%cmat, ndets, eigv, rwork, info)
            else
                call pheev_wrapper(job, 'U', ndets, hamil%cmat, 1, 1,               &
                            proc_blacs_info%desc_m, eigv, ceigvec, 1, 1, &
                            proc_blacs_info%desc_m, rwork, size(rwork), info)
            end if
            deallocate(rwork, stat=ierr)
            call check_deallocate('rwork',ierr)
            if (fci_in%analyse_fci_wfn /= 0) then
                write(6,'(1x,a36)') "Complex wavefunction analysis and printing not implemented. Skipping."
            end if
        else
            if (nprocs == 1) then
                call syev_wrapper(job, 'U', ndets, hamil%rmat, ndets, eigv, info)
            else
                call psyev_wrapper(job, 'U', ndets, hamil%rmat, 1, 1,           &
                            proc_blacs_info%desc_m, eigv, eigvec, 1, 1, &
                            proc_blacs_info%desc_m, info)
            end if
            nwfn = fci_in%analyse_fci_wfn
            if (nwfn < 0) nwfn = ndets
            do i = 1, nwfn
                if (nprocs == 1) then
                    call analyse_wavefunction(sys, hamil%rmat(:,i), dets, proc_blacs_info)
                else
                    call analyse_wavefunction(sys, eigvec(:,i), dets, proc_blacs_info)
                end if
            end do
            nwfn = fci_in%print_fci_wfn
            if (nwfn < 0) nwfn = ndets
            do i = 1, nwfn
                if (nprocs == 1) then
                    call print_wavefunction(fci_in%print_fci_wfn_file, dets, proc_blacs_info, hamil%rmat(:,i))
                else
                    call print_wavefunction(fci_in%print_fci_wfn_file, dets, proc_blacs_info, eigvec(:,i))
                end if
            end do
        end if

        if (nprocs > 1) then
            if (allocated(eigvec)) then
                deallocate(eigvec, stat=ierr)
                call check_deallocate('eigvec',ierr)
            end if
            if (allocated(ceigvec)) then
                deallocate(ceigvec, stat=ierr)
                call check_deallocate('ceigvec',ierr)
            end if
        end if

    end subroutine lapack_diagonalisation

    subroutine get_rdm_eigenvalues(basis, subsys_info, ndets, dets, eigvec, rdm, rdm_eigv)

        ! Take an FCI ground-state wave function (eigvec) and use it to
        ! calculate and return the RDM eigenvalues, for the RDM and subsystem
        ! specified by subsys_info.

        ! In:
        !    basis: type containing information about the basis set.
        !    ndets: number of determinants in the Hilbert space.
        !    dets: list of determinants in the Hilbert space (bit-string
        !        representation).
        !    eigvec: the eigenvector from which the RDM is to be calculated.
        ! In/Out:
        !    subsys_info: information about the subsystem for the RDM to be
        !        calculated.
        ! Out:
        !    rdm: space for the RDM to be accumulated. The array passed in
        !        must be the correct size. On output the RDM will have been
        !        destroyed and so should *not* be used.
        !    rdm_eigv: The eigenvalues of the calculated RDM.

        use basis_types, only: basis_t
        use linalg, only: syev_wrapper
        use dmqmc_data, only: subsys_t
        use dmqmc_procedures, only: decode_dm_bitstring

        type(basis_t), intent(in) :: basis
        type(subsys_t), intent(inout) :: subsys_info(:)
        integer, intent(in) :: ndets
        integer(i0), intent(in) :: dets(:,:)
        real(p), intent(in) :: eigvec(:,:)
        real(p), intent(out) :: rdm(:,:)
        real(p), intent(out) :: rdm_eigv(size(rdm,1))

        integer(i0) :: f1(basis%string_len), f2(basis%string_len)
        integer(i0) :: f3(2*basis%string_len)
        integer :: i, j, rdm_size, info
        integer(i0) :: rdm_f1(subsys_info(1)%string_len), rdm_f2(subsys_info(1)%string_len)
        real(p) :: rdm_element

        write(6,'(1x,a36)') "Setting up reduced density matrix..."

        ! Loop over all elements of the density matrix and add all contributing elements to the RDM.
        do i = 1, ndets
            do j = 1, ndets
                f1 = iand(subsys_info(1)%B_masks(:,1),dets(:,i))
                f2 = iand(subsys_info(1)%B_masks(:,1),dets(:,j))
                ! If the two bitstrings are the same after bits corresponding to subsystem B have
                ! been unset, then these two bitstrings contribute to the RDM.
                if (sum(abs(f1-f2)) == 0) then
                    ! In f3, concatenate the two bitstrings.
                    f3(1:basis%string_len) = dets(:,i)
                    f3(basis%string_len+1:basis%string_len*2) = dets(:,j)

                    ! Get the position in the RDM of this density matrix element.
                    call decode_dm_bitstring(basis, f3, 1, subsys_info(1), rdm_f1, rdm_f2)
                    rdm_f1 = rdm_f1 + 1
                    rdm_f2 = rdm_f2 + 1

                    ! The ground state wave function is stored in eigvec(:,1).
                    rdm_element = eigvec(i,1)*eigvec(j,1)
                    ! Finally add in the contribution from this density matrix element.
                    rdm(rdm_f1(1),rdm_f2(1)) = rdm(rdm_f1(1),rdm_f2(1)) + rdm_element
                end if
            end do
        end do

        ! Now the RDM is completley calculated.

        ! Calculate the eigenvalues:
        write(6,'(1x,a39,/)') "Diagonalising reduced density matrix..."

        rdm_size = size(rdm, 1)
        call syev_wrapper('N', 'U', rdm_size, rdm, rdm_size, rdm_eigv, info)

    end subroutine get_rdm_eigenvalues

end module fci_lapack

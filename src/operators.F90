module operators

! Module for operators other than the Hamiltonian operator.

implicit none

contains

!--- Kinetic energy ---

    pure function kinetic0_hub_k(f) result(kin)

        ! In:
        !    f: bit string representation of the Slater determinant.
        ! Returns:
        !    < D_i | T | D_i >, the diagonal Kinetic matrix elements, for
        !        the Hubbard model in momentum space.

        use basis, only: basis_length, basis_fns
        use determinants, only: decode_det
        use system, only: nel

        use const, only: p, i0

        real(p) :: kin
        integer(i0), intent(in) :: f(basis_length)

        integer :: i, occ(nel)

        ! Kinetic operator is diagonal in the Hubbard model in the Bloch basis.

        ! <D|T|D> = \sum_k \epsilon_k
        kin = 0.0_p
        call decode_det(f, occ)
        do i = 1, nel
            kin = kin + basis_fns(occ(i))%sp_eigv
        end do

    end function kinetic0_hub_k

!-- Debug/test routines for operating on exact wavefunction ---

    subroutine analyse_wavefunction(wfn)

        ! Analyse an exact wavefunction using the desired operator(s).

        ! In:
        !    wfn: exact wavefunction to be analysed.  wfn(i) = c_i, where
        !    |\Psi> = \sum_i c_i|D_i>.

        use const, only: i0, p
        use basis, only: nbasis, basis_fns, basis_length
        use calc, only: proc_blacs_info, distribute, distribute_off
        use determinants, only: dets_list, ndets
        use parallel

        real(p), intent(in) :: wfn(proc_blacs_info%nrows)

        integer :: idet, iorb, i, ii, ilocal

#ifdef PARALLEL
        integer :: ierr
#endif

        ! TODO: tidy and generalise.  Currently not functional!

        do i = 1, proc_blacs_info%nrows, block_size
            do ii = 1, min(block_size, proc_blacs_info%nrows - i + 1)
                ilocal = i - 1 + ii
                idet =  (i-1)*nproc_rows + proc_blacs_info%procx* block_size + ii
                ! TODO: Evaluate operator.
            end do
        end do

    end subroutine analyse_wavefunction

    subroutine print_wavefunction(filename, wfn)

        ! Print out an exact wavefunction.

        ! In:
        !    filename: file to be printed to.
        !    wfn: exact wavefunction to be printed out.  wfn(i) = c_i, where
        !    |\Psi> = \sum_i c_i|D_i>.

        use const, only: i0, p
        use basis, only: nbasis, basis_fns, basis_length
        use calc, only: proc_blacs_info, distribute, distribute_off
        use determinants, only: dets_list, ndets

        use checking, only: check_allocate, check_deallocate
        use parallel

        character(*), intent(in) :: filename
        real(p), intent(in) :: wfn(proc_blacs_info%nrows)

        integer :: idet, i, ii, ilocal

#ifdef PARALLEL
        integer :: ierr
        integer :: proc_info(2, 0:nprocs-1), info(2), proc
        real(p), allocatable :: wfn_recv(:)
        integer :: stat(MPI_STATUS_SIZE)
        integer, parameter :: comm_tag = 123
#endif

        if (parent) open(12,file=filename,status='unknown')

        if (nprocs == 1) then
            do idet = 1, size(wfn)
                write (12,*) dets_list(:,idet), wfn(idet)
            end do
        else
#ifdef PARALLEL
            ! Set up for receiving parts of the wavefunction from other
            ! processors.
            info = (/proc_blacs_info%nrows, proc_blacs_info%procx/)
            call mpi_gather(info, 2, mpi_integer, proc_info, 2, &
                             mpi_integer, root, mpi_comm_world, ierr)

            if (parent) then
                allocate(wfn_recv(maxval(proc_info(1,:))), stat=ierr)
                call check_allocate('wfn_local', maxval(proc_info(1,:)), ierr)
                ! Write out from root.
                call write_wavefunction_parallel(proc_blacs_info%nrows, proc_blacs_info%procx, wfn)
                ! Write out from other processors.
                do proc = 1, nprocs-1
                    call mpi_recv(wfn_recv, proc_info(1,proc), mpi_preal, proc, comm_tag, mpi_comm_world, stat, ierr)
                    call write_wavefunction_parallel(proc_info(1,proc), proc_info(2,proc), wfn_recv)
                end do
                deallocate(wfn_recv, stat=ierr)
                call check_deallocate('wfn_local', ierr)
            else
                ! Send data from from other processors.
                call mpi_send(wfn, proc_blacs_info%nrows, mpi_preal, root, comm_tag, mpi_comm_world, ierr)
            end if
#endif
        end if

        if (parent) close(12,status='keep')

        contains

            subroutine write_wavefunction_parallel(nrows, procx, wfn_curr)

                integer, intent(in) :: nrows, procx
                real(p), intent(in) :: wfn_curr(nrows)

                do i = 1, nrows, block_size
                    do ii = 1, min(block_size, nrows - i + 1)
                        ilocal = i - 1 + ii
                        idet =  (i-1)*nproc_rows + procx* block_size + ii
                        write (12,*) dets_list(:,idet), wfn_curr(ilocal)
                    end do
                end do

            end subroutine write_wavefunction_parallel

    end subroutine print_wavefunction

end module operators

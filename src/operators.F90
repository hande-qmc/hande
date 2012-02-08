module operators

! Module for operators other than the Hamiltonian operator.

implicit none

contains

    function calc_orb_occ(f, orb_mask) result(nocc)

        ! Evaluate n_i | D > = f_i | D >
        ! where n_i is the number operator of the i-th orbital, D is
        ! a determinant and f_i is the occupancy of the i-th orbital in the
        ! determinant.

        ! In:
        !    f: bit-string representation of the determinant.
        !    orb_mask: bit-string mask with the bits corresponding to the
        !    desired orbitals set.
        ! Returns:
        !    nocc: number of electrons in the orbitals specified in orb_mask.

        use basis, only: basis_length
        use bit_utils, only: count_set_bits
        use const, only: i0

        integer :: nocc
        integer(i0), intent(in) :: f(basis_length), orb_mask(basis_length)

        nocc = sum(count_set_bits(iand(f,orb_mask)))

    end function calc_orb_occ

    subroutine analyse_wavefunction(wfn)

        ! Analyse an exact wavefunction using the desired operator(s).

        ! In:
        !    wfn: exact wavefunction to be analysed.  wfn(i) = c_i, where
        !    |\Psi> = \sum_i c_i|D_i>.

        use const, only: i0, p
        use basis, only: nbasis, basis_fns, basis_length, set_orb_mask
        use calc, only: proc_blacs_info, distribute, distribute_off
        use determinants, only: dets_list, ndets
        use parallel

        real(p), intent(in) :: wfn(proc_blacs_info%nrows)

        integer(i0) :: orb_mask(basis_length)

        integer :: idet, iorb, i, ii, ilocal
        integer :: ldone(nbasis), l2

        real(p) :: orb_occ(nbasis)
#ifdef PARALLEL
        integer :: ierr
        real(p) :: orb_occ_recv(nbasis)
#endif


        ! Find the number of electrons in each group of symmetry-related orbitals.
        ldone = -1
        orb_occ = 0.0_p

        write (6,*) 'analyse', iproc, size(wfn), ndets, proc_blacs_info%nrows, 'shape',shape(dets_list)

        do iorb = 1, nbasis
            l2 = dot_product(basis_fns(iorb)%l, basis_fns(iorb)%l)
            if (any(ldone == l2)) cycle
            call set_orb_mask(l2, orb_mask)
            if (nprocs == 1) then
                do idet = 1, size(wfn)
                    orb_occ(iorb) = orb_occ(iorb) + wfn(idet)**2*calc_orb_occ(dets_list(:,idet), orb_mask)
                end do
            else
                do i = 1, proc_blacs_info%nrows, block_size
                    do ii = 1, min(block_size, proc_blacs_info%nrows - i + 1)
                        ilocal = i - 1 + ii
                        idet =  (i-1)*nproc_rows + proc_blacs_info%procx* block_size + ii
                        orb_occ(iorb) = orb_occ(iorb) + wfn(ilocal)**2*calc_orb_occ(dets_list(:,idet), orb_mask)
                        if (iorb == 3) write (26+iproc,*) dets_list(:,idet),wfn(ilocal), calc_orb_occ(dets_list(:,idet), orb_mask)
                    end do
                end do
            end if
            ldone(iorb) = l2
            write (6,*) 'results', iproc, iorb, orb_occ(iorb)
        end do

#ifdef PARALLEL
        call mpi_reduce(orb_occ, orb_occ_recv, nbasis, mpi_preal, mpi_sum, root, mpi_comm_world, ierr)
        orb_occ = orb_occ_recv
#endif

        if (parent) then
            do iorb = 1, nbasis
                if (ldone(iorb) /= -1) then
                    write (6,'(1X,i6,f16.8)') iorb, orb_occ(iorb)
                end if
            end do
        end if

    end subroutine analyse_wavefunction

    subroutine print_wavefunction(filename, wfn)

        ! Print out an exact wavefunction.

        ! In:
        !    filename: file to be printed to.
        !    wfn: exact wavefunction to be printed out.  wfn(i) = c_i, where
        !    |\Psi> = \sum_i c_i|D_i>.

        use const, only: i0, p
        use basis, only: nbasis, basis_fns, basis_length, set_orb_mask
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

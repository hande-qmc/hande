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
                    if (iorb == 3) write (16,*) dets_list(:,idet),wfn(idet), calc_orb_occ(dets_list(:,idet), orb_mask)
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

end module operators

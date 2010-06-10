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
        use determinants, only: dets_list, ndets

        real(p), intent(in) :: wfn(ndets)

        integer(i0) :: orb_mask(basis_length)

        integer :: idet, iorb
        integer :: ldone(nbasis), l2

        real(p) :: orb_occ(nbasis)


        ! Find the number of electrons in each group of symmetry-related orbitals.
        ldone = -1
        orb_occ = 0.0_p

        do iorb = 1, nbasis
            l2 = dot_product(basis_fns(iorb)%l, basis_fns(iorb)%l)
            if (any(ldone == l2)) cycle
            call set_orb_mask(l2, orb_mask)
            do idet = 1, ndets
                orb_occ(iorb) = orb_occ(iorb) + wfn(idet)**2*calc_orb_occ(dets_list(:,idet), orb_mask)
            end do
            ldone(iorb) = l2
        end do

    end subroutine analyse_wavefunction

end module operators

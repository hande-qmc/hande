module hamiltonian_molecular

! Module for evaluating Hamiltonian matrix elements for molecular systems (ie
! systems where the integrals have been read in from file).

! TODO:
! * optimise integral cache access.

use const, only: p, i0

implicit none

contains

    pure function get_hmatel_mol(f1, f2) result(hmatel)

        ! In:
        !    f1, f2: bit string representation of the Slater
        !        determinants D1 and D2 respectively.
        ! Returns:
        !    Hamiltonian matrix element between the two determinants, 
        !    < D1 | H | D2 >, where the determinants are formed from
        !    molecular orbitals read in from an FCIDUMP file.
        !
        ! Used in the read_in system only.

        use basis, only: basis_length
        use determinants, only: decode_det
        use excitations, only: excit, get_excitation
        use system, only: nel

        real(p) :: hmatel
        integer(i0), intent(in) :: f1(basis_length), f2(basis_length)

        type(excit) :: excitation
        integer :: occ_list(nel)

        hmatel = 0.0_p

        ! Test to see if matrix element is non-zero.
        excitation = get_excitation(f1, f2)

        ! Connected determinants can differ by (at most) 2 spin orbitals.
        if (excitation%nexcit <= 2) then

            select case(excitation%nexcit)
            ! Apply Slater--Condon rules.

            case(0)

                ! < D | H | D > = \sum_i < i | h(i) | i > + \sum_i \sum_{j>i} < ij || ij >
                hmatel = slater_condon0_mol(f1)

            case(1)

                ! < D | H | D_i^a > = < i | h(a) | a > + \sum_j < ij || aj >
                call decode_det(f1, occ_list)
                hmatel = slater_condon1_mol(occ_list, excitation%from_orb(1), &
                                            excitation%to_orb(1), excitation%perm)

            case(2)

                ! < D | H | D_{ij}^{ab} > = < ij || ab >

                ! Two electron operator
                hmatel = slater_condon2_mol(excitation%from_orb(1), excitation%from_orb(2), &
                                            & excitation%to_orb(1), excitation%to_orb(2), excitation%perm)
            end select

        end if

    end function get_hmatel_mol

    pure function slater_condon0_mol(f) result(hmatel)

        ! In:
        !    f: bit string representation of a Slater determinant, |D>.
        ! Returns:
        !    <D|H|D>, the diagonal Hamiltonian matrix element involving D for
        !    systems defined by integrals read in from an FCIDUMP file.

        use basis, only: basis_length
        use determinants, only: decode_det
        use molecular_integrals, only: get_one_body_int_mol, get_two_body_int_mol, &
                                       get_two_body_int_mol_nonzero, one_e_h_integrals, &
                                       coulomb_integrals
        use system, only: nel, Ecore

        real(p) :: hmatel
        integer(i0), intent(in) :: f(basis_length)

        integer :: occ_list(nel)
        integer :: iel, jel, i, j

        ! < D | H | D > = Ecore + \sum_i < i | h(i) | i > + \sum_i \sum_{j>i} < ij || ij >
        call decode_det(f, occ_list)

        hmatel = Ecore
        do iel = 1, nel
            i = occ_list(iel)
            hmatel = hmatel + get_one_body_int_mol(one_e_h_integrals, i, i)
            do jel = iel, nel
                j = occ_list(jel)
                hmatel = hmatel + get_two_body_int_mol_nonzero(coulomb_integrals, i, j, i, j) &
                                - get_two_body_int_mol(coulomb_integrals, i, j, j, i)
            end do
        end do

    end function slater_condon0_mol

    pure function slater_condon1_mol(occ_list, i, a, perm) result(hmatel)

        ! In:
        !    occ_list: list of orbitals occupied in a Slater Determinant, D.
        !    i: the spin-orbital from which an electron is excited in
        !       the reference determinant.
        !    a: the spin-orbital into which an electron is excited in
        !       the excited determinant.
        !    perm: true if an odd number of permutations are required for
        !       D and D_i^a to be maximally coincident.
        ! Returns:
        !    < D | H | D_i^a >, the Hamiltonian matrix element between a 
        !        determinant and a single excitation of it for systems defined
        !        by integrals read in from an FCIDUMP file.

        use molecular_integrals, only: get_one_body_int_mol, get_two_body_int_mol, &
                                       one_e_h_integrals, coulomb_integrals
        use system, only: nel

        real(p) :: hmatel
        integer, intent(in) :: occ_list(nel), i, a
        logical, intent(in) :: perm

        integer :: iel

        ! < D | H | D_i^a > = < i | h(a) | a > + \sum_j < ij || aj >

        hmatel = get_one_body_int_mol(one_e_h_integrals, i, a)

        do iel = 1, nel
            if (occ_list(iel) /= i) &
                hmatel = hmatel + get_two_body_int_mol(coulomb_integrals, i, occ_list(iel), a, occ_list(iel)) &
                                - get_two_body_int_mol(coulomb_integrals, i, occ_list(iel), occ_list(iel), a)
        end do

        if (perm) hmatel = -hmatel

    end function slater_condon1_mol

    elemental function slater_condon2_mol(i, j, a, b, perm) result(hmatel)

        ! In:
        !    i: the spin-orbital from which an electron is excited in
        !       the reference determinant.
        !    j: the spin-orbital from which an electron is excited in
        !       the reference determinant.
        !    a: the spin-orbital into which an electron is excited in
        !       the excited determinant.
        !    b: the spin-orbital into which an electron is excited in
        !       the excited determinant.
        !    perm: true if an odd number of permutations are required for
        !       D and D_i^a to be maximally coincident.
        ! Returns:
        !    < D | H | D_i^a >, the Hamiltonian matrix element between a 
        !        determinant and a single excitation of it for systems defined
        !        by integrals read in from an FCIDUMP file.
        !
        ! Note that we don't actually need to directly refer to |D>.

        use molecular_integrals, only: get_two_body_int_mol, coulomb_integrals
        use system, only: nel

        real(p) :: hmatel
        integer, intent(in) :: i, j, a, b
        logical, intent(in) :: perm

        ! < D | H | D_{ij}^{ab} > = < ij || ab >

        hmatel = get_two_body_int_mol(coulomb_integrals, i, j, a, b) - get_two_body_int_mol(coulomb_integrals, i, j, b, a)

        if (perm) hmatel = -hmatel

    end function slater_condon2_mol

endmodule hamiltonian_molecular

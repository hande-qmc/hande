module hamiltonian_trotter

! Module for evaluating matrix elements for commutators of VTV type in
! Trotterization studies.

use const

implicit none

contains

    pure function get_hmatel_trotter(sys, f1, f2) result(hmatel)

        ! In:
        !    sys: system to be studied.
        !    f1, f2: bit string representation of the Slater determinants
        !        D1 and D2 respectively.
        ! Returns:
        !    Hamiltonian matrix element between the two determinants,
        !    < D1 | H | D2 >, where the determinants are formed from
        !    molecular orbitals read in from an FCIDUMP file.

        ! Used in the Trotter VTV-type bound calculation only.

        use determinants, only: decode_det
        use excitations, only: excit_t, get_excitation
        use system, only: sys_t
        use hamiltonian_data

        type(hmatel_t) :: hmatel
        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f1(sys%basis%tot_string_len), f2(sys%basis%tot_string_len)

        type(excit_t) :: excitation
        integer :: occ_list(sys%nel)

        excitation = get_excitation(sys%nel, sys%basis, f1, f2)

        select case(excitation%nexcit)
        case(1)

            call decode_det(sys%basis, f1, occ_list)
            hmatel%r = offdiagonal_element_trotter(sys, occ_list, excitation%from_orb(1), excitation%to_orb(1))

        case default

            hmatel%r = 0.0_p

        end select

    end function get_hmatel_trotter

    pure function offdiagonal_element_trotter(sys, occ_list, i, a) result(hmatel)

        ! In:
        !    sys: system to be studied.
        !    occ_list: list of orbitals occupied in a Slater Determinant, D.
        !    i: the spin-orbital from which an electron is excited in
        !       the reference determinant.
        !    a: the spin-orbital into which an electron is excited in
        !       the excited determinant.
        ! Returns:
        !    < D | H | D_i^a >, the commutator matrix element between a
        !        determinant and a single excitation of it.

        use molecular_integrals, only: get_one_body_int_mol_real, get_two_body_int_mol_real
        use system, only: sys_t

        real(p) :: hmatel, one_body_elem
        type(sys_t), intent(in) :: sys
        integer, intent(in) :: occ_list(sys%nel), i, a

        integer :: iel

        associate(one_e_ints=>sys%read_in%one_e_h_integrals, coulomb_ints=>sys%read_in%coulomb_integrals)
            one_body_elem = get_one_body_int_mol_real(one_e_ints, i, a, sys)
            hmatel = 0.0_p

            do iel = 1, sys%nel
                if (occ_list(iel) /= i) then
                    hmatel = hmatel &
                        + get_two_body_int_mol_real(coulomb_ints, i, occ_list(iel), i, occ_list(iel), sys)
                    hmatel = hmatel - get_two_body_int_mol_real(coulomb_ints, a, occ_list(iel), &
                                        a, occ_list(iel), sys)
                end if
            end do
        end associate

        ! Note that we take the absolute value. This is because we are
        ! implementing the sign problem free version of the operator.
        hmatel = - abs(one_body_elem) * (hmatel**2)

    end function offdiagonal_element_trotter

    pure function slater_condon0_trotter(sys, f) result(hmatel)

        ! In:
        !    sys: system to be studied.
        !    f: bit string representation of a Slater determinant, |D>.
        ! Returns:
        !    <D|H|D>, the diagonal Hamiltonian matrix element involving D for
        !    systems defined by integrals read in from an FCIDUMP file.

        use system, only: sys_t

        real(p) :: hmatel
        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f(sys%basis%tot_string_len)

        ! This commutator as zero diagonal elements.
        hmatel = 0.0_p

    end function slater_condon0_trotter

    pure function slater_condon1_trotter_excit(sys, occ_list, i, a, perm) result(hmatel)

        ! In:
        !    sys: system to be studied.
        !    occ_list: list of orbitals occupied in a Slater Determinant, D.
        !    i: the spin-orbital from which an electron is excited in
        !       the reference determinant.
        !    a: the spin-orbital into which an electron is excited in
        !       the excited determinant.
        !    perm: true if an odd number of permutations are required for
        !       D and D_i^a to be maximally coincident.
        ! Returns:
        !    < D | H | D_i^a >, the commutator matrix element between a
        !        determinant and a single excitation of it.

        use molecular_integrals, only: get_one_body_int_mol_nonzero, get_two_body_int_mol_nonzero
        use system, only: sys_t
        use hamiltonian_data, only: hmatel_t

        type(hmatel_t) :: hmatel
        real(p) :: one_body_elem
        type(sys_t), intent(in) :: sys
        integer, intent(in) :: occ_list(sys%nel), i, a
        logical, intent(in) :: perm

        integer :: iel

        associate(basis_fns=>sys%basis%basis_fns, &
                  one_e_ints=>sys%read_in%one_e_h_integrals, &
                  coulomb_ints=>sys%read_in%coulomb_integrals)
            one_body_elem = get_one_body_int_mol_nonzero(one_e_ints, i, a, basis_fns)

            hmatel%r = 0.0_p
            do iel = 1, sys%nel
                if (occ_list(iel) /= i) then
                    hmatel%r = hmatel%r &
                                + get_two_body_int_mol_nonzero(coulomb_ints, i, occ_list(iel), i, occ_list(iel), basis_fns)
                    hmatel%r = hmatel%r &
                            - get_two_body_int_mol_nonzero(coulomb_ints, a, occ_list(iel), a, occ_list(iel), basis_fns)
                end if
            end do
        end associate

        ! Note that we take the absolute value. This is because we are
        ! implementing the sign problem free version of the operator.
        hmatel%r = - abs(one_body_elem) * (hmatel%r**2)

    end function slater_condon1_trotter_excit

    pure function slater_condon2_trotter_excit(sys, i, j, a, b, perm) result(hmatel)

        ! In:
        !    sys: system to be studied.
        !    i: the spin-orbital from which an electron is excited in
        !       the reference determinant.
        !    j: the spin-orbital from which an electron is excited in
        !       the reference determinant.
        !    a: the spin-orbital into which an electron is excited in
        !       the excited determinant.
        !    b: the spin-orbital into which an electron is excited in
        !       the excited determinant.
        !    perm: true if an odd number of permutations are required for
        !       D and D_{ij}^{ab} to be maximally coincident.
        ! Returns:
        !    < D | H | D_{ij}^{ab} >, the commutator matrix element of
        !        the commutator, which is always equal to 0.

        use molecular_integrals, only: get_two_body_int_mol_nonzero
        use system, only: sys_t
        use hamiltonian_data, only: hmatel_t

        type(hmatel_t) :: hmatel
        type(sys_t), intent(in) :: sys
        integer, intent(in) :: i, j, a, b
        logical, intent(in) :: perm

        ! There are no double excitations for VTV commutators, so we
        ! can just return 0.
        hmatel%r = 0.0_p

    end function slater_condon2_trotter_excit

end module hamiltonian_trotter

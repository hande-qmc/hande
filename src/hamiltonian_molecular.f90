module hamiltonian_molecular

! Module for evaluating Hamiltonian matrix elements for molecular systems (ie
! systems where the integrals have been read in from file).

use const, only: p, i0

implicit none

contains

    pure function get_hmatel_mol(sys, f1, f2) result(hmatel)

        ! In:
        !    sys: system to be studied.
        !    f1, f2: bit string representation of the Slater
        !        determinants D1 and D2 respectively.
        ! Returns:
        !    Hamiltonian matrix element between the two determinants,
        !    < D1 | H | D2 >, where the determinants are formed from
        !    molecular orbitals read in from an FCIDUMP file.
        !
        ! Used in the read_in system only.

        use determinants, only: decode_det
        use excitations, only: excit_t, get_excitation
        use system, only: sys_t
        use hamiltonian_data

        type(hmatel_t) :: hmatel
        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f1(sys%basis%tot_string_len), f2(sys%basis%tot_string_len)

        type(excit_t) :: excitation
        integer :: occ_list(sys%nel)

        hmatel%r = 0.0_p

        ! Test to see if matrix element is non-zero.
        excitation = get_excitation(sys%nel, sys%basis, f1, f2)

        ! Connected determinants can differ by (at most) 2 spin orbitals.
        if (excitation%nexcit <= 2) then

            select case(excitation%nexcit)
            ! Apply Slater--Condon rules.

            case(0)

                ! < D | H | D > = \sum_i < i | h(i) | i > + \sum_i \sum_{j>i} < ij || ij >
                hmatel%r = slater_condon0_mol(sys, f1)

            case(1)

                ! < D | H | D_i^a > = < i | h(a) | a > + \sum_j < ij || aj >
                call decode_det(sys%basis, f1, occ_list)
                hmatel%r = slater_condon1_mol(sys, occ_list, excitation%from_orb(1), &
                                            excitation%to_orb(1), excitation%perm)

            case(2)

                ! < D | H | D_{ij}^{ab} > = < ij || ab >

                ! Two electron operator
                hmatel%r = slater_condon2_mol(sys, excitation%from_orb(1), excitation%from_orb(2), &
                                            & excitation%to_orb(1), excitation%to_orb(2), excitation%perm)
            end select

        end if

    end function get_hmatel_mol

    pure function slater_condon0_mol(sys, f) result(hmatel)

        ! In:
        !    sys: system to be studied.
        !    f: bit string representation of a Slater determinant, |D>.
        ! Returns:
        !    <D|H|D>, the diagonal Hamiltonian matrix element involving D for
        !    systems defined by integrals read in from an FCIDUMP file.

        use determinants, only: decode_det
        use molecular_integrals, only: get_one_body_int_mol_nonzero, get_two_body_int_mol_nonzero
        use system, only: sys_t

        real(p) :: hmatel
        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f(sys%basis%tot_string_len)

        integer :: occ_list(sys%nel)

        ! < D | H | D > = Ecore + \sum_i < i | h(i) | i > + \sum_i \sum_{j>i} < ij || ij >
        call decode_det(sys%basis, f, occ_list)

        hmatel = slater_condon0_mol_orb_list(sys, occ_list)

    end function slater_condon0_mol

    pure function slater_condon0_mol_orb_list(sys, occ_list) result(hmatel)

        ! In:
        !    sys: system to be studied.
        !    occ_list: List of occpied orbitals in Slater determinant |D>.
        ! Returns:
        !    <D|H|D>, the diagonal Hamiltonian matrix element involving D for
        !    systems defined by integrals read in from an FCIDUMP file.

        use molecular_integrals, only: get_one_body_int_mol_nonzero, get_two_body_int_mol_nonzero
        use system, only: sys_t

        real(p) :: hmatel
        type(sys_t), intent(in) :: sys
        integer, intent(in) :: occ_list(:)

        integer :: iel, jel, i, j

        associate(one_e_ints=>sys%read_in%one_e_h_integrals, coulomb_ints=>sys%read_in%coulomb_integrals)
            hmatel = sys%read_in%Ecore
            do iel = 1, sys%nel
                i = occ_list(iel)
                hmatel = hmatel + get_one_body_int_mol_nonzero(one_e_ints, i, i, sys%basis%basis_fns)
                do jel = iel+1, sys%nel
                    j = occ_list(jel)
                    hmatel = hmatel + get_two_body_int_mol_nonzero(coulomb_ints, i, j, i, j, sys%basis%basis_fns)
                    if (sys%basis%basis_fns(i)%Ms == sys%basis%basis_fns(j)%Ms) &
                                  hmatel = hmatel - get_two_body_int_mol_nonzero(coulomb_ints, i, j, j, i, sys%basis%basis_fns)
                end do
            end do
        end associate

    end function slater_condon0_mol_orb_list

    pure function slater_condon1_mol(sys, occ_list, i, a, perm) result(hmatel)

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
        !    < D | H | D_i^a >, the Hamiltonian matrix element between a
        !        determinant and a single excitation of it for systems defined
        !        by integrals read in from an FCIDUMP file.

        use molecular_integrals, only: get_one_body_int_mol_real, get_two_body_int_mol_real
        use system, only: sys_t

        real(p) :: hmatel
        type(sys_t), intent(in) :: sys
        integer, intent(in) :: occ_list(sys%nel), i, a
        logical, intent(in) :: perm

        integer :: iel

        ! Check that this excitation is symmetry allowed.
        if (sys%basis%basis_fns(i)%sym/=sys%basis%basis_fns(a)%sym) then
            hmatel = 0.0_p
        else
            ! < D | H | D_i^a > = < i | h(a) | a > + \sum_j < ij || aj >

            associate(one_e_ints=>sys%read_in%one_e_h_integrals, coulomb_ints=>sys%read_in%coulomb_integrals)
                hmatel = get_one_body_int_mol_real(one_e_ints, i, a, sys)

                do iel = 1, sys%nel
                    if (occ_list(iel) /= i) &
                        hmatel = hmatel &
                            + get_two_body_int_mol_real(coulomb_ints, i, occ_list(iel), a, occ_list(iel), &
                                                    sys) &
                            - get_two_body_int_mol_real(coulomb_ints, i, occ_list(iel), occ_list(iel), a, &
                                                    sys)
                end do
            end associate

            if (perm) hmatel = -hmatel
        end if

    end function slater_condon1_mol

    pure function slater_condon1_mol_excit(sys, occ_list, i, a, perm) result(hmatel)

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
        !    < D | H | D_i^a >, the Hamiltonian matrix element between a
        !        determinant and a single excitation of it for systems defined
        !        by integrals read in from an FCIDUMP file.

        ! WARNING: This function assumes that the D_i^a is a symmetry allowed
        ! excitation from D (and so the matrix element is *not* zero by
        ! symmetry).  This is less safe that slater_condon1_mol but much faster
        ! as it allows symmetry checking to be skipped in the integral lookups.

        use molecular_integrals, only: get_one_body_int_mol_nonzero, get_two_body_int_mol_nonzero
        use system, only: sys_t
        use hamiltonian_data, only: hmatel_t

        type(hmatel_t) :: hmatel
        type(sys_t), intent(in) :: sys
        integer, intent(in) :: occ_list(sys%nel), i, a
        logical, intent(in) :: perm

        integer :: iel

        ! < D | H | D_i^a > = < i | h(a) | a > + \sum_j < ij || aj >

        associate(basis_fns=>sys%basis%basis_fns, &
                  one_e_ints=>sys%read_in%one_e_h_integrals, &
                  coulomb_ints=>sys%read_in%coulomb_integrals)
            hmatel%r = get_one_body_int_mol_nonzero(one_e_ints, i, a, basis_fns)

            do iel = 1, sys%nel
                if (occ_list(iel) /= i) then
                    hmatel%r = hmatel%r &
                                + get_two_body_int_mol_nonzero(coulomb_ints, i, occ_list(iel), a, occ_list(iel), basis_fns)
                    if (basis_fns(occ_list(iel))%Ms == basis_fns(i)%Ms) &
                        hmatel%r = hmatel%r &
                                    - get_two_body_int_mol_nonzero(coulomb_ints, i, occ_list(iel), occ_list(iel), a, basis_fns)
                end if
            end do
        end associate

        if (perm) hmatel%r = -hmatel%r

    end function slater_condon1_mol_excit

    elemental function slater_condon2_mol(sys, i, j, a, b, perm) result(hmatel)

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
        !    < D | H | D_{ij}^{ab} >, the Hamiltonian matrix element between a
        !        determinant and a single excitation of it for systems defined
        !        by integrals read in from an FCIDUMP file.
        !
        ! Note that we don't actually need to directly refer to |D>.

        use molecular_integrals, only: get_two_body_int_mol_real
        use system, only: sys_t
        use point_group_symmetry, only: cross_product_pg_sym

        real(p) :: hmatel
        type(sys_t), intent(in) :: sys
        integer, intent(in) :: i, j, a, b
        logical, intent(in) :: perm

        ! < D | H | D_{ij}^{ab} > = < ij || ab >

        hmatel = get_two_body_int_mol_real(sys%read_in%coulomb_integrals, i, j, a, b, sys) &
                 - get_two_body_int_mol_real(sys%read_in%coulomb_integrals, i, j, b, a, sys)

        if (perm) hmatel = -hmatel

    end function slater_condon2_mol

    pure function slater_condon2_mol_excit(sys, i, j, a, b, perm) result(hmatel)

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
        !    < D | H | D_{ij}^{ab} >, the Hamiltonian matrix element between a
        !        determinant and a single excitation of it for systems defined
        !        by integrals read in from an FCIDUMP file.
        !
        ! Note that we don't actually need to directly refer to |D>.
        !
        ! WARNING: This function assumes that the D_{ij}^{ab} is a symmetry allowed
        ! excitation from D (and so the matrix element is *not* zero by
        ! symmetry).  This is less safe that slater_condon2_mol but much faster
        ! as it allows symmetry checking to be skipped in the integral lookups.

        use molecular_integrals, only: get_two_body_int_mol_nonzero
        use system, only: sys_t
        use hamiltonian_data, only: hmatel_t

        type(hmatel_t) :: hmatel
        type(sys_t), intent(in) :: sys
        integer, intent(in) :: i, j, a, b
        logical, intent(in) :: perm

        ! < D | H | D_{ij}^{ab} > = < ij || ab >

        hmatel%r = 0.0_p

        if (sys%basis%basis_fns(i)%Ms == sys%basis%basis_fns(a)%Ms) &
            hmatel%r = get_two_body_int_mol_nonzero(sys%read_in%coulomb_integrals, i, j, a, b, sys%basis%basis_fns)
        if (sys%basis%basis_fns(i)%Ms == sys%basis%basis_fns(b)%Ms) &
            hmatel%r = hmatel%r - get_two_body_int_mol_nonzero(sys%read_in%coulomb_integrals, i, j, b, a, sys%basis%basis_fns)

        if (perm) hmatel%r = -hmatel%r

    end function slater_condon2_mol_excit

    pure function hf_hamiltonian_energy_mol(sys, f) result(hf_energy)

        ! Work out expectation value of Hartree-Fock Hamiltonian between
        ! two determinants:
        !    i.e. <f|H^HF|f> = \sum_{i_occ} epsilon^{HF}_i.

        ! In:
        !    sys: system being studied.
        !    occ_list: list of occupied orbitals in determinant.
        ! Returns:
        !    H0_energy: sum of hf eigenvalues

        use system, only: sys_t
        use determinants, only: decode_det, sum_sp_eigenvalues_occ_list

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f(sys%basis%tot_string_len)

        integer :: occ_list(sys%nel)
        real(p) :: hf_energy

        call decode_det(sys%basis, f, occ_list)

        hf_energy = sum_sp_eigenvalues_occ_list(sys, occ_list) + sys%read_in%Ecore

    end function hf_hamiltonian_energy_mol

    pure function double_counting_correction_mol(sys, occ_list) result(double_count)

        ! Work out the double counting correction between the Hartree-Fock
        ! energy and the sum of the Hartree-Fock orbitals.

        ! In:
        !    sys: system being studied.
        !    occ_list: list of occupied orbitals in determinant.
        ! Returns:
        !    double_count: double counting correction,
        !    i.e. <D|H|D> - \sum_{i_occ} epsilon_i.

        use system, only: sys_t
        use determinants, only: sum_sp_eigenvalues_occ_list

        type(sys_t), intent(in) :: sys
        integer, intent(in) :: occ_list(:)

        real(p) :: sc0, sum_hf, double_count

        sc0 = slater_condon0_mol_orb_list(sys, occ_list)
        sum_hf = sum_sp_eigenvalues_occ_list(sys, occ_list) + sys%read_in%Ecore

        double_count = sc0 - sum_hf

    end function double_counting_correction_mol

    pure function get_two_e_int_mol(sys, i, j, a, b) result(intgrl)

        ! In:
        !   sys: system being studied
        !   i,j,a,b: spin orbital indices
        ! Returns:
        !   The antisymmetrised integral < ij || ab >

        use system, only: sys_t
        use molecular_integrals, only: get_two_body_int_mol

        type(sys_t), intent(in) :: sys
        integer, intent(in) :: i, j, a, b
        real(p) :: intgrl

        intgrl = get_two_body_int_mol(sys%read_in%coulomb_integrals, i, j, a, b, sys) - &
                 get_two_body_int_mol(sys%read_in%coulomb_integrals, i, j, b, a, sys)

    end function get_two_e_int_mol

    pure function get_one_e_int_mol(sys, i, j) result(intgrl)

        ! In:
        !    sys: system being studied.
        !    i: index of a real-space basis function.
        !    j: index of a real-space basis function.
        ! Returns:
        !    h_ij

        use system, only: sys_t
        use molecular_integrals, only: get_one_body_int_mol

        type(sys_t), intent(in) :: sys
        integer, intent(in) :: i, j
        real(p) :: intgrl

        intgrl = get_one_body_int_mol(sys%read_in%one_e_h_integrals, i, j, sys)

    end function get_one_e_int_mol

end module hamiltonian_molecular

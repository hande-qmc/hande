module hamiltonian_molecular_complex

! Module for evaluating Hamiltonian matrix elements for molecular systems (ie
! systems where the integrals have been read in from file).

use const, only: p, i0

implicit none

private
public :: get_hmatel_mol_comp, slater_condon1_mol_excit_complex, slater_condon2_mol_excit_complex
public :: slater_condon0_mol_complex
contains

    pure function get_hmatel_mol_comp(sys, f1, f2) result(hmatel)

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

        complex(p) :: hmatel
        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f1(sys%basis%string_len), f2(sys%basis%string_len)

        type(excit_t) :: excitation
        integer :: occ_list(sys%nel)

        hmatel = cmplx(0.0, 0.0, p)

        ! Test to see if matrix element is non-zero.
        excitation = get_excitation(sys%nel, sys%basis, f1, f2)

        ! Connected determinants can differ by (at most) 2 spin orbitals.
        if (excitation%nexcit <= 2) then

            select case(excitation%nexcit)
            ! Apply Slater--Condon rules.

            case(0)

                ! < D | H | D > = \sum_i < i | h(i) | i > + \sum_i \sum_{j>i} < ij || ij >
                hmatel = cmplx(slater_condon0_mol_complex(sys, f1), 0.0_p, p)

            case(1)

                ! < D | H | D_i^a > = < i | h(a) | a > + \sum_j < ij || aj >
                call decode_det(sys%basis, f1, occ_list)
                hmatel = slater_condon1_mol_complex(sys, occ_list, excitation%from_orb(1), &
                                            excitation%to_orb(1), excitation%perm)

            case(2)

                ! < D | H | D_{ij}^{ab} > = < ij || ab >

                ! Two electron operator
                hmatel = slater_condon2_mol_complex(sys, excitation%from_orb(1), excitation%from_orb(2), &
                                            & excitation%to_orb(1), excitation%to_orb(2), excitation%perm)
            case default
                ! If f1 & f2 differ by more than two spin orbitals integral value is zero.

            end select

        end if

    end function get_hmatel_mol_comp

    pure function slater_condon0_mol_complex(sys, f) result(hmatel)

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
        complex(p) :: hmatel_comp
        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f(sys%basis%string_len)

        integer :: occ_list(sys%nel)

        ! < D | H | D > = Ecore + \sum_i < i | h(i) | i > + \sum_i \sum_{j>i} < ij || ij >
        call decode_det(sys%basis, f, occ_list)

        hmatel_comp = slater_condon0_mol_orb_list_complex(sys, occ_list)

        ! As Hamiltonian should be Hermitian < D | H | D > must be real, so for easier compatability..
        hmatel = real(hmatel_comp, p)
        ! If wanted to check this could try to write test; as inside pure procedure more
        ! trouble than worth. Further optimisation could remove accretion of imaginary
        ! component in slater_condon_mol_orb_list_complex.

    end function slater_condon0_mol_complex

    pure function slater_condon0_mol_orb_list_complex(sys, occ_list) result(hmatel)

        ! In:
        !    sys: system to be studied.
        !    occ_list: List of occpied orbitals in Slater determinant |D>.
        ! Returns:
        !    <D|H|D>, the diagonal Hamiltonian matrix element involving D for
        !    systems defined by integrals read in from an FCIDUMP file.

        use molecular_integrals, only: get_one_body_int_mol_nonzero, get_two_body_int_mol_nonzero
        use system, only: sys_t

        complex(p) :: hmatel
        type(sys_t), intent(in) :: sys
        integer, intent(in) :: occ_list(:)
        real(p) :: re, im

        integer :: iel, jel, i, j

        associate(one_e_ints=>sys%read_in%one_e_h_integrals, coulomb_ints=>sys%read_in%coulomb_integrals, &
                    one_e_ints_im=>sys%read_in%one_e_h_integrals_imag, coulomb_ints_im=>sys%read_in%coulomb_integrals_imag)
            hmatel = sys%read_in%Ecore
            do iel = 1, sys%nel
                i = occ_list(iel)
                re = get_one_body_int_mol_nonzero(one_e_ints, i, i, sys%basis%basis_fns)
                im = get_one_body_int_mol_nonzero(one_e_ints_im, i, i, sys%basis%basis_fns)
                hmatel = hmatel + cmplx(re, im, p)
                do jel = iel+1, sys%nel
                    j = occ_list(jel)
                    re =  get_two_body_int_mol_nonzero(coulomb_ints, i, j, i, j, sys%basis%basis_fns)
                    im =  get_two_body_int_mol_nonzero(coulomb_ints_im, i, j, i, j, sys%basis%basis_fns)
                    hmatel = hmatel + cmplx(re, im, p)
                    if (sys%basis%basis_fns(i)%Ms == sys%basis%basis_fns(j)%Ms) then
                        re =  get_two_body_int_mol_nonzero(coulomb_ints, i, j, j, i, sys%basis%basis_fns)
                        im =  get_two_body_int_mol_nonzero(coulomb_ints_im, i, j, j, i, sys%basis%basis_fns)
                        hmatel = hmatel - cmplx(re, im, p)
                    end if
                end do
            end do
        end associate

    end function slater_condon0_mol_orb_list_complex


    pure function slater_condon1_mol_complex(sys, occ_list, i, a, perm) result(hmatel)

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

        use molecular_integrals, only: get_one_body_int_mol, get_two_body_int_mol
        use system, only: sys_t

        complex(p) :: hmatel
        type(sys_t), intent(in) :: sys
        integer, intent(in) :: occ_list(sys%nel), i, a
        logical, intent(in) :: perm

        integer :: iel

        ! Check that this excitation is symmetry allowed.
        if (sys%basis%basis_fns(i)%sym/=sys%basis%basis_fns(a)%sym) then
            hmatel = cmplx(0.0_p, 0.0_p, p)
        else
            ! < D | H | D_i^a > = < i | h(a) | a > + \sum_j < ij || aj >

            associate(one_e_ints=>sys%read_in%one_e_h_integrals, coulomb_ints=>sys%read_in%coulomb_integrals,&
                     one_e_ints_im=>sys%read_in%one_e_h_integrals_imag,&
                     coulomb_ints_im=>sys%read_in%coulomb_integrals_imag)

                hmatel = get_one_body_int_mol(one_e_ints, i, a, sys%basis%basis_fns, sys%read_in%pg_sym, one_e_ints_im)

                do iel = 1, sys%nel
                    if (occ_list(iel) /= i) &
                        hmatel = hmatel &
                            + get_two_body_int_mol(coulomb_ints, i, occ_list(iel), a, occ_list(iel), &
                                                    sys%basis%basis_fns, sys%read_in%pg_sym, coulomb_ints_im) &
                            - get_two_body_int_mol(coulomb_ints, i, occ_list(iel), occ_list(iel), a, &
                                                    sys%basis%basis_fns, sys%read_in%pg_sym, coulomb_ints_im)
                end do
            end associate

            if (perm) hmatel = -hmatel
        end if

    end function slater_condon1_mol_complex

    pure function slater_condon1_mol_excit_complex(sys, occ_list, i, a, perm) result(hmatel)

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

        complex(p) :: hmatel
        type(sys_t), intent(in) :: sys
        integer, intent(in) :: occ_list(sys%nel), i, a
        logical, intent(in) :: perm
        real(p) :: re, im

        integer :: iel

        ! < D | H | D_i^a > = < i | h(a) | a > + \sum_j < ij || aj >

        associate(basis_fns=>sys%basis%basis_fns, &
                  one_e_ints=>sys%read_in%one_e_h_integrals, &
                  coulomb_ints=>sys%read_in%coulomb_integrals,&
                  one_e_ints_im=>sys%read_in%one_e_h_integrals_imag,&
                  coulomb_ints_im=>sys%read_in%coulomb_integrals_imag)

            re = get_one_body_int_mol_nonzero(one_e_ints, i, a, basis_fns)
            im = get_one_body_int_mol_nonzero(one_e_ints_im, i, a, basis_fns)
            hmatel = cmplx(re, im, p)
            do iel = 1, sys%nel
                if (occ_list(iel) /= i) then
                    re = get_two_body_int_mol_nonzero(coulomb_ints, i, occ_list(iel), a, occ_list(iel), basis_fns)
                    im = get_two_body_int_mol_nonzero(coulomb_ints_im, i, occ_list(iel), a, occ_list(iel), basis_fns)
                    hmatel = hmatel + cmplx(re, im, p)
                    if (basis_fns(occ_list(iel))%Ms == basis_fns(i)%Ms) then
                        re = get_two_body_int_mol_nonzero(coulomb_ints, i, occ_list(iel), occ_list(iel), a, basis_fns)
                        im = get_two_body_int_mol_nonzero(coulomb_ints_im, i, occ_list(iel), occ_list(iel), a, basis_fns)
                        hmatel = hmatel - cmplx(re, im, p)
                    end if
                end if
            end do
        end associate

        if (perm) hmatel = -hmatel

    end function slater_condon1_mol_excit_complex

    pure function slater_condon2_mol_complex(sys, i, j, a, b, perm) result(hmatel)

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

        use molecular_integrals, only: get_two_body_int_mol
        use system, only: sys_t
        use point_group_symmetry, only: cross_product_pg_sym

        complex(p) :: hmatel
        type(sys_t), intent(in) :: sys
        integer, intent(in) :: i, j, a, b
        logical, intent(in) :: perm

        ! < D | H | D_{ij}^{ab} > = < ij || ab >

        hmatel = get_two_body_int_mol(sys%read_in%coulomb_integrals, i, j, a, b, sys%basis%basis_fns, sys%read_in%pg_sym, &
                                        sys%read_in%coulomb_integrals_imag) &
                 - get_two_body_int_mol(sys%read_in%coulomb_integrals, i, j, b, a, sys%basis%basis_fns, sys%read_in%pg_sym, &
                                        sys%read_in%coulomb_integrals_imag)

        if (perm) hmatel = -hmatel

    end function slater_condon2_mol_complex

    elemental function slater_condon2_mol_excit_complex(sys, i, j, a, b, perm) result(hmatel)

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

        complex(p) :: hmatel
        real(p) :: re, im
        type(sys_t), intent(in) :: sys
        integer, intent(in) :: i, j, a, b
        logical, intent(in) :: perm

        ! < D | H | D_{ij}^{ab} > = < ij || ab >

        hmatel = cmplx(0.0, 0.0)

        if (sys%basis%basis_fns(i)%Ms == sys%basis%basis_fns(a)%Ms) then
            re = get_two_body_int_mol_nonzero(sys%read_in%coulomb_integrals, i, j, a, b, sys%basis%basis_fns)
            im = get_two_body_int_mol_nonzero(sys%read_in%coulomb_integrals_imag, i, j, a, b, sys%basis%basis_fns)
            hmatel = cmplx(re, im, p)
        end if
        if (sys%basis%basis_fns(i)%Ms == sys%basis%basis_fns(b)%Ms) then
            re =  get_two_body_int_mol_nonzero(sys%read_in%coulomb_integrals, i, j, b, a, sys%basis%basis_fns)
            im =  get_two_body_int_mol_nonzero(sys%read_in%coulomb_integrals_imag, i, j, b, a, sys%basis%basis_fns)
            hmatel = hmatel - cmplx(re, im, p)
        end if
        if (perm) hmatel = -hmatel

    end function slater_condon2_mol_excit_complex

end module hamiltonian_molecular_complex

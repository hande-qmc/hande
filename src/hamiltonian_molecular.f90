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

        use molecular_integrals, only: get_one_body_int_mol_nonzero, get_two_body_int_mol_nonzero, &
                                        get_two_body_exchange_pbc_int_nonzero
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
! [todo] - AJWT: replace with a pointer initialized at the start of the function to speed this up?
                    if (sys%basis%basis_fns(i)%Ms == sys%basis%basis_fns(j)%Ms) then
                        if (sys%read_in%extra_exchange_integrals) then
                            hmatel = hmatel - get_two_body_exchange_pbc_int_nonzero(sys%read_in%additional_exchange_ints, &
                                                    i, j, j, i, sys%basis%basis_fns)
                        else
                            hmatel = hmatel - get_two_body_int_mol_nonzero(coulomb_ints, i, j, j, i, sys%basis%basis_fns)
                        end if
                    end if
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

        use molecular_integrals, only: get_one_body_int_mol_real, get_two_body_int_mol_real, &
                                       get_two_body_exchange_pbc_int_real
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
                    if (occ_list(iel) /= i) then
                        hmatel = hmatel &
                            + get_two_body_int_mol_real(coulomb_ints, i, occ_list(iel), a, occ_list(iel), &
                                                    sys)
! [todo] - AJWT: replace with a pointer initialized at the start of the function to speed this up?
                        if (sys%read_in%extra_exchange_integrals) then
                            hmatel = hmatel - get_two_body_exchange_pbc_int_real(sys%read_in%additional_exchange_ints, &
                                                    i, occ_list(iel), occ_list(iel), a, sys)
                        else
                            hmatel = hmatel - get_two_body_int_mol_real(coulomb_ints, i, occ_list(iel), &
                                                occ_list(iel), a, sys)
                        end if
                    end if
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

        use molecular_integrals, only: get_one_body_int_mol_nonzero, get_two_body_int_mol_nonzero, &
                                       get_two_body_exchange_pbc_int_nonzero
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
                    if (basis_fns(occ_list(iel))%Ms == basis_fns(i)%Ms) then
! [todo] - AJWT: replace with a pointer initialized at the start of the function to speed this up?
                        if (sys%read_in%extra_exchange_integrals) then
                            hmatel%r = hmatel%r - get_two_body_exchange_pbc_int_nonzero(sys%read_in%additional_exchange_ints, &
                                                    i, occ_list(iel), occ_list(iel), a, basis_fns)
                        else
                            hmatel%r = hmatel%r &
                                    - get_two_body_int_mol_nonzero(coulomb_ints, i, occ_list(iel), occ_list(iel), a, basis_fns)
                        end if
                    end if
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

    pure subroutine create_weighted_excitation_list_mol(sys, i, b, a_list, a_list_len, weights, weight_tot)

        ! Generate a list of allowed excitations from i to one a of a_list with their weights based on
        ! sqrt(|<ia|ai>|) (Power-Pitzer) or sqrt(|<ia|ia>|) (Cauchy Schwarz)
        !
        ! In:
        !    sys:   The system in which the orbitals live
        !    i:  integer specifying the from orbital
        !    b:  integer specifying the other to orbital if found aready. If not found already, a 0 was passed.
        !    a_list:   list of integers specifying the basis functions we're allowed to excite to
        !    a_list_len:   The length of a_list
        ! Out:
        !    weights:   A list of reals (length a_list_len), with the weight of each of the to_list orbitals
        !    weight_tot: The sum of all the weights.

        use system, only: sys_t
        use proc_pointers, only: get_two_body_int_cou_ex_mol_real_ptr

        type(sys_t), intent(in) :: sys
        integer, intent(in) :: i, b, a_list_len, a_list(:)
        real(p), intent(out) :: weight_tot, weights(:)

        integer :: k
        real(p) :: weight

        weight_tot = 0
        do k = 1, a_list_len
            ! check whether a and b are identical (if one of them has already been found). 
            ! If one of these checks is true, assign zero weight.
            if (a_list(k) /= b) then
                ! [review] - JSS: could avoid the abs with an assumption or a one-off O(N2) check during initialisation?
                ! [review] - VAN: where in initialisation would that check be?
                ! If exchange integral, then should be +ve, but best abs below in case!
                weight = get_two_body_int_cou_ex_mol_real_ptr(sys, i, a_list(k))
                weights(k) = sqrt(abs(weight))
                weight_tot = weight_tot + weights(k)
            else
                ! a == b which is forbidden, assign zero weight.
                weights(k) = 0.0_p
            end if
        end do

    end subroutine create_weighted_excitation_list_mol
    
    pure function get_two_body_int_ex_mol_real(sys, i, a) result(integral)
        ! Get integral <ia|ai> (exchange) for real system.
        !
        ! In:
        !   sys: system object
        !   i: orbital to excite from
        !   a: orbital to excite to
        ! Out:
        !   integral: resulting exchange integral
        
        use system, only: sys_t
        use molecular_integrals, only: get_two_body_int_mol_real

        type(sys_t), intent(in) :: sys
        integer, intent(in) :: i, a
        real(p) :: integral

        integral = get_two_body_int_mol_real(sys%read_in%coulomb_integrals, i, a, a, i, sys)

    end function get_two_body_int_ex_mol_real

    pure function get_two_body_int_cou_mol_real(sys, i, a) result(integral)
        
        ! Get integral <ia|ia> (Coulomb) for real system.
        !
        ! In:
        !   sys: system object
        !   i: orbital to excite from
        !   a: orbital to excite to
        ! Out:
        !   integral: resulting Coulomb integral
        
        use system, only: sys_t
        use molecular_integrals, only: get_two_body_int_mol_real

        type(sys_t), intent(in) :: sys
        integer, intent(in) :: i, a
        real(p) :: integral

        integral = get_two_body_int_mol_real(sys%read_in%coulomb_integrals, i, a, i, a, sys)

    end function get_two_body_int_cou_mol_real

    pure function abs_hmatel_mol(hmatel) result(abs_hmatel)
        ! Return absolute value of hmatel. As we do not have a complex function here, hmatel is real.
        use hamiltonian_data, only: hmatel_t
        type(hmatel_t), intent(in) :: hmatel
        real(p) :: abs_hmatel

        abs_hmatel = abs(hmatel%r)
    end function abs_hmatel_mol

    pure function single_excitation_weight_mol(sys, ref, i, a) result(weight)
        ! Note that this function is very similar to single_excitation_weight_periodic in hamiltonian_periodic_complex.
        
        ! In:
        !    sys: system to be studied.
        !    ref: reference determinant information
        !    i: the spin-orbital from which an electron is excited in
        !       the reference determinant.
        !    a: the spin-orbital into which an electron is excited in
        !       the excited determinant.
        ! Returns:
        !    weight: used in the Power Pitzer Order N excitation generator as
        !            pre-calculated weight for the single excitations i -> a.

        ! [todo] - get_two_body_int_mol_complex imported here as single_excitation_weight_periodic_complex uses that and
        ! [todo] - they share a proc pointer.
        use molecular_integrals, only: get_two_body_int_mol_real, get_two_body_int_mol_complex
        use system, only: sys_t
        use reference_determinant, only: reference_t

        type(sys_t), intent(in) :: sys
        integer, intent(in) :: i, a
        type(reference_t), intent(in) :: ref
        real(p) :: weight
        integer :: n_jb, pos, j, b, virt_pos, occ_pos, nvirt, virt_list(sys%basis%nbasis - sys%nel)

        ! Note that we do not check whether Brillouin's theorem is obeyed or not.
        
        ! Assuming that we often have the RHF reference which obeys Brillouin's theorem
        ! (no single excitations from the reference), the main contribution to single
        ! excitations can be the single excitation from a double excitation a step
        ! back towards the reference (D_{ij}^{ab} -> D_j^b). For a given single excitation
        ! i -> a, we sum over j and b and look at the weights of D_j^b -> D_{ij}^{ab} to
        ! determine that weight of i -> a. 
        ! < D_j^b | H | D_{ij}^{ab} > = h1_i^a + \sum_k{occ. in ref.} (<ik|ak> - <ik|ka>)
        ! - < ij|aj> + <ij|ja> + <ib|ab> - <ib|ba>
        ! where h1_i^a + \sum_k{occ. in ref.} (<ik|ak> - <ik|ka>), called ref_term here, is zero in the
        ! case of an SCF reference (todo: check).
    
        associate(basis_fns=>sys%basis%basis_fns, &
                  one_e_ints=>sys%read_in%one_e_h_integrals, &
                  coulomb_ints=>sys%read_in%coulomb_integrals)
            virt_pos = 1
            occ_pos = 1
            nvirt = sys%basis%nbasis - sys%nel
            do pos = 1, sys%basis%nbasis
                if (occ_pos > sys%nel) then
                    virt_list(virt_pos) = pos
                    virt_pos = virt_pos + 1
                else
                    if (ref%occ_list0(occ_pos) == pos) then
                        occ_pos = occ_pos + 1
                    else
                        virt_list(virt_pos) = pos
                        virt_pos = virt_pos + 1
                    end if
                end if
            end do

            n_jb = 0
            weight = 0.0_p
            ! Let j be an occupied spinorbital in reference (consistent with equation above).
            ! [todo] - should j and b be more general? - might need to worry about double counting then
            ! take mean of abs values. [todo] - what if some j b combinations are not valid? count them too?
            ! Is there a chance of a bias? Can the weight ever be zero?
            do j = 1, sys%nel
                do b = 1, nvirt
                    n_jb = n_jb + 1
                    weight = weight + &
                        abs(get_two_body_int_mol_real(coulomb_ints, i, ref%occ_list0(j), ref%occ_list0(j), a, sys) - &
                        get_two_body_int_mol_real(coulomb_ints, i, ref%occ_list0(j), a, ref%occ_list0(j), sys) + &
                        get_two_body_int_mol_real(coulomb_ints, i, virt_list(b), a, virt_list(b), sys) - &
                        get_two_body_int_mol_real(coulomb_ints, i, virt_list(b), virt_list(b), a, sys))
                end do
            end do
            weight = weight/real(n_jb,p)

        end associate
    
    end function single_excitation_weight_mol

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

module hamiltonian_periodic_complex

! Module for evaluating Hamiltonian matrix elements for periodic
! systems where the integrals have been read in from file.

use const, only: p, i0, depsilon

implicit none

private
public :: get_hmatel_periodic_complex, slater_condon1_periodic_excit_complex
public :: slater_condon1_periodic_complex
public :: slater_condon2_periodic_excit_complex, slater_condon0_periodic_complex
public :: create_weighted_excitation_list_periodic_complex, abs_hmatel_periodic_complex
public :: single_excitation_weight_periodic
contains

    pure function get_hmatel_periodic_complex(sys, f1, f2) result(hmatel)

        ! In:
        !    sys: system to be studied.
        !    f1, f2: bit string representation of the Slater
        !        determinants D1 and D2 respectively.
        ! Returns:
        !    Hamiltonian matrix element between the two determinants,
        !    < D1 | H | D2 >, where the determinants are formed from
        !    orbitals read in from an FCIDUMP file.
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

        hmatel%c = cmplx(0.0, 0.0, p)

        ! Test to see if matrix element is non-zero.
        excitation = get_excitation(sys%nel, sys%basis, f1, f2)

        ! Connected determinants can differ by (at most) 2 spin orbitals.
        if (excitation%nexcit <= 2) then

            select case(excitation%nexcit)
            ! Apply Slater--Condon rules.

            case(0)

                ! < D | H | D > = \sum_i < i | h(i) | i > + \sum_i \sum_{j>i} < ij || ij >
                hmatel%c = cmplx(slater_condon0_periodic_complex(sys, f1), 0.0_p, p)

            case(1)

                ! < D | H | D_i^a > = < i | h(a) | a > + \sum_j < ij || aj >
                call decode_det(sys%basis, f1, occ_list)
                hmatel%c = slater_condon1_periodic_complex(sys, occ_list, excitation%from_orb(1), &
                                            excitation%to_orb(1), excitation%perm)

            case(2)

                ! < D | H | D_{ij}^{ab} > = < ij || ab >

                ! Two electron operator
                hmatel%c = slater_condon2_periodic_complex(sys, excitation%from_orb(1), excitation%from_orb(2), &
                                            & excitation%to_orb(1), excitation%to_orb(2), excitation%perm)
            case default
                ! If f1 & f2 differ by more than two spin orbitals integral value is zero.

            end select

        end if

    end function get_hmatel_periodic_complex

    pure function slater_condon0_periodic_complex(sys, f) result(hmatel)

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
        integer(i0), intent(in) :: f(sys%basis%tot_string_len)

        integer :: occ_list(sys%nel)

        ! < D | H | D > = Ecore + \sum_i < i | h(i) | i > + \sum_i \sum_{j>i} < ij || ij >
        call decode_det(sys%basis, f, occ_list)

        hmatel_comp = slater_condon0_periodic_orb_list_complex(sys, occ_list)

        ! As Hamiltonian should be Hermitian < D | H | D > must be real, so for easier compatability..
        hmatel = real(hmatel_comp, p)
        ! If wanted to check this could try to write test; as inside pure procedure more
        ! trouble than worth. Further optimisation could remove accretion of imaginary
        ! component in slater_condon_periodic_orb_list_complex.

    end function slater_condon0_periodic_complex

    pure function slater_condon0_periodic_orb_list_complex(sys, occ_list) result(hmatel)

        ! In:
        !    sys: system to be studied.
        !    occ_list: List of occpied orbitals in Slater determinant |D>.
        ! Returns:
        !    <D|H|D>, the diagonal Hamiltonian matrix element involving D for
        !    systems defined by integrals read in from an FCIDUMP file.

        use molecular_integrals, only: get_one_body_int_mol_nonzero, get_two_body_int_mol_nonzero, &
                                        get_two_body_exchange_pbc_int_nonzero
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
                        if (sys%read_in%extra_exchange_integrals) then
                            re = get_two_body_exchange_pbc_int_nonzero(sys%read_in%additional_exchange_ints, &
                                        i, j, j, i, sys%basis%basis_fns)
                            im = get_two_body_exchange_pbc_int_nonzero(sys%read_in%additional_exchange_ints_imag, &
                                        i, j, j, i, sys%basis%basis_fns)
                        else
                            re =  get_two_body_int_mol_nonzero(coulomb_ints, i, j, j, i, sys%basis%basis_fns)
                            im =  get_two_body_int_mol_nonzero(coulomb_ints_im, i, j, j, i, sys%basis%basis_fns)
                        end if
                        hmatel = hmatel - cmplx(re, im, p)
                    end if
                end do
            end do
        end associate

    end function slater_condon0_periodic_orb_list_complex

    pure function slater_condon1_periodic_complex(sys, occ_list, i, a, perm) result(hmatel)

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

        use molecular_integrals, only: get_one_body_int_mol, get_two_body_int_mol, &
                                        get_two_body_exchange_pbc_int_complex
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

                hmatel = get_one_body_int_mol(one_e_ints, one_e_ints_im, i, a, sys)

                do iel = 1, sys%nel
                    if (occ_list(iel) /= i) then
                        hmatel = hmatel &
                            + get_two_body_int_mol(coulomb_ints, coulomb_ints_im, i, occ_list(iel), a, occ_list(iel), &
                                                    sys)
                        if (sys%read_in%extra_exchange_integrals) then
                            hmatel = hmatel &
                                        - get_two_body_exchange_pbc_int_complex(sys%read_in%additional_exchange_ints, &
                                        sys%read_in%additional_exchange_ints_imag, i, occ_list(iel), occ_list(iel), a, sys)

                        else
                            hmatel = hmatel &
                                - get_two_body_int_mol(coulomb_ints, coulomb_ints_im, i, occ_list(iel), occ_list(iel), a, &
                                                    sys)
                        end if
                    end if
                end do
            end associate

            if (perm) hmatel = -hmatel
        end if

    end function slater_condon1_periodic_complex

    pure function slater_condon1_periodic_excit_complex(sys, occ_list, i, a, perm) result(hmatel)

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
            hmatel%c = cmplx(re, im, p)
            do iel = 1, sys%nel
                if (occ_list(iel) /= i) then
                    re = get_two_body_int_mol_nonzero(coulomb_ints, i, occ_list(iel), a, occ_list(iel), basis_fns)
                    im = get_two_body_int_mol_nonzero(coulomb_ints_im, i, occ_list(iel), a, occ_list(iel), basis_fns)
                    hmatel%c = hmatel%c + cmplx(re, im, p)
                    if (basis_fns(occ_list(iel))%Ms == basis_fns(i)%Ms) then
                        if (sys%read_in%extra_exchange_integrals) then
                            re = get_two_body_exchange_pbc_int_nonzero(sys%read_in%additional_exchange_ints, &
                                        i, occ_list(iel), occ_list(iel), a, sys%basis%basis_fns)
                            im = get_two_body_exchange_pbc_int_nonzero(sys%read_in%additional_exchange_ints_imag, &
                                        i, occ_list(iel), occ_list(iel), a, sys%basis%basis_fns)
                        else
                            re =  get_two_body_int_mol_nonzero(coulomb_ints, i, occ_list(iel), occ_list(iel), a, &
                                                                sys%basis%basis_fns)
                            im =  get_two_body_int_mol_nonzero(coulomb_ints_im, i, occ_list(iel), occ_list(iel), &
                                                                a, sys%basis%basis_fns)
                        end if
                        hmatel%c = hmatel%c - cmplx(re, im, p)
                    end if
                end if
            end do
        end associate

        if (perm) hmatel%c = -hmatel%c

    end function slater_condon1_periodic_excit_complex

    pure function slater_condon2_periodic_complex(sys, i, j, a, b, perm) result(hmatel)

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

        hmatel = get_two_body_int_mol(sys%read_in%coulomb_integrals, sys%read_in%coulomb_integrals_imag, &
                                        i, j, a, b, sys) &
                 - get_two_body_int_mol(sys%read_in%coulomb_integrals, sys%read_in%coulomb_integrals_imag, &
                                        i, j, b, a, sys)

        if (perm) hmatel = -hmatel

    end function slater_condon2_periodic_complex

    pure function slater_condon2_periodic_excit_complex(sys, i, j, a, b, perm) result(hmatel)

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
        ! symmetry).  This is less safe that slater_condon2_periodic but much faster
        ! as it allows symmetry checking to be skipped in the integral lookups.

        use molecular_integrals, only: get_two_body_int_mol_nonzero
        use system, only: sys_t
        use hamiltonian_data, only: hmatel_t

        type(hmatel_t) :: hmatel
        real(p) :: re, im
        type(sys_t), intent(in) :: sys
        integer, intent(in) :: i, j, a, b
        logical, intent(in) :: perm

        ! < D | H | D_{ij}^{ab} > = < ij || ab >

        hmatel%c = cmplx(0.0, 0.0)

        if (sys%basis%basis_fns(i)%Ms == sys%basis%basis_fns(a)%Ms) then
            re = get_two_body_int_mol_nonzero(sys%read_in%coulomb_integrals, i, j, a, b, sys%basis%basis_fns)
            im = get_two_body_int_mol_nonzero(sys%read_in%coulomb_integrals_imag, i, j, a, b, sys%basis%basis_fns)
            hmatel%c = cmplx(re, im, p)
        end if
        if (sys%basis%basis_fns(i)%Ms == sys%basis%basis_fns(b)%Ms) then
            re =  get_two_body_int_mol_nonzero(sys%read_in%coulomb_integrals, i, j, b, a, sys%basis%basis_fns)
            im =  get_two_body_int_mol_nonzero(sys%read_in%coulomb_integrals_imag, i, j, b, a, sys%basis%basis_fns)
            hmatel%c = hmatel%c - cmplx(re, im, p)
        end if
        if (perm) hmatel%c = -hmatel%c

    end function slater_condon2_periodic_excit_complex

    pure subroutine create_weighted_excitation_list_periodic_complex(sys, i, b, a_list, a_list_len, weights, weight_tot)

        ! Generate a list of allowed excitations from i to one a of a_list with their weights based on
        ! sqrt(|<ia|ai>|)
        !
        ! In:
        !    sys:   The system in which the orbitals live
        !    i:  integer specifying the from orbital
        !    b:  integer specifying the other to orbital if found aready. If not found already, a 0 is passed.
        !    a_list:   list of integers specifying the basis functions we're allowed to excite to
        !    a_list_len:   The length of a_list
        ! Out:
        !    weights:   A list of reals (length a_list_len), with the weight of each of the to_list orbit
        !    weight_tot: The sum of all the weights.

        use system, only: sys_t
        use molecular_integrals, only: get_two_body_int_mol_complex

        type(sys_t), intent(in) :: sys
        integer, intent(in) :: i, b, a_list_len, a_list(:)
        real(p), intent(out) :: weight_tot, weights(:)

        integer :: k
        complex(p) :: weight

        weight_tot = 0
        do k = 1, a_list_len
            ! check whether a and b are identical (if one of them has already been found). 
            ! If one of these checks is true, assign zero weight.
            if (a_list(k) /= b) then
                ! [review] - JSS: could avoid the abs with an assumption or a one-off O(N2) check during init
                ! [review] - VAN: need abs here for complex
                weight = get_two_body_int_mol_complex(sys%read_in%coulomb_integrals, sys%read_in%coulomb_integrals_imag, i, &
                                                      a_list(k), a_list(k), i, sys)
                weights(k) = sqrt(abs(weight))
                weight_tot = weight_tot + weights(k)
            else
                ! a == b which is forbidden, assign zero weight.
                weights(k) = 0.0_p
            end if
        end do

    end subroutine create_weighted_excitation_list_periodic_complex

    pure function abs_hmatel_periodic_complex(hmatel) result(abs_hmatel)
        ! Return absolute value of hmatel. As we have a complex function here, hmatel is complex.
        ! This version is more numerically stable than some versions in fortran libraries.      
        ! https://math.stackexchange.com/questions/750648/absolute-value-of-complex-number-in-numerical-recipes
        
        ! Following Cabs() for C as described on page 949 "Numerical Recipes in C - the art of scientific computing.
        ! Second Edition" by W. H. Press, S. A. Teukolsky, W. T. Vetterling, B. P. Flannery.
        ! Cambridge University Press, 1992.
        use hamiltonian_data, only: hmatel_t
        type(hmatel_t), intent(in) :: hmatel
        real(p) :: abs_hmatel
        real(p) :: a, b

        a = real(hmatel%c)
        b = aimag(hmatel%c)

        if (abs(a) < depsilon) then
            abs_hmatel = abs(b)
        else if (abs(b) < depsilon) then
            abs_hmatel = abs(a)
        else if (abs(a) > abs(b)) then
            abs_hmatel = abs(a) * sqrt(1.0_p + (b/a)**2)
        else
            abs_hmatel = abs(b) * sqrt(1.0_p + (a/b)**2)
        end if

    end function abs_hmatel_periodic_complex

    pure function single_excitation_weight_periodic(sys, ref, i, a) result(weight)
        ! Note that this function is very similar to single_excitation_weight_mol in hamiltonian_molecular.

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

        ! [todo] - get_two_body_int_mol_real is imported here as single_excitation_weight_mol uses that and
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
                  one_e_ints=>sys%read_in%one_e_h_integrals, & ![todo] delete when not needed
                  c_ints=>sys%read_in%coulomb_integrals, &
                  c_ints_im=>sys%read_in%coulomb_integrals_imag)
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
                        abs(get_two_body_int_mol_complex(c_ints, c_ints_im, i, ref%occ_list0(j), ref%occ_list0(j), a, sys) - &
                        get_two_body_int_mol_complex(c_ints, c_ints_im, i, ref%occ_list0(j), a, ref%occ_list0(j), sys) + &
                        get_two_body_int_mol_complex(c_ints, c_ints_im, i, virt_list(b), a, virt_list(b), sys) - &
                        get_two_body_int_mol_complex(c_ints, c_ints_im, i, virt_list(b), virt_list(b), a, sys))
                end do
            end do
            weight = weight/real(n_jb,p)

        end associate

    end function single_excitation_weight_periodic

end module hamiltonian_periodic_complex

module hamiltonian_ringium

! Module for evaluating Hamiltonian matrix elements for ringium (J Chem Phys 138 164124 (2013)).

! Note that ms = 1 for all occupied ringium orbitals.

use const

implicit none

contains

    pure function get_hmatel_ringium(sys, f1, f2) result(hmatel)

        ! In:
        !    sys: system to be studied.
        !    f1, f2: bit string representation of the Slater
        !        determinants D1 and D2 respectively.
        ! Returns:
        !    Hamiltonian matrix element between the two determinants,
        !    < D1 | H | D2 >, where the determinants are formed from
        !    complex exponential functions.

        ! Used in ringium only.

        use excitations, only: excit_t, get_excitation
        use system, only: sys_t

        real(p) :: hmatel
        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f1(sys%basis%string_len), f2(sys%basis%string_len)
        type(excit_t) :: excitation

        ! Test to see if Hamiltonian matrix element is non-zero.

        excitation = get_excitation(sys%nel, sys%basis, f1,f2)

        ! Connected determinants can differ by (at most) 2 spin orbitals.
        ! Ringium has only double excitations due to angular momentum conservation
        select case(excitation%nexcit)
        ! Apply Slater--Condon rules.
        case(0)

            ! < D | H | D > = \sum_i < i | h(i) | i > + \sum_i \sum_{j>i} < ij || ij >
            hmatel = slater_condon0_ringium(sys, f1)

        case(2)

            ! < D | H | D_{ij}^{ab} > = < ij || ab >

            ! Two electron operator
            hmatel = slater_condon2_ringium(sys, excitation%from_orb(1), excitation%from_orb(2), &
                                      & excitation%to_orb(1), excitation%to_orb(2), excitation%perm)
        case default

            hmatel = 0.0_p

        end select

    end function get_hmatel_ringium

    pure function slater_condon0_ringium(sys, f) result(hmatel)

        ! In:
        !    sys: system to be studied.
        !    f: bit string representation of the Slater determinant.
        ! Returns:
        !    < D_i | H | D_i >, the diagonal Hamiltonian matrix elements, for
        !        ringium in the RHF basis

        use determinants, only: decode_det
        use system, only: sys_t
        use ringium_system, only: get_two_e_int_ringium

        real(p) :: hmatel
        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f(sys%basis%string_len)
        integer :: occ_list(sys%nel)

        integer :: i, j

        call decode_det(sys%basis, f, occ_list)

        ! < D | H | D > = \sum_i < i | h(i) | i > + \sum_i \sum_{j>i} < ij || ij >
        hmatel = 0.0_p

        ! One electron operator: kinetic term
        do i = 1, sys%nel
            hmatel = hmatel + sys%basis%basis_fns(occ_list(i))%sp_eigv
        end do

        ! Two electron operator: Coulomb term.
        do i = 1, sys%nel
            do j = i+1, sys%nel
                hmatel = hmatel + get_two_e_int_ringium(sys, occ_list(i), occ_list(j), occ_list(i), occ_list(j))
            end do
        end do

    end function slater_condon0_ringium

    pure function slater_condon2_ringium(sys, i, j, a, b, perm) result(hmatel)

        ! In:
        !    sys: system to be studied.
        !    i,j:  index of the spin-orbital from which an electron is excited in
        !          the reference determinant.
        !    a,b:  index of the spin-orbital into which an electron is excited in
        !          the excited determinant.
        !    perm: true if D and D_i^a are connected by an odd number of
        !          permutations.
        ! Returns:
        !    < D | H | D_ij^ab >, the Hamiltonian matrix element between a
        !    determinant and a double excitation of it in the UEG.

        use ringium_system, only: get_two_e_int_ringium
        use system, only: sys_t

        real(p) :: hmatel
        type(sys_t), intent(in) :: sys
        integer, intent(in) :: i, j, a, b
        logical, intent(in) :: perm

        hmatel = get_two_e_int_ringium(sys, i, j, a, b)

        if (perm) hmatel = -hmatel

    end function slater_condon2_ringium

    end module hamiltonian_ringium

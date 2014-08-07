module hamiltonian_ueg

! Module for evaluating Hamiltonian matrix elements for the uniform electron
! gas.

use const

implicit none

contains

    pure function get_hmatel_ueg(sys, f1, f2) result(hmatel)

        ! In:
        !    sys: system to be studied.
        !    f1, f2: bit string representation of the Slater
        !        determinants D1 and D2 respectively.
        ! Returns:
        !    Hamiltonian matrix element between the two determinants,
        !    < D1 | H | D2 >, where the determinants are formed from
        !    real space basis functions.

        ! Used in the UEG only.

        use excitations, only: excit, get_excitation
        use system, only: sys_t

        real(p) :: hmatel
        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f1(sys%basis%string_len), f2(sys%basis%string_len)
        type(excit) :: excitation

        ! Test to see if Hamiltonian matrix element is non-zero.

        ! Assume D1 and D2 are of the same symmetry.  Namely:

        !     We assume Ms is conserved (ie has already been checked for).

        !     The overall crystal momentum must be conserved (i.e. satisfy
        !     translational symmetry).  We assume this is also already checked.

        excitation = get_excitation(sys%nel, sys%basis, f1,f2)

        ! Connected determinants can differ by (at most) 2 spin orbitals.
        ! UEG (at least in the RHF basis of plane waves) has only double
        ! excitations, as the kinetic operator is diagonal in a plane wave
        ! basis.
        select case(excitation%nexcit)
        ! Apply Slater--Condon rules.
        case(0)

            ! < D | H | D > = \sum_i < i | h(i) | i > + \sum_i \sum_{j>i} < ij || ij >
            hmatel = slater_condon0_ueg(sys, f1)

        case(2)

            ! < D | H | D_{ij}^{ab} > = < ij || ab >

            ! Two electron operator
            hmatel = slater_condon2_ueg(sys, excitation%from_orb(1), excitation%from_orb(2), &
                                      & excitation%to_orb(1), excitation%to_orb(2),excitation%perm)
        case default

            hmatel = 0.0_p

        end select

    end function get_hmatel_ueg

    pure function slater_condon0_ueg(sys, f) result(hmatel)

        ! In:
        !    sys: system to be studied.
        !    f: bit string representation of the Slater determinant.
        ! Returns:
        !    < D_i | H | D_i >, the diagonal Hamiltonian matrix elements, for
        !        the Hubbard model in momentum space.

        use determinants, only: decode_det
        use system, only: sys_t
        use ueg_system, only: exchange_int_ueg

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
                ! Coulomb term is infinite but cancels exactly with the
                ! infinities in the electron-background and
                ! background-background interactions.
                if (mod(occ_list(i),2) == mod(occ_list(j),2)) then
                    ! Have an exchange term
                    hmatel = hmatel - exchange_int_ueg(sys, occ_list(i), occ_list(j))
                end if
            end do
        end do

    end function slater_condon0_ueg

    pure function slater_condon2_ueg(sys, i, j, a, b, perm) result(hmatel)

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

        use ueg_system, only: get_two_e_int_ueg
        use system, only: sys_t

        real(p) :: hmatel
        type(sys_t), intent(in) :: sys
        integer, intent(in) :: i, j, a, b
        logical, intent(in) :: perm

        hmatel = get_two_e_int_ueg(sys, i, j, a, b)

        if (perm) hmatel = -hmatel

    end function slater_condon2_ueg

    pure function slater_condon2_ueg_excit(sys, i, j, a, b, perm) result(hmatel)

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

        ! WARNING: This function assumes that the D_{ij}^{ab} is a symmetry allowed
        ! excitation from D (and so the matrix element is *not* zero by
        ! symmetry).  This is less safe that slater_condon2_ueg but much faster
        ! as it allows symmetry checking to be skipped in the integral
        ! calculation.

        use ueg_system, only: coulomb_int_ueg
        use system, only: sys_t

        real(p) :: hmatel
        type(sys_t), intent(in) :: sys
        integer, intent(in) :: i, j, a, b
        logical, intent(in) :: perm

        hmatel = 0.0_p

        if (sys%basis%basis_fns(i)%Ms == sys%basis%basis_fns(a)%Ms) hmatel = coulomb_int_ueg(sys, i, a)
        if (sys%basis%basis_fns(i)%Ms == sys%basis%basis_fns(b)%Ms) hmatel = hmatel - coulomb_int_ueg(sys, i, b)

        if (perm) hmatel = -hmatel

    end function slater_condon2_ueg_excit

end module hamiltonian_ueg

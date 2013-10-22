module hamiltonian_hub_k

! Module for evaluating Hamiltonian matrix elements for the Hubbard model in
! momentum space (i.e.  using a Bloch orbital basis set).

use const

implicit none

contains

    pure function get_hmatel_hub_k(f1, f2) result(hmatel)

        ! In:
        !    f1, f2: bit string representation of the Slater
        !        determinants D1 and D2 respectively.
        ! Returns:
        !    Hamiltonian matrix element between the two determinants,
        !    < D1 | H | D2 >, where the determinants are formed from
        !    momentum space basis functions.

        ! Used in the momentum space formulation of the Hubbard model only.

        use determinants, only: basis_length
        use excitations, only: excit, get_excitation

        real(p) :: hmatel
        integer(i0), intent(in) :: f1(basis_length), f2(basis_length)
        logical :: non_zero
        type(excit) :: excitation

        hmatel = 0.0_p
        non_zero = .false.

        ! Test to see if Hamiltonian matrix element is non-zero.

        ! Assume D1 and D2 are of the same symmetry.  Namely:

        !     We assume Ms is conserved (ie has already been checked for).

        !     In the momentum space description the overall crystal
        !     momentum must be conserved up to a reciprocal sys_global%lattice%lattice
        !     vector (i.e. satisfy translational symmetry).
        !     We assume this is also already checked.

        excitation = get_excitation(f1,f2)
        ! Connected determinants can differ by (at most) 2 spin orbitals.
        if (excitation%nexcit <= 2) then
            non_zero = .true.
        end if

        if (non_zero) then
            select case(excitation%nexcit)
            ! Apply Slater--Condon rules.
            case(0)

                ! < D | H | D > = \sum_i < i | h(i) | i > + \sum_i \sum_{j>i} < ij || ij >
                hmatel = slater_condon0_hub_k(f1)

!            case(1)

                ! < D | H | D_i^a > = < i | h(a) | a > + \sum_j < ij || aj >

                ! Single excitations are not connected in the momentum space
                ! basis.

                ! One electron operator
                ! The kinetic operator is diagonal in the momentum space basis.

                ! Two electron operator
                ! < ij | aj > = 0 only if crystal momentum is conserved up to
                ! a reciprocal sys_global%lattice%lattice vector.
                ! As k_i /= k_j, this cannot be met.

            case(2)

                ! < D | H | D_{ij}^{ab} > = < ij || ab >

                ! Two electron operator
                hmatel = slater_condon2_hub_k(excitation%from_orb(1), excitation%from_orb(2), &
                                            & excitation%to_orb(1), excitation%to_orb(2),excitation%perm)

            end select
        end if

    end function get_hmatel_hub_k

    pure function slater_condon0_hub_k(f) result(hmatel)

        ! In:
        !    f: bit string representation of the Slater determinant.
        ! Returns:
        !    < D_i | H | D_i >, the diagonal Hamiltonian matrix elements, for
        !        the Hubbard model in momentum space.

        use determinants, only: decode_det, basis_fns, basis_length
        use system

        real(p) :: hmatel
        integer(i0), intent(in) :: f(basis_length)
        integer :: occ_list(sys_global%nel)
        integer :: i

        call decode_det(f, occ_list)

        ! < D | H | D > = \sum_i < i | h(i) | i > + \sum_i \sum_{j>i} < ij || ij >

        ! Two electron operator
        ! 1/2 \sum_i \sum_j < ij || ij >
        ! Some points to note:
        !   a) Crystal momentum is conserved for all < ij | ij > and < ij | ji > integrals
        !      by defintion.
        !   b) If i,j are of the same spin, then < ij | ij > = < ij | ji > and
        !      so < ij || ij > = 0.
        !   c) < ij | ij > = U/sys_global%lattice%nsites for all i,j.
        !   d) The double sum has 2*sys_global%nalpha*sys_global%nbeta terms corresponding to i,j of
        !      different spins.
        !   e) Thus  1/2 \sum_i \sum_j < ij || ij > = sys_global%nalpha*sys_global%nbeta*U/sys_global%lattice%nsites.
        hmatel = sys_global%nalpha*sys_global%nbeta*sys_global%hubbard%coulomb_k

        ! One electron operator
        ! Get directly rather than incur the cost of the if test in get_one_e_int_k.
        do i = 1, sys_global%nel
            hmatel = hmatel + basis_fns(occ_list(i))%sp_eigv
        end do

    end function slater_condon0_hub_k

    pure function slater_condon2_hub_k(i, j, a, b, perm) result(hmatel)

        ! In:
        !    i,j:  index of the spin-orbital from which an electron is excited in
        !          the reference determinant.
        !    a,b:  index of the spin-orbital into which an electron is excited in
        !          the excited determinant.
        !    perm: true if D and D_i^a are connected by an odd number of
        !          permutations.
        ! Returns:
        !    < D | H | D_ij^ab >, the Hamiltonian matrix element between a
        !    determinant and a double excitation of it in the momemtum space
        !    formulation of the Hubbard model.

        use hubbard_k, only: get_two_e_int_k

        real(p) :: hmatel
        integer, intent(in) :: i, j, a, b
        logical, intent(in) :: perm

        hmatel = get_two_e_int_k(i, j, a, b)

        if (perm) hmatel = -hmatel

    end function slater_condon2_hub_k

    pure subroutine slater_condon2_hub_k_excit(f, connection, hmatel)

        ! Generate the matrix element between a determinant and a double
        ! excitation in the momentum space formulation of the Hubbard model.
        ! WARNING: this routine assumes that the excitation is allowed (i.e.
        ! conserves crystal momentum).  It is, however, faster as symmetry
        ! checking is skipped.

        ! In:
        !    f: bit string representation of the Slater determinant, D.
        ! In/Out:
        !    connection: excit type describing the excitation between |D> and
        !    |D_ij^ab>.  On entry, only the from_orb and to_orb fields must be
        !    set.  On exit the from_orb and to_orb fields will be ordered
        !    and the perm field will be set.
        ! Out:
        !    hmatel: < D | H | D_ij^ab >, the Hamiltonian matrix element between a
        !    determinant and a double excitation of it in the momemtum space
        !    formulation of the Hubbard model.

        use excitations, only: excit, find_excitation_permutation2
        use basis, only: basis_length
        use system

        integer(i0), intent(in) :: f(basis_length)
        type(excit), intent(inout) :: connection
        real(p), intent(out) :: hmatel

        integer :: tmp

        ! The permuting algorithm works by lining up the min(i,j) with
        ! min(a,b) and max(i,j) with max(a,b) and hence we can find out
        ! whether the Coulomb or exchange integral is non-zero.
        ! Thus (i,j) and (a,b) must be ordered.
        if (connection%from_orb(1) > connection%from_orb(2)) then
            ! Swap.
            tmp = connection%from_orb(1)
            connection%from_orb(1) = connection%from_orb(2)
            connection%from_orb(2) = tmp
        end if
        if (connection%to_orb(1) > connection%to_orb(2)) then
            ! Swap
            tmp = connection%to_orb(1)
            connection%to_orb(1) = connection%to_orb(2)
            connection%to_orb(2) = tmp
        end if

        ! a) Sign from value of U as well---U can be negative!
        hmatel = sys_global%hubbard%coulomb_k

        ! b) Negative sign from permuting the determinants so that they line
        ! up?
        call find_excitation_permutation2(f, connection)
        if (connection%perm) then
            ! Matrix element gets a -sign from rearranging determinants so
            ! that they maximally line up.
            hmatel = -hmatel
        end if

        ! c) Because the only non-zero excitations are when i is alpha and
        ! j is beta or vice-versa, only the Coulomb integral or the exchange
        ! integral is non-zero.  If it's the exchange
        ! integral, then we obtain an additional minus sign.
        if (mod(connection%from_orb(1)-connection%to_orb(1),2) /= 0) then
            ! (i',a') are (alpha,beta) or (beta,alpha).
            ! Thus it is the exchange integral which contributes to the
            ! connecting matrix element.
            hmatel = -hmatel
        end if

    end subroutine slater_condon2_hub_k_excit

end module hamiltonian_hub_k

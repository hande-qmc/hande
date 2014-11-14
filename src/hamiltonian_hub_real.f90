module hamiltonian_hub_real

! Module for evaluating Hamiltonian matrix elements Hubbard model in real space
! (i.e.  using a local orbital basis set).

use const

implicit none

contains

    pure function get_hmatel_hub_real(sys, f1, f2) result(hmatel)

        ! In:
        !    sys: system to be studied.
        !    f1, f2: bit string representation of the Slater
        !        determinants D1 and D2 respectively.
        ! Returns:
        !    Hamiltonian matrix element between the two determinants,
        !    < D1 | H | D2 >, where the determinants are formed from
        !    real space basis functions.

        ! Used in the real space formulation of the Hubbard model only.

        use excitations, only: excit_t, get_excitation
        use system, only: sys_t

        real(p) :: hmatel
        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f1(sys%basis%string_len), f2(sys%basis%string_len)
        logical :: non_zero
        type(excit_t) :: excitation

        hmatel = 0.0_p
        non_zero = .false.

        ! Test to see if Hamiltonian matrix element is non-zero.

        ! We assume Ms is conserved (ie has already been checked for).
        excitation = get_excitation(sys%nel, sys%basis, f1, f2)
        ! Connected determinants can differ by (at most) 2 spin orbitals.
        if (excitation%nexcit <= 2) then
            ! Space group symmetry not currently implemented.
            non_zero = .true.
        end if

        ! Matrix elements in the real space formulation are quite simple.

        ! 1. < i | T | i > = 0
        !    Thus the one-electron terms only occur between single excitation
        !    matrix elements.
        ! 2. < m,s1 n,s2 | U | p,s1 q,s2 > = U \delta_{m,n} \delta_{m,p} \delta_{m,q} , s1/=s2
        !    Thus the Coulomb integrals that occur in < D | H | D_i^a > and
        !    < D | H | D_{ij}^{ab} > are zero.

        if (non_zero) then
            select case(excitation%nexcit)
            ! Apply Slater--Condon rules.
            case(0)

                ! < D | H | D > = \sum_i < i | h(i) | i > + \sum_i \sum_{j>i} < ij || ij >
                hmatel = slater_condon0_hub_real(sys, f1)

            case(1)

                hmatel = slater_condon1_hub_real(sys, excitation%from_orb(1), excitation%to_orb(1), excitation%perm)

!            case(2)

                ! < D | H | D_{ij}^{ab} > = < ij || ab >
                !                         = 0 within the real space formulation
                !                             of the Hubbard model.

            end select
        end if

    end function get_hmatel_hub_real

    pure function slater_condon0_hub_real(sys, f) result(hmatel)

        ! In:
        !    sys: system to be studied.
        !    f: bit string representation of the Slater determinant.
        ! Returns:
        !    < D_i | H | D_i >, the diagonal Hamiltonian matrix elements, for
        !        the Hubbard model in real space.

        use determinants, only: decode_det
        use real_lattice, only: get_one_e_int_real, get_coulomb_matel_real
        use system, only: sys_t

        real(p) :: hmatel
        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f(sys%basis%string_len)
        integer :: root_det(sys%nel)
        integer :: i

        ! < D | H | D > = \sum_i < i | h(i) | i > + \sum_i \sum_{j>i} < ij || ij >
        hmatel = 0.0_p

        ! < i | T | i > = 0 within the real space formulation of the
        ! Hubbard model, unless site i is its own periodic image, in
        ! which case it has a kinetic interaction with its self-image.
        ! This only arises if there is at least one crystal cell vector
        ! which is a unit cell vector.
        if (sys%real_lattice%t_self_images) then
            call decode_det(sys%basis, f, root_det)
            do i = 1, sys%nel
                hmatel = hmatel + get_one_e_int_real(sys, root_det(i), root_det(i))
            end do
        end if

        ! Two electron operator
        hmatel = hmatel + get_coulomb_matel_real(sys, f)

    end function slater_condon0_hub_real

    pure function slater_condon1_hub_real(sys, i, a, perm) result(hmatel)

        ! In:
        !    sys: system to be studied.
        !    i: index of the spin-orbital from which an electron is excited in
        !        the reference determinant.
        !    a: index of the spin-orbital into which an electron is excited in
        !        the excited determinant.
        !    perm: true if D and D_i^a are connected by an odd number of
        !        permutations.
        ! Returns:
        !    < D | H | D_i^a >, the Hamiltonian matrix element between a
        !        determinant and a single excitation of it.

        use real_lattice, only: get_one_e_int_real
        use system, only: sys_t

        real(p) :: hmatel
        type(sys_t), intent(in) :: sys
        integer, intent(in) :: i, a
        logical, intent(in) :: perm

        ! < D | H | D_i^a > = < i | h(a) | a > + \sum_j < ij || aj >

        ! One electron operator
         hmatel = get_one_e_int_real(sys, i, a)

        ! Two electron operator
        ! < D | U | D_i^a > = 0 within the real space formulation of the
        ! Hubbard model.

        if (perm) hmatel = -hmatel

    end function slater_condon1_hub_real

    pure subroutine slater_condon1_hub_real_excit(sys, f, connection, hmatel)

        ! Generate the matrix element between a determinant and a single
        ! excitation in the real space formulation of the Hubbard model.
        ! WARNING: this routine assumes that the excitation is allowed (i.e.
        ! the excitation is from an occupied orbital, i, to an unoccupied
        ! orbital, a, which is connected to i).  By skipping such checks, it is,
        ! however, faster.

        ! In:
        !    sys: system to be studied.
        !    f: bit string representation of the Slater determinant, D.
        ! In/Out:
        !    connection: excit_t type describing the excitation between |D> and
        !    |D_i^a>.  On entry, only the from_orb and to_orb fields must be
        !    set.  On exit the perm field will also be set.
        ! Out:
        !    hmatel: < D | H | D_i^a >, the Hamiltonian matrix element between a
        !    determinant and a single excitation of it in the real space
        !    formulation of the Hubbard model.

        use excitations, only: excit_t, find_excitation_permutation1
        use system, only: sys_t

        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f(sys%basis%string_len)
        type(excit_t), intent(inout) :: connection
        real(p), intent(out) :: hmatel

        ! a) Find out permutation required to line up determinants.
        call find_excitation_permutation1(f, connection)

        ! b) The matrix element connected |D> and |D_i^a> is <i|h|a> = -t.
        if (connection%perm) then
            hmatel = sys%hubbard%t
        else
            hmatel = -sys%hubbard%t
        end if

    end subroutine slater_condon1_hub_real_excit

end module hamiltonian_hub_real

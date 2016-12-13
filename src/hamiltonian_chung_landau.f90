module hamiltonian_chung_landau

! Module for evaluating Hamiltonian matrix elements in the Chung--Landau
! Hamiltonian (Phys Rev B 85 (2012) 115115).

! Due to the similarities with the Hubbard model, we can essentially piggy-back
! on the Hubbard model (in a local orbital basis) for all but the diagonal
! matrix elements.

! Note: to achieve equivalence with the Chung--Landau results, the
! sys%real_lattice%finite_cluster option must be set.

use const

implicit none

contains

    pure function get_hmatel_chung_landau(sys, f1, f2) result(hmatel)

        ! In:
        !    sys: system to be studied.
        !    f1, f2: bit string representation of the Slater
        !        determinants D1 and D2 respectively.
        ! Returns:
        !    Hamiltonian matrix element between the two determinants,
        !    < D1 | H | D2 >, where the determinants are formed from
        !    real space basis functions.

        ! Used in the Chung--Landau model only.

        use excitations, only: excit_t, get_excitation
        use hamiltonian_hub_real, only: slater_condon1_hub_real
        use system, only: sys_t
        use hamiltonian_data

        type(hmatel_t) :: hmatel
        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f1(sys%basis%tot_string_len), f2(sys%basis%tot_string_len)
        logical :: non_zero
        type(excit_t) :: excitation

        hmatel%r = 0.0_p
        non_zero = .false.

        ! Test to see if Hamiltonian matrix element is non-zero.

        ! We assume Ms is conserved (ie has already been checked for).
        excitation = get_excitation(sys%nel, sys%basis, f1, f2)

        ! Matrix elements in this model are quite simple.
        ! < D | H | D > = \sum_i < i | h(i) | i > + \sum_i \sum_{j>i} < ij || ij >
        !               = -t \sum_{<i,j>} ( c^{\dagger}_{i} c_{j} + c^{\dagger}_{j} c_{i}  )
        !                 + U \sum_{<i,j>} n_{i} n_{j}
        ! where <i,j> indicates i and j are nearest neighbours and j>i.
        ! This is not too dissimilar to the Hubbard model!

        ! Connected determinants can differ by (at most) 2 spin orbitals, in
        ! general.
        select case(excitation%nexcit)
        ! Apply Slater--Condon rules.
        case(0)

            ! < D | H | D > = \sum_i < i | h(i) | i > + \sum_i \sum_{j>i} < ij || ij >
            hmatel%r = slater_condon0_chung_landau(sys, f1)

        case(1)

            ! Identical to the Hubbard model in a local orbital basis.
            hmatel%r = slater_condon1_hub_real(sys, excitation%from_orb(1), excitation%to_orb(1), excitation%perm)

!        case(2)

            ! < D | H | D_{ij}^{ab} > = < ij || ab >
            !                         = 0 within the Chung--Landau Hamiltonian.

        end select

    end function get_hmatel_chung_landau

    pure function slater_condon0_chung_landau(sys, f) result(hmatel)

        ! In:
        !    sys: system to be studied.
        !    f: bit string representation of the Slater determinant.
        ! Returns:
        !    < D_i | H | D_i >, the diagonal Hamiltonian matrix elements, for
        !        the Chung--Landau model.

        use determinants, only: decode_det
        use real_lattice, only: get_one_e_int_real
        use system, only: sys_t

        real(p) :: hmatel
        type(sys_t), intent(in) :: sys
        integer(i0), intent(in) :: f(sys%basis%tot_string_len)
        integer :: root_det(sys%nel)
        integer :: i, j, indi, indj, posi, posj

        hmatel = 0.0_p
        call decode_det(sys%basis, f, root_det)

        ! < i | T | i > = 0 unless site i is its own periodic image, in
        ! which case it has a kinetic interaction with its self-image.
        ! This only arises if there is at least one crystal cell vector
        ! which is a unit cell vector.
        if (sys%real_lattice%t_self_images) then
            do i = 1, sys%nel
                hmatel = hmatel + get_one_e_int_real(sys, root_det(i), root_det(i))
            end do
        end if

        ! Two electron operator
        ! This can be done efficiently with bit string operations (see
        ! get_coulomb_matel_real) in 1D.
        do i = 1, sys%nel
            ! Connected to self-image?
            posi = sys%basis%bit_lookup(1,root_det(i))
            indi = sys%basis%bit_lookup(2,root_det(i))
            ! Diagonal term is non-zero if i is connected to its own periodic
            ! image.
            if (btest(sys%real_lattice%tmat(indi, root_det(i)), posi)) then
                hmatel = hmatel + sys%hubbard%u
            end if
            do j = i+1, sys%nel
                ! i <-> j via periodic boundary conditions.
                if (btest(sys%real_lattice%tmat(indi, root_det(j)), posi)) then
                    hmatel = hmatel + sys%hubbard%u
                end if
                ! i <-> j directly.
                posj = sys%basis%bit_lookup(1,root_det(j))
                indj = sys%basis%bit_lookup(2,root_det(j))
                if (btest(sys%real_lattice%tmat(indj, root_det(i)), posj)) then
                    hmatel = hmatel + sys%hubbard%u
                end if
            end do
        end do

    end function slater_condon0_chung_landau

end module hamiltonian_chung_landau

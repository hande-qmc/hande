module hamiltonian_chung_landau

! Module for evaluating Hamiltonian matrix elements in the Chung--Landau
! Hamiltonian (Phys Rev B 85 (2012) 115115).

! Due to the similarities with the Hubbard model, we can essentially piggy-back
! on the Hubbard model (in a local orbital basis) for all but the diagonal
! matrix elements.

! Note: to achieve equivalence with the Chung--Landau results, the
! finite_cluster option must be set.

use const

implicit none

contains

    pure function get_hmatel_chung_landau(f1, f2) result(hmatel)

        ! In:
        !    f1, f2: bit string representation of the Slater
        !        determinants D1 and D2 respectively.
        ! Returns:
        !    Hamiltonian matrix element between the two determinants,
        !    < D1 | H | D2 >, where the determinants are formed from
        !    real space basis functions.

        ! Used in the Chung--Landau model only.

        use determinants, only: basis_length
        use excitations, only: excit, get_excitation
        use hamiltonian_hub_real, only: slater_condon1_hub_real

        real(p) :: hmatel
        integer(i0), intent(in) :: f1(basis_length), f2(basis_length)
        logical :: non_zero
        type(excit) :: excitation

        hmatel = 0.0_p
        non_zero = .false.

        ! Test to see if Hamiltonian matrix element is non-zero.

        ! We assume Ms is conserved (ie has already been checked for).
        excitation = get_excitation(f1, f2)

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
            hmatel = slater_condon0_chung_landau(f1)

        case(1)

            ! Identical to the Hubbard model in a local orbital basis.
            hmatel = slater_condon1_hub_real(excitation%from_orb(1), excitation%to_orb(1), excitation%perm)

!        case(2)

            ! < D | H | D_{ij}^{ab} > = < ij || ab >
            !                         = 0 within the Chung--Landau Hamiltonian.

        end select

    end function get_hmatel_chung_landau

    pure function slater_condon0_chung_landau(f) result(hmatel)

        ! In:
        !    f: bit string representation of the Slater determinant.
        ! Returns:
        !    < D_i | H | D_i >, the diagonal Hamiltonian matrix elements, for
        !        the Chung--Landau model.

        use basis, only: bit_lookup
        use determinants, only: decode_det, basis_length
        use hubbard_real, only: t_self_images, tmat, get_one_e_int_real
        use system, only: nel, hubu

        real(p) :: hmatel
        integer(i0), intent(in) :: f(basis_length)
        integer :: root_det(nel)
        integer :: i, j, indi, indj, posi, posj

        hmatel = 0.0_p
        call decode_det(f, root_det)

        ! < i | T | i > = 0 unless site i is its own periodic image, in
        ! which case it has a kinetic interaction with its self-image.
        ! This only arises if there is at least one crystal cell vector
        ! which is a unit cell vector.
        if (t_self_images) then
            do i = 1, nel
                hmatel = hmatel + get_one_e_int_real(root_det(i), root_det(i))
            end do
        end if

        ! Two electron operator
        ! This can be done efficiently with bit string operations (see
        ! get_coulomb_matel_real) in 1D.
        do i = 1, nel
            ! Connected to self-image?
            posi = bit_lookup(1,root_det(i))
            indi = bit_lookup(2,root_det(i))
            ! Diagonal term is non-zero if i is connected to its own periodic
            ! image.
            if (btest(tmat(indi, root_det(i)), posi)) then
                hmatel = hmatel + hubu
            end if
            do j = i+1, nel
                ! i <-> j via periodic boundary conditions.
                if (btest(tmat(indi, root_det(j)), posi)) then
                    hmatel = hmatel + hubu
                end if
                ! i <-> j directly.
                posj = bit_lookup(1,root_det(j))
                indj = bit_lookup(2,root_det(j))
                if (btest(tmat(indj, root_det(i)), posj)) then
                    hmatel = hmatel + hubu
                end if
            end do
        end do

    end function slater_condon0_chung_landau

end module hamiltonian_chung_landau

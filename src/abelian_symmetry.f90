module abelian_symmetry

! Module for symmetry subroutines/fucntions shared between pg and translational
! symmetry.

use const

implicit none

contains

    pure function cross_product_basis_abelian(sys, i, j) result(sym_ij)

        ! In:
        !    sys: information on the system under consideration.
        !    i,j: (indices of) spin-orbitals.
        ! Returns:
        !    The symmetry representation
        !    The representation of the irreducible representation
        !    formed from the direct product i%sym \cross j%sym, where i%sym is
        !    the irreducible representation spanned by the i-th spin-orbital
        !    which can also potentially be point group symmetry with an Lz
        !    symmetry or a translational symmetry.

        use system, only: sys_t

        integer :: sym_ij
        type(sys_t), intent(in) :: sys
        integer, intent(in) :: i, j

        sym_ij = sys%read_in%cross_product_sym_ptr(sys%read_in, &
                    sys%basis%basis_fns(i)%sym, sys%basis%basis_fns(j)%sym)

    end function cross_product_basis_abelian

    elemental function is_basis_abelian_sym(sys,sym) result(valid)

        ! In:
        !    sys: system being studied
        !    sym: representation of an irreducible representation.
        !
        ! Returns:
        !   True if sym could have been one of the basis function symmetries
        !   False otherwise

        use system, only: sys_t

        logical :: valid
        type(sys_t), intent(in) :: sys
        integer, intent(in) :: sym

        valid = sym>=sys%sym0 .and. sym<=sys%sym_max

    end function is_basis_abelian_sym

    elemental function is_gamma_irrep_abelian(pg_sym, sym) result(is_gamma)

        ! In:
        !    pg_sym: information on the symmetries of the basis functions.
        !    sym: representation of an irreducible representation.
        ! Returns:
        !    True if sym represents the Gamma_1 (totally symmetric) irreducible
        !    representation.

        use symmetry_types, only: pg_sym_t

        logical :: is_gamma
        type(pg_sym_t), intent(in) :: pg_sym
        integer, intent(in) :: sym

        ! Value is set appropriately for both pg and translational symmetry.
        is_gamma = pg_sym%gamma_sym == sym

    end function is_gamma_irrep_abelian

    pure function symmetry_orb_list_abelian(read_in, basis, orb_list) result(isym)

        ! Function to obtain symmetry of a given orb list. Uses cross product
        ! pointer to enable use for both molecular and periodic systems.

        ! In:
        !    read_in: information on the system being studied. We use the
        !       symmetry information.
        !    basis: info about the single particle basis set.
        !    orb_list: list of orbitals (e.g. determinant).
        ! Returns:
        !    symmetry index of list (i.e. direct product of the representations
        !    of all the orbitals in the list).

        use basis_types, only: basis_t
        use system, only: sys_read_in_t

        integer :: isym
        type(sys_read_in_t), intent(in) :: read_in
        type(basis_t), intent(in) :: basis
        integer, intent(in) :: orb_list(:)

        integer :: i

        ! This value is set appropriately in both pg_sym and mom_sym
        ! initialisation.
        isym = read_in%pg_sym%gamma_sym

        do i = lbound(orb_list, dim=1), ubound(orb_list, dim=1)
            isym = read_in%cross_product_sym_ptr(read_in, isym, basis%basis_fns(orb_list(i))%sym)
        end do

    end function symmetry_orb_list_abelian

end module abelian_symmetry

module read_in_symmetry

! Module for symmetry subroutines/functions shared between pg and translational
! symmetry. For more specific comments on either symmetry see the comments at
! the start of pg_symmetry.f90 and momentum_sym_read_in.f90 respectively.
!
! Here we will briefly detail the logical structure of the symmetry
! implementations contained.
!
! Implementation Considerations
! -----------------------------
!
! These implementations are designed to enable easy implementation of symmetry
! constraints simultaneously for point group and translation symmetries without
! significant overhead or complexity, while avoiding compile-time specification.
!
! This is (hopefully) achieved by utilising standardised interfaces for basic
! symmetry operations (cross products and complex conjugates). These functions
! are accessible by pointers (sys%read_in%cross_product_sym_ptr and
! sys%read_in%sym_conj_ptr) set during either init_pg_symmetry or
! init_read_in_momentum_symmetry
!
! Using a standard interface requires both symmetries be represented by a
! single 32-bit integer value. For pg_sym this corresponds to the bit string
! representation of the irreducible representation, while for translational
! symmetry this is an index obtained by mapping all possible kpoint values
! to a cuboid.
!
! Code Structure
! --------------
!
! This module contains various utility functions/subroutines for
! commonly required operations (eg. the cross product of two basis function
! symmetries or obtaining the symmetry of an orbital list). These functions
! will function correctly for point group symmetry (with or without Lz
! symmetry) and translational symmetry.
!
! pg_symmetry.f90 contains functions/subroutines to:
! - initialise pg_symmetry for a non-model system and print required information to
!       report.
! - perform pg_sym cross product and conjugation (best accessed through pointers
!       for non-specific uses).
! - obtain the Lz value of a given symmetry representation.
!
! momentum_sym_read_in.f90 contains functions/subroutines to:
! - initialise translational symmetry for non-model system and print required
!       information to report. This includes functions converting between the
!       various symmetry represenations used (see notes at start or
!       momentum_sym_read_in.f90 or in symmetry_types.f90) and functions performing
!       symmetry operations within the kpoint vector format.
! - perform translational sym cross product and conjugation using symmetry indexes
!       (best accessed through pointers for non-specific uses).
!
! symmetry.f90 contains more general functions to obtain the cross product of two
! symmetries and the overall symmetry of an orbital list, also allowing use with
! model periodic systems. As this requires a select statement the use of these
! functions is by no means optimal.
!
! symmetry_types.f90 contains all derived types relevant to symmetries of any
! description, as well as functions to fully deallocate these derived types.
!
! Various other modules contain functions for model system symmetries. These are
! not detailed here.

use const

implicit none

contains

    pure function cross_product_basis_read_in(sys, i, j) result(sym_ij)

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

    end function cross_product_basis_read_in

    elemental function is_in_read_in_basis_sym(sys,sym) result(valid)

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

        ! Check that the symmetry index is within the allowed range.
        valid = sym>=sys%sym0 .and. sym<=sys%sym_max

    end function is_in_read_in_basis_sym

    elemental function is_gamma_irrep_read_in(pg_sym, sym) result(is_gamma)

        ! In:
        !    pg_sym: information on the symmetries of the basis functions.
        !    sym: representation of a symmetry.
        ! Returns:
        !    True if sym represents the Gamma_1 (totally symmetric) irreducible
        !    representation.

        use symmetry_types, only: pg_sym_t

        logical :: is_gamma
        type(pg_sym_t), intent(in) :: pg_sym
        integer, intent(in) :: sym

        ! Value is set appropriately for both pg and translational symmetry.
        is_gamma = pg_sym%gamma_sym == sym

    end function is_gamma_irrep_read_in

    pure function symmetry_orb_list_read_in(read_in, basis, orb_list) result(isym)

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

    end function symmetry_orb_list_read_in

end module read_in_symmetry

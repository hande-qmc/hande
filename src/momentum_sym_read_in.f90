module momentum_sym_read_in

! Module for handing crystal momentum symmetry routines unique to real, periodic systems.

use system

implicit none

contains

    pure function symmetry_orb_list_read_in(sys, orb_list) result(isym)
        use const, only: int_64, p
        type(sys_t), intent(in) :: sys
        integer, intent(in) :: orb_list(:)
        integer :: isym(3)

        integer :: i

        call decompose_abelian_sym(sys%read_in%mom_sym%gamma_sym, &
                sys%read_in%mom_sym%propbitlen, isym)
        do i = lbound(orb_list, dim = 1), ubound(orb_list, dim = 1)
            call cross_product_read_in_abelian(sys%read_in%mom_sym%nprop, &
                sys%basis%basis_fns(orb_list(i))%l, isym, isym)
        end do

    end function symmetry_orb_list_read_in

    pure subroutine get_kpoint_inverse(s1, nprop, inv)
        ! Takes quantum numbers corresponding to a single kpoint and finds inverse reciprocal
        ! lattice vector.
        integer, intent(in) :: s1, nprop
        integer, intent(out) :: inv

        inv = nprop - s1
    end subroutine get_kpoint_inverse

    pure function is_gamma_sym_periodic_read_in(mom_sym, sym) result(is_gamma_sym)
        use symmetry_types, only: mom_sym_t
        type(mom_sym_t), intent(in) :: mom_sym
        integer, intent(in) :: sym(3)
        logical :: is_gamma_sym

        is_gamma_sym = all(sym == mom_sym%gamma_point)
    end function is_gamma_sym_periodic_read_in

! Various possible cross products to be cut down later when decide what we actually need.

    pure function cross_product_translational_sym(mom_sym, s1, s2) result(prod)

        use const, only: int_64
        use symmetry_types, only: mom_sym_t

        integer(int_64), intent(in) :: s1, s2
        type(mom_sym_t), intent(in) :: mom_sym
        integer(int_64) :: prod
        integer :: abelian1(3), abelian2(3), res(3)

        call decompose_abelian_sym(s1, mom_sym%propbitlen, abelian1)
        call decompose_abelian_sym(s2, mom_sym%propbitlen, abelian2)
        call cross_product_read_in_abelian(mom_sym%nprop, abelian1, abelian2, res)
        call compose_abelian_sym(res, mom_sym%propbitlen, prod)

    end function cross_product_translational_sym

    pure subroutine cross_product_abelian_basis(mom_sym, b1, b2, basis_fns, prod)
        use basis_types, only: basis_fn_t
        use symmetry_types, only: mom_sym_t
        type(mom_sym_t), intent(in) :: mom_sym
        integer, intent(in) :: b1, b2
        type(basis_fn_t), intent(in) :: basis_fns(:)
        integer, intent(out) :: prod(3)

        call cross_product_read_in_abelian(mom_sym%nprop, basis_fns(b1)%l, basis_fns(b2)%l, prod)

    end subroutine cross_product_abelian_basis

    pure subroutine cross_product_read_in_abelian(nprop, a1, a2, prod)

        use const, only: int_64

        integer, intent(in) :: a1(3), a2(3), nprop(3)
        integer, intent(out) :: prod(3)

        prod = modulo(a1 + a2, nprop)

    end subroutine cross_product_read_in_abelian

! Indexing conversion routines:

    pure subroutine decompose_abelian_sym(isym, propbitlen, abel_sym)
        ! Takes symmetry index for translationally symmetric wavefunction and
        ! returns abelian representation of three "quantum numbers". In accordance
        ! with approach used in NECI, values stored according to:
        !   isym = 1 + \sum_i sym(i) * 2 ** (propbitlen * (i-1))

        use const, only: int_32, int_64

        integer(int_64), intent(in) :: isym
        integer, intent(in) :: propbitlen
        integer, intent(out) :: abel_sym(3)

        ! Use Iand and mask to select only bits in first propbitlen bits of
        ! isym.
        abel_sym(1) = int(Iand(isym, 2_int_64**(propbitlen-1)), int_32)
        ! Bit shift to access correct bits of isym.
        abel_sym(2) = int(Iand(Ishft(isym, -propbitlen), &
                            2_int_64 ** (propbitlen - 1)), int_32)
        abel_sym(3) = int(Ishft(isym, -(propbitlen * 2)), int_32)

    end subroutine decompose_abelian_sym

    pure subroutine compose_abelian_sym(abel_sym, propbitlen, isym)
        ! Takes abelian symmetry "quantum numbers" and combines into single value.
        ! Opposite effect to decompose_abelian_sym.
        use const, only: int_32, int_64

        integer, intent(in) :: abel_sym(3)
        integer, intent(in) :: propbitlen
        integer(int_64), intent(out) :: isym

        isym = abel_sym(1) + Ishft(abel_sym(2), propbitlen) + &
                    Ishft(abel_sym(3), propbitlen*2)
    end subroutine compose_abelian_sym

    pure function get_kpoint_index(a, nprop) result(ind)
        ! Converts from abelian symmetry quantum numbers into unique index.
        ! If we know size of unit cell, can calculate unique index by tiling
        ! first along axis 1, then 2, then 3.
        integer, intent(in) :: a(3), nprop(3)
        integer :: ind
        ! Want to start from index 1 at gamma point (0,0,0)
        ind = 1 + a(1) + (nprop(1) - 1) * a(2) + (nprop(1) - 1) * (nprop(2) - 1) * a(3)
    end function get_kpoint_index

    pure subroutine get_kpoint_numbers(ind, nprop, a)

        integer, intent(in) :: ind, nprop(3)
        integer, intent(out) :: a(3)

        a(1) = modulo(ind, nprop(2)-1)
        a(2) = modulo(ind - a(1), (nprop(1)-1))
        a(3) = (ind - a(1) - a(2)*(nprop(1)-1)) / ((nprop(1)-1)*(nprop(2)-1))

    end subroutine get_kpoint_numbers

end module momentum_sym_read_in

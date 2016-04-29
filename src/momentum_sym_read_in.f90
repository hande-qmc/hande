module momentum_sym_read_in

! Module for handing crystal momentum symmetry routines unique to real, periodic systems.

use system

implicit none

contains

    subroutine init_basis_momentum_symmetry_info(sys)

        use checking, only: check_allocate, check_deallocate
        use errors, only: stop_all

        type(sys_t), intent(inout) :: sys
        integer, allocatable :: nbasis_sym(:), current_index(:)
        integer :: i, ierr

        allocate(nbasis_sym(1:sys%nsym), stat=ierr)
        call check_allocate('nbasis_sym',sys%nsym,ierr)

        nbasis_sym = 0

        do i = 1, sys%basis%nbasis/2
            nbasis_sym(sys%basis%basis_fns(2*i-1)%sym) = &
                nbasis_sym(sys%basis%basis_fns(2*i-1)%sym) + 1
        end do
        if (.not. all(nbasis_sym == nbasis_sym(1))) call stop_all('init_basis_momentum_sym_info', 'Read in system has different &
                    &number of bands per kpoint.')
        sys%read_in%mom_sym%nbands = nbasis_sym(1)

        deallocate(nbasis_sym, stat=ierr)
        call check_deallocate('nbasis_sym',ierr)

        allocate(sys%read_in%mom_sym%basis_sym(sys%nsym,sys%read_in%mom_sym%nbands), stat=ierr)
        call check_allocate('basis_sym',sys%nsym*sys%read_in%mom_sym%nbands,ierr)

        allocate(current_index(sys%nsym), stat=ierr)
        call check_allocate('current_index',sys%nsym,ierr)
        current_index = 1

        associate(basis_sym=>sys%read_in%mom_sym%basis_sym, basis_fns=>sys%basis%basis_fns)
            do i = 1, sys%basis%nbasis/2
                basis_sym(basis_fns(2*i-1)%sym, current_index(basis_fns(2*i-1)%sym)) = i
                basis_fns(2*i-1:2*i)%sym_spin_index = current_index(basis_fns(2*i-1)%sym)
                current_index(basis_fns(2*i-1)%sym) = current_index(basis_fns(2*i-1)%sym) + 1
            end do
        end associate

        deallocate(current_index, stat=ierr)
        call check_deallocate('current_index',ierr)
    end subroutine init_basis_momentum_symmetry_info

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

        is_gamma_sym = all(modulo(sym, mom_sym%nprop) == mom_sym%gamma_point)
    end function is_gamma_sym_periodic_read_in

! Various possible cross products to be cut down later when decide what we actually need.

    pure function cross_product_abelian_basis(mom_sym, b1, b2, basis_fns) result(prod)
        use basis_types, only: basis_fn_t
        use symmetry_types, only: mom_sym_t
        type(mom_sym_t), intent(in) :: mom_sym
        integer, intent(in) :: b1, b2
        type(basis_fn_t), intent(in) :: basis_fns(:)
        integer :: prod

        prod = mom_sym%sym_table(basis_fns(b1)%sym, basis_fns(b2)%sym)

    end function cross_product_abelian_basis

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
        abel_sym(1) = int(Iand(isym, 2_int_64 ** propbitlen - 1), int_32)
        ! Bit shift to access correct bits of isym.
        abel_sym(2) = int(Iand(Ishft(isym, -propbitlen), &
                            2_int_64 ** (propbitlen) - 1), int_32)
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
        ind = 1 + a(1) + nprop(1) * a(2) + nprop(1) * nprop(2) * a(3)
    end function get_kpoint_index

    pure subroutine get_kpoint_numbers(ind, nprop, a)

        integer, intent(in) :: ind, nprop(3)
        integer, intent(out) :: a(3)
        integer :: scratch

        scratch = real(ind-1)/real(nprop(1)*nprop(2))
        a(3) = int(scratch)
        scratch = ind - 1 - a(3) * nprop(1) * nprop(2)
        a(2) = int(real(scratch)/real(nprop(1)))
        a(1) = scratch - a(2) * nprop(1)

    end subroutine get_kpoint_numbers

end module momentum_sym_read_in

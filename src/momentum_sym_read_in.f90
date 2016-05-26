module momentum_sym_read_in

! Module for handing crystal momentum symmetry routines unique to real, periodic systems.

use system

implicit none

contains

    subroutine init_basis_momentum_symmetry_info(sys)
        ! Initialises all required information for use of basis kpoint symmetry.

        ! In/Out:
        !   sys: system being studied. On output all information about basis function
        !       symmetry in read_in%mom_sym set as required.


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

        allocate(sys%read_in%mom_sym%sym_spin_basis(sys%read_in%mom_sym%nbands, 2, sys%nsym), stat=ierr)
        call check_allocate('basis_sym',sys%nsym*sys%read_in%mom_sym%nbands*2, ierr)

        allocate(current_index(sys%nsym), stat=ierr)
        call check_allocate('current_index',sys%nsym,ierr)
        current_index = 1

        associate(sym_spin_basis=>sys%read_in%mom_sym%sym_spin_basis, basis_fns=>sys%basis%basis_fns)
            do i = 1, sys%basis%nbasis/2
                sym_spin_basis(current_index(basis_fns(2*i-1)%sym), 2, basis_fns(2*i-1)%sym) = 2*i-1
                sym_spin_basis(current_index(basis_fns(2*i)%sym), 1, basis_fns(2*i)%sym) = 2*i
                basis_fns(2*i-1:2*i)%sym_spin_index = current_index(basis_fns(2*i-1)%sym)
                current_index(basis_fns(2*i-1)%sym) = current_index(basis_fns(2*i-1)%sym) + 1
            end do
        end associate

        deallocate(current_index, stat=ierr)
        call check_deallocate('current_index',ierr)

    end subroutine init_basis_momentum_symmetry_info

    pure function symmetry_orb_list_periodic_read_in(mom_sym, basis, orb_list) result(isym)

        ! In:
        !   mom_sym: basis function symmetry information.
        !   orb_list: list of orbitals (eg. in determinant).
        ! Returns:
        !   symmetry index of list (direct product of representation of all
        !       substituent orbitals.

        ! For momentum symmetry in real (read in from an FCIDUMP), periodic
        ! systems.

        use const, only: int_64, p
        use symmetry_types, only: mom_sym_t
        type(mom_sym_t), intent(in) :: mom_sym
        type(basis_t), intent(in) :: basis
        integer, intent(in) :: orb_list(:)
        integer :: isym

        integer :: i

        isym = int(mom_sym%gamma_sym)
        do i = lbound(orb_list, dim = 1), ubound(orb_list, dim = 1)
            isym = cross_product_periodic_read_in(mom_sym, &
                basis%basis_fns(orb_list(i))%sym, isym)
        end do

    end function symmetry_orb_list_periodic_read_in

    pure function is_gamma_sym_periodic_read_in(mom_sym, sym) result(is_gamma_sym)
        ! Checks if symmetry given is the gamma point symmetry.
        ! In:
        !   mom_sym: basis function symmetry information.
        !   sym: kpoint to compare, expressed via 3 integers.
        ! Returns:
        !   true: if symmetry provided is gamma sym.
        !   false: otherwise.

        ! For momentum symmetry in real (read in from an FCIDUMP), periodic
        ! systems.

        use symmetry_types, only: mom_sym_t
        type(mom_sym_t), intent(in) :: mom_sym
        integer, intent(in) :: sym(3)
        logical :: is_gamma_sym

        is_gamma_sym = all(modulo(sym, mom_sym%nprop) == mom_sym%gamma_point)

    end function is_gamma_sym_periodic_read_in

    pure function mom_sym_conj(mom_sym, sym) result(rsym)
        ! Returns symmetry of complex conjugate of provided sym.
        ! Since using complex plane waves, e^(ik.r), this will in
        ! general be e^(-ik.r) so just return inverse symmetry.

        ! In:
        !   mom_sym: basis function symmetry information.
        !   sym: symmetry to compare.
        ! Returns:
        !   Symmetry of complex conjugate of sym.

        ! For momentum symmetry in real (read in from an FCIDUMP), periodic
        ! systems.

        use symmetry_types, only: mom_sym_t
        type(mom_sym_t), intent(in) :: mom_sym
        integer, intent(in) :: sym
        integer :: rsym

        rsym = mom_sym%inv_sym(sym)

    end function mom_sym_conj

!--- Cross products ---

    pure function cross_product_periodic_read_in(mom_sym, a1, a2) result(prod)
        ! In:
        !   mom_sym: basis function symmetry information.
        !   a1, a2: symmetries to return cross product of.
        ! Returns:
        !   Direct product of symmetries a1 & a2.

        use symmetry_types, only: mom_sym_t
        type(mom_sym_t), intent(in) :: mom_sym
        integer, intent(in) :: a1, a2
        integer :: prod

        prod = mom_sym%sym_table(a1, a2)

    end function cross_product_periodic_read_in

    ! [review] - FDM: Is this function just a duplicate of the one above?
    pure function cross_product_periodic_basis(mom_sym, b1, b2, basis_fns) result(prod)
        ! In:
        !   mom_sym: basis function symmetry information.
        !   b1, b2: indicies of spinorbitals.
        ! Returns:
        !   Symmetry index corresponding to the product of
        !       symmetries of b1 & b2.

        use basis_types, only: basis_fn_t
        use symmetry_types, only: mom_sym_t
        type(mom_sym_t), intent(in) :: mom_sym
        integer, intent(in) :: b1, b2
        type(basis_fn_t), intent(in) :: basis_fns(:)
        integer :: prod

        prod = cross_product_periodic_read_in(mom_sym, basis_fns(b1)%sym, basis_fns(b2)%sym)

    end function cross_product_periodic_basis

!--- Indexing conversion routines: ---

    ! [review] - FDM: I think this (and the following two) function(s) could do with an expanded explanation.
    pure subroutine decompose_abelian_sym(isym, propbitlen, abel_sym)
        ! Takes symmetry index for translationally symmetric wavefunction and
        ! returns abelian representation of three "quantum numbers". In accordance
        ! with approach used in NECI, values stored according to:
        !   isym = 1 + \sum_i sym(i) * 2 ** (propbitlen * (i-1))

        ! In:
        !   isym: symmetry index provided in FCIDUMP file.
        !   propbitlen: length of bit representation of each "quantum number"
        !       within isym.
        ! Out:
        !   abel_sym: array of 3 quantum numbers stored in isym.

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

    ! [review] - FDM: I think this module could do with some more linebreaks to
    ! [review] - FDM: keep formatting consistent (between function names and docstrings and
    ! [review] - FDM: module use statements and variable declarations for instance.
    pure function get_kpoint_index(a, nprop) result(ind)
        ! Converts from abelian symmetry quantum numbers into unique index.
        ! If we know size of unit cell, can calculate unique index by tiling
        ! first along axis 1, then 2, then 3 in 3 dimensions.
        ! In:
        !   a: array containing representation of kpoint wavevectors.
        !   nprop: condition of periodic bounary conditions used, ie the
        !       supercell dimension.
        ! Returns:
        !   Index of given kpoint within indexing scheme.
        integer, intent(in) :: a(3), nprop(3)
        integer :: ind
        ! Want to start from index 1 at gamma point (0,0,0)
        ind = 1 + a(1) + nprop(1) * a(2) + nprop(1) * nprop(2) * a(3)

    end function get_kpoint_index

    pure subroutine get_kpoint_numbers(ind, nprop, a)
        ! Get kpoint index, in 3D array, from index defined in get_kpoint_index.
        ! In:
        !   ind: index to decode.
        !   nprop: condition of periodic bounary conditions used, ie the
        !       supercell dimension.
        ! Out:
        !   a: array containing kpoint identifier in terms of 3 "quantum numbers".

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

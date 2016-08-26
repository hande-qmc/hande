module momentum_sym_read_in

! Module for handing crystal momentum symmetry routines unique to real, periodic systems.

use system

implicit none

contains

    subroutine init_basis_momentum_symmetry_info(sys)

        ! Initialises all required information for use of basis kpoint symmetry.

        ! Specifically, allocates and sets pg_sym%nbasis_sym,
        ! pg_sym%nbasis_sym_spin and pg_sym%sym_spin_basis_fns appropriately.

        ! In/Out:
        !   sys: system being studied. On output all information about basis function
        !       symmetry in read_in%mom_sym set as required.

        use checking, only: check_allocate, check_deallocate
        use errors, only: stop_all

        type(sys_t), intent(inout) :: sys
        integer, allocatable :: current_index(:,:)
        integer :: i, ierr, ims

        allocate(sys%read_in%pg_sym%nbasis_sym_spin(2,1:sys%nsym), stat=ierr)
        call check_allocate('nbasis_sym_spin',sys%nsym,ierr)

        allocate(sys%read_in%pg_sym%nbasis_sym(1:sys%nsym), stat=ierr)
        call check_allocate('nbasis_sym',sys%nsym,ierr)

        sys%read_in%pg_sym%nbasis_sym = 0
        sys%read_in%pg_sym%nbasis_sym_spin = 0

        associate(basis_fns=>sys%basis%basis_fns, &
                    nbasis_sym_spin=>sys%read_in%pg_sym%nbasis_sym_spin, &
                    nbasis_sym=>sys%read_in%pg_sym%nbasis_sym)
            do i = 1, sys%basis%nbasis
                ims = (basis_fns(i)%Ms+3)/2
                nbasis_sym_spin(ims, basis_fns(i)%sym) = &
                            nbasis_sym_spin(ims, basis_fns(i)%sym) + 1
                nbasis_sym(basis_fns(i)%sym) = nbasis_sym(basis_fns(i)%sym) + 1
            end do
        end associate

        allocate(sys%read_in%pg_sym%sym_spin_basis_fns(maxval(sys%read_in%pg_sym%nbasis_sym_spin), &
                        2, sys%sym0_tot:sys%sym_max_tot), stat=ierr)
        call check_allocate('nbasis_sym_spin', &
                sys%nsym*maxval(sys%read_in%pg_sym%nbasis_sym_spin)*2, ierr)


        allocate(current_index(sys%sym0_tot:sys%sym_max_tot, 2), stat=ierr)
        call check_allocate('current_index',sys%nsym,ierr)
        current_index = 1

        associate(sym_spin_basis_fns=>sys%read_in%pg_sym%sym_spin_basis_fns, basis_fns=>sys%basis%basis_fns)
            do i = 1, sys%basis%nbasis
                ims = (basis_fns(i)%Ms+3)/2
                sym_spin_basis_fns(current_index(basis_fns(i)%sym, ims), ims, basis_fns(i)%sym) = i
                basis_fns(i)%sym_spin_index = current_index(basis_fns(i)%sym, ims)
                current_index(basis_fns(i)%sym, ims) = current_index(basis_fns(i)%sym, ims) + 1
            end do
        end associate

        deallocate(current_index, stat=ierr)
        call check_deallocate('current_index',ierr)

        sys%read_in%cross_product_sym_ptr => cross_product_periodic_read_in
        sys%read_in%sym_conj_ptr => mom_sym_conj

    end subroutine init_basis_momentum_symmetry_info

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

    pure function mom_sym_conj(read_in, sym) result(rsym)

        ! Returns symmetry index of complex conjugate of provided
        ! sym index.
        ! Since using complex plane waves, e^(ik.r), this will in
        ! general be e^(-ik.r) so just return inverse symmetry.

        ! In:
        !   sys: information about system being studied. We use the
        !       momentum symmetry information.
        !   sym: symmetry to compare.
        ! Returns:
        !   Symmetry of complex conjugate of sym.

        ! For momentum symmetry in real (read in from an FCIDUMP), periodic
        ! systems.

        use system, only: sys_read_in_t

        type(sys_read_in_t), intent(in) :: read_in
        integer, intent(in) :: sym
        integer :: rsym

        rsym = read_in%mom_sym%inv_sym(sym)

    end function mom_sym_conj

!--- Cross products ---

    pure function cross_product_periodic_read_in(read_in, a1, a2) result(prod)

        ! Return the index of the symmetry resulting from the cross
        ! product of the symmetries of indexes a1 and a2

        ! In:
        !   mom_sym: basis function symmetry information.
        !   a1, a2: symmetry indexes to return cross product of.
        ! Returns:
        !   Index of direct product of symmetries a1 & a2.

        use system, only: sys_read_in_t

        type(sys_read_in_t), intent(in) :: read_in
        integer, intent(in) :: a1, a2
        integer :: prod

        prod = read_in%mom_sym%sym_table(a1, a2)

    end function cross_product_periodic_read_in

!--- Indexing conversion routines: ---

    pure subroutine decompose_trans_sym(isym, propbitlen, abel_sym)

        ! Takes symmetry index for translationally symmetric wavefunction and
        ! returns representation of three "quantum numbers". In accordance
        ! with approach used in NECI, values stored according to:
        !   isym = 1 + \sum_i sym(i) * 2 ** (propbitlen * (i-1))

        ! For an example of this in practice see symmetry_types.f90

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

    end subroutine decompose_trans_sym

    pure function get_kpoint_index(a, nprop) result(ind)

        ! Converts from abelian symmetry quantum numbers into unique index.
        ! If we know size of unit cell, can calculate unique index by tiling
        ! first along axis 1, then 2, then 3 in 3 dimensions.

        ! By defining a unique mapping between the kpoint grid and a contiguous
        ! indexes we can refer to symmetry via only this index and arrays
        ! specifying the resultant products and conjugates/inverses. This
        ! enables easy compatibility with much of the pg_sym functionality
        ! once these arrays are constructed.

        ! As such this functonality is mainly for use when initialising
        ! symmetry information within read_in%mom_sym.

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

    pure subroutine get_kpoint_vector(ind, nprop, a)

        ! Get kpoint index, in 3D array, from index defined in get_kpoint_index.

        ! this functonality is mainly for use when initialising symmetry
        ! information within read_in%mom_sym.

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

    end subroutine get_kpoint_vector

end module momentum_sym_read_in

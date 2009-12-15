module basis

! Basis function information.

use kpoints

implicit none

! The kpoint type is used to specify a spin orbital in momentum space.
type basis_fn
    ! l is used in two different contexts depending upon whether a momentum
    ! space description or real space description is being used.
    ! Momentum space:
    !     l is the wavevector in terms of the reciprocal lattice vectors of the crystal cell.
    ! Real space:
    !     l is the position of the basis function within the crystal cell in
    !     units of the lattice vectors of the primitive unit cell.
    ! Obviously we should not convert between the two descritions within one
    ! calculation! ;-)
    integer, pointer :: l(:) => NULL()
    ! Spin of the electron (1 or -1).
    integer :: ms
    ! Kinetic energy.
    real(dp) :: kinetic
end type basis_fn

! Store of information about the (spin) basis functions of the system.
type(basis_fn), allocatable :: basis_fns(:) ! (nbasis)

! number of basis functions.  Equal to 2*number of sites as there are
! 2 spin orbitals per site.
integer :: nbasis

! The determinants are stored as a bit string.  Each element of an array is
! an integer of kind i0 (containing i0_length bits).
! (The bit type has just been deleted from the forthcoming F2008 standard, so we
! won't hold out breath until we can use bits directly......)
! basis_length is the length of the byte array necessary to contain a bit for
! each basis function, i.e. ceiling(nbasis/i0_length).
integer :: basis_length

! A determinant is stored in the array f(nbasis).  A basis function is occupied
! in the determinant if the relevant bit is set.  The relevant bit is given by
! bit_element, the element of the array which contains the bit corresponding to
! the basis function, and bit_position, which contains the position of the bit
! within the given element.  bit_lookup(:,i) gives the (/ bit_position,
! bit_element /) of the i-th basis function.
! Note fortran numbers bits starting from 0.
integer, allocatable :: bit_lookup(:,:) ! (2, nbasis)

! The reverse lookup to bit_lookup.
! basis_lookup(i,j) gives the basis function corresponding to
! the i-th bit in the j-th element of a determinant array.
integer, allocatable :: basis_lookup(:,:) ! (i0_length, basis_length)

contains

    pure subroutine init_basis_fn(b, l, ms)

        ! Initialise a variable of type basis_fn.
        ! In:
        !    k (optional): quantum numbers of the basis function.
        !                  Momentum space formulation:
        !                      wavevector in units of the reciprocal lattice vectors
        !                      of the crystal cell.
        !                  Real space formulation:
        !                      position of basis function within the crystal cell
        !                      in units of the primitive lattice vectors.
        !    ms (optional): set spin of an electron occupying the basis function.
        ! Out:
        !    b: initialsed basis function.  The wavevector and kinetic energy
        !      components are set if the k arguments is given and the ms component
        !      is set if the ms argument is given.  If no optional arguments are
        !      specified then a completely blank variable is returned.
        !
        ! This should be called even if l and ms are not specified so that the
        ! l component can be correctly allocated.

        use system, only: system_type, hub_real

        type(basis_fn), intent(out) :: b
        integer, intent(in), optional  :: l(ndim)
        integer, intent(in), optional  :: ms
        integer :: ierr

        if (.not.associated(b%l)) then
            allocate(b%l(ndim),stat=ierr)
        end if

        if (present(l)) then
            b%l = l
            if (system_type /= hub_real) then
                b%kinetic = calc_kinetic(l)
            else
                b%kinetic = 0.0_dp
            end if
        end if

        if (present(ms)) b%ms = ms

    end subroutine init_basis_fn

    subroutine write_basis_fn(b, iunit, new_line)

        use system, only: system_type, hub_real

        ! Print out information stored in b.
        ! In:
        !    b: basis_fn variable.
        !    iunit (optional): io unit to which the output is written.
        !        Default: 6 (stdout).
        !    new_line (optional): if true, then a new line is written at
        !        the end of the list of occupied orbitals.  Default: no
        !        new line.

        type(basis_fn), intent(in) :: b
        integer, intent(in), optional :: iunit
        logical, intent(in), optional :: new_line
        integer :: i, io

        if (present(iunit)) then
            io = iunit
        else
            io = 6
        end if

        write (io,'(1X,"(")', advance='no')
        write (io,'(i2)',advance='no') b%l(1)
        do i = 2,ndim
            write (io,'(",",i2)',advance='no') b%l(i)
        end do
        write (io,'(")")', advance='no')
        write (io,'(5X,i2)', advance='no') b%ms
        if (system_type /= hub_real) write (io,'(4X,f12.8)', advance='no') b%kinetic
        if (present(new_line)) then
            if (new_line) write (io,'()')
        end if

    end subroutine write_basis_fn

    pure function spin_symmetry(i, j) result(spin_match)

        ! In:
        !    i: index of a basis function
        !    j: index of a basis function
        ! Returns:
        !    true if i and j refer to basis functions of the same spin.

        logical :: spin_match
        integer, intent(in) :: i, j

        spin_match = basis_fns(i)%ms == basis_fns(j)%ms

    end function spin_symmetry

end module basis

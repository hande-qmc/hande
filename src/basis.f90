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
    real(p) :: kinetic
end type basis_fn

! Store of information about the (spin) basis functions of the system.
! The *odd* indices contain the alpha (spin up) functions.  This is in
! contrast to the bit strings used to refer to determinants where the *even*
! bits refer to alpha (spin up) functions.  This difference arises because 
! fortran numbers bits from 0...
type(basis_fn), allocatable :: basis_fns(:) ! (nbasis)

! number of basis functions.  Equal to 2*number of sites as there are
! 2 spin orbitals per site.
integer :: nbasis

! The determinants are stored as a bit string.  Each element of an array is
! an integer of kind i0 (containing i0_length bits).
! (The bit type has just been deleted from the forthcoming F2008 standard, so we
! won't hold our breath until we can use bits directly......)
! basis_length is the length of the byte array necessary to contain a bit for
! each basis function, i.e. ceiling(nbasis/i0_length).
! If separate_strings is true, then we actually store the alpha and beta
! strings separately, and so basis_length is 2*ceiling(nbasis/(2*i0_length)).
integer :: basis_length

! All bits in the determinant bit array correspond to a basis function apart
! from the last element in the bit array (which can contain some excess).
! last_basis_ind is the index of the last basis function in the last element of
! the bit array. 
integer :: last_basis_ind

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

    subroutine init_basis_fn(b, l, ms)

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

        use checking, only: check_allocate
        use system, only: system_type, hub_real

        type(basis_fn), intent(out) :: b
        integer, intent(in), optional  :: l(ndim)
        integer, intent(in), optional  :: ms
        integer :: ierr

        if (.not.associated(b%l)) then
            allocate(b%l(ndim),stat=ierr)
            call check_allocate('b%l',ndim,ierr)
        end if

        if (present(l)) then
            b%l = l
            if (system_type /= hub_real) then
                b%kinetic = calc_kinetic(l)
            else
                b%kinetic = 0.0_p
            end if
        end if

        if (present(ms)) b%ms = ms

    end subroutine init_basis_fn

    subroutine write_basis_fn(b, iunit, new_line, print_full)

        ! Print out information stored in b.
        ! In:
        !    b: basis_fn variable.
        !    iunit (optional): io unit to which the output is written.
        !        Default: 6 (stdout).
        !    new_line (optional): if true, then a new line is written at
        !        the end of the list of occupied orbitals.  Default: no
        !        new line.
        !    print_full (optional): if true (default) then the quantum numbers,
        !         spin and (for momentum space models) kinetic energy associated
        !         with the basis function are printed.  If false, only the
        !         quantum numbers are printed.

        use system, only: system_type, hub_real

        type(basis_fn), intent(in) :: b
        integer, intent(in), optional :: iunit
        logical, intent(in), optional :: new_line
        logical, intent(in), optional :: print_full
        logical :: print_all
        integer :: i, io

        if (present(iunit)) then
            io = iunit
        else
            io = 6
        end if

        if (present(print_full)) then
            print_all = print_full
        else
            print_all = .true.
        end if

        write (io,'(1X,"(")', advance='no')
        write (io,'(i2)',advance='no') b%l(1)
        do i = 2,ndim
            write (io,'(",",i2)',advance='no') b%l(i)
        end do
        write (io,'(")")', advance='no')
        if (print_all) then
            write (io,'(5X,i2)', advance='no') b%ms
            if (system_type /= hub_real) write (io,'(4X,f12.8)', advance='no') b%kinetic
        end if
        if (present(new_line)) then
            if (new_line) write (io,'()')
        end if

    end subroutine write_basis_fn

    subroutine init_basis_fns()

        ! Produce the basis functions.  The number of wavevectors is
        ! equal to the number of sites in the crystal cell (ie the number
        ! of k-points used to sample the FBZ of the primitive cell).
        ! From the cell parameters and the "tilt" used (if any) generate
        ! the list of wavevectors and hence the kinetic energy associated
        ! with each basis function (two per wavevector to account for spin).

        use checking, only: check_allocate, check_deallocate
        use m_mrgref, only: mrgref
        use errors, only: stop_all
        use parallel, only: parent

        integer :: limits(3,3), nmax(3), kp(3) ! Support a maximum of 3 dimensions.
        integer :: i, j, k, ibasis, ierr
        type(basis_fn), allocatable, target :: tmp_basis_fns(:)
        type(basis_fn), pointer :: basis_fn_p
        integer, allocatable :: basis_fns_ranking(:)

        nbasis = 2*nsites

        ! Find basis functions.

        ! We use a minimal basis: the hubbard model consisting of two
        ! spin-orbitals per lattice site.

        ! In the momentum space formulation the basis functions consist of a
        ! set of wavevectors/k-points that lie within the first Brillouin zone.

        ! In the real space formulation the basis functions used are those
        ! residing at the lattice points: we just need to find which lattice
        ! points fall within the crystal cell.

        ! Momentum space:

        ! Fold the crystal cell into the FBZ.
        ! The k-points must be integer multiples of the reciprocal lattice
        ! vectors of the crystal cell (so that the wavefunction is periodic in
        ! the crystal cell) and fall within the first Brillouin zone of the
        ! primitive unit cell (so that a unique set of k-points are chosen).
        ! The volume of the FBZ is inversely proportional to the volume of the
        ! cell, and so the number of sites in the crystal cell is equal to the
        ! number of reciprocal crystal cells in the FBZ of the unit cell and 
        ! hence this gives the required number of wavevectors.

        ! Real space:

        ! The same procedure applies as for the momentum space: we find which
        ! lattice sites lie within the Wigner--Seitz cell.  In fact, due to the
        ! relationship between reciprocal space and real space (and due to how
        ! we store the wavevectors in terms of the reciprocal lattice vectors of
        ! the crystal cell), *exactly* the same approach is needed, so we're
        ! just going to abuse the same code. Shocking, I know.

        ! Maximum limits...
        ! [Does it show that I've been writing a lot of python recently?]
        nmax = 0 ! Set nmax(i) to be 0 for unused higher dimensions.
        limits = 0
        ! forall is a poor substitute for list comprehension. ;-)
        forall (i=1:ndim)
            forall (j=1:ndim, lattice(i,j) /= 0) 
                limits(i,j) = abs(nint(box_length(i)**2/(2*lattice(i,j))))
            end forall
            nmax(i) = maxval(limits(:,i))
        end forall

        allocate(basis_fns(nbasis), stat=ierr)
        call check_allocate('basis_fns',nbasis,ierr)
        allocate(tmp_basis_fns(nbasis/2), stat=ierr)
        call check_allocate('tmp_basis_fns',nbasis/2,ierr)
        allocate(basis_fns_ranking(nbasis/2), stat=ierr)
        call check_allocate('basis_fns_ranking',nbasis/2,ierr)

        ! Find all alpha spin orbitals.
        ibasis = 0
        do k = -nmax(3), nmax(3)
            do j = -nmax(2), nmax(2)
                do i = -nmax(1), nmax(1)
                    ! kp is the Wavevector in terms of the reciprocal lattice vectors of
                    ! the crystal cell.
                    kp = (/ i, j, k /)
                    if (in_FBZ(kp(1:ndim))) then
                        if (ibasis==nbasis) then
                            call stop_all('init_basis_fns','Too many basis functions found.')
                        else
                            ! Have found an allowed wavevector/site.
                            ! Add 2 spin orbitals to the set of the basis functions.
                            ibasis = ibasis + 1
                            call init_basis_fn(tmp_basis_fns(ibasis), kp(1:ndim), 1)
                        end if
                    end if
                end do
            end do
        end do

        if (ibasis /= nbasis/2) call stop_all('init_basis_fns','Not enough basis functions found.')

        ! Rank by kinetic energy (applies to momentum space formulation only).
        select case(system_type)
        case(hub_k)
            call mrgref(tmp_basis_fns(:)%kinetic, basis_fns_ranking)
        case(hub_real)
            forall (i=1:nbasis/2) basis_fns_ranking(i) = i
        end select

        ! Form the list of sorted basis functions with both alpha and beta
        ! spins.
        do i = 1, nbasis/2
            ! Can't set a kpoint equal to another kpoint as then the k pointers
            ! can be assigned whereas we want to *copy* the values.
            basis_fn_p => tmp_basis_fns(basis_fns_ranking(i))
            call init_basis_fn(basis_fns(2*i-1), basis_fn_p%l, basis_fn_p%ms)
            call init_basis_fn(basis_fns(2*i), basis_fn_p%l, -basis_fn_p%ms)
            deallocate(tmp_basis_fns(basis_fns_ranking(i))%l, stat=ierr)
            call check_deallocate('tmp_basis_fns(basis_fns_ranking(i',ierr)
        end do
        deallocate(tmp_basis_fns, stat=ierr)
        call check_deallocate('tmp_basis_fns',ierr)
        deallocate(basis_fns_ranking, stat=ierr)
        call check_deallocate('basis_fns_ranking',ierr)

        if (parent) then
            write (6,'(1X,a15,/,1X,15("-"),/)') 'Basis functions'
            write (6,'(1X,a27)') 'Spin given in units of 1/2.'
            if (system_type == hub_real) then
                write (6,'(1X,a63,/)') 'Site positions given in terms of the primitive lattice vectors.'
                write (6,'(1X,a5,3X,a4,3X)', advance='no') 'index','site'
            else
                write (6,'(1X,a78)') 'k-points given in terms of the reciprocal lattice vectors of the crystal cell.'
                if (any(abs(ktwist) > 0.0_p)) then
                    write (6,'(1X,a26)', advance='no') 'Applying a twist angle of:'
                    write (6,'(1X,"(",f6.4)', advance='no') ktwist(1)
                    do i = 2, ndim
                        write (6,'(",",f6.4)', advance='no') ktwist(i)
                    end do
                    write (6,'(").")')
                end if
                write (6,'()')
                write (6,'(1X,a5,3X,a7)', advance='no') 'index','k-point'
            end if
            do i = 1, ndim
                write (6,'(3X)', advance='no')
            end do
            write (6,'(a2)', advance='no') 'ms'
            if (system_type == hub_real) then
                write(6,'()')
            else
                write(6,'(5X,a14)') 'kinetic energy'
            end if
            do i = 1, nbasis
                write (6,'(1X,i5,2X)',advance='no') i
                call write_basis_fn(basis_fns(i), new_line=.true.)
            end do
            write (6,'()')
        end if

    end subroutine init_basis_fns

    subroutine end_basis_fns()

        ! Clean up basis functions.

        use checking, only: check_deallocate

        integer :: ierr, i

        do i = 1, nbasis
            deallocate(basis_fns(i)%l, stat=ierr)
            call check_deallocate('basis_fns(i',ierr)
        end do
        deallocate(basis_fns, stat=ierr)
        call check_deallocate('basis_fns',ierr)

    end subroutine end_basis_fns

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

    subroutine set_orb(f,iorb)

        ! In:
        !    f: bit string of orbitals.
        !    iorb: orbital index.
        ! Out:
        !    f: bit string of orbitals with the bit corresponding to iorb set. 
        
        ! Note that f must be zerod before first using this procedure.

        integer, intent(in) :: iorb
        integer(i0), intent(inout) :: f(basis_length)
        integer :: pos, ind

        pos = bit_lookup(1,iorb)
        ind = bit_lookup(2,iorb)
        f(ind) = ibset(f(ind),pos)

    end subroutine set_orb

    subroutine set_orb_mask(lmag2, orb_mask)

        ! Set a mask with bits set for symmetry-related orbitals.

        ! In:
        !    lmag2: magnitude squared of the l quantum vector (component of the
        !      basis_fn type) which corresponds to the desired set of
        !      symmetry-related orbitals.
        ! Out:
        !    orb_mask: bit-string where only bits are set that correspond to the
        !      set of symmetry-related orbitals with input value of lmag2.

        integer, intent(in) :: lmag2
        integer(i0), intent(out) :: orb_mask(basis_length)

        integer :: i, ipos, iel

        orb_mask = 0

        do i = 1, nbasis
            if (dot_product(basis_fns(i)%l,basis_fns(i)%l) == lmag2) then
                ipos = bit_lookup(1,i)
                iel = bit_lookup(2,i)
                orb_mask(iel) = ibset(orb_mask(iel), ipos)
            end if
        end do

    end subroutine set_orb_mask

end module basis

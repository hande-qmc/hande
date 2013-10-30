module basis

! Basis function information.

use kpoints

implicit none

! Define a spin orbital.
type basis_fn
    ! Set of quantum numbers describing the basis function.
    ! l is used in two different contexts depending upon whether the orbitals
    ! are defined in momentum space or in real space.  Applies only to model
    ! Hamiltonians (e.g. Hubbard model).
    ! Momentum space:
    !     l is the wavevector in terms of the reciprocal lattice vectors of the crystal cell.
    ! Real space:
    !     l is the position of the basis function within the crystal cell in
    !     units of the lattice vectors of the primitive unit cell.
    ! Obviously we should not convert between descriptions within one
    ! calculation! ;-)
    integer, pointer :: l(:) => NULL()
    integer :: spatial_index
    ! Index of the irreducible representation spanned by the orbital.  Used only
    ! in systems where point group symmetry is used (e.g.  molecules).  See
    ! notes in pg_symmetry.
    integer :: sym = 0
    ! Index of basis function within the symmetry block.  sym_index = n
    ! indicates the basis function is the fifth in basis_fns array to have the
    ! symmetry given by sym.
    ! Used only with point_group symmetry.
    integer :: sym_index = 0
    ! Index of basis function within the symmetry block.  sym_spin_index = n
    ! indicates the basis function is the fifth in basis_fns array to have the
    ! symmetry given by sym *and* with the spin given by ms.
    ! Used only with point_group symmetry.
    integer :: sym_spin_index = 0
    ! Spin of the electron (1 or -1).
    integer :: ms
    ! single-particle energy of basis function.
    ! model Hamiltonians in momentum space:
    !     sp_eigv is the kinetic energy of the basis function.
    ! model Hamiltonians in real space:
    !     sp_eigv is not set/used.
    ! molecular systems:
    !     sp_eigv is the single-particle energy read in from the FCIDUMP file
    !     (e.g. Hartree--Fock or Kohn--Sham eigenvalue).
    real(p) :: sp_eigv
end type basis_fn

! Store of information about the (spin) basis functions of the system.
! The *odd* indices contain the alpha (spin up) functions.  This is in
! contrast to the bit strings used to refer to determinants where the *even*
! bits refer to alpha (spin up) functions.  This difference arises because
! fortran numbers bits from 0...
type(basis_fn), allocatable :: basis_fns(:) ! (nbasis)

! number of basis functions.
! For the Hubbard model is equal to twice the number of sites as there are
! 2 spin orbitals per site, for the Heisenberg model to the number of sites,
! for UEG equal to twice the number of k-points within the energy cutoff and for
! read in (e.g. molecular) systems the number of single-particle states read in.
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

! DMQMC uses two determinants for each psip to refer to the two components
! of the relevant matrix element. Hence, the bitstring which is stored in DMQMC has
! 2*basis_length components. There are some procedures which required basis_length
! when used for stahndard FCIQMC but 2*basis_length when used for DMQMC. It is
! therefore useful to have a quantity which equal to 2*basis_length for DMQMC and
! equal to basis_length for other methods. Then a procedure can use this quantity
! and will work for both methods, making it general. This quantity is total_basis_length.
! total_basis_length can then be used when we want to refer to *both* determinants
! in DMQMC, and hence the entire bitstring.
integer :: total_basis_length

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
integer, allocatable :: basis_lookup(:,:) ! (0:i0_end, basis_length)

contains

    subroutine init_basis_fn(b, l, sym, ms)

        ! Initialise a variable of type basis_fn.
        ! In:
        !    l (optional): quantum numbers of the basis function.  Used only in
        !        model Hamiltonians.
        !        Momentum space formulation:
        !            wavevector in units of the reciprocal lattice vectors of
        !            the crystal cell.
        !        Real space formulation:
        !            position of basis function within the crystal cell in units
        !            of the primitive lattice vectors.
        !    sym (optional): symmetry label of the basis function.  Used only in
        !        systems with point group symmetry (i.e. read in from an FCIDUMP
        !        file).
        !    ms (optional): set spin of an electron occupying the basis function.
        ! Out:
        !    b: initialsed basis function.  The wavevector and (if appropriate
        !      to the system) single-particle eigenvalue components are set if
        !      the l arguments is given and the ms component is set if the ms
        !      argument is given.  If no optional arguments are specified then
        !      a completely blank variable is returned.
        !
        ! This should be called even if l and ms are not specified so that the
        ! l component can be correctly allocated.

        use checking, only: check_allocate
        use system, only: system_type, hub_real

        type(basis_fn), intent(out) :: b
        integer, intent(in), optional  :: l(ndim)
        integer, intent(in), optional  :: sym, ms
        integer :: ierr

        if (.not.associated(b%l)) then
            allocate(b%l(ndim),stat=ierr)
            call check_allocate('b%l',ndim,ierr)
        end if

        if (present(l)) then
            b%l = l
            if (system_type == hub_k .or. system_type == ueg) then
                b%sp_eigv = calc_kinetic(l)
            else
                b%sp_eigv = 0.0_p
            end if
        end if

        if (present(sym)) b%sym = sym

        if (present(ms)) b%ms = ms

    end subroutine init_basis_fn

    subroutine write_basis_fn_header(iunit, print_full)

        ! Print out header for a table of basis functions.
        ! Format in line with write_basis_fn.
        !
        ! In:
        !    iunit (optional): io unit to which the output is written.
        !        Default: 6 (stdout).
        !    print_full (optional): if true (default) then print out header info
        !        for the symmetry and spin quantum numbers and (if appropriate)
        !        single-particle energy associated with the basis function.
        !        If false, only information about the quantum numbers is
        !        printed.

        use system, only: system_type, read_in, hub_k, hub_real, heisenberg

        integer, intent(in), optional :: iunit
        logical, intent(in), optional :: print_full

        integer :: io, i
        logical :: print_long

        if (present(iunit)) then
            io = iunit
        else
            io = 6
        end if

        ! If print_full is false, then the spin and single-particle eigenvalues
        ! are also printed out.
        if (present(print_full)) then
            print_long = print_full
        else
            print_long = .true.
        end if

        ! Title
        write (6,'(1X,a15,/,1X,15("-"),/)') 'Basis functions'

        ! Describe information.
        if (system_type /= heisenberg) write (6,'(1X,a27)') 'Spin given in units of 1/2.'

        select case(system_type)
        case(hub_real,heisenberg, chung_landau)
            write (6,'(1X,a63,/)') 'Site positions given in terms of the primitive lattice vectors.'
            write (6,'(1X,a5,3X,a4,3X)', advance='no') 'index','site'
        case(hub_k,ueg)
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
        case(read_in)
            write (6,'(/,1X,a5,2X,a7,X,a8,X,a9,2X)', advance='no') 'index','spatial','symmetry','sym_index'
        end select

        if (system_type /= read_in) then
            do i = 1, ndim
                write (6,'(4X)', advance='no')
            end do
        end if

        if (print_long) then
            if (system_type /= heisenberg .and. system_type /= chung_landau) &
                write (6,'(a2)', advance='no') 'ms'

            select case(system_type)
            case(hub_real, heisenberg, chung_landau)
                write(6,'()')
            case default
                write(6,'(5X,a7)') '<i|h|i>'
            end select
        else
            write (6,'()')
        end if

    end subroutine write_basis_fn_header

    subroutine write_basis_fn(b, ind, iunit, new_line, print_full)

        ! Print out information stored in b.
        ! Format in line with write_basis_fn_header.
        ! Please ensure formats are changed in both write_basis_fn and
        ! write_basis_fn_header.
        !
        ! In:
        !    b: basis_fn variable.
        !    ind: index of basis function.  Only printed out if present and
        !        positive.
        !    iunit (optional): io unit to which the output is written.
        !        Default: 6 (stdout).
        !    new_line (optional): if true, then a new line is written at
        !        the end of the list of occupied orbitals.  Default: no
        !        new line.
        !    print_full (optional): if true (default) then the symmetry and spin
        !        quantum numbers and (if appropriate) single-particle energy
        !        associated with the basis function are printed.  If false, only
        !        the quantum numbers are printed.

        use system, only: system_type, read_in, hub_k, ueg, heisenberg, chung_landau

        type(basis_fn), intent(in) :: b
        integer, intent(in), optional :: ind
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

        if (present(ind)) then
            if (ind >= 0) write (6,'(1X,i5,2X)',advance='no') ind
        end if

        if (system_type == read_in) then
            write (io, '(i5,2(3X,i5),X)',advance='no') b%spatial_index, b%sym, b%sym_index
        else
            write (io,'(1X,"(")', advance='no')
            write (io,'(i3)',advance='no') b%l(1)
            do i = 2,ndim
                write (io,'(",",i3)',advance='no') b%l(i)
            end do
            write (io,'(")")', advance='no')
        end if
        if (print_all) then
            select case (system_type)
            case(heisenberg, chung_landau)
            case default
                write (io,'(5X,i2)', advance='no') b%ms
            end select
            select case (system_type)
            case(heisenberg, chung_landau, hub_real)
            case default
                write (io,'(4X,f12.8)', advance='no') b%sp_eigv
            end select
        end if
        if (present(new_line)) then
            if (new_line) write (io,'()')
        end if

    end subroutine write_basis_fn

    subroutine init_model_basis_fns(store_info)

        ! Produce the basis functions for model Hamiltonian systems.
        !
        ! In:
        !    store_info (optional): if true (default) then store the data read
        !    in.  Otherwise the basis set is simply printed out.
        !
        ! The number of wavevectors is equal to the number of sites in the
        ! crystal cell (ie the number of k-points used to sample the FBZ of the
        ! primitive cell).  From the cell parameters and the "tilt" used (if
        ! any) generate the list of wavevectors and hence the kinetic energy
        ! associated with each basis function (two per wavevector to account for
        ! spin).

        use checking, only: check_allocate, check_deallocate
        use system
        use ranking, only: insertion_rank_rp
        use errors, only: stop_all
        use parallel, only: parent

        logical, intent(in), optional :: store_info

        logical :: t_store

        integer :: limits(3,3), nmax(3), kp(3) ! Support a maximum of 3 dimensions.
        integer :: i, j, k, ibasis, ierr, nspatial
        type(basis_fn), allocatable, target :: tmp_basis_fns(:)
        type(basis_fn), pointer :: basis_fn_p
        integer, allocatable :: basis_fns_ranking(:)

        if (present(store_info)) then
            t_store = store_info
        else
            t_store = .true.
        end if

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

        ! UEG:

        ! This is identical again to the real space formulation, except the FBZ
        ! is essentially infinite (as there is no underlying crystal lattice).

        ! For the Heisenberg model, each site has a single spin which must be either
        ! up or down, so only need 1 bit for each site => nbasis = nsites
        ! For the fermionic systems, each spatial has an alpha and beta spin orbital,
        ! each of which could be occupied or unoccupied.
        ! Initially we will only consider the spatial orbitals when constructing
        ! the basis and will expand it out to spin orbitals later.
        select case(system_type)
        case(hub_k, hub_real, heisenberg, chung_landau)
            nspatial = nsites
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
        case(ueg)
            ! UEG has 2 spin-orbitals per wavevector.  We include all
            ! wavevectors up to an energy cutoff, so we won't precisely count
            ! the number of basis functions right now.
            ! We actually generate all wavevectors in the smallest
            ! line/square/cube which encloses all wavevectors within the cutoff
            ! and discard those outside the cutoff.
            nmax = 0
            forall (i=1:ndim) nmax(i) = ceiling(sqrt(2*ueg_ecutoff))
            nspatial = (2*nmax(1)+1)**ndim
        end select

        allocate(tmp_basis_fns(nspatial), stat=ierr)
        call check_allocate('tmp_basis_fns',nspatial,ierr)

        ! Find all alpha spin orbitals.
        ibasis = 0
        do k = -nmax(3), nmax(3)
            do j = -nmax(2), nmax(2)
                do i = -nmax(1), nmax(1)
                    ! kp is the Wavevector in terms of the reciprocal lattice vectors of
                    ! the crystal cell.
                    kp = (/ i, j, k /)
                    if (in_FBZ(kp(1:ndim))) then
                        if (ibasis == nspatial) then
                            call stop_all('init_basis_fns','Too many basis functions found.')
                        else
                            ! Have found an allowed wavevector/site.
                            ! Add 2 spin orbitals to the set of the basis functions.
                            ibasis = ibasis + 1
                            call init_basis_fn(tmp_basis_fns(ibasis), l=kp(1:ndim), ms=1)
                            if (system_type==ueg .and. real(dot_product(kp,kp),p)/2 > ueg_ecutoff) then
                                ! Have found a wavevector with too large KE.
                                ! Discard.
                                ! Note that we don't use the calculated kinetic
                                ! energy as it's in a.u. (ueg_ecutoff is in
                                ! scaled units) and includes any twist.
                                ! Avoid the chance of having allocated
                                ! additional %l elements (eg if rejecting the
                                ! final basis function tested).
                                deallocate(tmp_basis_fns(ibasis)%l, stat=ierr)
                                call check_deallocate('tmp_basis_fns(basis_fns_ranking(i',ierr)
                                ibasis = ibasis - 1
                            end if
                        end if
                    end if
                end do
            end do
        end do

        select case(system_type)
        case(hub_k, hub_real, heisenberg, chung_landau)
            if (ibasis /= nspatial) call stop_all('init_basis_fns','Not enough basis functions found.')
        case(ueg)
        end select

        ! Convert nbasis to being in terms of spin-orbitals if applicable.
        select case(system_type)
        case(ueg)
            ! Set nbasis to be the number of basis functions found to be within
            ! the energy cutoff.
            ! Yes, I know this could be evaluated as one knows the 'volume'
            ! occupied by each wavevector, but I just cba.
            nbasis = 2*ibasis
            nspatial = ibasis
            nvirt = nbasis - nel
        case(hub_k, hub_real)
            ! nvirt set in init_system
            nbasis = 2*nspatial
        case(heisenberg, chung_landau)
            ! nvirt set in init_system
            nbasis = nspatial
        end select

        allocate(basis_fns_ranking(nspatial), stat=ierr)
        call check_allocate('basis_fns_ranking',nspatial,ierr)

        ! Rank by kinetic energy (applies to momentum space basis sets only).
        select case(system_type)
        case(hub_k, ueg)
            call insertion_rank_rp(tmp_basis_fns(:nspatial)%sp_eigv, basis_fns_ranking, tolerance=depsilon)
        case(hub_real, heisenberg, chung_landau)
            forall (i=1:nsites) basis_fns_ranking(i) = i
        end select

        allocate(basis_fns(nbasis), stat=ierr)
        call check_allocate('basis_fns',nbasis,ierr)

        ! Form the list of sorted basis functions with both alpha and beta
        ! spins.
        do i = 1, nspatial
            ! Can't set a kpoint equal to another kpoint as then the k pointers
            ! can be assigned whereas we want to *copy* the values.
            basis_fn_p => tmp_basis_fns(basis_fns_ranking(i))
            select case(system_type)
            case(heisenberg, chung_landau)
                call init_basis_fn(basis_fns(i), l=basis_fn_p%l)
            case default
                call init_basis_fn(basis_fns(2*i-1), l=basis_fn_p%l, ms=basis_fn_p%ms)
                call init_basis_fn(basis_fns(2*i), l=basis_fn_p%l, ms=-basis_fn_p%ms)
            end select
            deallocate(tmp_basis_fns(basis_fns_ranking(i))%l, stat=ierr)
            call check_deallocate('tmp_basis_fns(basis_fns_ranking(i',ierr)
        end do
        deallocate(tmp_basis_fns, stat=ierr)
        call check_deallocate('tmp_basis_fns',ierr)
        deallocate(basis_fns_ranking, stat=ierr)
        call check_deallocate('basis_fns_ranking',ierr)

        if (parent) then
            call write_basis_fn_header()
            do i = 1, nbasis
                call write_basis_fn(basis_fns(i), ind=i, new_line=.true.)
            end do
            write (6,'()')
        end if

        if (.not.t_store) then
            ! Should tidy up and deallocate everything we allocated.
            call end_basis_fns()
        end if

    end subroutine init_model_basis_fns

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

    subroutine end_basis_fns()

        ! Clean up basis functions.

        use checking, only: check_deallocate

        integer :: ierr, i

        if (allocated(basis_fns)) then
            do i = 1, nbasis
                deallocate(basis_fns(i)%l, stat=ierr)
                call check_deallocate('basis_fns(i',ierr)
            end do
            deallocate(basis_fns, stat=ierr)
            call check_deallocate('basis_fns',ierr)
        end if

    end subroutine end_basis_fns

end module basis

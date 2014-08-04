module basis

! Basis function information.

use const, only: p, i0
use basis_types

implicit none

! Temporary global object to aid removal of global state.
type(basis_t) :: basis_global
 
contains

    subroutine init_basis_fn(sys, b, l, sym, lz, ms)

        ! Initialise a variable of type basis_fn_t.
        ! In:
        !    sys: system being studied.
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
        !    lz  (optional): lz of the basis function.  Used only in
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

        use kpoints, only: calc_kinetic
        use system

        type(sys_t), intent(in) :: sys
        type(basis_fn_t), intent(out) :: b
        integer, intent(in), optional  :: l(sys%lattice%ndim)
        integer, intent(in), optional  :: sym, ms, lz
        integer :: ierr

        if (.not.allocated(b%l)) then
            allocate(b%l(sys%lattice%ndim),stat=ierr)
            call check_allocate('b%l',sys%lattice%ndim,ierr)
        end if

        if (present(l)) then
            b%l = l
            if (sys%system == hub_k .or. sys%system == ueg) then
                b%sp_eigv = calc_kinetic(sys, l)
            else
                b%sp_eigv = 0.0_p
            end if
        end if

        if (present(sym)) b%sym = sym

        if (present(lz)) b%lz = lz
        if (present(ms)) b%ms = ms

    end subroutine init_basis_fn

    subroutine write_basis_fn_header(sys, iunit, print_full)

        ! Print out header for a table of basis functions.
        ! Format in line with write_basis_fn.
        !
        ! In:
        !    sys: system being studied.
        !    iunit (optional): io unit to which the output is written.
        !        Default: 6 (stdout).
        !    print_full (optional): if true (default) then print out header info
        !        for the symmetry and spin quantum numbers and (if appropriate)
        !        single-particle energy associated with the basis function.
        !        If false, only information about the quantum numbers is
        !        printed.

        use system

        type(sys_t), intent(in) :: sys
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
        if (sys%system /= heisenberg) write (6,'(1X,a27)') 'Spin given in units of 1/2.'

        select case(sys%system)
        case(hub_real,heisenberg, chung_landau)
            write (6,'(1X,a63,/)') 'Site positions given in terms of the primitive lattice vectors.'
            write (6,'(1X,a5,3X,a4,3X)', advance='no') 'index','site'
        case(hub_k,ueg)
            write (6,'(1X,a78)') 'k-points given in terms of the reciprocal lattice vectors of the crystal cell.'
            if (any(abs(sys%lattice%ktwist) > 0.0_p)) then
                write (6,'(1X,a26)', advance='no') 'Applying a twist angle of:'
                write (6,'(1X,"(",f6.4)', advance='no') sys%lattice%ktwist(1)
                do i = 2, sys%lattice%ndim
                    write (6,'(",",f6.4)', advance='no') sys%lattice%ktwist(i)
                end do
                write (6,'(").")')
            end if
            write (6,'()')
            write (6,'(1X,a5,3X,a7)', advance='no') 'index','k-point'
        case(read_in)
            write (6,'(/,1X,a5,2X,a7,X,a8,X,a9,X,a2,5X)', advance='no') 'index','spatial','symmetry','sym_index','lz'
        end select

        if (sys%system /= read_in) then
            do i = 1, sys%lattice%ndim
                write (6,'(4X)', advance='no')
            end do
        end if

        if (print_long) then
            if (sys%system /= heisenberg .and. sys%system /= chung_landau) &
                write (6,'(a2)', advance='no') 'ms'

            select case(sys%system)
            case(hub_real, heisenberg, chung_landau)
                write(6,'()')
            case default
                write(6,'(5X,a7)') '<i|h|i>'
            end select
        else
            write (6,'()')
        end if

    end subroutine write_basis_fn_header

    subroutine write_basis_fn(sys, b, ind, iunit, new_line, print_full)

        ! Print out information stored in b.
        ! Format in line with write_basis_fn_header.
        ! Please ensure formats are changed in both write_basis_fn and
        ! write_basis_fn_header.
        !
        ! In:
        !    sys: system being studied.
        !    b: basis_fn_t variable.
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

        use system

        type(sys_t), intent(in) :: sys
        type(basis_fn_t), intent(in) :: b
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

        if (sys%system == read_in) then
            write (io, '(i5,3(3X,i5),X)',advance='no') b%spatial_index, b%sym, b%sym_index,b%lz
        else
            write (io,'(1X,"(")', advance='no')
            write (io,'(i3)',advance='no') b%l(1)
            do i = 2,sys%lattice%ndim
                write (io,'(",",i3)',advance='no') b%l(i)
            end do
            write (io,'(")")', advance='no')
        end if
        if (print_all) then
            select case (sys%system)
            case(heisenberg, chung_landau)
            case default
                write (io,'(5X,i2)', advance='no') b%ms
            end select
            select case (sys%system)
            case(heisenberg, chung_landau, hub_real)
            case default
                write (io,'(4X,f12.8)', advance='no') b%sp_eigv
            end select
        end if
        if (present(new_line)) then
            if (new_line) write (io,'()')
        end if

    end subroutine write_basis_fn

    subroutine init_model_basis_fns(sys, store_info)

        ! Produce the basis functions for model Hamiltonian systems.
        !
        ! In/Out:
        !    sys: system being studied.  On output the nvirt component is set
        !         (UEG *only*; unmodified for all other systems).
        ! In:
        !    store_info (optional): if true (default) then store the data read
        !         in.  Otherwise the basis set is simply printed out.
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

        type(sys_t), intent(inout) :: sys
        logical, intent(in), optional :: store_info

        logical :: t_store

        integer :: limits(3,3), nmax(3), kp(3) ! Support a maximum of 3 dimensions.
        integer :: i, j, k, ibasis, ierr, nspatial
        type(basis_fn_t), allocatable, target :: tmp_basis_fns(:)
        type(basis_fn_t), pointer :: basis_fn_p
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
        select case(sys%system)
        case(hub_k, hub_real, heisenberg, chung_landau)
            nspatial = sys%lattice%nsites
            ! Maximum limits...
            ! [Does it show that I've been writing a lot of python recently?]
            nmax = 0 ! Set nmax(i) to be 0 for unused higher dimensions.
            limits = 0
            ! forall is a poor substitute for list comprehension. ;-)
            forall (i=1:sys%lattice%ndim)
                forall (j=1:sys%lattice%ndim, sys%lattice%lattice(i,j) /= 0)
                    limits(i,j) = abs(nint(sys%lattice%box_length(i)**2/(2*sys%lattice%lattice(i,j))))
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
            forall (i=1:sys%lattice%ndim) nmax(i) = ceiling(sqrt(2*sys%ueg%ecutoff))
            nspatial = (2*nmax(1)+1)**sys%lattice%ndim
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
                    if (in_FBZ(sys%system, sys%lattice, kp(1:sys%lattice%ndim))) then
                        if (ibasis == nspatial) then
                            call stop_all('init_basis_fns','Too many basis functions found.')
                        else
                            ! Have found an allowed wavevector/site.
                            ! Add 2 spin orbitals to the set of the basis functions.
                            ibasis = ibasis + 1
                            call init_basis_fn(sys, tmp_basis_fns(ibasis), l=kp(1:sys%lattice%ndim), ms=1)
                            if (sys%system==ueg .and. real(dot_product(kp,kp),p)/2 > sys%ueg%ecutoff) then
                                ! Have found a wavevector with too large KE.
                                ! Discard.
                                ! Note that we don't use the calculated kinetic
                                ! energy as it's in a.u. (sys%ueg%ecutoff is in
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

        select case(sys%system)
        case(hub_k, hub_real, heisenberg, chung_landau)
            if (ibasis /= nspatial) call stop_all('init_basis_fns','Not enough basis functions found.')
        case(ueg)
        end select

        ! Convert nbasis to being in terms of spin-orbitals if applicable.
        select case(sys%system)
        case(ueg)
            ! Set nbasis to be the number of basis functions found to be within
            ! the energy cutoff.
            ! Yes, I know this could be evaluated as one knows the 'volume'
            ! occupied by each wavevector, but I just cba.
            basis_global%nbasis = 2*ibasis
            nspatial = ibasis
            sys%nvirt = basis_global%nbasis - sys%nel
        case(hub_k, hub_real)
            ! sys%nvirt set in init_system
            basis_global%nbasis = 2*nspatial
        case(heisenberg, chung_landau)
            ! sys%nvirt set in init_system
            basis_global%nbasis = nspatial
        end select

        allocate(basis_fns_ranking(nspatial), stat=ierr)
        call check_allocate('basis_fns_ranking',nspatial,ierr)

        ! Rank by kinetic energy (applies to momentum space basis sets only).
        select case(sys%system)
        case(hub_k, ueg)
            call insertion_rank_rp(tmp_basis_fns(:nspatial)%sp_eigv, basis_fns_ranking, tolerance=depsilon)
        case(hub_real, heisenberg, chung_landau)
            forall (i=1:sys%lattice%nsites) basis_fns_ranking(i) = i
        end select

        allocate(basis_global%basis_fns(basis_global%nbasis), stat=ierr)
        call check_allocate('basis_global%basis_fns',basis_global%nbasis,ierr)

        ! Form the list of sorted basis functions with both alpha and beta
        ! spins.
        do i = 1, nspatial
            ! Can't set a kpoint equal to another kpoint as then the k pointers
            ! can be assigned whereas we want to *copy* the values.
            basis_fn_p => tmp_basis_fns(basis_fns_ranking(i))
            select case(sys%system)
            case(heisenberg, chung_landau)
                call init_basis_fn(sys, basis_global%basis_fns(i), l=basis_fn_p%l)
            case default
                call init_basis_fn(sys, basis_global%basis_fns(2*i-1), l=basis_fn_p%l, ms=basis_fn_p%ms)
                call init_basis_fn(sys, basis_global%basis_fns(2*i), l=basis_fn_p%l, ms=-basis_fn_p%ms)
            end select
            deallocate(tmp_basis_fns(basis_fns_ranking(i))%l, stat=ierr)
            call check_deallocate('tmp_basis_fns(basis_fns_ranking(i',ierr)
        end do
        deallocate(tmp_basis_fns, stat=ierr)
        call check_deallocate('tmp_basis_fns',ierr)
        deallocate(basis_fns_ranking, stat=ierr)
        call check_deallocate('basis_fns_ranking',ierr)

        if (parent) then
            call write_basis_fn_header(sys)
            do i = 1, basis_global%nbasis
                call write_basis_fn(sys, basis_global%basis_fns(i), ind=i, new_line=.true.)
            end do
            write (6,'()')
        end if

        if (.not.t_store) then
            ! Should tidy up and deallocate everything we allocated.
            call end_basis_fns()
        end if

    end subroutine init_model_basis_fns

    pure function spin_symmetry(basis_fns, i, j) result(spin_match)

        ! In:
        !    basis_fns: list of single-particle basis functions.
        !    i: index of a basis function
        !    j: index of a basis function
        ! Returns:
        !    true if i and j refer to basis functions of the same spin.

        logical :: spin_match
        type(basis_fn_t), intent(in) :: basis_fns(:)
        integer, intent(in) :: i, j

        spin_match = basis_fns(i)%ms == basis_fns(j)%ms

    end function spin_symmetry

    subroutine set_orb(bit_lookup, f,iorb)

        ! In:
        !    bit_lookup: bit lookup table for basis functions (see basis_t).
        !    f: bit string of orbitals.
        !    iorb: orbital index.
        ! Out:
        !    f: bit string of orbitals with the bit corresponding to iorb set.

        ! Note that f must be zerod before first using this procedure.

        integer, intent(in) :: bit_lookup(:,:)
        integer, intent(in) :: iorb
        integer(i0), intent(inout) :: f(:)
        integer :: pos, ind

        pos = bit_lookup(1,iorb)
        ind = bit_lookup(2,iorb)
        f(ind) = ibset(f(ind),pos)

    end subroutine set_orb

    subroutine end_basis_fns()

        ! Clean up basis functions.

        use checking, only: check_deallocate

        integer :: ierr, i

        if (allocated(basis_global%basis_fns)) then
            do i = 1, basis_global%nbasis
                if (allocated(basis_global%basis_fns(i)%l)) deallocate(basis_global%basis_fns(i)%l, stat=ierr)
                call check_deallocate('basis_global%basis_fns(i',ierr)
            end do
            deallocate(basis_global%basis_fns, stat=ierr)
            call check_deallocate('basis_global%basis_fns',ierr)
        end if

    end subroutine end_basis_fns

end module basis

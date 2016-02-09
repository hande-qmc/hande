module system

! Module to store system information.  See also the basis module.

! Hubbard, Chung--Landau and Heisenberg systems:
!
! We assume that the lattice is defined by unit vectors and sites are located at
! all integer combinations of the primitive vectors.
! We never actually store the lattice vectors or reciprocal lattice vectors of
! the primitive unit cell.

! UEG:
!
! The lattice defines the simulation cell, to which periodic boundary conditions
! are applied.  Note that there is no underlying primitive lattice and so we
! only consider the case where the simulation cell is aligned with the
! coordinate system.

use const
use basis_types, only: basis_t
use molecular_integral_types
use ueg_types, only: ueg_basis_t
use symmetry_types, only: pg_sym_t, mom_sym_t

implicit none

! --- Interfaces to system-specific procedure pointers ---

abstract interface
    ! UEG-specific integral procedure pointers.
    ! The integral routines are different for 2D and UEG.  Abstract them using
    ! procedure pointers.
    pure function i_int_ueg(cell_param, b, i, a) result(intgrl)
        import :: p, basis_t
        real(p) :: intgrl
        real(p), intent(in) :: cell_param
        type(basis_t), intent(in) :: b
        integer, intent(in) :: i, a
    end function i_int_ueg
end interface

! --- System type constants ---

! Parameters to used to specify the system type.
enum, bind(c)
    enumerator :: hub_k
    enumerator :: hub_real
    enumerator :: read_in
    enumerator :: heisenberg
    enumerator :: ueg
    enumerator :: chung_landau
    enumerator :: ringium
end enum

! --- System-specific data structures ---

type sys_lattice_t

    ! Model (lattice/UEG) Hamiltonians
    ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    ! See also specific structures.

    ! 1, 2 or 3 dimensions.
    ! Set to a nonsensical value for error checking in input parser.
    ! It is useful to default to 0 so we can still safely allocate arrays using it even
    ! if we're not doing a lattice model.
    integer :: ndim = 0

    ! Is a triangular lattice is being used? (not UEG)
    logical :: triangular_lattice = .false.
    ! Is the lattice is bipartite or is it geometrically frustrated? (not UEG)
    logical :: bipartite_lattice = .false.

    ! Number of sites in crystal cell (Hubbard; Heisenberg).
    integer :: nsites

    ! Lattice vectors of crystal cell. (:,i) is the i-th vector.
    ! Not used for the UEG (see box_length).
    integer, allocatable :: lattice(:,:)  ! ndim, ndim.

    ! As we are working in an orthogonal space, the reciprocal lattice vectors are
    ! easily obtained:
    ! b_i = 2\pi/|a_i|^2 a_i in general,
    !     = 2\pi/ L_i for UEG, where L_i is box_length(i)
    real(p), allocatable :: rlattice(:,:) ! ndim, ndim. (:,i) is 1/(2pi)*b_i.

    ! Lengths of lattice vectors.
    ! This also defines the cubic simulation cell used in the UEG.
    real(p), allocatable :: box_length(:) ! ndim.

end type sys_lattice_t

type sys_k_lattice_t

    ! Model Hamiltonians in reciprocal space
    ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    ! See also specific (i.e. Hubbard and UEG) structures.

    ! Twist applied to wavevectors (systems in Bloch basis).
    real(p), allocatable :: ktwist(:)

end type sys_k_lattice_t

type sys_real_lattice_t

    ! Model Hamiltonians in real space
    ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    ! See also specific (i.e. Hubbard and Heisenberg) structures.

    ! The kinetic term is constant in the real space formulation of the Hubbard and Chung--Landau Hamiltonians:
    ! only the connectivity of the lattice matters.
    ! tmat(:,i) is a bit string.  The j-th bit corresponding to a basis function
    ! (as given by bit_lookup) is set if i and j are connected.
    ! We need to distinguish between connections within the cell and those due to
    ! periodic boundaries.  We do this by the following strategy:
    !   a) j>i.
    !          If the j-th bit is set then i and j are connected within the crystal
    !          cell.
    !   b) j<=i.
    !          If the i-th bit of tmat(:,j) is set, then i and j are connected due
    !          to periodic boundary conditions.
    ! This may seem like a somewhat arbitrary choice, but it enables for the
    ! correct evaluation of the kinetic energy using bit operations.
    ! Further it enables us to pick up cases such as the 2x2 (non-tilted) system,
    ! where a site is connected to a different site and that site's periodic image.
    integer(i0), allocatable :: tmat(:,:) ! (string_len, nbasis)

    ! Orbitals i and j are connected if the j-th bit of connected_orbs(:,i) is
    ! set.  This is a bit like tmat but without a bit set for a site being its own
    ! periodic image.
    integer(i0), allocatable :: connected_orbs(:,:) ! (string_len, nbasis)

    ! connected_sites(0,i) contains the number of unique sites connected to i.
    ! connected_sites(1:,i) contains the list of sites connected to site i (ie is the
    ! decoded/non-bit list form of connected_orbs).
    ! If connected_orbs(j,i) is 0 then it means there are fewer than 2*ndim unique sites
    ! that are connected to i that are not a periodic image of i (or connected to
    ! i both directly and via periodic boundary conditions).
    ! For the triangular lattice there are 3*ndim bonds and ndim must equal 2,
    ! so each site is connected to 6 other sites.
    integer, allocatable :: connected_sites(:,:) ! (0:2ndim, nbasis) or (0:3dim, nbasis)

    ! next_nearest_orbs(i,j) gives the number of paths by which sites i and j are
    ! are next nearest neighbors. For example, on a square lattice in the
    ! Heisenberg model, if we consider a spin, we can get to a next-nearest
    ! neighbor spin by going one right then one up, or to the same spin by going
    ! one up and then one right - there are two different paths, so the correpsonding
    ! value of next_nearest_orbs would be 2 for these spins. This is an important
    ! number to know when calculating the thermal energy squared in DMQMC.
    ! If two spins are not next-nearest neighbors by any path then this quantity is 0.
    ! By next nearest neighbors, it is meant sites which can be joined by exactly two
    ! bonds - any notion one may have of where the spins are located spatially is unimportant.
    integer(i0), allocatable :: next_nearest_orbs(:,:) ! (nbasis, nbasis)

    ! True if any site is its own periodic image.
    ! This is the case if one dimension (or more) has only one site per simulation
    ! cell.  If so then the an orbital can incur a kinetic interaction with itself.
    ! This is the only way that the integral < i | T | i >, where i is a basis
    ! function centred on a lattice site, can be non-zero.
    logical :: t_self_images

    ! True if we are actually only modelling a finite system (e.g. a H_2 molecule)
    ! False if we are modelling an infinite lattice
    ! The code is set up to model inifinite lattices by default, however in order
    ! to model only a finite "cluster" of sites, all one need do is set the
    ! connection matrix elements corresponding to connections across cell
    ! boundaries (i.e. periodic boundary conditions) to 0
    logical :: finite_cluster = .false. ! default to infinite crystals

end type sys_real_lattice_t

type sys_hubbard_t

    ! Hubbard model
    ! ^^^^^^^^^^^^^

    ! Hubbard T and U parameters specifying the kinetic energy and Coulomb
    ! interaction respectively.
    real(p) :: u = 1, t = 1
    ! The Coulomb integral in the momentum space formulation of the Hubbard model
    ! is constant, so it's convenient to store it.
    real(p) :: coulomb_k
    ! Crystal momentum symmetry information
    type(mom_sym_t) :: mom_sym

end type sys_hubbard_t

type sys_heisenberg_t

    ! Heisenberg model
    ! ^^^^^^^^^^ ^^^^^

    ! Coupling constant J.
    real(p) :: J = 1.0_p
    ! External magnetic field h, in the z direction (the z direction is defined
    ! the same direction as the external field).
    real(p) :: magnetic_field = 0.0_p
    ! Staggered magnetisation operator.
    ! staggered_magnetic_field gives the constant of proportionality:
    ! \hat{H} = -J \sum_{i,j} \sigma_i \sigma_j -
    !                      staggered_magnetic_field \sum_{i}(-1)^{\zeta(i)}\sigma_{i,z}
    ! where \zeta(i) gives \pm 1 depending upon which sublattice site i is on.
    ! NOTE: applicable only to bipartite lattices.
    real(p) :: staggered_magnetic_field = 0.0_p
    ! Number of bonds in the crystal cell.
    integer :: nbonds

    ! For the Heisenberg model, certain lattices can be split into two sublattices such
    ! that all the sites on one sublattice only have neighbors on the other sublattice.
    ! This is, for example, when finding the staggered magnetisation.  lattice_mask masks
    ! out one sublattice from a bit string representation of a spin product basis
    ! function.
    integer(i0), allocatable :: lattice_mask(:)

end type sys_heisenberg_t

type sys_ueg_t

    ! UEG
    ! ^^^

    ! Electron density.
    real(p) :: r_s = 1.0_p
    ! Energy cutoff for basis.
    ! This is in provided in scaled units of (2*pi/L)^2.
    real(p) :: ecutoff = 3.0_p

    ! UEG-specific basis lookup tables, etc.
    type(ueg_basis_t) :: basis

    ! Fermi-wavevector (Infinite system).
    real(p) :: kf = 1.0_p
    ! Fermi-Energy (Infinite system).
    real(p) :: ef = 1.0_p

    ! System-specific (e.g. dimensionality, potential, etc) integral procedures.
    procedure(i_int_ueg), pointer, nopass :: coulomb_int
    procedure(i_int_ueg), pointer, nopass :: exchange_int

    ! Symmetry information: we only store the gamma point
    integer :: gamma_sym

end type sys_ueg_t

type sys_read_in_t

    ! Read-in system
    ! ^^^^^^^^^^^^^^

    ! FCIDUMP filename (contains integrals which make up the Hamiltonian matrix).
    character(255) :: fcidump = 'FCIDUMP'
    ! Dipole integral file, contains integrals <i|O|a>, where O=x, y or z.
    character(255) :: dipole_int_file = ''
    ! UHF or RHF orbitals?
    logical :: uhf = .false.
    ! Core energy (e.g. nuclear-nuclear terms, contributions from frozen orbitals...
    real(p) :: Ecore
    ! Contribution from frozen core orbitals and nucleii terms to dipole moment
    real(p) :: dipole_core

    ! Single-particle orbitals transform according to Lz symmetry
    logical :: useLz = .false.

    ! Store for <i|h|j>, where h is the one-electron Hamiltonian operator.
    type(one_body_t) :: one_e_h_integrals

    ! Store for <i|o|j>, where o is a one-electron operator.
    type(one_body_t) :: one_body_op_integrals

    ! Store for the two-body integrals, <ij|1/r_12|ab>, where i,j,a,b are spin basis
    ! functions and 1/r_12 is the Coulomb operator.
    type(two_body_t) :: coulomb_integrals

    ! Data about the orbital symmetries
    type(pg_sym_t) :: pg_sym

end type sys_read_in_t

type sys_ringium_t

    ! Ringium
    ! ^^^^^^^

    ! Radius of the ring the electrons are confined to
    real(p) :: radius
    ! Lz cutoff for basis
    integer :: maxlz

end type sys_ringium_t

! --- System data container structure ---

type sys_t

    ! General System variables
    ! ^^^^^^^^^^^^^^^^^^^^^^^^

    ! Which system are we examining?  Hubbard (real space)? Hubbard (k space)? ...?
    integer :: system

    ! # of electrons
    ! *NOTE*: For the Heisenberg model nel refers to the number of spins up in the
    ! basis functions. This is to reuse code for the Hubbard model, where it simply
    ! refers to the number of electrons in the system
    integer :: nel = 0
    ! # number of virtual orbitals
    ! *NOTE*: For the Heisenberg model nvirt refers to the number of spins down, or
    ! the number of 0's in the basis functions
    integer :: nvirt

    ! Ms spin quantum number of states.
    integer :: Ms = huge(1)

    ! Spin polarisation is set in set_spin_polarisation in the system module.
    ! Only used in Fermionic systems; see note for nel and nvirt for how the spin
    ! polarisation is handled in the Heisenberg system.
    ! Note: Chung-Landau Hamiltonian is for spinless fermions; hence nalpha is
    ! (without loss of generality) set to nel and nbeta to 0.
    ! # number of alpha, beta electrons
    integer :: nalpha, nbeta
    ! # number of virtual alpha, beta spin-orbitals
    integer :: nvirt_alpha, nvirt_beta

    ! The specific symmetry sector to be used. A negative value implies that no symmetry
    ! is to be used. This is necessary when dealing the the uniform electron gas in very
    ! large basis sets and is only useful for canonical energy calculations.
    integer :: symmetry = huge(1)

    ! Number of symmetries.  These are specified for orbitals.
    integer :: nsym = 1
    ! Index of lowest symmetry (normally 0 or 1).
    integer :: sym0 = 1
    ! Index of highest symmetry (i.e. nsym+(sym0-1))
    integer :: sym_max

    ! Because Lz symmetry is not closed, we store some additional irreps and information about them:
    ! These functions allow one to loop over those syms.
    ! We go up to cross products involving three orbitals.  See comments in point_group_symmetry
    ! and init_pg_symmetry.
    integer :: nsym_tot = 1
    integer :: sym0_tot = 1
    integer :: sym_max_tot

    ! System uses Bloch states for basis functions (i.e. hub_k and UEG)?
    logical :: momentum_space = .false.

    ! Complete active space of basis.  Valid only in systems where the
    ! single-particle basis can be ordered by the single-particle eigenvalues (e.g.
    ! not in the real-space formulation of the Hubbard model or Heisenberg model).
    integer :: CAS(2) = (/ -1, -1 /)

    ! The maximum number of excitations which a system can have.
    ! Fermion systems: number of active elections.
    ! Heisenberg model: number of antiparallel pairs of spins which can be
    ! flipped.  Used only in DMQMC.
    integer :: max_number_excitations

    ! Basis set information
    ! ^^^^^^^^^^^^^^^^^^^^^

    type(basis_t) :: basis

    ! System-specific structure handles
    ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    type(sys_lattice_t) :: lattice
    type(sys_k_lattice_t) :: k_lattice
    type(sys_real_lattice_t) :: real_lattice
    type(sys_hubbard_t) :: hubbard
    type(sys_heisenberg_t) :: heisenberg
    type(sys_ueg_t) :: ueg
    type(sys_read_in_t) :: read_in
    type(sys_ringium_t) :: ringium

end type sys_t

contains

    subroutine init_system(sys)

        ! Initialise system based upon input parameters.
        ! In/Out:
        !    sys: system being studied.

        use checking, only: check_allocate, check_deallocate
        use errors, only: stop_all

        type(sys_t), intent(inout) :: sys

        integer :: i, ivec, ierr, counter

        associate(sl=>sys%lattice, sr=>sys%real_lattice, sk=>sys%k_lattice, su=>sys%ueg, sh=>sys%heisenberg)

            if (.not. sys%system == read_in) then

                allocate(sl%box_length(sl%ndim), stat=ierr)
                call check_allocate('sys%lattice%box_length', sl%ndim, ierr)
                allocate(sl%rlattice(sl%ndim,sl%ndim), stat=ierr)
                call check_allocate('sys%lattice%rlattice', sl%ndim**2, ierr)

                select case(sys%system)
                case(ueg)

                    ! UEG specific information.
                    if (sl%ndim /= 2 .and. sl%ndim /= 3) call stop_all('init_system','2D or 3D UEG not specified in input.')

                    ! Lattice vectors (integers) are not used in the UEG as the
                    ! simulation cell might well be defined by non-integers.
                    if (allocated(sl%lattice)) then
                        deallocate(sl%lattice, stat=ierr)
                        call check_deallocate('sys%lattice%lattice',ierr)
                        write (6,'(1X,a)') 'Ignoring lattice input for the UEG.'
                    end if

                    ! Use a cubic simulation cell.
                    ! The system is uniquely defined by two out of the number of
                    ! electrons, the density and the simulation cell lattice parameter.
                    ! It is most convenient to have the first two as input parameters.
                    select case(sl%ndim)
                    case(2)
                        sl%box_length = sys%ueg%r_s*sqrt(pi*sys%nel)
                    case(3)
                        sl%box_length = sys%ueg%r_s*((4*pi*sys%nel)/3)**(1.0_p/3.0_p)
                    end select
                    sl%rlattice = 0.0_p
                    forall (ivec=1:sl%ndim) sl%rlattice(ivec,ivec) = 1.0_p/sl%box_length(ivec)

                case(ringium)
                case default

                    ! Lattice-model specific information.

                    forall (ivec=1:sl%ndim) sl%box_length(ivec) = sqrt(real(dot_product(sl%lattice(:,ivec),sl%lattice(:,ivec)),p))
                    sl%nsites = nint(product(sl%box_length))

                    forall (ivec=1:sl%ndim) sl%rlattice(:,ivec) = sl%lattice(:,ivec)/sl%box_length(ivec)**2

                    ! If performing a calculation in all symmetry sectors, set Ms to be its maximum value
                    ! so that the necessary arrays will be allocated to their maximum size.

                end select

                if (allocated(sk%ktwist)) then
                    if (size(sk%ktwist) /= sl%ndim) &
                        call stop_all('init_system', 'Twist vector dimensions not consistent with lattice dimenstions.')
                else
                    allocate(sk%ktwist(sl%ndim), stat=ierr)
                    call check_allocate('sys%lattice%ktwist',sl%ndim,ierr)
                    sk%ktwist = 0.0_p
                end if

                ! For the Heisenberg model, we have Ms and nsites defined, but nel not.
                ! Here nel means the number of up spins, nvirt means the number of down spins.
                ! For other models, both Ms and nel have already been set, and nel refers
                ! to the number of electrons in the system.
                select case(sys%system)
                case(heisenberg)
                    sys%nel = (sl%nsites+sys%Ms)/2
                    sys%nvirt = (sl%nsites-sys%Ms)/2
                case(hub_k, hub_real)
                    sys%nvirt = 2*sl%nsites - sys%nel
                case(chung_landau)
                    sys%nvirt = sl%nsites - sys%nel
                    ! Spinless fermions.  Set all fermions to be spin-up.
                    sys%Ms = sys%nel
                case(ueg)
                    ! set nvirt in basis once the basis set has been generated.
                case(ringium)
                    ! Spin manifolds are degenerate in 1D - set all electrons spin up
                    sys%nvirt = 2*(sys%ringium%maxlz + 1) - sys%nel
                    sys%ms = sys%nel
                end select

                if (sys%system /= ueg .and. sys%system /= ringium) then
                    ! This checks if the lattice is the correct shape and correct size to be bipartite. If so it
                    ! sets the logical variable bipartite_lattice to be true, which allows staggered magnetizations
                    ! to be calculated.
                    counter = 0
                    do i = 1,sl%ndim
                        if ( abs(sum(sl%lattice(:,i)) - sl%box_length(i)) < depsilon.and. &
                             mod(ceiling(sl%box_length(i)), 2) == 0) counter = counter + 1
                    end do
                    if (counter == sl%ndim) sl%bipartite_lattice = .true.
                end if

                ! This logical variable is set true if the system being used has basis functions
                ! in momentum space - the Heisenberg and real Hubbard models are in real space.
                sys%momentum_space = .not.(      sys%system == hub_real   &
                                        .or. sys%system == heisenberg &
                                        .or. sys%system == chung_landau &
                                        .or. sys%system == read_in    &
                                      )

                sys%hubbard%coulomb_k = sys%hubbard%u/sl%nsites

            end if

        end associate

    end subroutine init_system

    subroutine end_lattice_system(sl, sk, sr)

        ! Clean up system allocations.

        ! In/Out:
        !    sl, sk, sr: sys_*lattice_t objects.  On output the components of all provided
        !       arguments are deallocated.

        use checking, only: check_deallocate

        type(sys_lattice_t), intent(inout), optional :: sl
        type(sys_k_lattice_t), intent(inout), optional :: sk
        type(sys_real_lattice_t), intent(inout), optional :: sr
        integer :: ierr

        if (present(sl)) then
            if (allocated(sl%box_length)) then
                deallocate(sl%box_length, stat=ierr)
                call check_deallocate('box_length',ierr)
            end if
            if (allocated(sl%rlattice)) then
                deallocate(sl%rlattice, stat=ierr)
                call check_deallocate('rlattice',ierr)
            end if
            if (allocated(sl%lattice)) then
                deallocate(sl%lattice, stat=ierr)
                call check_deallocate('lattice',ierr)
            end if
        end if
        if (present(sk)) then
            if (allocated(sk%ktwist)) then
                deallocate(sk%ktwist, stat=ierr)
                call check_deallocate('ktwist',ierr)
            end if
        end if
        if (present(sr)) then
            if (allocated(sr%tmat)) then
                deallocate(sr%tmat, stat=ierr)
                call check_deallocate('sr%tmat',ierr)
            end if
            if (allocated(sr%connected_orbs)) then
                deallocate(sr%connected_orbs, stat=ierr)
                call check_deallocate('sr%connected_orbs',ierr)
            end if
            if (allocated(sr%connected_sites)) then
                deallocate(sr%connected_sites, stat=ierr)
                call check_deallocate('sr%connected_sites',ierr)
            end if
        end if

    end subroutine end_lattice_system

    pure function in_FBZ(system_type,sl,k)

        ! Test if k is in the FBZ of the primitive unit cell.
        ! In:
        !    system_type: parameter defining system (equal to parameter at top
        !       of module).
        !    sl: lattice system
        !    k: wavevector in units of the reciprocal lattice of the crystal
        !       cell.

        logical :: in_FBZ
        integer, intent(in) :: system_type
        type(sys_lattice_t), intent(in) :: sl
        integer, intent(in) :: k(:) ! sl%ndim

        integer :: i
        real(p) :: kc(sl%ndim)

        select case(system_type)
        case(UEG)
            ! FBZ is infinite.
            in_FBZ = .true.
        case default
            ! Convert to cartesian units.
            forall (i=1:sl%ndim) kc(i) = sum(k*sl%rlattice(i,:))

            ! This test only works because the underlying lattice is orthogonal.
            ! The asymmetry of the boundary conditions prevent the acceptance of
            ! all wavevectors on the boundaries...
            in_FBZ = all(kc<=(0.50_p+depsilon)).and.all(kc>(-0.50_p+1.e-5_p))
        end select

    end function in_FBZ

    subroutine set_spin_polarisation(nbasis, sys)

        ! Set the spin polarisation information stored in components of sys.
        !    nalpha, nbeta: number of alpha, beta electrons.
        !    nvirt_alpha, nvirt_beta: number of alpha, beta virtual spin-orbitals.

        ! In:
        !    nbasis: number of single-particle basis functions.
        ! In/Out:
        !    sys: initialised system object describing system. On output
        !       components related to spin-polarisation are set.

        use errors, only: stop_all

        integer, intent(in) :: nbasis
        type(sys_t), intent(inout) :: sys

        character(len=*), parameter :: proc_name = 'set_spin_polarisation'
        character(len=*), parameter :: err_fmt = '("Required Ms not possible: ",i11, ".")'
        character(len=40) :: err

        select case(sys%system)

        case(heisenberg)

            ! Spin polarization is different (see comments in system) as the
            ! Heisenberg model is a collection of spins rather than electrons.
            ! See comments in init_system and at module-level.
            if (abs(sys%Ms) > sys%lattice%nsites) then
                write (err, err_fmt) sys%Ms
                call stop_all(proc_name, err)
            end if
            sys%nel = (sys%lattice%nsites + sys%Ms)/2
            sys%nvirt = (sys%lattice%nsites - sys%Ms)/2
            ! The Heisenberg model doesn't use these values, but they need to be
            ! initialized to something sensible as we allocate memory using them in
            ! alloc_det_info_t.
            sys%nalpha = 0
            sys%nbeta = 0
            sys%nvirt_alpha = 0
            sys%nvirt_beta = 0
            sys%max_number_excitations = min(sys%nel, (sys%lattice%nsites-sys%nel))

        case(chung_landau)

            ! Spinless system.  Similarly for the Heisenberg model but treat all fermions as alpha electrons.
            sys%nalpha = sys%nel
            sys%nvirt_alpha = sys%lattice%nsites - sys%nalpha
            sys%nbeta = 0
            sys%nvirt_beta = 0
            sys%max_number_excitations = min(sys%nel, (sys%lattice%nsites-sys%nel))

        case (ringium)

            if (sys%ms /= sys%nel) call stop_all(proc_name, "Ringium must be completely spin polarised.")

            sys%nalpha = sys%nel
            sys%nbeta = 0

            sys%nvirt_alpha = nbasis/2 - sys%nel
            sys%nvirt_beta = 0
            sys%max_number_excitations = min(sys%nalpha, sys%nvirt_alpha)

        case default

            ! Find the number of determinants with the required spin.
            if (abs(mod(sys%Ms,2)) /= mod(sys%nel,2) .or. abs(sys%Ms) > sys%nel) then
                write (err, err_fmt) sys%Ms
                call stop_all(proc_name, err)
            end if

            sys%nbeta = (sys%nel - sys%Ms)/2
            sys%nalpha = (sys%nel + sys%Ms)/2

            sys%nvirt_alpha = nbasis/2 - sys%nalpha
            sys%nvirt_beta = nbasis/2 - sys%nbeta
            sys%max_number_excitations = min(sys%nel, sys%nvirt)

        end select

        call set_fermi_energy(sys)

    end subroutine set_spin_polarisation

    elemental subroutine copy_sys_spin_info(sys1, sys2)

        ! Copy all spin information used in *any* system type (irrespective of
        ! which system type is actually being studied...)

        ! In:
        !    sys1: reference system object.
        ! In/Out:
        !    sys2: target system object.  On output, its spin-related components
        !       have the identical values as those in sys1.

        type(sys_t), intent(in) :: sys1
        type(sys_t), intent(inout) :: sys2

        sys2%nel = sys1%nel
        sys2%nvirt = sys1%nvirt
        sys2%nalpha = sys1%nalpha
        sys2%nbeta = sys1%nbeta
        sys2%nvirt_alpha = sys1%nvirt_alpha
        sys2%nvirt_beta = sys1%nvirt_beta

    end subroutine copy_sys_spin_info

    subroutine set_fermi_energy(sys)

        ! Set the Fermi energy and wavevector.
        ! Currently only implemented for the fully polarised/unpolarised
        ! electron gas.

        ! In/Out:
        !    sys: sys_t object, on output Fermi energy and wavevector are set.

        use errors, only: warning
        use parallel, only: parent

        type(sys_t), intent(inout) :: sys

        real(p):: pol_factor
        real(p) :: dim_factor

        select case(sys%system)
        case (ueg)
            ! Polarisation factor = 2 for polarised system, 1 for unpolarised.
            ! Only deal with fully spin (un)polarised system, so that kf^{up} =
            ! kf^{down}.
            pol_factor = 1.0_p + abs(real(sys%nalpha-sys%nbeta,p)/sys%nel)
            if (pol_factor-int(pol_factor) > depsilon .and. parent) &
                call warning('set_fermi_energy','Fermi energy not calculated correctly for given &
                             &spin polarisation. Please implement.')
            select case(sys%lattice%ndim)
            case (3)
                dim_factor = (9.0_p*pol_factor*pi/4.0_p)**(1.0_p/3.0_p)
            case (2)
                dim_factor = sqrt(2.0_p*pol_factor)
            case (1)
                dim_factor = pi*pol_factor / 4.0_p
            end select
            ! Fermi wavevector.
            sys%ueg%kf = dim_factor / sys%ueg%r_s
            ! Fermi energy.
            sys%ueg%ef = 0.5_p * sys%ueg%kf**2
        end select

    end subroutine set_fermi_energy

    subroutine sys_t_json(js, sys, terminal)

        ! Serialise a sys_t object in JSON format.

        ! In/Out:
        !   js: json_out_t controlling the output unit and handling JSON internal state.  Unchanged on output.
        ! In:
        !   sys: sys_t object describing system being studied.
        !   terminal (optional): if true, this is the last entry in the enclosing JSON object.  Default: false.

        use json_out

        type(json_out_t), intent(inout) :: js
        type(sys_t), intent(in) :: sys
        logical, intent(in), optional :: terminal
        character(10) :: isite
        integer :: i
        logical :: lattice_system

        lattice_system = sys%system == heisenberg .or. sys%system == hub_k .or. &
                         sys%system == hub_real .or. sys%system == chung_landau

        call json_object_init(js, 'system')

        call json_write_key(js, 'nbasis', sys%basis%nbasis)
        call json_write_key(js, 'nel', sys%nel)
        call json_write_key(js, 'nvirt', sys%nvirt)
        call json_write_key(js, 'Ms', sys%Ms)
        call json_write_key(js, 'nalpha', sys%nalpha)
        call json_write_key(js, 'nbeta', sys%nbeta)
        call json_write_key(js, 'nvirt_alpha', sys%nvirt_alpha)
        call json_write_key(js, 'nvirt_beta', sys%nvirt_beta)
        call json_write_key(js, 'nsym', sys%nsym)
        call json_write_key(js, 'sym0', sys%sym0)
        call json_write_key(js, 'sym_max', sys%sym_max)
        call json_write_key(js, 'nsym_tot', sys%nsym_tot)
        call json_write_key(js, 'sym0_tot', sys%sym0_tot)
        call json_write_key(js, 'sym_max_tot', sys%sym_max_tot)
        call json_write_key(js, 'symmetry', sys%symmetry)
        call json_write_key(js, 'max_number_excitations', sys%max_number_excitations)

        if (lattice_system) then
            call json_object_init(js, 'lattice')
            call json_write_key(js, 'ndim', sys%lattice%ndim)
            call json_write_key(js, 'nsites', sys%lattice%nsites)
            call json_write_key(js, 'lattice', sys%lattice%lattice, terminal=sys%system == hub_k)
            if (sys%system /= hub_k) then
                call json_write_key(js, 'triangular_lattice', sys%lattice%triangular_lattice)
                call json_write_key(js, 'bipartite_lattice', sys%lattice%bipartite_lattice)
                call json_object_init(js, 'connected_sites')
                do i = 1, sys%basis%nbasis
                    write (isite,'(i0)') i
                    call json_write_key(js, trim(isite), sys%real_lattice%connected_sites(1:,i), i==sys%basis%nbasis)
                end do
                call json_object_end(js)
                call json_write_key(js, 'self_image', sys%real_lattice%t_self_images)
                call json_write_key(js, 'finite_cluster', sys%real_lattice%finite_cluster, terminal=.true.)
            end if
            call json_object_end(js)
        end if

        if (sys%system == hub_real .or. sys%system == hub_k .or. sys%system == chung_landau) then
            call json_object_init(js, 'hubbard')
            call json_write_key(js, 'U', sys%hubbard%u)
            call json_write_key(js, 't', sys%hubbard%t, terminal=sys%system/=hub_k)
            if (sys%system == hub_k) call json_write_key(js, 'ktwist', sys%k_lattice%ktwist, terminal=.true.)
            call json_object_end(js, terminal=.true.)
        end if

        if (sys%system == heisenberg) then
            call json_object_init(js, 'heisenberg')
            call json_write_key(js, 'J', sys%heisenberg%J)
            call json_write_key(js, 'magnetic_field', sys%heisenberg%magnetic_field)
            call json_write_key(js, 'staggered_magnetic_field', sys%heisenberg%staggered_magnetic_field)
            call json_write_key(js, 'nbonds', sys%heisenberg%nbonds, terminal=.true.)
            call json_object_end(js, terminal=.true.)
        end if

        if (sys%system == ueg) then
            call json_object_init(js, 'ueg')
            call json_write_key(js, 'r_s', sys%ueg%r_s)
            call json_write_key(js, 'ecutoff', sys%ueg%ecutoff)
            call json_write_key(js, 'k_fermi', sys%ueg%kf)
            call json_write_key(js, 'E_fermi', sys%ueg%ef)
            call json_write_key(js, 'ktwist', sys%k_lattice%ktwist)
            call json_write_key(js, 'L', sys%lattice%box_length, terminal=.true.)
            call json_object_end(js, terminal=.true.)
        end if

        if (sys%system == read_in) then
            call json_object_init(js, 'read_in')
            call json_write_key(js, 'int_file', trim(sys%read_in%fcidump))
            call json_write_key(js, 'uhf', sys%read_in%uhf)
            call json_write_key(js, 'Ecore', sys%read_in%Ecore)
            if (trim(sys%read_in%dipole_int_file) /= '') then
                call json_write_key(js, 'dipole_int_file', trim(sys%read_in%dipole_int_file))
                call json_write_key(js, 'dipole_core', sys%read_in%dipole_core)
            end if
            call json_write_key(js, 'CAS', sys%CAS)
            call json_write_key(js, 'useLz', sys%read_in%useLz, terminal=.true.)
            call json_object_end(js, terminal=.true.)
        end if

        if (sys%system == ringium) then
            call json_object_init(js, 'ringium')
            call json_write_key(js, 'radius', sys%ringium%radius)
            call json_write_key(js, 'maxlz', sys%ringium%maxlz, terminal=.true.)
            call json_object_end(js, terminal=.true.)
        end if

        call json_object_end(js, terminal)

    end subroutine sys_t_json

end module system

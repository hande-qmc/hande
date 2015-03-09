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
    ! Chemical potential.
    real(p) :: chem_pot = 1.0_p

    ! When creating an arbitrary excitation, k_i,k_j->k_a,k_b, we must conserve
    ! crystal momentum, k_i+k_j-k_a-k_b=0.  Hence once we've chosen k_i, k_j and
    ! k_a, k_b is uniquely defined.  Further, once we've chosen k_i and k_j and if
    ! we require k_b to exist in the basis, then only certain values of k_a are
    ! permitted.  sys%ueg%ternary_conserve(0,k1,k2,k3) gives how many k_a are permitted
    ! for k_i+k_j = (k1,k2,k3) and sys%ueg%ternary_conserve(1:,k1,k2,k3) gives a bit
    ! string with only bytes set corresponding to permitted k_a values.  Note only
    ! basis functions corresponding to *alpha* orbitals are set.
    ! For systems with dimensionality lower than 3, the higher ki values are set to
    ! 0, i.e. dimensions:
    ! (0:string_len,-N:N,0,0) (1D)
    ! (0:string_len,-N:N,-N:N,0) (2D)
    ! (0:string_len,-N:N,-N:N,-N:N) (3D)
    ! NOTE: this contains values of k_i+k_j which cannot be formed by the basis with
    ! the energy cutoff.  Memory can be saved by not using a cubic array for
    ! k_i+k_j...
    integer(i0), allocatable :: ternary_conserve(:,:,:,:)

    ! System-specific (e.g. dimensionality, potential, etc) integral procedures.
    procedure(i_int_ueg), pointer, nopass :: coulomb_int
    procedure(i_int_ueg), pointer, nopass :: exchange_int

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

end type sys_read_in_t

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

    ! Spin polarisation is set in set_spin_polarisation in the system module.
    ! Only used in Fermionic systems; see note for nel and nvirt for how the spin
    ! polarisation is handled in the Heisenberg system.
    ! Note: Chung-Landau Hamiltonian is for spinless fermions; hence nalpha is
    ! (without loss of generality) set to nel and nbeta to 0.
    ! # number of alpha, beta electrons
    integer :: nalpha, nbeta
    ! # number of virtual alpha, beta spin-orbitals
    integer :: nvirt_alpha, nvirt_beta

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

end type sys_t

contains

    subroutine init_system(sys)

        ! Initialise system based upon input parameters.

        use calc, only: ms_in
        use fciqmc_data, only: all_spin_sectors

        use checking, only: check_allocate, check_deallocate
        use errors, only: stop_all

        type(sys_t), intent(inout) :: sys

        integer :: i, ivec, ierr, counter

        associate(sl=>sys%lattice, sr=>sys%real_lattice, sk=>sys%k_lattice, su=>sys%ueg, sh=>sys%heisenberg)

            if (.not. sys%system == read_in) then

                allocate(sl%box_length(sl%ndim), stat=ierr)
                call check_allocate('sys%lattice%box_length',size(sl%box_length),ierr)
                allocate(sl%rlattice(sl%ndim,sl%ndim), stat=ierr)
                call check_allocate('sys%lattice%rlattice',size(sl%rlattice),ierr)

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

                case default

                    ! Lattice-model specific information.

                    forall (ivec=1:sl%ndim) sl%box_length(ivec) = sqrt(real(dot_product(sl%lattice(:,ivec),sl%lattice(:,ivec)),p))
                    sl%nsites = nint(product(sl%box_length))

                    forall (ivec=1:sl%ndim) sl%rlattice(:,ivec) = sl%lattice(:,ivec)/sl%box_length(ivec)**2

                    ! If performing a calculation in all symmetry sectors, set ms_in to be its maximum value
                    ! so that the necessary arrays will be allocated to their maximum size.
                    if (all_spin_sectors) ms_in = sl%nsites

                end select

                if (.not.allocated(sk%ktwist)) then
                    allocate(sk%ktwist(sl%ndim), stat=ierr)
                    call check_allocate('sys%lattice%ktwist',sl%ndim,ierr)
                    sk%ktwist = 0.0_p
                end if

                ! For the Heisenberg model, we have ms_in and nsites defined, but nel not.
                ! Here nel means the number of up spins, nvirt means the number of down spins.
                ! For other models, both ms_in and nel have already been set, and nel refers
                ! to the number of electrons in the system.
                select case(sys%system)
                case(heisenberg)
                    sys%nel = (sl%nsites+ms_in)/2
                    sys%nvirt = (sl%nsites-ms_in)/2
                case(hub_k, hub_real)
                    sys%nvirt = 2*sl%nsites - sys%nel
                case(chung_landau)
                    sys%nvirt = sl%nsites - sys%nel
                    ! Spinless fermions.  Set all fermions to be spin-up.
                    ms_in = sys%nel
                case(ueg)
                    ! set nvirt in basis once the basis set has been generated.
                end select

                if (sys%system /= ueg) then
                    ! This checks if the lattice is the correct shape and correct size to be bipartite. If so it
                    ! sets the logical variable bipartite_lattice to be true, which allows staggered magnetizations
                    ! to be calculated.
                    counter = 0
                    do i = 1,sl%ndim
                        if ( sum(sl%lattice(:,i)) == sl%box_length(i) .and. &
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

            select case(sys%system)
            case(heisenberg)
                if (all_spin_sectors) then
                    sys%max_number_excitations = sl%nsites/2
                else
                    sys%max_number_excitations = min(sys%nel, (sl%nsites-sys%nel))
                end if
            case(chung_landau)
                sys%max_number_excitations = min(sys%nel, (sl%nsites-sys%nel))
            case default
                sys%max_number_excitations = sys%nel
            end select

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

    subroutine set_spin_polarisation(nbasis, Ms, sys)

        ! Set the spin polarisation information stored in components of sys.
        !    nalpha, nbeta: number of alpha, beta electrons.
        !    nvirt_alpha, nvirt_beta: number of alpha, beta virtual spin-orbitals.

        ! In:
        !    nbasis: number of single-particle basis functions.
        !    Ms: spin of determinants that are being considered.  Ignored for the
        !       Chung--Landau model.
        ! In/Out:
        !    sys: initialised system object describing system. On output
        !       components related to spin-polarisation are set.

        use errors, only: stop_all

        integer, intent(in) :: nbasis, Ms
        type(sys_t), intent(inout) :: sys

        character(len=*), parameter :: proc_name = 'set_spin_polarisation'
        character(len=*), parameter :: err_fmt = '("Required Ms not possible: ",i11, ".")'
        character(len=40) :: err

        select case(sys%system)

        case(heisenberg)

            ! Spin polarization is different (see comments in system) as the
            ! Heisenberg model is a collection of spins rather than electrons.
            ! See comments in init_system and at module-level.
            if (abs(Ms) > sys%lattice%nsites) then
                write (err, err_fmt) Ms
                call stop_all(proc_name, err)
            end if
            sys%nel = (sys%lattice%nsites + Ms)/2
            sys%nvirt = (sys%lattice%nsites - Ms)/2
            ! The Heisenberg model doesn't use these values, but they need to be
            ! initialized to something sensible as we allocate memory using them in 
            ! alloc_det_info_t.
            sys%nalpha = 0
            sys%nbeta = 0
            sys%nvirt_alpha = 0
            sys%nvirt_beta = 0

        case(chung_landau)

            ! Spinless system.  Similarly for the Heisenberg model but treat all fermions as alpha electrons.
            sys%nalpha = sys%nel
            sys%nvirt_alpha = sys%lattice%nsites - sys%nalpha
            sys%nbeta = 0
            sys%nvirt_beta = 0

        case default

            ! Find the number of determinants with the required spin.
            if (abs(mod(Ms,2)) /= mod(sys%nel,2) .or. abs(Ms) > sys%nel) then
                write (err, err_fmt) Ms
                call stop_all(proc_name, err)
            end if

            sys%nbeta = (sys%nel - Ms)/2
            sys%nalpha = (sys%nel + Ms)/2

            sys%nvirt_alpha = nbasis/2 - sys%nalpha
            sys%nvirt_beta = nbasis/2 - sys%nbeta

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
            if (pol_factor-int(pol_factor) > depsilon) &
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
            if (parent) write (6,'(1X,a13,1X,f10.8,/)') 'Fermi Energy:', sys%ueg%ef
        end select

    end subroutine set_fermi_energy

end module system

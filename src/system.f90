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

implicit none

! --- System type constants ---

! Parameters to used to specify the system type.
integer, parameter :: hub_k = 0
integer, parameter :: hub_real = 1
integer, parameter :: read_in = 2
integer, parameter :: heisenberg = 3
integer, parameter :: ueg = 4
integer, parameter :: chung_landau = 5

! --- System-specific data structures ---

type sys_lattice_t

    ! Model (lattice/UEG) Hamiltonians
    ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    ! See also specific structures.

    ! 1, 2 or 3 dimensions.
    integer :: ndim = -1 ! set to a nonsensical value for error checking in input parser.

    ! Is a triangular lattice is being used? (not UEG)
    logical :: triangular_lattice = .false.
    ! Is the lattice is bipartite or is it geometrically frustrated? (not UEG)
    logical :: bipartite_lattice = .false.

    ! Number of sites in crystal cell (Hubbard; Heisenberg).
    integer :: nsites

    ! Twist applied to wavevectors (systems in Bloch basis).
    real(p), allocatable :: ktwist(:)

    ! Lattice vectors of crystal cell. (:,i) is the i-th vector.
    ! Not used for the UEG (see box_length).
    integer, allocatable :: lattice(:,:)  ! ndim, ndim.

    ! As we are working in an orthogonal space, the reciprocal lattice vectors are
    ! easily obtained:
    ! b_i = 2\pi/|a_i|^2 a_i
    !     = 2\pi/ L_i for UEG, where L_i is box_length(i)
    real(p), allocatable :: rlattice(:,:) ! ndim, ndim. (:,i) is 1/(2pi)*b_i.

    ! Lengths of lattice vectors.
    ! This defines the cubic simulation cell used in the UEG.
    real(p), allocatable :: box_length(:) ! ndim.

    ! lvecs contains all combinations of the above lattice vectors, where the
    ! amplitude for each lattice vector can be either -1, 0 or +1. lvec(:,i)
    ! stores the i'th such combination.
    ! TODO: move to hubbard_real (ie sole place where it is actually used).
    integer, allocatable :: lvecs(:,:) ! ndim, 3**ndim

    ! Contains integer lattice lengths. If less than 3 dimensions are used
    ! then the corresponding unused components are set to 1.
    ! This is useful for making loops over all dimension general.
    integer :: lattice_size(3)

end type sys_lattice_t

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
    real(p) :: J = 1
    ! External magnetic field h, in the z direction (the z direction is defined
    ! the same direction as the external field).
    real(p) :: magnetic_field = 0
    ! Staggered magnetisation operator.
    ! staggered_magnetic_field gives the constant of proportionality:
    ! \hat{H} = -J \sum_{i,j} \sigma_i \sigma_j -
    !                      staggered_magnetic_field \sum_{i}(-1)^{\zeta(i)}\sigma_{i,z}
    ! where \zeta(i) gives \pm 1 depending upon which sublattice site i is on.
    ! NOTE: applicable only to bipartite lattices.
    real(p) :: staggered_magnetic_field = 0
    ! Number of bonds in the crystal cell.
    integer :: nbonds

end type sys_heisenberg_t

type sys_ueg_t

    ! UEG
    ! ^^^

    ! Electron density.
    real(p) :: r_s = 1.0_p
    ! Energy cutoff for basis.
    ! This is in provided in scaled units of (2*pi/L)^2.
    real(p) :: ecutoff = 3.0_p

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

    ! Spin polarisation is set in set_spin_polarisation in the determinants module.
    ! Only used in Fermionic systems; see note for nel and nvirt for how the spin
    ! polarisation is handled in the Heisenberg system.
    ! Note: Chung-Landau Hamiltonian is for spinless fermions; hence nalpha is
    ! (without loss of generality) set to nel and nbeta to 0.
    ! # number of alpha, beta electrons
    integer :: nalpha, nbeta
    ! # number of virtual alpha, beta spin-orbitals
    integer :: nvirt_alpha, nvirt_beta

    ! Number of symmetries.
    integer :: nsym = 1
    ! Index of lowest symmetry (normally 0 or 1).
    integer :: sym0 = 1
    ! Index of highest symmetry (i.e. nsym+(sym0-1))
    integer :: sym_max
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

    ! System-specific structure handles
    ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    type(sys_lattice_t) :: lattice
    type(sys_hubbard_t) :: hubbard
    type(sys_heisenberg_t) :: heisenberg
    type(sys_ueg_t) :: ueg
    type(sys_read_in_t) :: read_in

end type sys_t

contains

    subroutine init_system(sys)

        ! Initialise system based upon input parameters.

        use calc, only: ms_in

        use checking, only: check_allocate, check_deallocate
        use errors, only: stop_all

        type(sys_t), intent(inout) :: sys

        integer :: i, ivec, ierr, counter

        associate(sl=>sys%lattice, su=>sys%ueg, sh=>sys%heisenberg)

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

                end select

                if (.not.allocated(sl%ktwist)) then
                    allocate(sl%ktwist(sl%ndim), stat=ierr)
                    call check_allocate('sys%lattice%ktwist',sl%ndim,ierr)
                    sl%ktwist = 0.0_p
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

                ! lattice_size is useful for loops over a general number of dimensions. This
                ! variable is only concerned with simple lattices which could be bipartite,
                ! as it is used in init_determinants to split a bipartite lattice into its two parts.
                sl%lattice_size = 1
                sl%lattice_size(1) = ceiling(sl%box_length(1), 2)
                if (sl%ndim > 1) sl%lattice_size(2) = ceiling(sl%box_length(2), 2)
                if (sl%ndim > 2) sl%lattice_size(3) = ceiling(sl%box_length(3), 2)

                if (sys%system /= ueg) then
                    ! This checks if the lattice is the correct shape and correct size to be bipartite. If so it
                    ! sets the logical variable bipartite_lattice to be true, which allows staggered magnetizations
                    ! to be calculated.
                    counter = 0
                    do i = 1,sl%ndim
                        if ( sum(sl%lattice(:,i)) == sl%box_length(i) .and. mod(sl%lattice_size(i), 2) == 0) counter = counter + 1
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

                if (sl%triangular_lattice) then
                    ! Triangular lattice, only in 2d. Each site has 6 bonds, but each bond is
                    ! connected to 2 sites, so we divide by 2 to avoid counting twice. So there
                    ! are 6*nsites/2 = 3*nsites bonds in total.
                    sh%nbonds = 3*sl%nsites
                else
                    ! For a simple rectangular lattice, each site has 2*ndim bonds, so there are
                    ! (ndim*nsites) bonds in total.
                    sh%nbonds = sl%ndim*sl%nsites
                end if

                sys%hubbard%coulomb_k = sys%hubbard%u/sl%nsites

                select case(sys%system)
                case(heisenberg,chung_landau)
                    sys%max_number_excitations = min(sys%nel, (sl%nsites-sys%nel))
                case default
                    sys%max_number_excitations = sys%nel
                end select

            end if

        end associate

    end subroutine init_system

    subroutine end_lattice_system(sl)

        ! Clean up system allocations.

        use checking, only: check_deallocate

        type(sys_lattice_t), intent(inout) :: sl
        integer :: ierr

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
        if (allocated(sl%ktwist)) then
            deallocate(sl%ktwist, stat=ierr)
            call check_deallocate('ktwist',ierr)
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
        !    Ms: spin of determinants that are being considered.
        ! In/Out:
        !    sys: initialised system object describing system. On output
        !       components related to spin-polarisation are set.

        use errors, only: stop_all

        integer, intent(in) :: nbasis, Ms
        type(sys_t), intent(inout) :: sys

        select case(sys%system)

        case(heisenberg)

            ! Spin polarization is different (see comments in system) as the
            ! Heisenberg model is a collection of spins rather than electrons.
            ! See comments in init_system and at module-level.
            sys%nel = (sys%lattice%nsites + Ms)/2
            sys%nvirt = (sys%lattice%nsites - Ms)/2

        case(chung_landau)

            ! Spinless system.

        case default

            ! Find the number of determinants with the required spin.
            if (abs(mod(Ms,2)) /= mod(sys%nel,2)) call stop_all('set_spin_polarisation','Required Ms not possible.')

            sys%nbeta = (sys%nel - Ms)/2
            sys%nalpha = (sys%nel + Ms)/2

            sys%nvirt_alpha = nbasis/2 - sys%nalpha
            sys%nvirt_beta = nbasis/2 - sys%nbeta

        end select

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

end module system

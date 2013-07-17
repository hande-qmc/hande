module system

! Module to store system information.  See also the basis module.

! Hubbard and Heisenberg systems:
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

! --- System type ---

! Parameters to used to specify the system type.
integer, parameter :: hub_k = 0
integer, parameter :: hub_real = 1
integer, parameter :: read_in = 2
integer, parameter :: heisenberg = 3
integer, parameter :: ueg = 4

! Which system are we examining?  Hubbard (real space)? Hubbard (k space)? ...?
integer :: system_type = hub_k

! True for systems which use Bloch states for basis functions (hub_k and UEG)
logical :: momentum_space = .false.

! True if the lattice is bipartite. False if it is geometrically frustrated.
logical :: bipartite_lattice = .false.

! --- QMC reference state and trial (importance-sampling) functions ---

! For the Heisenberg model, several different trial functions can be used in the
! energy estimator. Only a single determinant can be used for the Hubbard model.
integer, parameter :: single_basis = 0
integer, parameter :: neel_singlet = 1

! Which trial function are we using? Only relevant to the Heisneberg model.
! trial_function will always be 0 for other models to represent a single determinant.
integer :: trial_function = 0

! For the Heisenberd model, a guiding function may be used,
! |psi_G> = \sum_{i} a_i |psi_i>, so that the new Hamiltonian matrix elements are
! H_ij^new = (a_i*H_ij)/a_j. This is just importance sampling. These functions
! represent the different types of functions whihc may be used.
integer, parameter :: no_guiding = 0
! Note that when we use the Neel singlet state as a guiding function, it must also
! be used as the trial function in calculating the projected energy.
integer, parameter :: neel_singlet_guiding = 1

! If we are not using importance sampling, this is set to 0. Else it is set to one
! of the above values to specify the corresponding guiding function being used.
integer :: guiding_function = 0

! --- General System variables ---

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

! Complete active space of basis.  Valid only in systems where the
! single-particle basis can be ordered by the single-particle eigenvalues (e.g.
! not in the real-space formulation of the Hubbard model or Heisenberg model).
integer :: CAS(2) = (/ -1, -1 /)

! --- Model Hamiltonians ---

! 1, 2 or 3 dimensions.
integer :: ndim = -1 ! set to a nonsensical value for error checking in input parser.

! Number of sites in crystal cell (Hubbard; Heisenberg).
integer :: nsites

! Number of bonds in the crystal cell (Heisenberg).
integer :: nbonds

! Electron density (UEG only)
real(p) :: r_s = 1.0_p
! Energy cutoff for basis (UEG only).
! This is in provided in scaled units of (2*pi/L)^2.
real(p) :: ueg_ecutoff = 3.0_p

! Lattice vectors of crystal cell. (:,i) is the i-th vector.
! Not used for the UEG (see box_length).
integer, allocatable :: lattice(:,:)  ! ndim, ndim.

! lvecs contains all combinations of the above lattice vectors, where the
! amplitude for each lattice vector can be either -1, 0 or +1. lvec(:,i)
! stores the i'th such combination.
integer, allocatable :: lvecs(:,:) ! ndim, 3**ndim

! If a triangular lattice is being used, this variable is true (Hubbard; Heisenberg).
logical :: triangular_lattice

! Lengths of lattice vectors.
! This defines the cubic simulation cell used in the UEG.
real(p), allocatable :: box_length(:) ! ndim.

! Contains integer lattice lengths. If less than 3 dimensions are used
! then the corresponding unused components are set to 1.
! This is useful for making loops over all dimension general.
integer :: lattice_size(3)

! As we are working in an orthogonal space, the reciprocal lattice vectors are
! easily obtained:
! b_i = 2\pi/|a_i|^2 a_i
!     = 2\pi/ L_i for UEG, where L_i is box_length(i)
real(p), allocatable :: rlattice(:,:) ! ndim, ndim. (:,i) is 1/(2pi)*b_i.

! Twist applied to wavevectors (Hubbard; UEG).
real(p), allocatable :: ktwist(:)

! The maximum number of excitations which a system can have.
integer :: max_number_excitations

! --- Hubbard model ---

! Hubbard T and U parameters specifying the kinetic energy and Coulomb
! interaction respectively.
real(p) :: hubu = 1, hubt = 1

! The Coulomb integral in the momentum space formulation of the Hubbard model
! is constant, so it's convenient to store it.
real(p) :: hub_k_coulomb

! --- Heisenberg model ---

! Coupling constant J In the Heisenberg model.
real(p) :: J_coupling = 1

! External magnetic field h in the Heisenberg model, in the z direction
! (the z direction is defined the same direction as the external field).
real(p) :: magnetic_field = 0

! This parameter allows a staggered magnetisation operator to be added
! to the Hamiltonian. staggered_magnetic_field gives the constant of proportionality:
! \hat{H} = -J \sum_{i,j} \sigma_i \sigma_j -
!                      staggered_magnetic_field \sum_{i}(-1)^{\zeta(i)}\sigma_{i,z}
! where \zeta(i) gives \pm 1 depending upon which sublattice site i is on.
! Applicable only to bipartite lattices.
real(p) :: staggered_magnetic_field = 0

! --- Read-in system ---

! FCIDUMP filename
character(255) :: fcidump = 'FCIDUMP'

! UHF or RHF orbitals?
logical :: uhf = .false.

! Core energy (e.g. nuclear-nuclear terms, contributions from frozen orbitals...
real(p) :: Ecore

contains

    subroutine init_system()

        ! Initialise system based upon input parameters.

        use calc, only: ms_in

        use checking, only: check_allocate, check_deallocate
        use errors, only: stop_all

        integer :: i, ivec, ierr, counter

        if (.not. system_type == read_in) then

            allocate(box_length(ndim), stat=ierr)
            call check_allocate('box_length',ndim,ierr)
            allocate(rlattice(ndim,ndim), stat=ierr)
            call check_allocate('rlattice',ndim*ndim,ierr)

            select case(system_type)
            case(ueg)

                ! UEG specific information.
                if (ndim /= 2 .and. ndim /= 3) call stop_all('init_system','2D or 3D UEG not specified in input.')

                ! Lattice vectors (integers) are not used in the UEG as the
                ! simulation cell might well be defined by non-integers.
                if (allocated(lattice)) then
                    deallocate(lattice, stat=ierr)
                    call check_deallocate('lattice',ierr)
                    write (6,'(1X,a)') 'Ignoring lattice input for the UEG.'
                end if

                ! Use a cubic simulation cell.
                ! The system is uniquely defined by two out of the number of
                ! electrons, the density and the simulation cell lattice parameter.
                ! It is most convenient to have the first two as input parameters.
                select case(ndim)
                case(2)
                    box_length = r_s*sqrt(pi*nel)
                case(3)
                    box_length = r_s*((4*pi*nel)/3)**(1.0_p/3.0_p)
                end select
                rlattice = 0.0_p
                forall (ivec=1:ndim) rlattice(ivec,ivec) = 1.0_p/box_length(ivec)

            case default

                ! Lattice-model specific information.

                forall (ivec=1:ndim) box_length(ivec) = sqrt(real(dot_product(lattice(:,ivec),lattice(:,ivec)),p))
                nsites = nint(product(box_length))

                forall (ivec=1:ndim) rlattice(:,ivec) = lattice(:,ivec)/box_length(ivec)**2

                ! This checks if the lattice is the correct shape and correct size to be bipartite. If so it
                ! sets the logical variable bipartite_lattice to be true, which allows staggered magnetizations
                ! to be calculated.
                counter = 0
                do i = 1,ndim
                    if ( sum(lattice(:,i)) == box_length(i) .and. mod(lattice_size(i), 2) == 0) counter = counter + 1
                end do
                if (counter == ndim) bipartite_lattice = .true.

            end select

            if (.not.allocated(ktwist)) then
                allocate(ktwist(ndim), stat=ierr)
                call check_allocate('ktwist',ndim,ierr)
                ktwist = 0.0_p
            end if

            ! For the Heisenberg model, we have ms_in and nsites defined, but nel not.
            ! Here nel means the number of up spins, nvirt means the number of down spins.
            ! For other models, both ms_in and nel have already been set, and nel refers
            ! to the number of electrons in the system.
            select case(system_type)
            case(heisenberg)
                nel = (nsites+ms_in)/2
                nvirt = (nsites-ms_in)/2
            case(hub_k, hub_real)
                nvirt = 2*nsites - nel
            case(ueg)
                ! set nvirt in basis once the basis set has been generated.
            end select

            ! lattice_size is useful for loops over a general number of dimensions. This
            ! variable is only concerned with simple lattices which could be bipartite,
            ! as it is used in init_determinants to split a bipartite lattice into its two parts.
            lattice_size = 1
            lattice_size(1) = ceiling(box_length(1), 2)
            if (ndim > 1) lattice_size(2) = ceiling(box_length(2), 2)
            if (ndim > 2) lattice_size(3) = ceiling(box_length(3), 2)


            ! This logical variable is set true if the system being used has basis functions
            ! in momentum space - the Heisenberg and real Hubbard models are in real space.
            momentum_space = .not.(system_type == hub_real .or. system_type == heisenberg .or. system_type == read_in)

            if (triangular_lattice) then
                ! Triangular lattice, only in 2d. Each site has 6 bonds, but each bond is
                ! connected to 2 sites, so we divide by 2 to avoid counting twice. So there
                ! are 6*nsites/2 = 3*nsites bonds in total.
                nbonds = 3*nsites
            else
                ! For a simple rectangular lattice, each site has 2*ndim bonds, so there are
                ! (ndim*nsites) bonds in total.
                nbonds = ndim*nsites
            end if

            hub_k_coulomb = hubu/nsites

            max_number_excitations = min(nel, (nsites-nel))

        end if

    end subroutine init_system

    subroutine end_system()

        ! Clean up system allocations.

        use checking, only: check_deallocate

        integer :: ierr

        if (allocated(box_length)) then
            deallocate(box_length, stat=ierr)
            call check_deallocate('box_length',ierr)
        end if
        if (allocated(rlattice)) then
            deallocate(rlattice, stat=ierr)
            call check_deallocate('rlattice',ierr)
        end if
        if (allocated(lattice)) then
            deallocate(lattice, stat=ierr)
            call check_deallocate('lattice',ierr)
        end if
        if (allocated(ktwist)) then
            deallocate(ktwist, stat=ierr)
            call check_deallocate('ktwist',ierr)
        end if

    end subroutine end_system

    pure function in_FBZ(k)

        ! Test if k is in the FBZ of the primitive unit cell.
        ! In:
        !    k: wavevector in units of the reciprocal lattice of the crystal
        !       cell.

        logical :: in_FBZ
        integer, intent(in) :: k(ndim)

        integer :: i
        real(p) :: kc(ndim)

        select case(system_type)
        case(UEG)
            ! FBZ is infinite.
            in_FBZ = .true.
        case default
            ! Convert to cartesian units.
            forall (i=1:ndim) kc(i) = sum(k*rlattice(i,:))

            ! This test only works because the underlying lattice is orthogonal.
            ! The asymmetry of the boundary conditions prevent the acceptance of
            ! all wavevectors on the boundaries...
            in_FBZ = all(kc<=(0.50_p+depsilon)).and.all(kc>(-0.50_p+1.e-5_p))
        end select

    end function in_FBZ

end module system

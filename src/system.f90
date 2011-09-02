module system

use const

! We assume that the lattice is defined by unit vectors and sites are located at
! all integer combinations of the primitive vectors.
! We never actually store the lattice vectors or reciprocal lattice vectors of
! the primitive unit cell.

implicit none

! Parameters to used to specify the system type.
integer, parameter :: hub_k = 0
integer, parameter :: hub_real = 1
integer, parameter :: heisenberg = 3

! Which system are we examining?  Hubbard (real space)? Hubbard (k space)? ...?
integer :: system_type = hub_k

! True for systems which use Bloch states for basis functions (hub_k and UEG)
logical :: momentum_space = .false.

! True if the lattice is bipartite. False if it is geometrically frustrated.
logical :: bipartite_lattice = .false.

! For the Heisenberg model, several different trial functions can be used in the
! energy estimator. Only a single determinant can be used for the Hubbard model.
integer, parameter :: single_basis = 0
integer, parameter :: uniform_combination = 1
integer, parameter :: neel_singlet = 2

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
integer, parameter :: gutzwiller_guiding = 2

! If we are not using importance sampling, this is set to 0. Else it is set to one
! of the above values to specify the corresponding guiding function being used.
integer :: guiding_function = 0

! 1, 2 or 3 dimensions.
integer :: ndim 

! Number of sites in crystal cell.
integer :: nsites

! Number of bonds in the crystal cell.
integer :: nbonds

! Lattice vectors of crystal cell. (:,i) is the i-th vector.
integer, allocatable :: lattice(:,:)  ! ndim, ndim.

! If a triangular lattice is being used, this variable is true.
logical :: triangular_lattice

! Lengths of lattice vectors.
real(p), allocatable :: box_length(:) ! ndim.

! Contains integer lattice lengths. If less than 3 dimensions are used
! then the corresponding unused components are set to 1.
! This is useful for making loops over all dimension general.
integer :: lattice_size(3)

! As we are working in an orthogonal space, the reciprocal lattice vectors are
! easily obtained:
! b_i = 2\pi/|a_i|^2 a_i
real(p), allocatable :: rlattice(:,:) ! ndim, ndim. (:,i) is 1/(2pi)*b_i.

! Twist applied to wavevectors.
real(p), allocatable :: ktwist(:)

! # of electrons
! *NOTE*: For the Heisenberg model nel refers to the number of spins up in the 
! basis functions. This is to reuse code for the Hubbard model, where it simply
! refers to the number of electrons in the system
integer :: nel = 0 
! # number of virtual orbitals
! *NOTE*: For the Heisenberg model nvirt refers to the number of spins down, or
! the number of 0's in the basis functions
integer :: nvirt

! Spin polarisation is set in set_spin_polarisation in the determinants modules
! Only used in Fermionic systems; see note for nel and nvirt for how the spin
! polarisation is handled in the Heisenberg system.
! # number of alpha, beta electrons
integer :: nalpha, nbeta
! # number of virtual alpha, beta spin-orbitals
integer :: nvirt_alpha, nvirt_beta

! Hubbard T and U parameters specifying the kinetic energy and Coulomb
! interaction respectively.
real(p) :: hubu = 1, hubt = 1

! Coupling constant J In the Heisenberg model.
real(p) :: J_coupling = 1

! External magnetic field h in the Heisenberg model, in the z direction
! (the z direction is defined the same direction as the external field).
real(p) :: h_field = 0

! This parameter allows a staggered magnetisation operator to be added
! to the Hamiltonian. staggered_field gives the constant of proportionality:
! \hat{H} = -J \sum_{i,j} \sigma_i \sigma_j - 
!                      staggered_field \sum_{i}(-1)^{\zeta}\sigma_{i}^{z}
real(p) :: staggered_field = 0

! For the Heisenberg model applied to bipartite lattices,
! a unitary transformation can be applied which rotates
! all spins on one of the sublattices by 180 degrees around the
! z axis. When working in a basis of spins in the z direction, this
! transformation does not affect the diagonal elements, but multiplies 
! all off-diagonal elements by -1 to make all elements either 0 or -2.
! unitary_factor is -1 if this transformation is applied or +1 if not.
integer :: unitary_factor = 1

! The Coulomb integral in the momentum space formulation of the Hubbard model
! is constant, so it's convenient to store it.
real(p) :: hub_k_coulomb

contains

    subroutine init_system()

        ! Initialise system based upon input parameters.

        use checking, only: check_allocate
        use calc, only: ms_in

        integer :: i, ivec, ierr, counter

        allocate(box_length(ndim), stat=ierr)
        call check_allocate('box_length',ndim,ierr)
        allocate(rlattice(ndim,ndim), stat=ierr)
        call check_allocate('rlattice',ndim*ndim,ierr)

        forall (ivec=1:ndim) box_length(ivec) = sqrt(real(dot_product(lattice(:,ivec),lattice(:,ivec)),p))
        nsites = nint(product(box_length))
        forall (ivec=1:ndim) rlattice(:,ivec) = lattice(:,ivec)/box_length(ivec)**2

        if (.not.allocated(ktwist)) then
            allocate(ktwist(ndim), stat=ierr)
            call check_allocate('ktwist',ndim,ierr)
            ktwist = 0.0_p
        end if

        hub_k_coulomb = hubu/nsites

        ! For the Heisenberg model, we have ms_in and nsites defined, but nel not.
        ! Here nel means the number of up spins, nvirt means the number of down spins.
        ! For other models, both ms_in and nel have already been set, and nel refers
        ! to the number of electrons in the system.
        if (system_type == heisenberg) then
            nel = (nsites+ms_in)/2
            nvirt = (nsites-ms_in)/2
        else
            nvirt = 2*nsites - nel
        end if
        
        ! lattice_size is useful for loops over a general number of dimensions. This
        ! variable is only concerned with simple lattices which could be bipartite,
        ! as it is used in init_determinants to split a bipartite lattice into its two parts.
        lattice_size = 1  
        lattice_size(1) = ceiling(box_length(1), 2)
        if (ndim > 1) lattice_size(2) = ceiling(box_length(2), 2)
        if (ndim > 2) lattice_size(3) = ceiling(box_length(3), 2)
        
        ! This checks if the lattice is the correct shape and correct size to be bipartite. If so it
        ! sets the logical variable bipartite_lattice to be true, which allows staggered magnetizations
        ! to be calculated.
        counter = 0
        do i = 1,ndim
            if ( sum(lattice(:,i)) == box_length(i) .and. mod(lattice_size(i), 2) == 0) counter = counter + 1 
        end do
        if (counter == ndim) bipartite_lattice = .true.
        
        ! This logical variable is set true if the system being used has basis functions
        ! in momentum space - the Heisenberg and real Hubbard models are in real space.
        momentum_space = .not.(system_type == hub_real .or. system_type == heisenberg)
        
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

    end subroutine init_system

    subroutine end_system()

        ! Clean up system allocations.

        use checking, only: check_deallocate

        integer :: ierr

        deallocate(box_length, stat=ierr)
        call check_deallocate('box_length',ierr)
        deallocate(rlattice, stat=ierr)
        call check_deallocate('rlattice',ierr)
        deallocate(lattice, stat=ierr)
        call check_deallocate('lattice',ierr)
        deallocate(ktwist, stat=ierr)
        call check_deallocate('ktwist',ierr)

    end subroutine end_system

    pure function in_FBZ(k)

        ! Test if k is in the FBZ of the primitive unit cell.
        ! In:
        !    k: wavevector in units of the reciprocal lattice of the crystal
        !       cell.

        logical :: in_FBZ
        integer, intent(in) :: k(ndim)

        integer :: i
        real :: kc(ndim)

        ! Convert to cartesian units.
        forall (i=1:ndim) kc(i) = sum(k*rlattice(i,:))

        ! This test only works because the underlying lattice is orthogonal.
        ! The asymmetry of the boundary conditions prevent the acceptance of
        ! all wavevectors on the boundaries...
        in_FBZ = all(kc<=(0.50_p+depsilon)).and.all(kc>(-0.50_p))

    end function in_FBZ

end module system

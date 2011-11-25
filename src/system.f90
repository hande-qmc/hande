module system

use const

! We assume that the lattice is defined by unit vectors and sites are located at
! all integer combinations of the primitive vectors.
! We never actually store the lattice vectors or reciprocal lattice vectors of
! the primitive unit cell.

implicit none

! --- System type ---

! Parameters to used to specify the system type.
integer, parameter :: hub_k = 0
integer, parameter :: hub_real = 1
integer, parameter :: read_in = 2

! Which system are we examining?  Hubbard (real space)? Hubbard (k space)? ...?
integer :: system_type = hub_k

! --- General System variables ---

! # of electrons
integer :: nel = 0 
! # number of virtual orbitals
integer :: nvirt

! Spin polarisation is set in set_spin_polarisation in the determinants module.
! # number of alpha, beta electrons
integer :: nalpha, nbeta
! # number of virtual alpha, beta spin-orbitals
integer :: nvirt_alpha, nvirt_beta

! Number of symmetries.
integer :: nsym

! Index of lowest symmetry (normally 0 or 1).
integer :: sym0

! Complete active space of basis.  Valid only in systems where the
! single-particle basis can be ordered by the single-particle eigenvalues (e.g.
! not in the real-space formulation of the Hubbard model)
integer :: CAS(2) = (/ -1, -1 /)

! --- Model Hamiltonians ---

! 1, 2 or 3 dimensions.
integer :: ndim 

! Number of sites in crystal cell.
integer :: nsites 

! Lattice vectors of crystal cell. (:,i) is the i-th vector.
integer, allocatable :: lattice(:,:)  ! ndim, ndim.

! Lengths of lattice vectors.
real(p), allocatable :: box_length(:) ! ndim.

! As we are working in an orthogonal space, the reciprocal lattice vectors are
! easily obtained:
! b_i = 2\pi/|a_i|^2 a_i
real(p), allocatable :: rlattice(:,:) ! ndim, ndim. (:,i) is 1/(2pi)*b_i.

! Twist applied to wavevectors.
real(p), allocatable :: ktwist(:)

! --- Hubbard model ---

! Hubbard T and U parameters specifying the kinetic energy and Coulomb
! interaction respectively.
real(p) :: hubu = 1, hubt = 1

! The Coulomb integral in the momentum space formulation of the Hubbard model
! is constant, so it's convenient to store it.
real(p) :: hub_k_coulomb

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

        use checking, only: check_allocate

        integer :: ivec, ierr

        if (.not. system_type == read_in) then

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
        real :: kc(ndim)

        ! Convert to cartesian units.
        forall (i=1:ndim) kc(i) = sum(k*rlattice(i,:))

        ! This test only works because the underlying lattice is orthogonal.
        ! The asymmetry of the boundary conditions prevent the acceptance of
        ! all wavevectors on the boundaries...
        in_FBZ = all(kc<=(0.50_p+depsilon)).and.all(kc>(-0.50_p))

    end function in_FBZ

end module system

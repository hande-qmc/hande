module system

use const

! We assume that the lattice is defined by unit vectors and sites are located at
! all integer combinations of the primitive vectors.
! We never actually store the lattice vectors or reciprocal lattice vectors of
! the primitive unit cell.

implicit none

integer :: ndim ! 1, 2 or 3 dimensions.
integer :: nsites ! Number of sites in crystal cell.
integer, allocatable :: lattice(:,:)  ! ndim, ndim.  Gives lattice vectors of crystal cell. (:,i) is the i-th vector.
real(dp), allocatable :: box_length(:) ! ndim.  Gives lengths of lattice vectors.

! As we are working in an orthogonal space, the reciprocal lattice vectors are
! easily obtained:
! b_i = 2\pi/|a_i|^2 a_i
real(dp), allocatable :: rlattice(:,:) ! ndim, ndim. (:,i) is 1/(2pi)*b_i.

integer :: nel = 0 ! # of electrons

real(dp) :: hubu = 1, hubt = 1

contains

    subroutine init_system()

        integer :: ivec, ierr

        allocate(box_length(ndim), stat=ierr)
        allocate(rlattice(ndim,ndim), stat=ierr)

        forall (ivec=1:ndim) box_length(ivec) = sqrt(real(dot_product(lattice(:,ivec),lattice(:,ivec)),dp))
        nsites = nint(product(box_length))
        forall (ivec=1:ndim) rlattice(:,ivec) = lattice(:,ivec)/box_length(ivec)**2

    end subroutine init_system

    pure function in_FBZ(k)

        logical :: in_FBZ
        integer, intent(in) :: k(ndim)

        integer :: i
        real :: kc(ndim)

        forall (i=1:ndim) kc(i) = sum(k*rlattice(i,:))

        in_FBZ = all(kc<=(0.50_dp+depsilon)).and.all(kc>(-0.50_dp))

    end function in_FBZ

end module system

module ueg_system

! Core routines for the uniform electron gas system.

use const

implicit none

! basis_fns stores the basis functions as a list ordered by the kinetic energy
! of the wavefunction.  This is inconvenient for the UEG where we need to find
! a given basis function from the sum of two other wavevectors.
! ueg_basis_lookup(ind) returns the index of the alpha spin-orbital with
! wavevector, k,  in the basis_fns array, where ind is defined by a function
! which considers a cube of wavevectors:
!   ind = (k_x + kmax) + (k_y + kmax)*N_kx + (k_z + kmax)*N_kx^2 + 1
!       = k_x + k_y*N_kx + k_z*N_kx*N_z + k_max*(1 + N_kx + N_kx^2) + 1
! (and analogously for the 2D case), where:
!    kmax is the maximum component of a wavevector in the smallest square/cube
!    which contains all the wavevectors in the basis.
!    N_kx is the number of k-points in each dimension
integer, allocatable :: ueg_basis_lookup(:) ! (N_kx^ndim)

! ueg_basis_dim = (1, N_kx, N_kx^2), so that
!    ind = ueg_basis_dim.k + ueg_basis_origin
! for ueg_basis_lookup.
integer, allocatable :: ueg_basis_dim(:) ! (ndim)

! ueg_basis_origin accounts for the fact that ueg_basis_lookup is a 1-indexed array.  
! ueg_basis_origin = k_max*(1 + N_x + N_x*N_y) + 1
integer :: ueg_basis_origin

! Max component of a wavevector in the UEG basis set, kmax.
! Note that N_x = 2*kmax+1
integer :: ueg_basis_max

contains

    subroutine init_ueg_indexing()

        ! Create arrays and data for index mapping needed for UEG.

        use basis, only: basis_fns, nbasis
        use system, only: box_length, ueg_ecutoff, ndim

        use checking, only: check_allocate

        integer :: ierr, i, N_kx

        ueg_basis_max = ceiling(sqrt(2*ueg_ecutoff))

        N_kx = 2*ueg_basis_max+1

        allocate(ueg_basis_dim(ndim), stat=ierr)
        call check_allocate('ueg_basis_dim', ndim, ierr)
        forall (i=1:ndim) ueg_basis_dim(i) = N_kx**(i-1)

        ueg_basis_origin = ueg_basis_max*(1 + N_kx + N_kx**2) + 1

        allocate(ueg_basis_lookup(N_kx**ndim), stat=ierr)
        call check_allocate('ueg_basis_lookup', N_kx**ndim, ierr)

        ! ueg_basis_lookup should be -1 for any wavevector that is in the
        ! square/cubic grid defined by ueg_basis_max but not in the actual basis
        ! set described by ueg_ecutoff.
        ueg_basis_lookup = -1

        ! Now fill in the values for the alpha orbitals which are in the basis.
        forall (i=1:nbasis:2) ueg_basis_lookup(dot_product(basis_fns(i)%l, ueg_basis_dim) + ueg_basis_origin) = i

    end subroutine init_ueg_indexing

    pure function ueg_basis_index(k, spin) result(indx)

        ! In:
        !    k: wavevector in units of 2\pi/L.
        !    spin: 1 for alpha orbital, -1 for beta orbital
        ! Returns:
        !    Index of basis function in the (energy-ordered) basis_fns array.
        !    Set to < 0 if the spin-orbital described by k and spin is not in the
        !    basis set.

        use system, only: ndim

        integer :: indx
        integer, intent(in) :: k(ndim), spin

        if (minval(k) < -ueg_basis_max .or. maxval(k) > ueg_basis_max) then
            indx = -1
        else
            ! ueg_basis_lookup contains the mapping between a wavevector in
            ! a given basis and its entry in the energy-ordered list of basis
            ! functions.
            indx = ueg_basis_lookup(dot_product(k,ueg_basis_dim) + ueg_basis_origin)
            ! ueg_basis_lookup only contains entries for the alpha spin-orbital.
            ! The corresponding beta orbital is the next entry in the basis_fns
            ! array.
            if (spin < 0) indx = indx + 1
        end if

    end function ueg_basis_index

    subroutine end_ueg_indexing()

        ! Clean up UEG index arrays.
        
        use checking, only: check_deallocate

        integer :: ierr

        if (allocated(ueg_basis_lookup)) then
            deallocate(ueg_basis_lookup, stat=ierr)
            call check_deallocate('ueg_basis_lookup', ierr)
        end if
        if (allocated(ueg_basis_dim)) then
            deallocate(ueg_basis_dim, stat=ierr)
            call check_deallocate('ueg_basis_dim', ierr)
        end if

    end subroutine end_ueg_indexing

end module ueg_system

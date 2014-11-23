module ueg_types

! Derived types and interfaces specific to the UEG.

! We only define the types here; all procedures are in ueg_system.

implicit none

type ueg_basis_t

    ! basis_fns stores the basis functions as a list ordered by the kinetic energy
    ! of the wavefunction.  This is inconvenient for the UEG where we need to find
    ! a given basis function from the sum of two other wavevectors.
    ! sys%ueg%basis%lookup(ind) returns the index of the alpha spin-orbital with
    ! wavevector, k,  in the basis_fns array, where ind is defined by a function
    ! which considers a cube of wavevectors:
    !   ind = (k_x + kmax) + (k_y + kmax)*N_kx + (k_z + kmax)*N_kx^2 + 1
    !       = k_x + k_y*N_kx + k_z*N_kx*N_z + k_max*(1 + N_kx + N_kx^2) + 1
    ! (and analogously for the 2D case), where:
    !    kmax is the maximum component of a wavevector in the smallest square/cube
    !    which contains all the wavevectors in the basis.
    !    N_kx is the number of k-points in each dimension
    integer, allocatable :: lookup(:) ! (N_kx^ndim)

    ! sys%ueg%basis%offset_inds = (1, N_kx, N_kx^2), so that
    !    ind = offset_inds.k + offset
    ! for sys%ueg%basis%lookup.
    integer, allocatable :: offset_inds(:) ! (ndim)

    ! offset accounts for the fact that lookup is a 1-indexed array.
    ! origin = k_max*(1 + N_x + N_x*N_y) + 1
    integer :: offset

    ! Max component of a wavevector in the UEG basis set.
    ! Note that N_kx = 2*kmax+1
    integer :: kmax

end type ueg_basis_t

end module ueg_types

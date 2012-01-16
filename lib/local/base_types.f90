module base_types

! Module containing derived types useful in various circumstances.

use const

implicit none

! Allocatable 1D array.  Useful for creating triangular arrays etc.

type alloc_int1d
    integer, allocatable :: v(:)
end type

type alloc_rp1d
    real(p), allocatable :: v(:)
end type

! Allocatable 2D array.  Useful for creating triangular spin arrays etc.

type alloc_int2d
    integer, allocatable :: v(:,:)
end type

type alloc_rp2d
    real(p), allocatable :: v(:,:)
end type

end module base_types

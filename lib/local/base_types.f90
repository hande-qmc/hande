module base_types

! Module containing derived types useful in various circumstances.

use const

implicit none

! Allocatable 1D array.  Useful for creating triangular arrays etc.
type alloc_r1d
    real(p), allocatable :: v(:)
end type

end module base_types

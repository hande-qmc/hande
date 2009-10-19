module system

implicit none

integer :: ndim ! 1, 2 or 3 dimensions.
integer, allocatable :: Length(:) ! ndim.  Gives box lengths.

end module system

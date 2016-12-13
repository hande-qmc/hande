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

type alloc_rdp1d
    real(dp), allocatable :: v(:)
end type

! Allocatable 2D array.  Useful for creating triangular spin arrays etc.

type alloc_int2d
    integer, allocatable :: v(:,:)
end type

type alloc_rp2d
    real(p), allocatable :: v(:,:)
end type

contains

    subroutine dealloc_rdp1d(array)

        ! Deallocates an allocatable array of 1D allocatable real
        ! arrays, all of arbitrary size.

        ! In/Out:
        !   array: array of type alloc_rdp1d. Fully deallocated on output.

        use checking, only: check_deallocate

        type(alloc_rdp1d), allocatable, intent(inout) :: array(:)
        integer :: i, ierr

        do i = lbound(array, dim=1), ubound(array, dim=1)
            deallocate(array(i)%v, stat=ierr)
            call check_deallocate('array%v', ierr)
        end do

        deallocate(array, stat=ierr)
        call check_deallocate('array', ierr)

    end subroutine dealloc_rdp1d

    subroutine dealloc_int2d(array)

        ! Deallocates an allocatable array of 2D allocatable integer
        ! arrays, all of arbitrary size.

        ! In/Out:
        !   array: array of type alloc_int2d. Fully deallocated on output.

        use checking, only: check_deallocate

        type(alloc_int2d), allocatable, intent(inout) :: array(:)
        integer :: i, ierr

        do i = lbound(array, dim=1), ubound(array, dim=1)
            deallocate(array(i)%v, stat=ierr)
            call check_deallocate('array%v', ierr)
        end do

        deallocate(array, stat=ierr)
        call check_deallocate('array', ierr)

    end subroutine dealloc_int2d

end module base_types

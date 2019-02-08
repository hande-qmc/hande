module base_types

! Module containing derived types useful in various circumstances.

use const
use shmem, only: shmem_handle_t

implicit none

! pointer 1D array.  Useful for creating triangular arrays etc.

type alloc_int1d
    integer, pointer :: v(:)
    ! shmem handle for safe allocation/deallocation.
    type(shmem_handle_t) :: shmem_handle
end type

type alloc_rp1d
    real(p), pointer :: v(:)
    ! shmem handle for safe allocation/deallocation.
    type(shmem_handle_t) :: shmem_handle
end type

type alloc_rdp1d
    real(dp), pointer :: v(:)
end type

! pointer 2D array.  Useful for creating triangular spin arrays etc.

type alloc_int2d
    integer, pointer :: v(:,:)
end type

type alloc_rp2d
    real(p), pointer :: v(:,:)
    type(shmem_handle_t) :: shmem_handle
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

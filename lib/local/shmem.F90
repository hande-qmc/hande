module shmem

! A wrapper around allocating shared memory.
! Some notes:
! * Using mmap may work in older systems but is rather fragile and depends upon a system giving
!  access to the shmem interface.  This may not always be true (e.g. not true for OS X).
! * Using MPI 3's shared memory interface is much cleaner and more portable but requires a fairly
!   recent MPI implementation.
! * Where neither are available, fall back to just doing standard allocate.

implicit none

private
public :: allocate_shared, deallocate_shared

interface allocate_shared
    !module procedure alloc_shared_real_sp_1D
    !module procedure alloc_shared_real_sp_2D
    module procedure alloc_shared_real_dp_1D
    module procedure alloc_shared_real_dp_2D
end interface

interface deallocate_shared
    !module procedure dealloc_shared_real_sp_1D
    !module procedure dealloc_shared_real_sp_2D
    module procedure dealloc_shared_real_dp_1D
    module procedure dealloc_shared_real_dp_2D
end interface deallocate_shared

interface
    subroutine alloc_shared_posix(A, handle, N) bind(C)
        ! Allocate a block of memory that is accessible from all cores on a node
        ! via mmap.
        ! In:
        !   handle: character string unique to shared memory block.
        !   N: number of bytes to allocate.
        ! Out:
        !   A: C pointer to memory.
        use, intrinsic :: iso_c_binding
        implicit none
        type(c_ptr), intent(out) :: A
        character(kind=c_char), intent(in) :: handle(*)
        integer(c_size_t), value :: N
    end subroutine alloc_shared_posix
    subroutine free_shared_posix(A,handle,N) bind(C)
        ! Deallocate a block of memory allocated using mmap.
        ! In:
        !   handle: character string unique to shared memory block.
        !   N: number of bytes allocated.
        ! Out:
        !   A: C pointer to memory.
        use, intrinsic :: iso_c_binding
        implicit none
        type(c_ptr), intent(out) :: A
        character(kind=c_char), intent(in) :: handle(*)
        integer(c_size_t), value :: N
    end subroutine free_shared_posix
end interface

contains

    ! alloc_shared_*

    subroutine alloc_shared_real_dp_1D(A, A_name, N1)

        use const, only: dp, int_64
        use checking, only: check_allocate
        use utils, only: fstring_to_carray
        use, intrinsic :: iso_c_binding, only: c_f_pointer, c_ptr, c_null_ptr

        real(dp), pointer, intent(out) :: A(:)
        character(*), intent(in) :: A_name
        integer(int_64), intent(in) :: N1
        integer :: ierr
        integer, parameter :: nbytes = 8
        type(c_ptr) :: ptr
        ptr = c_null_ptr

#if defined SHMEM_POSIX
        ierr = 0
        call alloc_shared_posix(ptr, fstring_to_carray(A_name), N1*nbytes)
        call c_f_pointer(ptr, A, [N1])
#else
        allocate(A(N1), stat=ierr)
        call check_allocate(A_name, N1, ierr)
#endif

    end subroutine alloc_shared_real_dp_1D

    subroutine alloc_shared_real_dp_2D(A, A_name, N1, N2)

        use const, only: dp, int_64
        use checking, only: check_allocate
        use, intrinsic :: iso_c_binding, only: c_f_pointer, c_ptr, c_null_ptr

        real(dp), pointer, intent(out) :: A(:,:)
        character(*), intent(in) :: A_name
        integer(int_64), intent(in) :: N1, N2
        integer :: ierr
        integer, parameter :: nbytes = 8
        type(c_ptr) :: ptr
        ptr = c_null_ptr

#if defined SHMEM_POSIX
        ierr = 0
        call alloc_shared_posix(ptr, A_name, N1*N2*nbytes)
        call c_f_pointer(ptr, A, [N1, N2])
#else
        allocate(A(N1, N2), stat=ierr)
        call check_allocate(A_name, N1*N2, ierr)
#endif

    end subroutine alloc_shared_real_dp_2D

    ! dealloc_shared_*

    subroutine dealloc_shared_real_dp_1D(A, A_name)

        use const, only: dp, int_64
        use checking, only: check_deallocate
        use, intrinsic :: iso_c_binding, only: c_loc, c_ptr, c_null_ptr

        real(dp), pointer, intent(inout) :: A(:)
        character(*), intent(in) :: A_name
        integer :: ierr
        integer, parameter :: nbytes = 8
        type(c_ptr) :: ptr
        ptr = c_null_ptr

#if defined SHMEM_POSIX
        ierr = 0
        ptr = c_loc(A)
        call free_shared_posix(ptr, A_name, size(A, kind=int_64)*nbytes)
        nullify(A)
#else
        deallocate(A, stat=ierr)
        call check_deallocate(A_name, ierr)
#endif

    end subroutine dealloc_shared_real_dp_1D

    subroutine dealloc_shared_real_dp_2D(A, A_name)

        use const, only: dp, int_64
        use checking, only: check_deallocate
        use, intrinsic :: iso_c_binding, only: c_loc, c_ptr, c_null_ptr

        real(dp), pointer, intent(inout) :: A(:,:)
        character(*), intent(in) :: A_name
        integer :: ierr
        integer, parameter :: nbytes = 8
        type(c_ptr) :: ptr
        ptr = c_null_ptr

#if defined SHMEM_POSIX
        ierr = 0
        ptr = c_loc(A)
        call free_shared_posix(ptr, A_name, size(A, kind=int_64)*nbytes)
        nullify(A)
#else
        deallocate(A, stat=ierr)
        call check_deallocate(A_name, ierr)
#endif

    end subroutine dealloc_shared_real_dp_2D

end module shmem

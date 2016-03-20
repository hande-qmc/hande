module shmem

! A wrapper around allocating shared memory.
! Some notes:
! * Using mmap may work in older systems but is rather fragile and depends upon a system giving
!  access to the shmem interface.  This may not always be true (e.g. not true for OS X).
! * Using MPI 3's shared memory interface is much cleaner and more portable but requires a fairly
!   recent MPI implementation.
! * Where neither are available, fall back to just doing standard allocate.

! If:
! * ENABLE_SHMEM_POSIX is defined: shm_open, mmap and other POSIX functions are used.
! * DISABLE_MPI3 is defined or PARALLEL is not defined, fall back to standard allocate (no shared memory).
! * otherwise MPI 3 is used (default).

! With MPI 3, we create a single array per node and provide a window to it on all processors.
! It is convenient to only do MPI operations (e.g. MPI_BCast) using the intra-node communicator.
! If the array is being repeatedly modified and read, then care must be taken to avoid data races.
! See the RMA section of MPI 3 for more details.  If, as in the integral use-case, this is just
! being used for a large lookup table, then a barrier across the node communicator after
! initialisation is sufficient.

! Avoiding such data races with POSIX is trickier(!) so assume that it's only
! used for large read-only lookup tables and each process can merrily overwrite
! the lookup table with the same information.

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

    ! helpers

    subroutine mpi3_shared_alloc(nbytes, ptr)

        use, intrinsic :: iso_c_binding, only: c_ptr
        use const, only: int_64
        use errors, only: stop_all
        use parallel

        integer(int_64), intent(in) :: nbytes
        type(c_ptr), intent(out) :: ptr
        integer :: win, win_info, ierr, disp
        integer(MPI_ADDRESS_KIND) :: n

#if defined PARALLEL && ! defined DISABLE_MPI3
        call MPI_Info_create(win_info, ierr)
        call MPI_Info_set(win_info, "alloc_shared_non_contig", "true", ierr)
        if (intra_node_comm == MPI_COMM_NULL) then
            ! Just get access.
            call MPI_Win_allocate_shared(0_MPI_ADDRESS_KIND, 1, win_info, inter_node_comm, ptr, win, ierr)
            call MPI_Win_shared_query(win, root, n, disp, ptr, ierr)
            if (disp /= 1 .or. n /= nbytes) call stop_all('mpi3_shared_alloc', 'Shared memory region not expected size.')
        else
            ! Allocate.
            call MPI_Win_allocate_shared(int(nbytes, MPI_ADDRESS_KIND), 1, win_info, inter_node_comm, ptr, win, ierr)
        end if
        call MPI_Info_free(win_info, ierr)
#else
        win = 0
        win_info = 0
        ierr = 0
        n = 0
        disp = 0
        call stop_all('mpi3_shared_alloc', 'MPI-3 is not available.')
#endif

    end subroutine mpi3_shared_alloc

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

#if defined ENABLE_SHMEM_POSIX
        ierr = 0
        call alloc_shared_posix(ptr, fstring_to_carray(A_name), N1*nbytes)
        call c_f_pointer(ptr, A, [N1])
#elif defined PARALLEL && ! defined DISABLE_MPI3
        call mpi3_shared_alloc(N1*nbytes, ptr)
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

#if defined ENABLE_SHMEM_POSIX
        ierr = 0
        call alloc_shared_posix(ptr, A_name, N1*N2*nbytes)
        call c_f_pointer(ptr, A, [N1, N2])
#elif defined PARALLEL && ! defined DISABLE_MPI3
        call mpi3_shared_alloc(N1*N2*nbytes, ptr)
        call c_f_pointer(ptr, A, [N1,N2])
#else
        allocate(A(N1, N2), stat=ierr)
        call check_allocate(A_name, N1*N2, ierr)
#endif

    end subroutine alloc_shared_real_dp_2D

    ! dealloc_shared_*

    subroutine dealloc_shared_real_dp_1D(A, A_name)

        use, intrinsic :: iso_c_binding, only: c_loc, c_ptr, c_null_ptr
        use const, only: dp, int_64
        use checking, only: check_deallocate
        use parallel

        real(dp), pointer, intent(inout) :: A(:)
        character(*), intent(in) :: A_name
        integer :: ierr
        integer, parameter :: nbytes = 8
        type(c_ptr) :: ptr
        ptr = c_null_ptr


#if defined ENABLE_SHMEM_POSIX
        ierr = 0
        ptr = c_loc(A)
        call free_shared_posix(ptr, A_name, size(A, kind=int_64)*nbytes)
        nullify(A)
#elif defined PARALLEL && ! defined DISABLE_MPI3
        ! [todo] - need to call MPI_Win_Free on the window object (which is not yet exposed outside of mpi3_shared_alloc).
        nullify(A)
#else
        deallocate(A, stat=ierr)
        call check_deallocate(A_name, ierr)
#endif

    end subroutine dealloc_shared_real_dp_1D

    subroutine dealloc_shared_real_dp_2D(A, A_name)

        use, intrinsic :: iso_c_binding, only: c_loc, c_ptr, c_null_ptr
        use const, only: dp, int_64
        use checking, only: check_deallocate
        use parallel

        real(dp), pointer, intent(inout) :: A(:,:)
        character(*), intent(in) :: A_name
        integer :: ierr
        integer, parameter :: nbytes = 8
        type(c_ptr) :: ptr
        ptr = c_null_ptr

#if defined ENABLE_SHMEM_POSIX
        ierr = 0
        ptr = c_loc(A)
        call free_shared_posix(ptr, A_name, size(A, kind=int_64)*nbytes)
        nullify(A)
#elif defined PARALLEL && ! defined DISABLE_MPI3
        ! [todo] - need to call MPI_Win_Free on the window object (which is not yet exposed outside of mpi3_shared_alloc).
        nullify(A)
#else
        deallocate(A, stat=ierr)
        call check_deallocate(A_name, ierr)
#endif

    end subroutine dealloc_shared_real_dp_2D

end module shmem

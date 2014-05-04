module hdf5_helper

    ! Helper routines (primarily for reading and writing) around the HDF5
    ! Fortran 2003 API.

    ! We are willing to do things repeatedly (e.g. opening and closing datasets)
    ! as the speed loss should be minimal and the gain in ease of coding quite
    ! nice...

#include "../../src/cdefs.h"

    use hdf5, only: hid_t

    implicit none

    private
    public :: hdf5_kinds_t, hdf5_kinds_init, hdf5_write, hdf5_read

    ! HDF5 kinds equivalent to the kinds defined in const.  Set in
    ! hdf5_init_kinds.
    type hdf5_kinds_t
        integer(hid_t) :: i32
        integer(hid_t) :: i64
        integer(hid_t) :: p
    end type hdf5_kinds_t

    interface hdf5_write
        module procedure write_string
        module procedure write_integer
        module procedure write_array_1d_int_32
        module procedure write_array_1d_int_64
        module procedure write_array_2d_int_32
        module procedure write_array_2d_int_64
        module procedure write_array_1d_real_p
        module procedure write_array_2d_real_p
    end interface hdf5_write

    interface hdf5_read
        module procedure read_integer
        module procedure read_array_1d_int_32
        module procedure read_array_1d_int_64
        module procedure read_array_2d_int_32
        module procedure read_array_2d_int_64
        module procedure read_array_1d_real_p
        module procedure read_array_2d_real_p
    end interface hdf5_read

    contains

        ! === Helper procedures: initialisation ===

        subroutine hdf5_kinds_init(kinds)

            use const, only: int_4, int_8, p
            use hdf5, only: H5_INTEGER_KIND, H5_REAL_KIND, h5kind_to_type

            type(hdf5_kinds_t), intent(out) :: kinds

            ! Convert our non-standard kinds to something HDF5 understands.
            ! [review] - JSS: should these kinds be parameters in const?
            ! [reply] - NSB: I tried this but you can't (unless I'm missing something?) do:
            ! [reply] - NSB: integer, parameter :: h5_int_4 = h5kind_to_type(int_4, H5_INTEGER_KIND)
            ! [reply] - NSB: so I had to remove the parameter statement and set h5_int_4 in a
            ! [reply] - NSB: routine like this one, so in the end I decided it wasn't worth it and
            ! [reply] - NSB: have changed it back.
            kinds%i32 = h5kind_to_type(int_4, H5_INTEGER_KIND)
            kinds%i64 = h5kind_to_type(int_8, H5_INTEGER_KIND)
            kinds%p = h5kind_to_type(p, H5_REAL_KIND)

        end subroutine hdf5_kinds_init

        ! === Helper procedures: writing ===

        subroutine write_string(id, dset, string)

            ! Write a string to an open HDF5 file/group.

            ! In:
            !    id: file or group HD5 identifier.
            !    dset: dataset name.
            !    string: string to write out.

            use hdf5

            integer(hid_t), intent(in) :: id
            character(*), intent(in) :: dset, string

            integer(hid_t) :: type_id, dspace_id, dset_id
            integer :: ierr

            ! Set up fortran string type of *this* length...
            call h5tcopy_f(H5T_FORTRAN_S1, type_id, ierr)
            call h5tset_size_f(type_id, len(string, SIZE_T), ierr)

            ! Create space and write string.
            call h5screate_f(H5S_SCALAR_F, dspace_id, ierr)
            call h5dcreate_f(id, dset, type_id, dspace_id, dset_id, ierr)
            call h5dwrite_f(dset_id, type_id, string, [0_HSIZE_T], ierr)
            call h5sclose_f(dspace_id, ierr)
            call h5dclose_f(dset_id, ierr)

            ! Release fortran string type.
            call h5tclose_f(type_id, ierr)

        end subroutine write_string

        subroutine write_integer(id, dset, val)

            ! Write an integer to an open HDF5 file/group.

            ! In:
            !    id: file or group HD5 identifier.
            !    dset: dataset name.
            !    val: integer to write out.

            use hdf5

            integer(hid_t), intent(in) :: id
            character(*), intent(in) :: dset
            integer, intent(in) :: val

            integer(hid_t) :: dspace_id, dset_id
            integer :: ierr

            call h5screate_f(H5S_SCALAR_F, dspace_id, ierr)
            call h5dcreate_f(id, dset, H5T_NATIVE_INTEGER, dspace_id, dset_id, ierr)

            call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, val, [0_HSIZE_T,0_HSIZE_T], ierr)

            call h5dclose_f(dset_id, ierr)
            call h5sclose_f(dspace_id, ierr)

        end subroutine write_integer

        ! I/O is not pure, so can't write elemental procedures...arse!  Please fill in with
        ! types and array dimensions as needed!

        subroutine write_array_1d_int_32(id, dset, kinds, arr_shape, arr)

            ! Write out 1D 32-bit integer array to an HDF5 file.

            ! In:
            !    id: file or group HD5 identifier.
            !    dset: dataset name.
            !    kinds: hdf5_kinds_t object containing the mapping between the non-default &
            !        kinds used in HANDE and HDF5 types.
            !    arr_shape: shape of array to be written (i.e. as given by shape(arr)).
            !    arr: array to be written.

            use, intrinsic :: iso_c_binding, only: c_ptr, c_loc
            use hdf5, only: hid_t, HSIZE_T
            use const, only: int_4

            integer(hid_t), intent(in) :: id
            character(*), intent(in) :: dset
            type(hdf5_kinds_t), intent(in) :: kinds
            integer, intent(in) :: arr_shape(:)
            ! Note assumed-shape arrays (e.g. arr(:)) are not C interoperable and hence
            ! cannot be passed to c_loc.
            integer(int_4), intent(in), target :: arr(arr_shape(1))

            type(c_ptr) :: ptr

            ptr = c_loc(arr)
            call write_ptr(id, dset, kinds%i32, size(arr_shape), int(arr_shape,HSIZE_T), ptr)

        end subroutine write_array_1d_int_32

        subroutine write_array_1d_int_64(id, dset, kinds, arr_shape, arr)

            ! Write out 1D 64-bit integer array to an HDF5 file.

            ! In:
            !    id: file or group HD5 identifier.
            !    dset: dataset name.
            !    kinds: hdf5_kinds_t object containing the mapping between the non-default &
            !        kinds used in HANDE and HDF5 types.
            !    arr_shape: shape of array to be written (i.e. as given by shape(arr)).
            !    arr: array to be written.

            use, intrinsic :: iso_c_binding, only: c_ptr, c_loc
            use hdf5, only: hid_t, HSIZE_T
            use const, only: int_8

            integer(hid_t), intent(in) :: id
            character(*), intent(in) :: dset
            type(hdf5_kinds_t), intent(in) :: kinds
            integer, intent(in) :: arr_shape(:)
            ! Note assumed-shape arrays (e.g. arr(:)) are not C interoperable and hence
            ! cannot be passed to c_loc.
            integer(int_8), intent(in), target :: arr(arr_shape(1))

            type(c_ptr) :: ptr

            ptr = c_loc(arr)
            call write_ptr(id, dset, kinds%i64, size(arr_shape), int(arr_shape,HSIZE_T), ptr)

        end subroutine write_array_1d_int_64

        subroutine write_array_2d_int_32(id, dset, kinds, arr_shape, arr)

            ! Write out 2D 32-bit integer array to an HDF5 file.

            ! In:
            !    id: file or group HD5 identifier.
            !    dset: dataset name.
            !    kinds: hdf5_kinds_t object containing the mapping between the non-default &
            !        kinds used in HANDE and HDF5 types.
            !    arr_shape: shape of array to be written (i.e. as given by shape(arr)).
            !    arr: array to be written.

            use, intrinsic :: iso_c_binding, only: c_ptr, c_loc
            use hdf5, only: hid_t, HSIZE_T
            use const, only: int_4

            integer(hid_t), intent(in) :: id
            character(*), intent(in) :: dset
            type(hdf5_kinds_t), intent(in) :: kinds
            integer, intent(in) :: arr_shape(:)
            ! Note assumed-shape arrays (e.g. arr(:)) are not C interoperable and hence
            ! cannot be passed to c_loc.
            integer(int_4), intent(in), target :: arr(arr_shape(1), arr_shape(2))

            type(c_ptr) :: ptr

            ptr = c_loc(arr)
            call write_ptr(id, dset, kinds%i32, size(arr_shape), int(arr_shape,HSIZE_T), ptr)

        end subroutine write_array_2d_int_32

        subroutine write_array_2d_int_64(id, dset, kinds, arr_shape, arr)

            ! Write out 2D 64-bit integer array to an HDF5 file.

            ! In:
            !    id: file or group HD5 identifier.
            !    dset: dataset name.
            !    kinds: hdf5_kinds_t object containing the mapping between the non-default &
            !        kinds used in HANDE and HDF5 types.
            !    arr_shape: shape of array to be written (i.e. as given by shape(arr)).
            !    arr: array to be written.

            use, intrinsic :: iso_c_binding, only: c_ptr, c_loc
            use hdf5, only: hid_t, HSIZE_T
            use const, only: int_8

            integer(hid_t), intent(in) :: id
            character(*), intent(in) :: dset
            type(hdf5_kinds_t), intent(in) :: kinds
            integer, intent(in) :: arr_shape(:)
            ! Note assumed-shape arrays (e.g. arr(:)) are not C interoperable and hence
            ! cannot be passed to c_loc.
            integer(int_8), intent(in), target :: arr(arr_shape(1), arr_shape(2))

            type(c_ptr) :: ptr

            ptr = c_loc(arr)
            call write_ptr(id, dset, kinds%i64, size(arr_shape), int(arr_shape,HSIZE_T), ptr)

        end subroutine write_array_2d_int_64

        subroutine write_array_1d_real_p(id, dset, kinds, arr_shape, arr)

            ! Write out 1D real(p) array to an HDF5 file.

            ! In:
            !    id: file or group HD5 identifier.
            !    dset: dataset name.
            !    kinds: hdf5_kinds_t object containing the mapping between the non-default &
            !        kinds used in HANDE and HDF5 types.
            !    arr_shape: shape of array to be written (i.e. as given by shape(arr)).
            !    arr: array to be written.

            use, intrinsic :: iso_c_binding, only: c_ptr, c_loc
            use hdf5, only: hid_t, HSIZE_T
            use const, only: p

            integer(hid_t), intent(in) :: id
            character(*), intent(in) :: dset
            type(hdf5_kinds_t), intent(in) :: kinds
            integer, intent(in) :: arr_shape(:)
            ! Note assumed-shape arrays (e.g. arr(:)) are not C interoperable and hence
            ! cannot be passed to c_loc.
            real(p), intent(in), target :: arr(arr_shape(1))

            type(c_ptr) :: ptr

            ptr = c_loc(arr)
            call write_ptr(id, dset, kinds%p, size(arr_shape), int(arr_shape, HSIZE_T), ptr)

        end subroutine write_array_1d_real_p

        subroutine write_array_2d_real_p(id, dset, kinds, arr_shape, arr)

            ! Write out 2D real(p) array to an HDF5 file.

            ! In:
            !    id: file or group HD5 identifier.
            !    dset: dataset name.
            !    kinds: hdf5_kinds_t object containing the mapping between the non-default &
            !        kinds used in HANDE and HDF5 types.
            !    arr_shape: shape of array to be written (i.e. as given by shape(arr)).
            !    arr: array to be written.

            use, intrinsic :: iso_c_binding, only: c_ptr, c_loc
            use hdf5, only: hid_t, HSIZE_T
            use const, only: p

            integer(hid_t), intent(in) :: id
            character(*), intent(in) :: dset
            type(hdf5_kinds_t), intent(in) :: kinds
            integer, intent(in) :: arr_shape(:)
            ! Note assumed-shape arrays (e.g. arr(:)) are not C interoperable and hence
            ! cannot be passed to c_loc.
            real(p), intent(in), target :: arr(arr_shape(1), arr_shape(2))

            type(c_ptr) :: ptr

            ptr = c_loc(arr)
            call write_ptr(id, dset, kinds%p, size(arr_shape), int(arr_shape, HSIZE_T), ptr)

        end subroutine write_array_2d_real_p

        subroutine write_ptr(id, dset, dtype, arr_rank, arr_dim, arr_ptr)

            ! Write an array to an open HDF5 file/group.

            ! In:
            !    id: file or group HD5 identifier.
            !    dset: dataset name.
            !    dtype: HDF5 data type of array.
            !    arr_rank: rank of array.
            !    arr_dim: size of array along each dimension.
            !    arr_ptr: C pointer to first element in array to be written out.

            ! NOTE: get dtype from h5kind_to_type if not using a native HDF5
            ! Fortran type.

            use hdf5
            use, intrinsic :: iso_c_binding

            integer(hid_t), intent(in) :: id
            character(*), intent(in) :: dset
            integer(hid_t), intent(in) :: dtype
            integer, intent(in) :: arr_rank
            integer(hsize_t), intent(in) :: arr_dim(:)
            type(c_ptr), intent(in) :: arr_ptr

            integer :: ierr
            integer(hid_t) :: dspace_id, dset_id

            call h5screate_simple_f(arr_rank, arr_dim, dspace_id, ierr)
            call h5dcreate_f(id, dset, dtype, dspace_id, dset_id, ierr)

            call h5dwrite_f(dset_id, dtype, arr_ptr, ierr)

            call h5dclose_f(dset_id, ierr)
            call h5sclose_f(dspace_id, ierr)

        end subroutine write_ptr

        ! === Helper procedures: reading ===

        subroutine read_integer(id, dset, val)

            ! Read an integer from an open HDF5 file/group.

            ! In:
            !    id: file or group HD5 identifier.
            !    dset: dataset name.
            ! Out:
            !    val: integer read from HDF5 file.

            use hdf5

            integer(hid_t), intent(in) :: id
            character(*), intent(in) :: dset
            integer, intent(out) :: val

            integer(hid_t) :: dset_id
            integer :: ierr

            call h5dopen_f(id, dset, dset_id, ierr)
            call h5dread_f(dset_id, H5T_NATIVE_INTEGER, val, [0_HSIZE_T,0_HSIZE_T], ierr)
            call h5dclose_f(dset_id, ierr)

        end subroutine read_integer

        ! I/O is not pure, so can't write elemental procedures...arse!  Please fill in with
        ! types and array dimensions as needed!

        subroutine read_array_1d_int_32(id, dset, kinds, arr_shape, arr)

            ! Read in 1D 32-bit integer array from an HDF5 file.

            ! In:
            !    id: file or group HD5 identifier.
            !    dset: dataset name.
            !    kinds: hdf5_kinds_t object containing the mapping between the non-default &
            !        kinds used in HANDE and HDF5 types.
            !    arr_shape: shape of array to be read (i.e. as given by shape(arr)).
            ! Out:
            !    arr: array to be read.

            use, intrinsic :: iso_c_binding, only: c_ptr, c_loc
            use hdf5, only: hid_t, HSIZE_T
            use const, only: int_4

            integer(hid_t), intent(in) :: id
            character(*), intent(in) :: dset
            type(hdf5_kinds_t), intent(in) :: kinds
            integer, intent(in) :: arr_shape(:)
            integer(int_4), intent(out), target :: arr(arr_shape(1))

            type(c_ptr) :: ptr

            ptr = c_loc(arr)
            call read_ptr(id, dset, kinds%i32, size(arr_shape), int(arr_shape,HSIZE_T), ptr)

        end subroutine read_array_1d_int_32

        subroutine read_array_1d_int_64(id, dset, kinds, arr_shape, arr)

            ! Read in 1D 64-bit integer array from an HDF5 file.

            ! In:
            !    id: file or group HD5 identifier.
            !    dset: dataset name.
            !    kinds: hdf5_kinds_t object containing the mapping between the non-default &
            !        kinds used in HANDE and HDF5 types.
            !    arr_shape: shape of array to be read (i.e. as given by shape(arr)).
            ! Out:
            !    arr: array to be read.

            use, intrinsic :: iso_c_binding, only: c_ptr, c_loc
            use hdf5, only: hid_t, HSIZE_T
            use const, only: int_8

            integer(hid_t), intent(in) :: id
            character(*), intent(in) :: dset
            type(hdf5_kinds_t), intent(in) :: kinds
            integer, intent(in) :: arr_shape(:)
            integer(int_8), intent(out), target :: arr(arr_shape(1))

            type(c_ptr) :: ptr

            ptr = c_loc(arr)
            call read_ptr(id, dset, kinds%i64, size(arr_shape), int(arr_shape,HSIZE_T), ptr)

        end subroutine read_array_1d_int_64

        subroutine read_array_2d_int_32(id, dset, kinds, arr_shape, arr)

            ! Read in 2D 32-bit integer array from an HDF5 file.

            ! In:
            !    id: file or group HD5 identifier.
            !    dset: dataset name.
            !    kinds: hdf5_kinds_t object containing the mapping between the non-default &
            !        kinds used in HANDE and HDF5 types.
            !    arr_shape: shape of array to be read (i.e. as given by shape(arr)).
            ! Out:
            !    arr: array to be read.

            use, intrinsic :: iso_c_binding, only: c_ptr, c_loc
            use hdf5, only: hid_t, HSIZE_T
            use const, only: int_4

            integer(hid_t), intent(in) :: id
            character(*), intent(in) :: dset
            type(hdf5_kinds_t), intent(in) :: kinds
            integer, intent(in) :: arr_shape(:)
            integer(int_4), intent(out), target :: arr(arr_shape(1), arr_shape(2))

            type(c_ptr) :: ptr

            ptr = c_loc(arr)
            call read_ptr(id, dset, kinds%i32, size(arr_shape), int(arr_shape,HSIZE_T), ptr)

        end subroutine read_array_2d_int_32

        subroutine read_array_2d_int_64(id, dset, kinds, arr_shape, arr)

            ! Read in 2D 64-bit integer array from an HDF5 file.

            ! In:
            !    id: file or group HD5 identifier.
            !    dset: dataset name.
            !    kinds: hdf5_kinds_t object containing the mapping between the non-default &
            !        kinds used in HANDE and HDF5 types.
            !    arr_shape: shape of array to be read (i.e. as given by shape(arr)).
            ! Out:
            !    arr: array to be read.

            use, intrinsic :: iso_c_binding, only: c_ptr, c_loc
            use hdf5, only: hid_t, HSIZE_T
            use const, only: int_8

            integer(hid_t), intent(in) :: id
            character(*), intent(in) :: dset
            type(hdf5_kinds_t), intent(in) :: kinds
            integer, intent(in) :: arr_shape(:)
            integer(int_8), intent(out), target :: arr(arr_shape(1), arr_shape(2))

            type(c_ptr) :: ptr

            ptr = c_loc(arr)
            call read_ptr(id, dset, kinds%i64, size(arr_shape), int(arr_shape,HSIZE_T), ptr)

        end subroutine read_array_2d_int_64

        subroutine read_array_1d_real_p(id, dset, kinds, arr_shape, arr)

            ! Read in 1D real(p) array from an HDF5 file.

            ! In:
            !    id: file or group HD5 identifier.
            !    dset: dataset name.
            !    kinds: hdf5_kinds_t object containing the mapping between the non-default &
            !        kinds used in HANDE and HDF5 types.
            !    arr_shape: shape of array to be read (i.e. as given by shape(arr)).
            ! Out:
            !    arr: array to be read.

            use, intrinsic :: iso_c_binding, only: c_ptr, c_loc
            use hdf5, only: hid_t, HSIZE_T
            use const, only: p

            integer(hid_t), intent(in) :: id
            character(*), intent(in) :: dset
            type(hdf5_kinds_t), intent(in) :: kinds
            integer, intent(in) :: arr_shape(:)
            real(p), intent(out), target :: arr(arr_shape(1))

            type(c_ptr) :: ptr

            ptr = c_loc(arr)
            call read_ptr(id, dset, kinds%p, size(arr_shape), int(arr_shape, HSIZE_T), ptr)

        end subroutine read_array_1d_real_p

        subroutine read_array_2d_real_p(id, dset, kinds, arr_shape, arr)

            ! Read in 2D real(p) array from an HDF5 file.

            ! In:
            !    id: file or group HD5 identifier.
            !    dset: dataset name.
            !    kinds: hdf5_kinds_t object containing the mapping between the non-default &
            !        kinds used in HANDE and HDF5 types.
            !    arr_shape: shape of array to be read (i.e. as given by shape(arr)).
            ! Out:
            !    arr: array to be read.

            use, intrinsic :: iso_c_binding, only: c_ptr, c_loc
            use hdf5, only: hid_t, HSIZE_T
            use const, only: p

            integer(hid_t), intent(in) :: id
            character(*), intent(in) :: dset
            type(hdf5_kinds_t), intent(in) :: kinds
            integer, intent(in) :: arr_shape(:)
            real(p), intent(out), target :: arr(arr_shape(1), arr_shape(2))

            type(c_ptr) :: ptr

            ptr = c_loc(arr)
            call read_ptr(id, dset, kinds%p, size(arr_shape), int(arr_shape, HSIZE_T), ptr)

        end subroutine read_array_2d_real_p

        subroutine read_ptr(id, dset, dtype, arr_rank, arr_dim, arr_ptr)

            ! Read an array from an open HDF5 file/group.

            ! In:
            !    id: file or group HD5 identifier.
            !    dset: dataset name.
            !    dtype: HDF5 data type of array.
            !    arr_rank: rank of array.
            !    arr_dim: size of array along each dimension.
            ! In/Out:
            !    arr_ptr: C pointer to first element in array to read.  On
            !        output, the dataset is store in the array pointed to by
            !        arr_ptr.

            ! NOTE: get dtype from h5kind_to_type if not using a native HDF5
            ! Fortran type.

            use, intrinsic :: iso_c_binding
            use hdf5

            use errors, only: stop_all

            integer(hid_t), intent(in) :: id
            character(*), intent(in) :: dset
            integer(hid_t), intent(in) :: dtype
            integer, intent(in) :: arr_rank
            integer(hsize_t), intent(in) :: arr_dim(:)
            type(c_ptr), intent(inout) :: arr_ptr

            integer :: ierr, arr_rank_stored
            integer(hsize_t) :: arr_dim_stored(size(arr_dim)), arr_max_dim_stored(size(arr_dim))
            integer(hid_t) :: dspace_id, dset_id, dtype_stored
            logical :: dtype_equal

            call h5dopen_f(id, dset, dset_id, ierr)

            call h5dget_type_f(dset_id, dtype_stored, ierr)
            call h5tequal_f(dtype, dtype_stored, dtype_equal, ierr)
            if (.not.dtype_equal) call stop_all('read_ptr', 'Reading mismatched data type.')

            call h5dget_space_f(dset_id, dspace_id, ierr)

            call h5sget_simple_extent_ndims_f(dspace_id, arr_rank_stored, ierr)
            if (arr_rank /= arr_rank_stored) call stop_all('read_ptr', 'Reading mismatched data rank.')

            call h5sget_simple_extent_dims_f(dspace_id, arr_dim_stored, arr_max_dim_stored, ierr)
            if (any(arr_dim(:arr_rank-1) - arr_dim_stored(:arr_rank-1) /= 0)) &
                call stop_all('read_ptr', 'Reading mismatched array bounds')
            if (arr_dim(arr_rank) - arr_dim_stored(arr_rank) < 0) &
                call stop_all('read_ptr', 'Reading mismatched array bounds')

            call h5sclose_f(dspace_id, ierr)

            call h5dread_f(dset_id, dtype, arr_ptr, ierr)
            call h5dclose_f(dset_id, ierr)

        end subroutine read_ptr

end module hdf5_helper

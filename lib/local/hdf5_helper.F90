module hdf5_helper

    ! Helper routines (primarily for reading and writing) around the HDF5
    ! Fortran 2003 API.

    !!! WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING !!!
    !!! The HDF5 library *must* be opened (using h5open_f) before any HDF5      !!!
    !!! procedures (including procedures below) are called.                     !!!
    !!! WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING !!!

    ! In general, you should do something like:

    ! call h5open_f(ierr)
    ! call hdf5_kinds_init(hdf5_kinds)
    ! <HDF5 calls to read/write/etc>
    ! call h5close_f(ierr)

    ! We are willing to do things repeatedly (e.g. opening and closing datasets)
    ! as the speed loss should be minimal and the gain in ease of coding quite
    ! nice...

    ! [todo] - Check error flags returned by HDF5 procedures.

#ifndef DISABLE_HDF5

    use hdf5, only: hid_t

    implicit none

    private
    public :: hdf5_kinds_t, hdf5_kinds_init, hdf5_write, hdf5_read, dtype_equal, dset_shape, hdf5_path


    ! HDF5 kinds equivalent to the kinds defined in const.  Set in
    ! hdf5_kinds_init.
    type hdf5_kinds_t
        integer(hid_t) :: i32
        integer(hid_t) :: i64
        integer(hid_t) :: int_p
        integer(hid_t) :: i0
        integer(hid_t) :: sp
        integer(hid_t) :: dp
    end type hdf5_kinds_t

    interface hdf5_write
        module procedure write_string
        module procedure write_integer
        module procedure write_boolean
        module procedure write_array_1d_int_32
        module procedure write_array_1d_int_64
        module procedure write_array_2d_int_32
        module procedure write_array_2d_int_64
        module procedure write_array_1d_real_sp
        module procedure write_array_1d_real_dp
        module procedure write_array_2d_real_sp
        module procedure write_array_2d_real_dp
    end interface hdf5_write

    interface hdf5_read
        module procedure read_string
        module procedure read_integer
        module procedure read_boolean
        module procedure read_array_1d_int_32
        module procedure read_array_1d_int_64
        module procedure read_array_2d_int_32
        module procedure read_array_2d_int_64
        module procedure read_array_1d_real_sp
        module procedure read_array_1d_real_dp
        module procedure read_array_2d_real_sp
        module procedure read_array_2d_real_dp
    end interface hdf5_read

    interface hdf5_path
        module procedure hdf5_path_2
        module procedure hdf5_path_3
        module procedure hdf5_path_4
    end interface hdf5_path

    integer, parameter :: hande_hdf5_false = 0, hande_hdf5_true = 1

    contains

        ! === Helper procedures: initialisation, properties ===

        subroutine hdf5_kinds_init(kinds)

            ! Out:
            !    kinds: HDF5 kinds corresponding to integer and real kinds used in the QMC algorithms.

            ! WARNING: Take care.  If called before h5open_f is called, then it will return garbage.
            !          Further, the kind values are *not* constant between closing and then re-opening
            !          the HDF5 library.  Learn from my (painful) experiences...

            use const, only: int_32, int_64, int_p, i0, sp, dp
            use hdf5, only: H5_INTEGER_KIND, H5_REAL_KIND, h5kind_to_type

            type(hdf5_kinds_t), intent(out) :: kinds

            ! Convert our non-standard kinds to something HDF5 understands.
            kinds%i32 = h5kind_to_type(int_32, H5_INTEGER_KIND)
            kinds%i64 = h5kind_to_type(int_64, H5_INTEGER_KIND)
            kinds%int_p = h5kind_to_type(int_p, H5_INTEGER_KIND)
            kinds%i0 = h5kind_to_type(i0, H5_INTEGER_KIND)
            kinds%sp = h5kind_to_type(sp, H5_REAL_KIND)
            kinds%dp = h5kind_to_type(dp, H5_REAL_KIND)

        end subroutine hdf5_kinds_init

        subroutine dataset_array_plist(arr_rank, arr_dim, plist_id, chunk_size, compress_lvl, chunk)

            ! Create a properties list for an array data set.

            ! In:
            !    arr_rank: rank of array.
            !    arr_dim: size of array along each dimension.
            !    chunk_size (optional): size of chunks to stored when using a chunked
            !        layout.  We currently chunk solely the last dimension.  Default: 100000.
            !    compress_lvl (optional): compression level (1-9).  Default: 6.
            !        The levels are those used by gzip (see man page) where 1 is fastest
            !        with the worst compression ratio and 9 is the slowest but
            !        gives (usually) the best compression ratio.  6 is the same default
            !        used by gzip.
            !    chunk (optional) :: if true, force chunking of the data set.  Default: chunk
            !        only arrays with more elements than chunk_size.
            ! Out:
            !    plist_id: properties list.  If the array has more than chunk_size entries,
            !        chunking and compression are enabled.  Currently we use the default
            !        compression (gzip) though other compression filters (usually via
            !        third-party plugins) are available.

            ! NOTE: the properties list should be closed after use with h5pclose_f in order
            ! to avoid a HDF5 resource leak.

            use hdf5, only: hsize_t, h5pcreate_f, h5pset_chunk_f, h5pset_deflate_f, H5P_DATASET_CREATE_F

            integer(hid_t), intent(in) :: arr_rank
            integer(hsize_t), intent(in) :: arr_dim(:)
            integer(hid_t), intent(out) :: plist_id
            integer, intent(in), optional :: chunk_size, compress_lvl
            logical, intent(in), optional :: chunk

            integer :: ierr
            integer :: chunk_size_loc
            integer :: compress_lvl_loc
            integer(hsize_t) :: chunk_dim(arr_rank)
            logical :: chunk_loc

            chunk_size_loc = 100000
            compress_lvl_loc = 6
            chunk_loc = .false.
            if (present(chunk_size)) chunk_size_loc = chunk_size
            if (present(compress_lvl)) compress_lvl_loc = compress_lvl
            if (present(chunk)) chunk_loc = chunk

            call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, ierr)

            if (product(arr_dim) > chunk_size_loc .or. chunk_loc) then
                chunk_dim(1:arr_rank-1) = arr_dim(1:arr_rank-1)
                chunk_dim(arr_rank) = chunk_size_loc / product(chunk_dim(1:arr_rank-1))
                call h5pset_chunk_f(plist_id, arr_rank, chunk_dim, ierr)
                call h5pset_deflate_f(plist_id, compress_lvl_loc, ierr)
            end if

        end subroutine dataset_array_plist

        function dtype_equal(id, dset, dtype)

            ! In:
            !    id: file or group HD5 identifier.
            !    dset: dataset name.
            !    dtype: HDF5 data type of array.
            ! Returns:
            !    True if dtype matches the data type of the dataset dset and false otherwise.

            use hdf5, only: hid_t, h5dopen_f, h5dget_type_f, h5tequal_f, h5dclose_f

            logical :: dtype_equal
            integer(hid_t), intent(in) :: id
            character(*), intent(in) :: dset
            integer(hid_t), intent(in) :: dtype

            integer :: ierr
            integer(hid_t) :: dset_id, dtype_stored

            call h5dopen_f(id, dset, dset_id, ierr)
            call h5dget_type_f(dset_id, dtype_stored, ierr)
            call h5tequal_f(dtype, dtype_stored, dtype_equal, ierr)
            call h5dclose_f(dset_id, ierr)

        end function dtype_equal

        subroutine dset_shape(id, dset, arr_shape)

            ! Find the shape of a dataset in an HDF5 file.

            ! In:
            !    id: file or group HD5 identifier.
            !    dset: dataset name.
            ! Out:
            !    arr_shape: shape of array stored in dataset.  Must have dimensions
            !        equal to the rank of the dataset (not checked).

            use hdf5

            integer(hid_t), intent(in) :: id
            character(*), intent(in) :: dset
            integer(hsize_t), intent(out) :: arr_shape(:)

            integer :: ierr
            integer(hid_t) :: dset_id, dspace_id
            integer(hsize_t) :: arr_max_shape(size(arr_shape))

            call h5dopen_f(id, dset, dset_id, ierr)
            call h5dget_space_f(dset_id, dspace_id, ierr)
            call h5sget_simple_extent_dims_f(dspace_id, arr_shape, arr_max_shape, ierr)
            call h5sclose_f(dspace_id, ierr)
            call h5dclose_f(dset_id, ierr)

        end subroutine dset_shape

        ! === Helper procedures: combining group/data/link names ===

        pure function hdf5_path_2(l1, l2) result(path)

            ! In:
            !    l1, l2: group or dataspace or link names, where l2 is inside l1.
            ! Returns:
            !    l1/l2, ie the complete HDF5 path to l2 from the group containing l1.

            character(*), intent(in) :: l1, l2
            character(len=len(l1)+len(l2)+1) :: path

            path = l1 // '/' // l2

        end function hdf5_path_2

        pure function hdf5_path_3(l1, l2, l3) result(path)

            ! In:
            !    l1, l2, l3: group or dataspace or link names, where l2 is inside l1 and l3 inside l2.
            ! Returns:
            !    l1/l2/l3, ie the complete HDF5 path to l3 from the group containing l1.

            character(*), intent(in) :: l1, l2, l3
            character(len=len(l1)+len(l2)+len(l3)+2) :: path

            path = hdf5_path(l1, l2)
            path = hdf5_path(trim(path), l3)

        end function hdf5_path_3

        pure function hdf5_path_4(l1, l2, l3, l4) result(path)

            ! In:
            !    l1, l2, l3, l4: group or dataspace or link names, where l2 is inside l1, l3 inside l2 and l4 inside l3.
            ! Returns:
            !    l1/l2/l3, ie the complete HDF5 path to l4 from the group containing l1.

            character(*), intent(in) :: l1, l2, l3, l4
            character(len=len(l1)+len(l2)+len(l3)+len(l4)+3) :: path

            path = hdf5_path(l1, l2, l3)
            path = hdf5_path(trim(path), l4)

        end function hdf5_path_4

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

        subroutine write_boolean(id, dset, val)

            ! Write a boolean to an open HDF5 file/group.

            ! In:
            !    id: file or group HD5 identifier.
            !    dset: dataset name.
            !    val: boolean to write out.

            ! NOTE: HDF5 can't handle boolean types, so instead we write out an integer
            ! which we interpret ourselves in a consistent fashion.

            integer(hid_t), intent(in) :: id
            character(*), intent(in) :: dset
            logical, intent(in) :: val

            if (val) then
                call hdf5_write(id, dset, hande_hdf5_true)
            else
                call hdf5_write(id, dset, hande_hdf5_false)
            end if

        end subroutine write_boolean

        ! I/O is not pure, so can't write elemental procedures...arse!  Please fill in with
        ! types and array dimensions as needed!

        subroutine write_array_1d_int_32(id, dset, kinds, arr_shape, arr, append)

            ! Write out 1D 32-bit integer array to an HDF5 file.

            ! In:
            !    id: file or group HD5 identifier.
            !    dset: dataset name.
            !    kinds: hdf5_kinds_t object containing the mapping between the non-default &
            !        kinds used in HANDE and HDF5 types.
            !    arr_shape: shape of array to be written (i.e. as given by shape(arr)).
            !    arr: array to be written.
            !    append (optional): if true, append to an existing dataset or create an
            !        extensible dataset if it doesn't already exist.

            use, intrinsic :: iso_c_binding, only: c_ptr, c_loc
            use hdf5, only: hid_t, HSIZE_T
            use const, only: int_32

            integer, parameter :: rank = 1

            integer(hid_t), intent(in) :: id
            character(*), intent(in) :: dset
            type(hdf5_kinds_t), intent(in) :: kinds
            integer, intent(in) :: arr_shape(:)
            ! Note assumed-shape arrays (e.g. arr(:)) are not C interoperable and hence
            ! cannot be passed to c_loc.
            integer(int_32), intent(in), target :: arr(arr_shape(1))
            logical, intent(in), optional :: append

            type(c_ptr) :: ptr

            ptr = c_loc(arr)
            call write_ptr(id, dset, kinds%i32, rank, int(arr_shape,HSIZE_T), ptr, append)

        end subroutine write_array_1d_int_32

        subroutine write_array_1d_int_64(id, dset, kinds, arr_shape, arr, append)

            ! Write out 1D 64-bit integer array to an HDF5 file.

            ! In:
            !    id: file or group HD5 identifier.
            !    dset: dataset name.
            !    kinds: hdf5_kinds_t object containing the mapping between the non-default &
            !        kinds used in HANDE and HDF5 types.
            !    arr_shape: shape of array to be written (i.e. as given by shape(arr)).
            !    arr: array to be written.
            !    append (optional): if true, append to an existing dataset or create an
            !        extensible dataset if it doesn't already exist.

            use, intrinsic :: iso_c_binding, only: c_ptr, c_loc
            use hdf5, only: hid_t, HSIZE_T
            use const, only: int_64

            integer, parameter :: rank = 1

            integer(hid_t), intent(in) :: id
            character(*), intent(in) :: dset
            type(hdf5_kinds_t), intent(in) :: kinds
            integer, intent(in) :: arr_shape(:)
            ! Note assumed-shape arrays (e.g. arr(:)) are not C interoperable and hence
            ! cannot be passed to c_loc.
            integer(int_64), intent(in), target :: arr(arr_shape(1))
            logical, intent(in), optional :: append

            type(c_ptr) :: ptr

            ptr = c_loc(arr)
            call write_ptr(id, dset, kinds%i64, rank, int(arr_shape,HSIZE_T), ptr, append)

        end subroutine write_array_1d_int_64

        subroutine write_array_2d_int_32(id, dset, kinds, arr_shape, arr, lim, append)

            ! Write out 2D 32-bit integer array to an HDF5 file.

            ! In:
            !    id: file or group HD5 identifier.
            !    dset: dataset name.
            !    kinds: hdf5_kinds_t object containing the mapping between the non-default &
            !        kinds used in HANDE and HDF5 types.
            !    arr_shape: shape of array to be written (i.e. as given by shape(arr)).
            !    arr: array to be written.
            !    lim (optional): Write out only arr(:,:lim).  Default: full array (as given
            !        in arr_shape) is written out.
            !    append (optional): if true, append to an existing dataset or create an
            !        extensible dataset if it doesn't already exist.

            ! NOTE: it is safer to use lim rather than to pass in an array slice (e.g. arr(:,:lim))
            ! which can trigger an expensive array temporary with some compilers.

            use, intrinsic :: iso_c_binding, only: c_ptr, c_loc
            use hdf5, only: hid_t, HSIZE_T
            use const, only: int_32

            integer, parameter :: rank = 2

            integer(hid_t), intent(in) :: id
            character(*), intent(in) :: dset
            type(hdf5_kinds_t), intent(in) :: kinds
            integer, intent(in) :: arr_shape(:)
            ! Note assumed-shape arrays (e.g. arr(:)) are not C interoperable and hence
            ! cannot be passed to c_loc.
            integer(int_32), intent(in), target :: arr(arr_shape(1), arr_shape(2))
            integer, intent(in), optional :: lim
            logical, intent(in), optional :: append

            integer(HSIZE_T) :: arr_dims(size(arr_shape))
            type(c_ptr) :: ptr

            ptr = c_loc(arr)
            arr_dims = int(arr_shape, HSIZE_T)
            if (present(lim)) arr_dims(rank) = lim
            call write_ptr(id, dset, kinds%i32, rank, arr_dims, ptr, append)

        end subroutine write_array_2d_int_32

        subroutine write_array_2d_int_64(id, dset, kinds, arr_shape, arr, lim, append)

            ! Write out 2D 64-bit integer array to an HDF5 file.

            ! In:
            !    id: file or group HD5 identifier.
            !    dset: dataset name.
            !    kinds: hdf5_kinds_t object containing the mapping between the non-default &
            !        kinds used in HANDE and HDF5 types.
            !    arr_shape: shape of array to be written (i.e. as given by shape(arr)).
            !    arr: array to be written.
            !    lim (optional): Write out only arr(:,:lim).  Default: full array (as given
            !        in arr_shape) is written out.
            !    append (optional): if true, append to an existing dataset or create an
            !        extensible dataset if it doesn't already exist.

            ! NOTE: it is safer to use lim rather than to pass in an array slice (e.g. arr(:,:lim))
            ! which can trigger an expensive array temporary with some compilers.

            use, intrinsic :: iso_c_binding, only: c_ptr, c_loc
            use hdf5, only: hid_t, HSIZE_T
            use const, only: int_64

            integer, parameter :: rank = 2

            integer(hid_t), intent(in) :: id
            character(*), intent(in) :: dset
            type(hdf5_kinds_t), intent(in) :: kinds
            integer, intent(in) :: arr_shape(:)
            ! Note assumed-shape arrays (e.g. arr(:)) are not C interoperable and hence
            ! cannot be passed to c_loc.
            integer(int_64), intent(in), target :: arr(arr_shape(1), arr_shape(2))
            integer, intent(in), optional :: lim
            logical, intent(in), optional :: append

            integer(HSIZE_T) :: arr_dims(size(arr_shape))
            type(c_ptr) :: ptr

            ptr = c_loc(arr)
            arr_dims = int(arr_shape, HSIZE_T)
            if (present(lim)) arr_dims(rank) = lim
            call write_ptr(id, dset, kinds%i64, rank, arr_dims, ptr, append)

        end subroutine write_array_2d_int_64

        subroutine write_array_1d_real_dp(id, dset, kinds, arr_shape, arr, append)

            ! Write out 1D real(dp) array to an HDF5 file.

            ! In:
            !    id: file or group HD5 identifier.
            !    dset: dataset name.
            !    kinds: hdf5_kinds_t object containing the mapping between the non-default &
            !        kinds used in HANDE and HDF5 types.
            !    arr_shape: shape of array to be written (i.e. as given by shape(arr)).
            !    arr: array to be written.
            !    append (optional): if true, append to an existing dataset or create an
            !        extensible dataset if it doesn't already exist.

            use, intrinsic :: iso_c_binding, only: c_ptr, c_loc
            use hdf5, only: hid_t, HSIZE_T
            use const, only: dp

            integer, parameter :: rank = 1

            integer(hid_t), intent(in) :: id
            character(*), intent(in) :: dset
            type(hdf5_kinds_t), intent(in) :: kinds
            integer, intent(in) :: arr_shape(:)
            ! Note assumed-shape arrays (e.g. arr(:)) are not C interoperable and hence
            ! cannot be passed to c_loc.
            real(dp), intent(in), target :: arr(arr_shape(1))
            logical, intent(in), optional :: append

            type(c_ptr) :: ptr

            ptr = c_loc(arr)
            call write_ptr(id, dset, kinds%dp, rank, int(arr_shape, HSIZE_T), ptr, append)

        end subroutine write_array_1d_real_dp

        subroutine write_array_2d_real_dp(id, dset, kinds, arr_shape, arr, lim, append)

            ! Write out 2D real(dp) array to an HDF5 file.

            ! In:
            !    id: file or group HD5 identifier.
            !    dset: dataset name.
            !    kinds: hdf5_kinds_t object containing the mapping between the non-default &
            !        kinds used in HANDE and HDF5 types.
            !    arr_shape: shape of array to be written (i.e. as given by shape(arr)).
            !    arr: array to be written.
            !    lim (optional): Write out only arr(:,:lim).  Default: full array (as given
            !        in arr_shape) is written out.
            !    append (optional): if true, append to an existing dataset or create an
            !        extensible dataset if it doesn't already exist.

            ! NOTE: it is safer to use lim rather than to pass in an array slice (e.g. arr(:,:lim))
            ! which can trigger an expensive array temporary with some compilers.

            use, intrinsic :: iso_c_binding, only: c_ptr, c_loc
            use hdf5, only: hid_t, HSIZE_T
            use const, only: dp

            integer, parameter :: rank = 2

            integer(hid_t), intent(in) :: id
            character(*), intent(in) :: dset
            type(hdf5_kinds_t), intent(in) :: kinds
            integer, intent(in) :: arr_shape(:)
            ! Note assumed-shape arrays (e.g. arr(:)) are not C interoperable and hence
            ! cannot be passed to c_loc.
            real(dp), intent(in), target :: arr(arr_shape(1), arr_shape(2))
            integer, intent(in), optional :: lim
            logical, intent(in), optional :: append

            integer(HSIZE_T) :: arr_dims(size(arr_shape))
            type(c_ptr) :: ptr

            ptr = c_loc(arr)
            arr_dims = int(arr_shape, HSIZE_T)
            if (present(lim)) arr_dims(rank) = lim
            call write_ptr(id, dset, kinds%dp, rank, arr_dims, ptr, append)

        end subroutine write_array_2d_real_dp

        subroutine write_array_1d_real_sp(id, dset, kinds, arr_shape, arr, append)

            ! Write out 1D real(sp) array to an HDF5 file.

            ! In:
            !    id: file or group HD5 identifier.
            !    dset: dataset name.
            !    kinds: hdf5_kinds_t object containing the mapping between the non-default &
            !        kinds used in HANDE and HDF5 types.
            !    arr_shape: shape of array to be written (i.e. as given by shape(arr)).
            !    arr: array to be written.
            !    append (optional): if true, append to an existing dataset or create an
            !        extensible dataset if it doesn't already exist.

            use, intrinsic :: iso_c_binding, only: c_ptr, c_loc
            use hdf5, only: hid_t, HSIZE_T
            use const, only: sp

            integer, parameter :: rank = 1

            integer(hid_t), intent(in) :: id
            character(*), intent(in) :: dset
            type(hdf5_kinds_t), intent(in) :: kinds
            integer, intent(in) :: arr_shape(:)
            ! Note assumed-shape arrays (e.g. arr(:)) are not C interoperable and hence
            ! cannot be passed to c_loc.
            real(sp), intent(in), target :: arr(arr_shape(1))
            logical, intent(in), optional :: append

            type(c_ptr) :: ptr

            ptr = c_loc(arr)
            call write_ptr(id, dset, kinds%sp, rank, int(arr_shape, HSIZE_T), ptr, append)

        end subroutine write_array_1d_real_sp

        subroutine write_array_2d_real_sp(id, dset, kinds, arr_shape, arr, lim, append)

            ! Write out 2D real(sp) array to an HDF5 file.

            ! In:
            !    id: file or group HD5 identifier.
            !    dset: dataset name.
            !    kinds: hdf5_kinds_t object containing the mapping between the non-default &
            !        kinds used in HANDE and HDF5 types.
            !    arr_shape: shape of array to be written (i.e. as given by shape(arr)).
            !    arr: array to be written.
            !    lim (optional): Write out only arr(:,:lim).  Default: full array (as given
            !        in arr_shape) is written out.
            !    append (optional): if true, append to an existing dataset or create an
            !        extensible dataset if it doesn't already exist.

            ! NOTE: it is safer to use lim rather than to pass in an array slice (e.g. arr(:,:lim))
            ! which can trigger an expensive array temporary with some compilers.

            use, intrinsic :: iso_c_binding, only: c_ptr, c_loc
            use hdf5, only: hid_t, HSIZE_T
            use const, only: sp

            integer, parameter :: rank = 2

            integer(hid_t), intent(in) :: id
            character(*), intent(in) :: dset
            type(hdf5_kinds_t), intent(in) :: kinds
            integer, intent(in) :: arr_shape(:)
            ! Note assumed-shape arrays (e.g. arr(:)) are not C interoperable and hence
            ! cannot be passed to c_loc.
            real(sp), intent(in), target :: arr(arr_shape(1), arr_shape(2))
            integer, intent(in), optional :: lim
            logical, intent(in), optional :: append

            integer(HSIZE_T) :: arr_dims(size(arr_shape))
            type(c_ptr) :: ptr

            ptr = c_loc(arr)
            arr_dims = int(arr_shape, HSIZE_T)
            if (present(lim)) arr_dims(rank) = lim
            call write_ptr(id, dset, kinds%sp, rank, arr_dims, ptr, append)

        end subroutine write_array_2d_real_sp

        subroutine write_ptr(id, dset, dtype, arr_rank, arr_dim, arr_ptr, arr_append)

            ! Write an array to an open HDF5 file/group.

            ! In:
            !    id: file or group HD5 identifier.
            !    dset: dataset name.
            !    dtype: HDF5 data type of array.
            !    arr_rank: rank of array.
            !    arr_dim: size of array along each dimension.
            !    arr_ptr: C pointer to first element in array to be written out.
            !    arr_append (optional): if true, append to an existing data set (created
            !        if necessary).  Appending is currently only supported along the
            !        outer-most index.

            ! NOTE: get dtype from h5kind_to_type if not using a native HDF5
            ! Fortran type.

            use hdf5
            use, intrinsic :: iso_c_binding
            use errors, only: stop_all

            integer(hid_t), intent(in) :: id
            character(*), intent(in) :: dset
            integer(hid_t), intent(in) :: dtype
            integer, intent(in) :: arr_rank
            integer(hsize_t), intent(in) :: arr_dim(:)
            type(c_ptr), intent(in) :: arr_ptr
            logical, intent(in), optional :: arr_append

            integer :: ierr
            integer(hid_t) :: dspace_id, dset_id, plist_id, memspace_id
            integer :: arr_rank_stored
            integer(hsize_t) :: arr_max_dim(size(arr_dim)), arr_dim_stored(size(arr_dim)), &
                                arr_max_dim_stored(size(arr_dim)), offset(size(arr_dim)), &
                                arr_tot_dim(size(arr_dim))
            logical :: append, exists, chunk

            append = .false.
            arr_max_dim = arr_dim

            if (present(arr_append)) append = arr_append
            if (append) then
                call h5lexists_f(id, dset, exists, ierr)
                if (.not.exists) then
                    ! Need to create a new extensible space.
                    append = .false.
                end if
                arr_max_dim(size(arr_dim)) = H5S_UNLIMITED_F
            end if

            if (append) then
                call h5dopen_f(id, dset, dset_id, ierr)
                call h5dget_space_f(dset_id, dspace_id, ierr)

                ! Get existing size.
                call h5sget_simple_extent_ndims_f(dspace_id, arr_rank_stored, ierr)
                if (arr_rank /= arr_rank_stored) call stop_all('write_ptr', 'Reading mismatched data rank.')
                call h5sget_simple_extent_dims_f(dspace_id, arr_dim_stored, arr_max_dim_stored, ierr)
                if (any(arr_dim(:arr_rank-1) - arr_dim_stored(:arr_rank-1) /= 0)) &
                    call stop_all('write_ptr', 'Reading mismatched array bounds')

                ! Extend.
                arr_tot_dim = arr_dim_stored
                arr_tot_dim(size(arr_dim)) = arr_dim_stored(size(arr_dim)) + arr_dim(size(arr_dim))
                call h5dset_extent_f(dset_id, arr_tot_dim, ierr)
                call h5screate_simple_f(arr_rank, arr_dim, memspace_id, ierr, arr_max_dim)
                ! Refresh dataspace following extension before accessing more data.
                call h5dget_space_f(dset_id, dspace_id, ierr)

                ! Write to extended part of dataset.
                offset = 0
                offset(size(arr_dim)) = arr_dim_stored(size(arr_dim))
                call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, arr_dim, ierr)
                call h5dwrite_f(dset_id, dtype, arr_ptr, ierr, memspace_id, dspace_id)

            else
                ! Create.
                call h5screate_simple_f(arr_rank, arr_dim, dspace_id, ierr, arr_max_dim)
                chunk = arr_max_dim(size(arr_dim)) == H5S_UNLIMITED_F
                call dataset_array_plist(arr_rank, arr_dim, plist_id, chunk=chunk)
                call h5dcreate_f(id, dset, dtype, dspace_id, dset_id, ierr, dcpl_id=plist_id)

                ! Write.
                call h5dwrite_f(dset_id, dtype, arr_ptr, ierr)

                call h5pclose_f(plist_id, ierr)
            end if

            call h5dclose_f(dset_id, ierr)
            call h5sclose_f(dspace_id, ierr)

        end subroutine write_ptr

        ! === Helper procedures: reading ===

        subroutine read_string(id, dset, length, string)

            ! Read a string from an open HDF5 file/group.

            ! In:
            !    id: file or group HD5 identifier.
            !    dset: dataset name.
            !    length: length of the string to read
            ! Out:
            !    string: string read from HDF5 file.

            use hdf5

            integer(hid_t), intent(in) :: id
            character(*), intent(in) :: dset
            integer, intent(in) :: length
            character(*), intent(out) :: string

            integer(hid_t) :: type_id, dset_id
            integer :: ierr

            ! Set up fortran string type of *this* length...
            call h5tcopy_f(H5T_FORTRAN_S1, type_id, ierr)
            call h5tset_size_f(type_id, int(length, SIZE_T), ierr)

            call h5dopen_f(id, dset, dset_id, ierr)
            call h5dread_f(dset_id, type_id, string, [0_HSIZE_T], ierr)
            call h5dclose_f(dset_id, ierr)

        end subroutine read_string

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

        subroutine read_boolean(id, dset, val)

            ! Read a boolean from an open HDF5 file/group.

            ! In:
            !    id: file or group HD5 identifier.
            !    dset: dataset name.
            ! Out:
            !    val: boolean read from HDF5 file.

            ! NOTE: HDF5 can't handle boolean types, so instead we read an integer
            ! which we interpret ourselves in a consistent fashion.

            use errors, only: stop_all

            integer(hid_t), intent(in) :: id
            character(*), intent(in) :: dset
            logical, intent(out) :: val

            integer :: val_int

            call hdf5_read(id, dset, val_int)
            select case(val_int)
            case(hande_hdf5_true)
                val = .true.
            case(hande_hdf5_false)
                val = .false.
            case default
                call stop_all('read_boolean', 'Illegal boolean flag detected.')
            end select

        end subroutine read_boolean

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
            use const, only: int_32

            integer(hid_t), intent(in) :: id
            character(*), intent(in) :: dset
            type(hdf5_kinds_t), intent(in) :: kinds
            integer, intent(in) :: arr_shape(:)
            integer(int_32), intent(out), target :: arr(arr_shape(1))

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
            use const, only: int_64

            integer(hid_t), intent(in) :: id
            character(*), intent(in) :: dset
            type(hdf5_kinds_t), intent(in) :: kinds
            integer, intent(in) :: arr_shape(:)
            integer(int_64), intent(out), target :: arr(arr_shape(1))

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
            use const, only: int_32

            integer(hid_t), intent(in) :: id
            character(*), intent(in) :: dset
            type(hdf5_kinds_t), intent(in) :: kinds
            integer, intent(in) :: arr_shape(:)
            integer(int_32), intent(out), target :: arr(arr_shape(1), arr_shape(2))

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
            use const, only: int_64

            integer(hid_t), intent(in) :: id
            character(*), intent(in) :: dset
            type(hdf5_kinds_t), intent(in) :: kinds
            integer, intent(in) :: arr_shape(:)
            integer(int_64), intent(out), target :: arr(arr_shape(1), arr_shape(2))

            type(c_ptr) :: ptr

            ptr = c_loc(arr)
            call read_ptr(id, dset, kinds%i64, size(arr_shape), int(arr_shape,HSIZE_T), ptr)

        end subroutine read_array_2d_int_64

        subroutine read_array_1d_real_sp(id, dset, kinds, arr_shape, arr)

            ! Read in 1D real(sp) array from an HDF5 file.

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
            use const, only: dp, sp, int_32, int_64
            use checking
            use errors
            use parallel, only: parent

            integer(hid_t), intent(in) :: id
            character(*), intent(in) :: dset
            type(hdf5_kinds_t), intent(in) :: kinds
            integer, intent(in) :: arr_shape(:)
            real(sp), intent(out), target :: arr(arr_shape(1))

            type(c_ptr) :: ptr
            integer(int_32), allocatable :: arr_32(:)
            integer(int_64), allocatable :: arr_64(:)
            real(dp), allocatable :: arr_dp(:)
            integer :: ierr
            character(*), parameter :: msg = 'Attempting conversion from integer to real for dataset: '

            if (dtype_equal(id, dset, kinds%sp)) then
                ptr = c_loc(arr)
                call read_ptr(id, dset, kinds%sp, size(arr_shape), int(arr_shape, HSIZE_T), ptr)
            else if (dtype_equal(id, dset, kinds%dp)) then
                if (parent) call warning('read_array_1d_real_sp', &
                                         'Converting from double to single precision for dataset: '//dset//'.', 2)
                allocate(arr_dp(arr_shape(1)), stat=ierr)
                call check_allocate('arr_dp', arr_shape(1), ierr)
                call hdf5_read(id, dset, kinds, arr_shape, arr_dp)
                arr = real(arr_dp,sp)
                deallocate(arr_dp, stat=ierr)
                call check_deallocate('arr_dp', ierr)
            else
                if (parent) call warning('read_array_1d_real_sp', msg//dset//'.', 2)
                if (dtype_equal(id, dset, kinds%i32)) then
                    allocate(arr_32(arr_shape(1)), stat=ierr)
                    call check_allocate('arr_32', size(arr_32), ierr)
                    call hdf5_read(id, dset, kinds, arr_shape, arr_32)
                    arr = arr_32
                    deallocate(arr_32, stat=ierr)
                    call check_deallocate('arr_32', ierr)
                else if (dtype_equal(id, dset, kinds%i64)) then
                    allocate(arr_64(arr_shape(1)), stat=ierr)
                    call check_allocate('arr_64', size(arr_64), ierr)
                    call hdf5_read(id, dset, kinds, arr_shape, arr_64)
                    arr = real(arr_64,sp)
                    deallocate(arr_64, stat=ierr)
                    call check_deallocate('arr_64', ierr)
                else
                    call stop_all('read_array_1d_real_sp', 'Reading mismatched data type and cannot convert.')
                end if
            end if

        end subroutine read_array_1d_real_sp

        subroutine read_array_2d_real_sp(id, dset, kinds, arr_shape, arr)

            ! Read in 2D real(sp) array from an HDF5 file.

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
            use const, only: dp, sp, int_32, int_64
            use checking
            use errors
            use parallel, only: parent

            integer(hid_t), intent(in) :: id
            character(*), intent(in) :: dset
            type(hdf5_kinds_t), intent(in) :: kinds
            integer, intent(in) :: arr_shape(:)
            real(sp), intent(out), target :: arr(arr_shape(1), arr_shape(2))

            type(c_ptr) :: ptr
            integer(int_32), allocatable :: arr_32(:,:)
            integer(int_64), allocatable :: arr_64(:,:)
            real(dp), allocatable :: arr_dp(:,:)
            integer :: ierr
            character(*), parameter :: msg = 'Attempting conversion from integer to real for dataset: '

            if (dtype_equal(id, dset, kinds%sp)) then
                ptr = c_loc(arr)
                call read_ptr(id, dset, kinds%sp, size(arr_shape), int(arr_shape, HSIZE_T), ptr)
            else if (dtype_equal(id, dset, kinds%dp)) then
                if (parent) call warning('read_array_2d_real_sp', &
                                         'Converting from double to single precision for dataset: '//dset//'.', 2)
                allocate(arr_dp(arr_shape(1),arr_shape(2)), stat=ierr)
                call check_allocate('arr_dp', size(arr_dp), ierr)
                call hdf5_read(id, dset, kinds, arr_shape, arr_dp)
                arr = real(arr_dp,sp)
                deallocate(arr_dp, stat=ierr)
                call check_deallocate('arr_dp', ierr)
            else
                if (parent) call warning('read_array_2d_real_sp', msg//dset//'.', 2)
                if (dtype_equal(id, dset, kinds%i32)) then
                    allocate(arr_32(arr_shape(1),arr_shape(2)), stat=ierr)
                    call check_allocate('arr_32', size(arr_32), ierr)
                    call hdf5_read(id, dset, kinds, arr_shape, arr_32)
                    arr = arr_32
                    deallocate(arr_32, stat=ierr)
                    call check_deallocate('arr_32', ierr)
                else if (dtype_equal(id, dset, kinds%i64)) then
                    allocate(arr_64(arr_shape(1),arr_shape(2)), stat=ierr)
                    call check_allocate('arr_64', size(arr_64), ierr)
                    call hdf5_read(id, dset, kinds, arr_shape, arr_64)
                    arr = real(arr_64,sp)
                    deallocate(arr_64, stat=ierr)
                    call check_deallocate('arr_64', ierr)
                else
                    call stop_all('read_array_1d_real_sp', 'Reading mismatched data type and cannot convert.')
                end if
            end if

        end subroutine read_array_2d_real_sp

        subroutine read_array_1d_real_dp(id, dset, kinds, arr_shape, arr)

            ! Read in 1D real(dp) array from an HDF5 file.

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
            use const, only: sp, dp, int_32, int_64
            use checking
            use errors
            use parallel, only: parent

            integer(hid_t), intent(in) :: id
            character(*), intent(in) :: dset
            type(hdf5_kinds_t), intent(in) :: kinds
            integer, intent(in) :: arr_shape(:)
            real(dp), intent(out), target :: arr(arr_shape(1))

            type(c_ptr) :: ptr
            integer(int_32), allocatable :: arr_32(:)
            integer(int_64), allocatable :: arr_64(:)
            real(sp), allocatable :: arr_sp(:)
            integer :: ierr
            character(*), parameter :: msg = 'Attempting conversion from integer to real for dataset: '

            if (dtype_equal(id, dset, kinds%dp)) then
                ptr = c_loc(arr)
                call read_ptr(id, dset, kinds%dp, size(arr_shape), int(arr_shape, HSIZE_T), ptr)
            else if (dtype_equal(id, dset, kinds%sp)) then
                if (parent) call warning('read_array_1d_real_dp', &
                                         'Converting from single to double precision for dataset: '//dset//'.', 2)
                allocate(arr_sp(arr_shape(1)), stat=ierr)
                call check_allocate('arr_sp', arr_shape(1), ierr)
                call hdf5_read(id, dset, kinds, arr_shape, arr_sp)
                arr = arr_sp
                deallocate(arr_sp, stat=ierr)
                call check_deallocate('arr_sp', ierr)
            else
                if (parent) call warning('read_array_1d_real_dp', msg//dset//'.', 2)
                if (dtype_equal(id, dset, kinds%i32)) then
                    allocate(arr_32(arr_shape(1)), stat=ierr)
                    call check_allocate('arr_32', size(arr_32), ierr)
                    call hdf5_read(id, dset, kinds, arr_shape, arr_32)
                    arr = arr_32
                    deallocate(arr_32, stat=ierr)
                    call check_deallocate('arr_32', ierr)
                else if (dtype_equal(id, dset, kinds%i64)) then
                    allocate(arr_64(arr_shape(1)), stat=ierr)
                    call check_allocate('arr_64', size(arr_64), ierr)
                    call hdf5_read(id, dset, kinds, arr_shape, arr_64)
                    arr = arr_64
                    deallocate(arr_64, stat=ierr)
                    call check_deallocate('arr_64', ierr)
                else
                    call stop_all('read_array_1d_real_dp', 'Reading mismatched data type and cannot convert.')
                end if
            end if

        end subroutine read_array_1d_real_dp

        subroutine read_array_2d_real_dp(id, dset, kinds, arr_shape, arr)

            ! Read in 2D real(dp) array from an HDF5 file.

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
            use const, only: sp, dp, int_32, int_64
            use checking
            use errors
            use parallel, only: parent

            integer(hid_t), intent(in) :: id
            character(*), intent(in) :: dset
            type(hdf5_kinds_t), intent(in) :: kinds
            integer, intent(in) :: arr_shape(:)
            real(dp), intent(out), target :: arr(arr_shape(1), arr_shape(2))

            type(c_ptr) :: ptr
            integer(int_32), allocatable :: arr_32(:,:)
            integer(int_64), allocatable :: arr_64(:,:)
            real(sp), allocatable :: arr_sp(:,:)
            integer :: ierr
            character(*), parameter :: msg = 'Attempting conversion from integer to real for dataset: '

            if (dtype_equal(id, dset, kinds%dp)) then
                ptr = c_loc(arr)
                call read_ptr(id, dset, kinds%dp, size(arr_shape), int(arr_shape, HSIZE_T), ptr)
            else if (dtype_equal(id, dset, kinds%sp)) then
                if (parent) call warning('read_array_2d_real_dp', &
                                         'Converting from single to double precision for dataset: '//dset//'.', 2)
                allocate(arr_sp(arr_shape(1),arr_shape(2)), stat=ierr)
                call check_allocate('arr_sp', size(arr_sp), ierr)
                call hdf5_read(id, dset, kinds, arr_shape, arr_sp)
                arr = arr_sp
                deallocate(arr_sp, stat=ierr)
                call check_deallocate('arr_sp', ierr)
            else
                if (parent) call warning('read_array_2d_real_dp', msg//dset//'.', 2)
                if (dtype_equal(id, dset, kinds%i32)) then
                    allocate(arr_32(arr_shape(1),arr_shape(2)), stat=ierr)
                    call check_allocate('arr_32', size(arr_32), ierr)
                    call hdf5_read(id, dset, kinds, arr_shape, arr_32)
                    arr = arr_32
                    deallocate(arr_32, stat=ierr)
                    call check_deallocate('arr_32', ierr)
                else if (dtype_equal(id, dset, kinds%i64)) then
                    allocate(arr_64(arr_shape(1),arr_shape(2)), stat=ierr)
                    call check_allocate('arr_64', size(arr_64), ierr)
                    call hdf5_read(id, dset, kinds, arr_shape, arr_64)
                    arr = arr_64
                    deallocate(arr_64, stat=ierr)
                    call check_deallocate('arr_64', ierr)
                else
                    call stop_all('read_array_1d_real_dp', 'Reading mismatched data type and cannot convert.')
                end if
            end if

        end subroutine read_array_2d_real_dp

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
            integer(hid_t) :: dspace_id, dset_id

            if (.not.dtype_equal(id, dset, dtype)) &
                call stop_all('read_ptr', 'Reading mismatched data type.')

            call h5dopen_f(id, dset, dset_id, ierr)

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

#endif

end module hdf5_helper

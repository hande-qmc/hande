module restart_utils

! Convert restart files for calculations with different parameters.

#ifndef DISABLE_HDF5

implicit none

private
public :: convert_dets, convert_ref

interface convert_dets
    module procedure convert_i32_to_i64
end interface convert_dets

interface convert_ref
    module procedure convert_ref_32_to_64
end interface convert_ref

contains

    subroutine convert_i32_to_i64(id, dset, kinds, dets)

        ! Convert determinants from 32 to 64 bit integers

        use const, only: int_64, int_32

        use hdf5
        use hdf5_helper, only: hdf5_kinds_t, hdf5_read, dset_shape

        integer(hid_t), intent(in) :: id
        character(*), intent(in) :: dset
        type(hdf5_kinds_t), intent(in) :: kinds
        integer(int_64), intent(out) :: dets(:,:)

        integer(int_32), allocatable :: dets_tmp(:,:)
        integer(hsize_t) :: dims(2)

        call dset_shape(id, dset, dims)
        allocate(dets_tmp(dims(1), dims(2)))
        call hdf5_read(id, dset, kinds, shape(dets_tmp), dets_tmp)

        dets = 0
        dets(:dims(1)/2,:dims(2)) = ishft(int(dets_tmp(2::2,:),int_64),32)
        dets(:,:dims(2)) = dets(:,:dims(2)) + dets_tmp(1::2,:)

        deallocate(dets_tmp)

    end subroutine convert_i32_to_i64

    subroutine convert_ref_32_to_64(id, dset, kinds, f0)
        
        use const, only: int_64, int_32

        use hdf5
        use hdf5_helper, only: hdf5_kinds_t, hdf5_read, dset_shape

        integer(hid_t), intent(in) :: id
        character(*), intent(in) :: dset
        type(hdf5_kinds_t), intent(in) :: kinds
        integer(int_64), intent(out) :: f0(:)

        integer(int_32), allocatable :: f0_tmp(:)
        integer(hsize_t) :: dims(1)

        call dset_shape(id, dset, dims)
        allocate(f0_tmp(dims(1)))
        call hdf5_read(id, dset, kinds, shape(f0_tmp), f0_tmp)

        f0 = 0
        f0(:dims(1)/2) = ishft(int(f0_tmp(2::2),int_64),32)
        f0 = f0 + f0_tmp(1::2)

        deallocate(f0_tmp)

    end subroutine convert_ref_32_to_64

#endif

end module restart_utils

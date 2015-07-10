module restart_utils

! Convert restart files for calculations with different parameters.

#ifndef DISABLE_HDF5

implicit none

private
public :: convert_dets

interface convert_dets
    module procedure convert_i32_to_i64
end interface convert_dets

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

        integer :: i

        call dset_shape(id, dset, dims)
        allocate(dets_tmp(dims(1), dims(2)))
        call hdf5_read(id, dset, kinds, shape(dets_tmp), dets_tmp)

        dets = 0
        dets(:dims(1)/2,:dims(2)) = ishft(int(dets_tmp(2::2,:),int_64),32)
        dets(:,:dims(2)) = dets(:,:dims(2)) + dets_tmp(1::2,:)

        deallocate(dets_tmp)

    end subroutine convert_i32_to_i64

#endif

end module restart_utils

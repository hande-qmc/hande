module restart_utils

! Convert restart files for calculations with different compile-time parameters.
! Note that we do not attempt to be memory efficient -- converting from one kind
! to another uses arrays of equal size of both kinds.  This could be reduced by
! using hyperslabs to read in part of the array at a time, so the array of the
! original kind (as in the HDF5 file) could be made smaller than the array on file.

#ifndef DISABLE_HDF5

implicit none

private
public :: convert_dets, convert_ref, convert_pops, change_pop_scaling, change_nbasis

interface convert_dets
    module procedure convert_dets_32_to_64
    module procedure convert_dets_64_to_32
end interface convert_dets

interface convert_ref
    module procedure convert_ref_32_to_64
    module procedure convert_ref_64_to_32
end interface convert_ref

interface convert_pops
    module procedure convert_pops_32_to_64
    module procedure convert_pops_64_to_32
end interface convert_pops

interface change_pop_scaling
    module procedure change_pop_scaling_32
    module procedure change_pop_scaling_64
end interface change_pop_scaling

contains

    subroutine convert_dets_32_to_64(id, dset, kinds, dets)

        ! Convert determinants from 32 to 64 bit integers from a restart file.

        ! In:
        !   id: file or group HD5 identifier,
        !   dset: dataset name.
        !   kinds: hdf5_kinds_t object containing the mapping between the non-default
        !       kinds used in HANDE and HDF5 datatypes.
        ! Out:
        !   dets: determinant list read from restart file converted to 64 bit.

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

        do i = 1, dims(2)
            ! Must convert each determinant separately as there is a different
            ! amount of padding if the number of 32 bit integers is odd
            dets(:,i) = transfer(dets_tmp(:,i), dets)
        end do

        deallocate(dets_tmp)

    end subroutine convert_dets_32_to_64

    subroutine change_nbasis(id, dset, kinds, dets)
        
        ! Read determinants in from restart file and change to a larger basis

        ! In:
        !   id: file or group HD5 identifier,
        !   dset: dataset name.
        !   kinds: hdf5_kinds_t object containing the mapping between the non-default
        !       kinds used in HANDE and HDF5 datatypes.
        ! Out:
        !   dets: determinant list read from restart file in larger basis.

        use const

        use hdf5
        use hdf5_helper

        integer(hid_t), intent(in) :: id
        character(*), intent(in) :: dset
        type(hdf5_kinds_t), intent(in) :: kinds
        integer(i0), intent(out) :: dets(:,:)

        integer(i0), allocatable :: dets_tmp(:,:)
        integer(hsize_t) :: dims(2)

        call dset_shape(id, dset, dims)
        allocate(dets_tmp(dims(1),dims(2)))
        call hdf5_read(id, dset, kinds, shape(dets_tmp), dets_tmp)

        ! Assume the old (small) basis corresponds to the first orbitals in the new basis
        dets = 0
        dets(:dims(1),:dims(2)) = dets_tmp

        deallocate(dets_tmp)

    end subroutine change_nbasis

    subroutine convert_ref_32_to_64(id, dset, kinds, f0)
 
        ! Convert a reference determinant from a restart file from 32 to 64 bit integers.

        ! In:
        !   id: file or group HD5 identifier,
        !   dset: dataset name.
        !   kinds: hdf5_kinds_t object containing the mapping between the non-default
        !       kinds used in HANDE and HDF5 datatypes.
        ! Out:
        !   f0: determinant read from restart file converted to 64 bit.
       
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
        f0 = transfer(f0_tmp, f0)

        deallocate(f0_tmp)

    end subroutine convert_ref_32_to_64

    subroutine convert_pops_32_to_64(id, dset, kinds, pops, scale_factor)

        ! Convert populations from 32 to 64 bit integers.  Needs to be a
        ! different routine from determinants because transfer will not
        ! handle negative numbers correctly.

        ! In:
        !   id: file or group HD5 identifier,
        !   dset: dataset name.
        !   kinds: hdf5_kinds_t object containing the mapping between the non-default
        !       kinds used in HANDE and HDF5 datatypes.
        ! In/Out:
        !   scale_factor: unaltered and unused.  For interface compatibility only.
        ! Out:
        !   pops: population list read from restart file converted to 64 bit.


        use const, only: int_32, int_64

        use hdf5
        use hdf5_helper, only: hdf5_kinds_t, hdf5_read, dset_shape

        integer(hid_t), intent(in) :: id
        character(*), intent(in) :: dset
        type(hdf5_kinds_t), intent(in) :: kinds
        integer(int_64), intent(out) :: pops(:,:)
        integer(int_64), intent(inout) :: scale_factor

        integer(int_32), allocatable :: pops_tmp(:,:)
        integer(hsize_t) :: dims(2)

        call dset_shape(id, dset, dims)
        allocate(pops_tmp(dims(1),dims(2)))
        call hdf5_read(id, dset, kinds, shape(pops_tmp), pops_tmp)

        pops(:,:dims(2)) = pops_tmp

        deallocate(pops_tmp)

    end subroutine convert_pops_32_to_64

    subroutine convert_dets_64_to_32(id, dset, kinds, dets)

        ! Convert determinants from 64 to 32 bit integers

        ! In:
        !   id: file or group HD5 identifier,
        !   dset: dataset name.
        !   kinds: hdf5_kinds_t object containing the mapping between the non-default
        !       kinds used in HANDE and HDF5 datatypes.
        ! Out:
        !   dets: determinant list read from restart file converted to 32 bit.

        use const, only: int_32, int_64

        use hdf5
        use hdf5_helper, only: hdf5_kinds_t, hdf5_read, dset_shape

        integer(hid_t), intent(in) :: id
        character(*), intent(in) :: dset
        type(hdf5_kinds_t), intent(in) :: kinds
        integer(int_32), intent(out) :: dets(:,:)

        integer(int_64), allocatable :: dets_tmp(:,:)
        integer(hsize_t) :: dims(2)

        integer :: i

        call dset_shape(id, dset, dims)
        allocate(dets_tmp(dims(1), dims(2)))
        call hdf5_read(id, dset, kinds, shape(dets_tmp), dets_tmp)

        dets = 0

        do i = 1, dims(2)
            ! Must convert each determinant separately as there is a different
            ! amount of padding if the number of 32 bit integers is odd
            dets(:,i) = transfer(dets_tmp(:,i), dets)
        end do

        deallocate(dets_tmp)

    end subroutine convert_dets_64_to_32

    subroutine convert_ref_64_to_32(id, dset, kinds, f0)

        ! Convert a reference determinant from a restart file from 64 to 32 bit integers.
 
        ! In:
        !   id: file or group HD5 identifier,
        !   dset: dataset name.
        !   kinds: hdf5_kinds_t object containing the mapping between the non-default
        !       kinds used in HANDE and HDF5 datatypes.
        ! Out:
        !   f0: determinant read from restart file converted to 32 bit.
       
        use const, only: int_32, int_64

        use hdf5
        use hdf5_helper, only: hdf5_kinds_t, hdf5_read, dset_shape

        integer(hid_t), intent(in) :: id
        character(*), intent(in) :: dset
        type(hdf5_kinds_t), intent(in) :: kinds
        integer(int_32), intent(out) :: f0(:)

        integer(int_64), allocatable :: f0_tmp(:)
        integer(hsize_t) :: dims(1)

        call dset_shape(id, dset, dims)
        allocate(f0_tmp(dims(1)))
        call hdf5_read(id, dset, kinds, shape(f0_tmp), f0_tmp)

        f0 = 0
        f0 = transfer(f0_tmp, f0)

        deallocate(f0_tmp)

    end subroutine convert_ref_64_to_32

    subroutine convert_pops_64_to_32(id, dset, kinds, pops, scale_factor)

        ! Convert populations from 64 to 32 bit integers.  Needs to be a
        ! different routine from determinants because transfer will not
        ! handle negative numbers correctly.

        ! In:
        !   id: file or group HD5 identifier,
        !   dset: dataset name.
        !   kinds: hdf5_kinds_t object containing the mapping between the non-default
        !       kinds used in HANDE and HDF5 datatypes.
        ! In/Out:
        !   scale_factor: On input population scaling factor to store populations.  On output, unaltered
        !       if integer weights were used or changed to the scaling for 32-bit fixed precision if real
        !       weights were used.
        ! Out:
        !   pops: population list read from restart file converted to 32 bit.

        use const, only: int_64, int_32

        use hdf5
        use hdf5_helper, only: hdf5_kinds_t, hdf5_read, dset_shape

        integer(hid_t), intent(in) :: id
        character(*), intent(in) :: dset
        type(hdf5_kinds_t), intent(in) :: kinds
        integer(int_32), intent(out) :: pops(:,:)
        integer(int_64), intent(inout) :: scale_factor

        integer(int_64), allocatable :: pops_tmp(:,:)
        integer(hsize_t) :: dims(2)

        call dset_shape(id, dset, dims)
        allocate(pops_tmp(dims(1),dims(2)))
        call hdf5_read(id, dset, kinds, shape(pops_tmp), pops_tmp)

        if (scale_factor == 2_int_64**31) then
            ! Need to change real population scaling factor
            call change_pop_scaling(pops_tmp, scale_factor, 2_int_64**11)
            scale_factor = 2_int_64**11
        end if
        pops(:,:dims(2)) = pops_tmp

        deallocate(pops_tmp)

    end subroutine convert_pops_64_to_32

    subroutine change_pop_scaling_32(pops, old_scaling, new_scaling)

        ! Change scaling factor for populations (for integer -> real or 32 -> 64 bit POP_SIZE)

        ! In/Out:
        !   pops: population list
        ! In:
        !   old_scaling: value of real_factor for pops on entry
        !   new_scaling: desired value of real_factor

        use const, only: int_64, int_32

        integer(int_32), intent(inout) :: pops(:,:)
        integer(int_64), intent(in) :: old_scaling, new_scaling

        integer(int_64) :: scaling_change

        if (old_scaling < new_scaling) then
            scaling_change = new_scaling/old_scaling
            pops = pops * scaling_change
        else if (old_scaling > new_scaling) then
            scaling_change = old_scaling/new_scaling
            pops = pops/scaling_change
        end if
    end subroutine change_pop_scaling_32

    subroutine change_pop_scaling_64(pops, old_scaling, new_scaling)

        ! Change scaling factor for populations (for integer -> real or 32 -> 64 bit POP_SIZE)

        ! In/Out:
        !   pops: population list
        ! In:
        !   old_scaling: value of real_factor for pops on entry
        !   new_scaling: desired value of real_factor

        use const, only: int_64

        integer(int_64), intent(inout) :: pops(:,:)
        integer(int_64), intent(in) :: old_scaling, new_scaling

        integer(int_64) :: scaling_change

        if (old_scaling < new_scaling) then
            scaling_change = new_scaling/old_scaling
            pops = pops * scaling_change
        else if (old_scaling > new_scaling) then
            scaling_change = old_scaling/new_scaling
            pops = pops/scaling_change
        end if
    end subroutine change_pop_scaling_64

#endif

end module restart_utils

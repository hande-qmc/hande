module restart_utils

! Convert restart files for calculations with different compile-time parameters.
! Note that we do not attempt to be memory efficient -- converting from one kind
! to another uses arrays of equal size of both kinds.  This could be reduced by
! using hyperslabs to read in part of the array at a time, so the array of the
! original kind (as in the HDF5 file) could be made smaller than the array on file.

#ifndef DISABLE_HDF5

implicit none

private
public :: convert_dets, convert_ref, convert_pops, change_pop_scaling, change_nbasis, change_ninfo

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

    subroutine convert_dets_32_to_64(id, dset, kinds, info_string_len, info_string_len_restart, dets)

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
        integer, intent(in) :: info_string_len, info_string_len_restart
        integer(int_64), intent(out) :: dets(:,:)

        integer(int_32), allocatable :: dets_tmp(:,:)
        integer(hsize_t) :: dims(2)

        integer(hsize_t) :: i

        call dset_shape(id, dset, dims)
        allocate(dets_tmp(dims(1), dims(2)))
        call hdf5_read(id, dset, kinds, shape(dets_tmp, kind=int_64), dets_tmp)

        dets = 0

        do i = 1, dims(2)
            ! Must convert each determinant separately as there is a different
            ! amount of padding if the number of 32 bit integers is odd
            call bit_string_32_to_64(dets_tmp(:dims(1)-info_string_len_restart,i), dets(:,i))
        end do
        ! First check there is information stored within information string
        ! and we want to store it again.
        if (info_string_len_restart > 0 .and. info_string_len > 0) then
            do i = 1, dims(2)
                ! NB previous checks will ensure info_string_len >= info_strin_len_restart
                dets(size(dets,dim=1)-info_string_len+1:size(dets,dim=1)-info_string_len+info_string_len_restart,i) = &
                        int(dets_tmp(dims(1)-info_string_len_restart+1:,i),kind=int_64)
            end do
        end if

        deallocate(dets_tmp)

    end subroutine convert_dets_32_to_64

    subroutine change_nbasis(id, dset, kinds, info_string_len, info_string_len_restart, dets)

        ! Read determinants in from restart file and change to a larger basis

        ! In:
        !   id: file or group HD5 identifier,
        !   dset: dataset name.
        !   kinds: hdf5_kinds_t object containing the mapping between the non-default
        !       kinds used in HANDE and HDF5 datatypes.
        !   info_string_len: number of integers we want to use to store additional
        !       information in bit strings once read in.
        !   info_string_len_restart: number of integers used to store additional
        !       information in bit strings in restart file.
        ! Out:
        !   dets: determinant list read from restart file in larger basis.

        use const

        use hdf5
        use hdf5_helper

        integer(hid_t), intent(in) :: id
        character(*), intent(in) :: dset
        type(hdf5_kinds_t), intent(in) :: kinds
        integer, intent(in) :: info_string_len, info_string_len_restart
        integer(i0), intent(out) :: dets(:,:)

        integer(i0), allocatable :: dets_tmp(:,:)
        integer(hsize_t) :: dims(2)

        call dset_shape(id, dset, dims)
        allocate(dets_tmp(dims(1),dims(2)))

        call hdf5_read(id, dset, kinds, shape(dets_tmp, kind=int_64), dets_tmp)

        ! Assume the old (small) basis corresponds to the first orbitals in the new basis
        dets = 0
        ! If we have a change in bitlength used for additional information we need to
        ! be careful to ensure any additional information is transferred correctly.
        ! NB previous checks will ensure info_string_len >= info_strin_len_restart unless
        ! info_string_len == 0.
        if (info_string_len > 0 .and. info_string_len_restart > 0) then
            dets(:dims(1)-info_string_len_restart,:dims(2)) = dets_tmp(:dims(1)-info_string_len_restart,:)
            dets(size(dets,dim=1)-info_string_len+1:,:dims(2)) = dets_tmp(dims(1)-info_string_len_restart+1:,:)
        else
            ! We know that dims(1)-info_string_len_restart <= size(dets,dim=1)-info_string_len
            dets(:dims(1)-info_string_len_restart,:dims(2)) = dets_tmp(:dims(1)-info_string_len_restart,:dims(2))
        end if
        deallocate(dets_tmp)

    end subroutine change_nbasis

    subroutine change_ninfo(id, dset, kinds, info_string_len, info_string_len_restart, ref)

        ! Read in reference and change number of bit string entries used for
        ! additional information. For determinants this is combned with
        ! change_nbasis, but for determinants this cannot be achieved.

        ! In:
        !   id: file or group HD5 identifier,
        !   dset: dataset name.
        !   kinds: hdf5_kinds_t object containing the mapping between the non-default
        !       kinds used in HANDE and HDF5 datatypes.
        !   info_string_len: number of integers we want to use to store additional
        !       information in bit strings once read in.
        !   info_string_len_restart: number of integers used to store additional
        !       information in bit strings in restart file.
        ! Out:
        !   ref: reference determinant bit string read from restart file with correct
        !       length of information string.

        use const

        use hdf5
        use hdf5_helper

        integer(hid_t), intent(in) :: id
        character(*), intent(in) :: dset
        type(hdf5_kinds_t), intent(in) :: kinds
        integer, intent(in) :: info_string_len, info_string_len_restart
        integer(i0), intent(out) :: ref(:)

        integer(i0), allocatable :: ref_tmp(:)
        integer(hsize_t) :: dims(1)

        call dset_shape(id, dset, dims)
        allocate(ref_tmp(dims(1)))

        call hdf5_read(id, dset, kinds, shape(ref_tmp, kind=int_64), ref_tmp)

        ref = 0
        ! NB previously know that info_string_len /= info_string_len_restart
        if (info_string_len_restart == 0 .or. info_string_len == 0) then
            ref(:size(ref,dim=1)-info_string_len) = ref_tmp(:dims(1)-info_string_len_restart)
        else if (info_string_len_restart < info_string_len) then
            ref(:size(ref,dim=1)-info_string_len) = ref_tmp(:dims(1)-info_string_len_restart)
            ref(size(ref,dim=1)-info_string_len+1:size(ref,dim=1)-info_string_len+info_string_len_restart) = &
                            ref_tmp(dims(1)-info_string_len_restart+1:)
        end if
        deallocate(ref_tmp)

    end subroutine change_ninfo

    subroutine convert_ref_32_to_64(id, dset, kinds, info_string_len, info_string_len_restart, f0)
 
        ! Convert a reference determinant from a restart file from 32 to 64 bit integers.

        ! In:
        !   id: file or group HD5 identifier,
        !   dset: dataset name.
        !   kinds: hdf5_kinds_t object containing the mapping between the non-default
        !       kinds used in HANDE and HDF5 datatypes.
        !   info_string_len: number of integers we want to use to store additional
        !       information in bit strings once read in.
        !   info_string_len_restart: number of integers used to store additional
        !       information in bit strings in restart file.
        ! Out:
        !   f0: determinant read from restart file converted to 64 bit.
       
        use const, only: int_64, int_32

        use hdf5
        use hdf5_helper, only: hdf5_kinds_t, hdf5_read, dset_shape

        integer(hid_t), intent(in) :: id
        character(*), intent(in) :: dset
        type(hdf5_kinds_t), intent(in) :: kinds
        integer, intent(in) :: info_string_len, info_string_len_restart
        integer(int_64), intent(out) :: f0(:)

        integer(int_32), allocatable :: f0_tmp(:)
        integer(hsize_t) :: dims(1)

        call dset_shape(id, dset, dims)
        allocate(f0_tmp(dims(1)))
        call hdf5_read(id, dset, kinds, shape(f0_tmp, kind=int_64), f0_tmp)

        f0 = 0
        if (info_string_len > 0 .and. info_string_len_restart > 0) then
            call bit_string_32_to_64(f0_tmp(:dims(1)-info_string_len_restart), f0)
            f0(size(f0,dim=1)-info_string_len+1:) = int(f0_tmp(dims(1)-info_string_len_restart+1:),kind=int_64)
        else
            call bit_string_32_to_64(f0_tmp, f0)
        end if
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
        call hdf5_read(id, dset, kinds, shape(pops_tmp, kind=int_64), pops_tmp)

        pops(:,:dims(2)) = pops_tmp

        deallocate(pops_tmp)

    end subroutine convert_pops_32_to_64

    subroutine convert_dets_64_to_32(id, dset, kinds, info_string_len, info_string_len_restart, dets)

        ! Convert determinants from 64 to 32 bit integers

        ! In:
        !   id: file or group HD5 identifier,
        !   dset: dataset name.
        !   kinds: hdf5_kinds_t object containing the mapping between the non-default
        !       kinds used in HANDE and HDF5 datatypes.
        !   info_string_len: number of integers we want to use to store additional
        !       information in bit strings once read in.
        !   info_string_len_restart: number of integers used to store additional
        !       information in bit strings in restart file.
        ! Out:
        !   dets: determinant list read from restart file converted to 32 bit.

        use const, only: int_32, int_64

        use hdf5
        use hdf5_helper, only: hdf5_kinds_t, hdf5_read, dset_shape

        integer(hid_t), intent(in) :: id
        character(*), intent(in) :: dset
        type(hdf5_kinds_t), intent(in) :: kinds
        integer, intent(in) :: info_string_len, info_string_len_restart
        integer(int_32), intent(out) :: dets(:,:)

        integer(int_64), allocatable :: dets_tmp(:,:)
        integer(hsize_t) :: dims(2)

        integer(hsize_t) :: i

        call dset_shape(id, dset, dims)
        allocate(dets_tmp(dims(1), dims(2)))
        call hdf5_read(id, dset, kinds, shape(dets_tmp, kind=int_64), dets_tmp)

        dets = 0

        do i = 1, dims(2)
            ! Must convert each determinant separately as there is a different
            ! amount of padding if the number of 32 bit integers is odd
            call bit_string_64_to_32(dets_tmp(:dims(1)-info_string_len_restart,i), dets(:,i))
        end do

        ! First check there is information stored within information string
        ! and we want to store it again.
        if (info_string_len_restart > 0 .and. info_string_len > 0) then
            do i = 1, dims(2)
                ! NB previous checks will ensure info_string_len >= info_strin_len_restart
                dets(size(dets,dim=1)-info_string_len+1:size(dets,dim=1)-info_string_len+info_string_len_restart,i) = &
                        int(dets_tmp(dims(1)-info_string_len_restart+1:,i),kind=int_32)
            end do
        end if

        deallocate(dets_tmp)

    end subroutine convert_dets_64_to_32

    subroutine convert_ref_64_to_32(id, dset, kinds, info_string_len, info_string_len_restart, f0)

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
        integer, intent(in) :: info_string_len, info_string_len_restart
        integer(int_32), intent(out) :: f0(:)

        integer(int_64), allocatable :: f0_tmp(:)
        integer(hsize_t) :: dims(1)

        call dset_shape(id, dset, dims)
        allocate(f0_tmp(dims(1)))
        call hdf5_read(id, dset, kinds, shape(f0_tmp, kind=int_64), f0_tmp)

        f0 = 0
        if (info_string_len > 0 .and. info_string_len_restart > 0) then
            call bit_string_64_to_32(f0_tmp(:dims(1)-info_string_len_restart), f0)
            f0(size(f0,dim=1)-info_string_len+1:) = int(f0_tmp(dims(1)-info_string_len_restart+1:),kind=int_32)
        else
            call bit_string_64_to_32(f0_tmp, f0)
        end if
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
        call hdf5_read(id, dset, kinds, shape(pops_tmp, kind=int_64), pops_tmp)

        if (scale_factor == 2_int_64**31) then
            ! Need to change real population scaling factor
            call change_pop_scaling(pops_tmp, scale_factor, 2_int_64**11)
            scale_factor = 2_int_64**11
        end if
        ! Assume populations now fit in a 32-bit integer.
        pops(:,:dims(2)) = int(pops_tmp, int_32)

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

        ! Note currently scaling values are 2^31 and 2^11, so scaling_change fits in a 32-bit integer comfortably.
        if (old_scaling < new_scaling) then
            scaling_change = new_scaling/old_scaling
            pops = pops * int(scaling_change,int_32)
        else if (old_scaling > new_scaling) then
            scaling_change = old_scaling/new_scaling
            pops = pops/int(scaling_change,int_32)
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

    subroutine bit_string_32_to_64(b32, b64)

        ! In:
        !   b32: a 32-bit string of size n containing up to n*32 bits
        ! Out:
        !   b64: a 64-bit string of size ceiling(n/2) containing the up to n*32 bits. All higher bits are set to 0.

        ! NOTE: we assume b64 is just big enough to hold all the information in b32.

        use const, only: int_32, int_64

        integer(int_32), intent(in) :: b32(:)
        integer(int_64), intent(out) :: b64(:)

        b64 = transfer(b32, b64, size(b64))
        if (mod(size(b32),2) == 1) then
            ! b64 string is 32 bits longer than the b32 string.
            ! Trailing bits are not set in the final integer and can be
            ! set to arbitary values: https://gcc.gnu.org/onlinedocs/gfortran/TRANSFER.html
            ! Easiest solution: simply mask out the trailing 32 bits.
            associate (b64_high => b64(size(b64)), mask32 => maskr(32,int_64))
                b64_high = iand(b64_high, mask32)
            end associate
        end if

    end subroutine bit_string_32_to_64

    subroutine bit_string_64_to_32(b64, b32)

        ! In:
        !   b64: a 64-bit string of size n containing up to n*64 bits
        ! Out:
        !   b32: a 64-bit string of size m. Set to hold the m*32 bits from b64 on output.

        ! NOTE: we assume b32 is big enough to hold all the information in b64 -- ie m*32 is either n*64 or n*64-32, in which case
        ! the top 32 bits of b64 are *not* set.

        use const, only: int_32, int_64

        integer(int_64), intent(in) :: b64(:)
        integer(int_32), intent(out) :: b32(:)

        ! Spec doesn't explicitly cover this but I believe transfer truncates the result from a 2*n array to a m array.
        b32 = transfer(b64, b32, size(b32))

    end subroutine bit_string_64_to_32

#endif

end module restart_utils

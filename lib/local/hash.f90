module hashing

! Provide interfaces and wrappers to hashing procedures that are written in C++.

! Warning: only working with i0 as a 32-bit or 64-bit integer.

use const

implicit none

interface
    ! Interfaces to hashing algorithms written in C++.
    !    key: array to be hashed.
    !    N: length of array to be hashed.
    !    seed: random(ish!) number to seed the hash (MurmurHash2 only).
    ! Note that MurmurHash2 algorithms destroy N so it's recommended to use the
    ! wrapper function below.
    function MurmurHash2(key, N, seed) result(hash) bind(c)
        use, intrinsic:: iso_c_binding
        use const
        integer(c_i0) :: hash
        type(c_ptr), value :: key
        integer(c_int), intent(in) :: N
        integer(c_int), intent(in) :: seed
    end function MurmurHash2
end interface

contains

    function murmurhash_bit_string(f, N, seed) result(hash)

        ! Wrapper around MurmurHash2.

        ! In:
        !    f: bit string.
        !    N: length of bit string.
        ! Returns:
        !    Hash of f using the MurmurHash2 algorithm.

        use, intrinsic:: iso_c_binding

        integer(i0) :: hash
        integer(i0), intent(in), target :: f(N)
        integer, intent(in) :: N
        integer(c_int), intent(in) :: seed
        type(c_ptr) :: key
        integer(c_int) :: tmp

        ! MurmurHash2 destroys the size paramter, so create a copy.
        ! The size parameter used in Murmurhash is the number of bytes...
        tmp = 4*N

        ! Unfortunately it seems c_loc is not required to be pure by the
        ! F2003 standards! :-(
        key = c_loc(f(1))

        hash = MurmurHash2(key, tmp, seed)

    end function murmurhash_bit_string

end module hashing

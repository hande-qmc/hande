module hashing

! Provide interfaces and wrappers to hashing procedures that are written in C++.

! Warning: only working with i0 as a 32-bit or 64-bit integer.

implicit none

interface
    ! Interfaces to hashing algorithms written in C.
    function MurmurHash2(key, N, seed) result(hash) bind(c, name='MurmurHash2')
        ! In:
        !    key: data to be hashed.
        !    seed: random(ish!) number to seed the hash.
        ! In/Out:
        !    N: number of bytes in data.
        ! Returns:
        !    32-bit hash of data.
        use, intrinsic:: iso_c_binding
        integer(c_int32_t) :: hash
        type(c_ptr), value :: key
        integer(c_int), value :: N
        integer(c_int32_t), value :: seed
    end function MurmurHash2
end interface

contains

    function murmurhash_bit_string(f, N, seed) result(hash)

        ! Wrapper around MurmurHash2.

        ! In:
        !    f: bit string.
        !    N: length of bits in bit string.  NOTE: we hash in multiples of 32 bits.
        ! Returns:
        !    Hash of f using the MurmurHash2 algorithm.

        use, intrinsic:: iso_c_binding
        use const, only: i0

        integer :: hash
        integer, intent(in) :: N
        integer(i0), intent(in), target :: f(N)
        integer(c_int32_t), intent(in) :: seed
        type(c_ptr) :: key
        integer(c_int) :: nbytes

        ! Pass MurmurHash2 a multiple of 32-bits.
        ! Note that DET_SIZE must currently be 32 or 64 bits, so this is safe!
        ! The size parameter used in Murmurhash is in bytes...
        nbytes = ceiling(real(N)/32)*4

        ! Unfortunately it seems c_loc is not required to be pure by the
        ! F2003 standards! :-(
        key = c_loc(f(1))

        hash = MurmurHash2(key, nbytes, seed)

    end function murmurhash_bit_string

end module hashing

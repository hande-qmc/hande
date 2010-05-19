module hashing

! Provide interfaces and wrappers to hashing procedures that are written in C++.

use const

implicit none

interface
    ! Interfaces to hashing algorithms written in C++.
    !    key: array to be hashed.
    !    N: length of array to be hashed.
    !    seed: random(ish!) number to seed the hash (MurmurHash2 only).
    ! Note that MurmurHash2 algorithms destroy N so it's recommended to use the
    ! wrapper function below.
    function MurmurHash2(key, N, seed) result(hash)
        use const
        integer(i0) :: hash
        integer(i0), intent(in) :: key(:)
        integer, intent(inout) :: N
        integer, intent(in) :: seed
    end function MurmurHash2
    function fnv1_hash(key, N) result(hash)
        use const
        integer(i0) :: hash
        integer(i0), intent(in) :: key(:)
        integer, intent(in) :: N
    end function fnv1_hash
    function fnv1a_hash(key, N) result(hash)
        use const
        integer(i0) :: hash
        integer(i0), intent(in) :: key(:)
        integer, intent(in) :: N
    end function fnv1a_hash
    function fnv1_hash32(key, N) result(hash)
        use const
        integer :: hash
        integer, intent(in) :: key(:)
        integer, intent(in) :: N
    end function fnv1_hash32
    function fnv1a_hash32(key, N) result(hash)
        use const
        integer :: hash
        integer, intent(in) :: key(:)
        integer, intent(in) :: N
    end function fnv1a_hash32
end interface

contains

    function murmurhash_bit_string(f, N) result(hash)

        ! Wrapper around MurmurHash2.

        ! In:
        !    f: bit string.
        !    N: length of bit string.
        ! Returns:
        !    Hash of f using the MurmurHash2 algorithm.

        integer(i0) :: hash
        integer(i0), intent(in) :: f(N)
        integer, intent(in) :: N
        integer, parameter :: seed = 7 ! doesn't really matter...
        integer :: tmp

        ! MurmurHash2 destroys the size paramter, so create a copy.
        tmp = N

        hash = MurmurHash2(f, tmp, seed)

    end function murmurhash_bit_string

end module hashing

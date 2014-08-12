module bit_utils

#include "cdefs.h"

! Module for bit utilities.

! Mainly mainipulate variables of type integer(i0) on a bit-wise basis
! as these are used in storing determinants.

use const

implicit none

interface operator(.bitstrgt.)
    module procedure bit_str_32_gt
    module procedure bit_str_64_gt
end interface

interface bit_str_cmp
    module procedure bit_str_32_cmp
    module procedure bit_str_64_cmp
end interface

contains

!--- Counting set bits ---

    elemental function naive_count_set_bits(b) result(nbits)

        ! In:
        !    b: bit string stored as an integer(i0)
        ! Returns:
        !    The number of bits set in  b.
        ! This is an exceptionally naive implementation and should only be used
        ! for testing more optimal versions.

        integer :: nbits
        integer(i0), intent(in) :: b
        integer :: i

        nbits = 0
        do i = 0, i0_end
            if (btest(b, i)) nbits = nbits + 1
        end do

    end function naive_count_set_bits

    elemental function count_set_bits(b) result(nbits)

        integer :: nbits
        integer(i0), intent(in) :: b
        integer(i0) :: tmp

        ! For 8 bit integers:
#if DET_SIZE == 8
        integer(i0), parameter :: m1 = Z'55'
        integer(i0), parameter :: m2 = Z'33'
        integer(i0), parameter :: m3 = Z'0F'

#elif DET_SIZE == 16
        ! For 16 bit integers:
        integer(i0), parameter :: m1 = Z'5555'
        integer(i0), parameter :: m2 = Z'3333'
        integer(i0), parameter :: m3 = Z'0F0F'
        integer(i0), parameter :: m4 = Z'0101'

#elif DET_SIZE == 32
        ! For 32 bit integers:
        integer(i0), parameter :: m1 = Z'55555555'
        integer(i0), parameter :: m2 = Z'33333333'
        integer(i0), parameter :: m3 = Z'0F0F0F0F'
        integer(i0), parameter :: m4 = Z'01010101'

#elif DET_SIZE == 64
        ! For 64 bit integers:
        integer(i0), parameter :: m1 = Z'5555555555555555'
        integer(i0), parameter :: m2 = Z'3333333333333333'
        integer(i0), parameter :: m3 = Z'0f0f0f0f0f0f0f0f'
        integer(i0), parameter :: m4 = Z'0101010101010101'
#endif

        ! This is quite cool.

        tmp = b

        ! For a more detailed explanation and discussion see:
        !   http://everything2.com/title/Counting+1+bits
        !   http://graphics.stanford.edu/~seander/bithacks.html
        !   http://gurmeetsingh.wordpress.com/2008/08/05/fast-bit-counting-routines/
        !   Chapter 5 of the excellent Hacker's Delight by Henry S. Warren.

        ! The general idea is to use a divide and conquer approach.
        ! * Set each 2 bit field to be the sum of the set bits in the two single
        !   bits originally in that field.
        ! * Set each 4 bit field to be the sum of the set bits in the two 2 bit
        !   fields originally in the 4 bit field.
        ! * Set each 8 bit field to be the sum of the set bits in the two 4 bit
        !   fields it contains.
        ! * etc.
        ! Thus we obtain an algorithm like:
        !     x = ( x & 01010101...) + ( (x>>1) & 01010101...)
        !     x = ( x & 00110011...) + ( (x>>2) & 00110011...)
        !     x = ( x & 00001111...) + ( (x>>4) & 00001111...)
        ! etc., where & indicates AND and >> is the shift right operator.
        ! Further optimisations are:
        ! * Any & operations can be omitted where there is no danger that
        ! a field's sum will carry over into the next field.
        ! * The first line can be replaced by:
        !     x = x - ( (x>>1) & 01010101...)
        !   thanks to the population (number of set bits) in an integer
        !   containing p bits being given by:
        !     pop(x) = \sum_{i=0}^{p-1} x/2^i
        ! * Summing 8 bit fields together can be performed via a multiplication
        !   followed by a right shift.
        ! Thus the following (extremely fast) algorithms.

#if DET_SIZE == 8
        ! For 8 bit integers:
        tmp = tmp - iand(ishft(tmp,-1), m1)
        tmp = iand(tmp, m2) + iand(ishft(tmp,-2), m2)
        nbits = iand(tmp + ishft(tmp,-4), m3)

#elif DET_SIZE == 16
        ! For 16 bit integers:
        tmp = tmp - iand(ishft(tmp,-1), m1)
        tmp = iand(tmp, m2) + iand(ishft(tmp,-2), m2)
        tmp = iand(tmp + ishft(tmp,-4), m3)*m4
        nbits = ishft(tmp, -8)

#elif DET_SIZE == 32
        ! For 32 bit integers:
        tmp = tmp - iand(ishft(tmp,-1), m1)
        tmp = iand(tmp, m2) + iand(ishft(tmp,-2), m2)
        tmp = iand((tmp + ishft(tmp,-4)), m3)*m4
        nbits = ishft(tmp, -24)

#elif DET_SIZE == 64
        ! For 64 bit integers:
        tmp = tmp - iand(ishft(tmp,-1), m1)
        tmp = iand(tmp, m2) + iand(ishft(tmp,-2), m2)
        tmp = iand(tmp, m3) + iand(ishft(tmp,-4), m3)
        nbits = ishft(tmp*m4, -56)
#endif

    end function count_set_bits

!--- I/O helpers ---

    elemental function bit_string(b) result(s)

        ! In:
        !    b: bit string stored as an integer(i0)
        ! Returns:
        !    A binary representation of the bit string as a character string.

        character(i0_length) :: s
        integer(i0), intent(in) :: b
        character(10) :: bit_fmt

        ! This is good for integers containing less than 1000 bits.
        ! Producing the format string each time is non-optimal, but this is only
        ! i/o.
        ! The format is something like (B8.8), which gives a bit string of
        ! length 8 and forces all 8 bits to be written (including leading 0s).
        write (bit_fmt,'("(B",I3,".",I3,")")') i0_length, i0_length

        write (s,bit_fmt) b

    end function bit_string

!--- Permutations of set bits in bit string ---

    function first_perm(n) result(p)

        ! In:
        !    n: number of bits to set.
        ! Returns:
        !    i0 bit string containing the lexicographically first permutation of n set bits.

        integer(i0) :: p
        integer, intent(in) :: n
        integer :: i

        p = 0
        do i = 0, n-1
            p = ibset(p,i)
        end do

    end function first_perm

    function bit_permutation(v) result(w)

        ! In:
        !    v: a bit string.
        ! Returns:
        !    The next permutation of the bit string in lexicographic order.
        !
        !    As we store the bit strings as i0 integers, overflow is possible,
        !    i.e. with 10 spin functions and 5 electrons, bit_permuation can
        !    return bits set in the 11th and higher sites.  Fortunately this
        !    only happens after all permutations involving just the first 10
        !    sites are exhausted (by design!), so only happens if bit_permuation
        !    is called too many times...

        integer(i0) :: w
        integer(i0), intent(in) :: v
        integer(i0) :: t1, t2

        ! From http://graphics.stanford.edu/~seander/bithacks.html.

        t1 = ior(v, v-1) + 1
        t2 = ishft(iand(t1,-t1)/iand(v,-v),-1) - 1
        w = ior(t1, t2)

    end function bit_permutation

!--- Converting bit strings ---

    pure subroutine decode_bit_string(b, d)

        ! In:
        !    b: bit string stored as an integer(i0)
        ! Out:
        !    d: list of bits set in b.  It is assumed that d is at least as
        !    large as the number of bits set in b.  If not, then all elements of
        !    d are set to -1.
        !    The bit string is 0-indexed.

        ! NOTE:
        !    This is a simple and rather naive algorithm and should not be used
        !    in performance-critical code.  It can be greatly improved by using
        !    a chunk-wise loop coupled a data table.

        integer(i0), intent(in) :: b
        integer, intent(out) :: d(:)

        integer :: ipos, i

        i = lbound(d, dim=1) - 1
        do ipos = 0, i0_end
            if (btest(b, ipos)) then
                i = i + 1
                if (i > ubound(d, dim=1)) then
                    d = -1
                    exit
                end if
                d(i) = ipos
            end if
        end do

    end subroutine decode_bit_string

!--- Comparison of bit strings---

    pure function bit_str_32_gt(b1, b2) result(gt)

        ! In:
        !    b1(:), b2(:) bit string.
        ! Returns:
        !    True if the first element of b1 which is not equal to the
        !    corresponding element of b2 is greater than the corresponding
        !    element in b2.

        logical :: gt
        integer(int_4), intent(in) :: b1(:), b2(:)

        integer :: i

        gt = .false.
        do i = 1, ubound(b1,dim=1)
            if (b1(i) > b2(i)) then
                gt = .true.
                exit
            else if (b1(i) < b2(i)) then
                gt = .false.
                exit
            end if
        end do

    end function bit_str_32_gt

    pure function bit_str_64_gt(b1, b2) result(gt)

        ! In:
        !    b1(:), b2(:) bit string.
        ! Returns:
        !    True if the first element of b1 which is not equal to the
        !    corresponding element of b2 is greater than the corresponding
        !    element in b2.

        logical :: gt
        integer(int_8), intent(in) :: b1(:), b2(:)

        integer :: i

        gt = .false.
        do i = 1, ubound(b1,dim=1)
            if (b1(i) > b2(i)) then
                gt = .true.
                exit
            else if (b1(i) < b2(i)) then
                gt = .false.
                exit
            end if
        end do

    end function bit_str_64_gt

    pure function bit_str_32_cmp(b1, b2) result(cmp)

        ! In:
        !    b1(:), b2(:): bit string.
        ! Returns:
        !    0 if b1 and b2 are identical;
        !    1 if the first non-identical element in b1 is smaller than the
        !    corresponding element in b2;
        !    -1 if the first non-identical element in b1 is greater than the
        !    corresponding element in b2;

        integer :: cmp
        integer(int_4), intent(in) :: b1(:), b2(:)

        integer :: i

        cmp = 0
        do i = 1, ubound(b1, dim=1)
            if (b1(i) < b2(i)) then
                cmp = 1
                exit
            else if (b1(i) > b2(i)) then
                cmp = -1
                exit
            end if
        end do

    end function bit_str_32_cmp

    pure function bit_str_64_cmp(b1, b2) result(cmp)

        ! In:
        !    b1(:), b2(:): bit string.
        ! Returns:
        !    0 if b1 and b2 are identical;
        !    1 if the first non-identical element in b1 is smaller than the
        !    corresponding element in b2;
        !    -1 if the first non-identical element in b1 is greater than the
        !    corresponding element in b2;

        integer :: cmp
        integer(int_8), intent(in) :: b1(:), b2(:)

        integer :: i

        cmp = 0
        do i = 1, ubound(b1, dim=1)
            if (b1(i) < b2(i)) then
                cmp = 1
                exit
            else if (b1(i) > b2(i)) then
                cmp = -1
                exit
            end if
        end do

    end function bit_str_64_cmp

end module bit_utils

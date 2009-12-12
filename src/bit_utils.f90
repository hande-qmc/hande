module bit_utils

! Module for bit utilities.

! Mainly mainipulate variables of type integer(i0) on a bit-wise basis
! as these are used in storing determinants.

use const

implicit none

contains

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
!        integer(i0), parameter :: m1 = Z'55'
!        integer(i0), parameter :: m2 = Z'33'
!        integer(i0), parameter :: m3 = Z'0F'

        ! For 16 bit integers:
!        integer(i0), parameter :: m1 = Z'5555'
!        integer(i0), parameter :: m2 = Z'3333'
!        integer(i0), parameter :: m3 = Z'0F0F'
!        integer(i0), parameter :: m4 = Z'0101'

        ! For 32 bit integers:
        integer(i0), parameter :: m1 = Z'55555555'
        integer(i0), parameter :: m2 = Z'33333333'
        integer(i0), parameter :: m3 = Z'0F0F0F0F'
        integer(i0), parameter :: m4 = Z'01010101'

        ! For 64 bit integers:
!        integer(i0), parameter :: m1 = Z'5555555555555555'
!        integer(i0), parameter :: m2 = Z'3333333333333333'
!        integer(i0), parameter :: m3 = Z'0f0f0f0f0f0f0f0f'
!        integer(i0), parameter :: m4 = Z'0101010101010101'

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
        !   containing p bits is given by:
        !     pop(x) = \sum_{i=0}^{p-1} x/2^i
        ! * Summing 8 bit fields together can be performed via a multiplication
        !   followed by a right shift.
        ! Thus the following (extremely fast) algorithms.

        ! For 8 bit integers:
!        tmp = tmp - iand(ishft(tmp,-1), m1)
!        tmp = iand(tmp, m2) + iand(ishft(tmp,-2), m2)
!        nbits = iand(tmp + ishft(tmp,-4), m3)

        ! For 16 bit integers:
!        tmp = tmp - iand(ishft(tmp,-1), m1)
!        tmp = iand(tmp, m2) + iand(ishft(tmp,-2), m2)
!        tmp = iand(tmp + ishft(tmp,-4), m3)*m4
!        nbits = ishft(tmp, -8)

        ! For 32 bit integers:
        tmp = tmp - iand(ishft(tmp,-1), m1)
        tmp = iand(tmp, m2) + iand(ishft(tmp,-2), m2)
        tmp = iand((tmp + ishft(tmp,-4)), m3)*m4
        nbits = ishft(tmp, -24)

        ! For 64 bit integers:
!        tmp = tmp - iand(ishft(tmp,-1), m1)
!        tmp = iand(tmp, m2) + iand(ishft(tmp,-2), m2)
!        tmp = iand(tmp, m3) + iand(ishft(tmp,-4), m3)
!        nbits = ishft(tmp*m4, -56)

    end function count_set_bits

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

end module bit_utils

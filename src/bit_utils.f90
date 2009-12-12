module bit_utils

! Module for bit utilities.

! Mainly mainipulate variables of type integer(i0) on a bit-wise basis
! as these are used in storing determinants.

use const

implicit none

contains

    elemental function count_set_bits(b) result(nbits)

        ! In:
        !    b: bit string stored as an integer(i0)
        ! Returns:
        !    The number of bits set in  b.

        integer :: nbits
        integer(i0), intent(in) :: b
        integer :: i

        nbits = 0
        do i = 0, i0_end
            if (btest(b, i)) nbits = nbits + 1
        end do

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

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
        do i = 0,7
            if (btest(b, i)) nbits = nbits + 1
        end do

    end function count_set_bits

    elemental function bit_string(b) result(s)
        
        ! In:
        !    b: bit string stored as an integer(i0)
        ! Returns:
        !    A binary representation of the bit string as a character string.

        character(8) :: s
        integer(i0), intent(in) :: b

        write (s,'(B8.8)') b

    end function bit_string

end module bit_utils

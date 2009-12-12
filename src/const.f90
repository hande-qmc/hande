module const
! Module containing precision of types and constant data.

implicit none

! i0 gives the equivalent of a byte type (8 bits)
! Allows data range of -128 to 127.
! Integers of kind i0 are used as a bit string to store determinants.
! Adjusting i0 varies the size of the integer and may improve performance.
! selected_int_kind(0): equivalent to a byte.  -128 <= int_i0 <= 127.
! selected_int_kind(3): equivalent to a 16 bit integer.  -32768 <= int_i0 <= 32767.
! selected_int_kind(6): equivalent to a 32 bit integer.  -2147483648 <= int_i0 <= 2147483647.
! selected_int_kind(10): equivalent to a 64 bit integer. -9223372036854775808 <= int_i0 <= 9223372036854775807.
! Note that the memory wasted (but not having the number of basis functions
! being a multiple of the number of bits in i0) can increase with the kind.
! However, the performance of intrinsic bit operations with 32 bit integers is
! far superior to that of 8 bit integers (tested on a 64-bit Xeon quad-core).
! If this is changed then the count_set_bits function in the bit_utils module
! must also be changed: code is provided (but commented out where appropriate)
! for 8, 16, 32 and 64 bit integers.
integer, parameter :: i0 = selected_int_kind(9)

! Number of bits in an integer of type i0.
integer, parameter :: i0_length = bit_size(int(0,i0))

! Index of the last bit in an integer of type i0.
! (Bit indexing in fortran ranges from 0 to bit_size-1.)
integer, parameter :: i0_end = bit_size(int(0,i0))-1

! Double precision kind.
! If this is changed then the lapack and scalapack calls must also be changed
! accordingly.
integer, parameter :: dp = selected_real_kind(15,307)

real(dp), parameter :: pi = 3.1415926535897931_dp

! depsilon is the precision used to compare floating point numbers.
real(dp), parameter :: depsilon = 1.e-8

end module const

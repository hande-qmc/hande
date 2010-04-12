module const

! Module containing precision of types and constant data.

#include "cdefs.h"

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
#if DET_SIZE == 8
integer, parameter :: i0 = selected_int_kind(0)
#elif DET_SIZE == 16
integer, parameter :: i0 = selected_int_kind(3)
#elif DET_SIZE == 32
integer, parameter :: i0 = selected_int_kind(6)
#elif DET_SIZE == 64
integer, parameter :: i0 = selected_int_kind(10)
#endif

! Number of bits in an integer of type i0.
integer, parameter :: i0_length = bit_size(int(0,i0))

! Index of the last bit in an integer of type i0.
! (Bit indexing in fortran ranges from 0 to bit_size-1.)
integer, parameter :: i0_end = bit_size(int(0,i0))-1

! Single precision kind.
integer, parameter :: sp = selected_real_kind(6,37)
! Double precision kind.
integer, parameter :: dp = selected_real_kind(15,307)

! Compile time choice of precision level.
! We use p for all real kinds unless double precision is *absolutely* required
! by the algorithm (e.g. the Mersenne Twister RNG).
#ifdef SINGLE_PRECISION
integer, parameter :: p = sp
# else
integer, parameter :: p = dp
#endif

real(p), parameter :: pi = 3.1415926535897931_p

! depsilon is the precision used to compare floating point numbers.
real(p), parameter :: depsilon = 1.e-8

end module const

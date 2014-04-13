module const

! Module containing precision of types and constant data.

#include "cdefs.h"

use, intrinsic :: iso_c_binding, only: c_int32_t, c_int64_t

implicit none

! i0 gives the equivalent of a byte type (8 bits)
! Allows data range of -128 to 127.
! Integers of kind i0 are used as a bit string to store determinants.
! Adjusting i0 varies the size of the integer and may improve performance.
! selected_int_kind(0): equivalent to a byte.  -128 <= int_i0 <= 127.
! selected_int_kind(3): equivalent to a 16 bit integer.  -32768 <= int_i0 <= 32767.
! selected_int_kind(6): equivalent to a 32 bit integer.  -2147483648 <= int_i0 <= 2147483647.
! selected_int_kind(15): equivalent to a 64 bit integer. -9223372036854775808 <= int_i0 <= 9223372036854775807.
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
! C int type which interoperates with i0.
integer, parameter :: c_i0 = c_int32_t
#elif DET_SIZE == 64
integer, parameter :: i0 = selected_int_kind(15)
integer, parameter :: c_i0 = c_int64_t
#endif

! int_p determines whether 32 or 64 integers are used for walker_population.
#if POP_SIZE == 32
integer, parameter :: int_p = selected_int_kind(6)
#elif POP_SIZE == 64
integer, parameter :: int_p = selected_int_kind(15)
#else
! Use 32-bit integers by default.
integer, parameter :: int_p = selected_int_kind(6)
#endif

! The sdata array holds both walker populations and determinants together.
! Therefore, if 64-bit integers are being used for either walker populations
! or determinants, int_s must be 64-bit. Otherwise it can be 32-bit.
#if POP_SIZE == 64 || DET_SIZE == 64
integer, parameter :: int_s = selected_int_kind(15)
#else
integer, parameter :: int_s = selected_int_kind(6)
#endif

! Number of bits in an integer of type i0.
! Note that pgi 10.3 has a bug are returns 32 if bit_size(int(0,i0)) is used.
integer, parameter :: i0_length = bit_size(0_i0)

! Index of the last bit in an integer of type i0.
! (Bit indexing in fortran ranges from 0 to bit_size-1.)
integer, parameter :: i0_end = i0_length - 1

! Single precision kind.
integer, parameter :: sp = selected_real_kind(6,37)
! Double precision kind.
integer, parameter :: dp = selected_real_kind(15,307)

! long integer; used for psip population (where we can exceed 2e9)
integer, parameter :: lint = selected_int_kind(15)

! Compile time choice of precision level.
! We use p for all real kinds unless double precision is *absolutely* required
! by the algorithm (e.g. the Mersenne Twister RNG).
#ifdef SINGLE_PRECISION
integer, parameter :: p = sp
# else
integer, parameter :: p = dp
#endif

real(p), parameter :: pi = 3.1415926535897931_p

! depsilon is the precision used to compare floating point numbers of kind=p.
#ifdef SINGLE_PRECISION
real(p), parameter :: depsilon = 1.e-8_p ! MACHEPS is ~10^-8
#else
real(p), parameter :: depsilon = 1.e-12_p ! MACHEPS is ~10^-16 but we won't go quite that far...
#endif

end module const

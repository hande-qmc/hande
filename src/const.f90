module const
! Module containing precision of types and constant data.

implicit none

! i0 gives the equivalent of a byte type (8 bits)
! Allows data range of -128 to 127.
integer, parameter :: i0 = selected_int_kind(0)

integer, parameter :: dp = kind(0.0d0)

real(dp), parameter :: pi = 3.1415926535897931_dp

! depsilon is the precision used to compare floating point numbers.
real(dp), parameter :: depsilon = 1.e-8

end module const

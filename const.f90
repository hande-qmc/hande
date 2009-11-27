module const
! Module containing precision of types and constant data.

implicit none

integer, parameter :: dp = kind(0.0d0)

real(dp), parameter :: pi = 3.1415926535897931_dp

! depsilon is the precision used to compare floating point numbers.
real(dp), parameter :: depsilon = 1.e-8

end module const

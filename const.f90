module const
! Module containing precision of types and data.

implicit none

integer, parameter :: dp = kind(0.0d0)

real(dp), parameter :: pi = 3.1415926535897931_dp

real(dp), parameter :: depsilon = 1.e-8

end module const

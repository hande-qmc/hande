module hamiltonian_data

! Small module to contain hamiltonian-related derived types to avoid cyclic
! dependencies until a better place for them is found.

use const, only: p

implicit none

! --- WARNING ---
! Only one component of hamil_t is in use at a time, depending upon the context.
! The other component is NOT initialised and so cannot be relied on to be zero.
! --- WARNING ---

type hmatel_t
    ! Derived type to contain all information on a particular hmatel, to
    ! allow compatibility between real and complex interfaces.

    ! Value for real hamiltonian elements.
    real(p) :: r

    ! Value for complex hamiltonian elements.
    complex(p) :: c
end type hmatel_t
   interface operator (*)
      module procedure hmatel_mul_rh
      module procedure hmatel_mul_hr
   end interface

contains

   elemental function hmatel_mul_rh(r1,e2) result(e3)
      real(p), intent(in) ::  r1
      type(hmatel_t), intent(in) ::  e2
      type(hmatel_t) :: e3
      e3%r = e2%r * r1
      e3%c = e2%c * r1
   end function hmatel_mul_rh

   elemental function hmatel_mul_hr(e1,r2) result(e3)
      type(hmatel_t), intent(in) ::  e1
      real(p), intent(in) ::  r2
      type(hmatel_t) :: e3
      e3%r = e1%r * r2
      e3%c = e1%c * r2
   end function hmatel_mul_hr

end module hamiltonian_data

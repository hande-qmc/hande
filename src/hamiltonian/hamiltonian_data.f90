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

end module hamiltonian_data

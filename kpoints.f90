module kpoints

use const
use system

implicit none

type kpoint
    integer, pointer :: k(:)
    real(dp) :: eigv
end type kpoint

contains

    subroutine init_kpoint(kp,k)

        type(kpoint) :: kp
        integer, intent(in), optional  :: k(ndim)
        integer :: ierr

        allocate(kp%k(ndim),stat=ierr)

        if (present(k)) then
            kp%k = k
            kp%eigv = calc_eigv(k)
        end if

    end subroutine init_kpoint

    real(dp) function calc_eigv(k)

        integer, intent(in) :: k(ndim)
        integer :: i

        calc_eigv = 0.0_dp
        do i = 1,ndim
            calc_eigv = calc_eigv + cos((2*pi*k(i))/Length(i))
        end do

        calc_eigv = -2*calc_eigv

    end function calc_eigv

end module kpoints

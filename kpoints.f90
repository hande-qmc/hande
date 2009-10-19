module kpoints

use const
use system

implicit none

type kpoint
    integer, pointer :: k(:)
    real(dp) :: kinetic
end type kpoint

contains

    subroutine init_kpoint(kp,k)

        type(kpoint) :: kp
        integer, intent(in), optional  :: k(ndim)
        integer :: ierr

        allocate(kp%k(ndim),stat=ierr)

        if (present(k)) then
            kp%k = k
            kp%kinetic = calc_kinetic(k)
        end if

    end subroutine init_kpoint

    real(dp) function calc_kinetic(k)

        integer, intent(in) :: k(ndim)
        integer :: i

        calc_kinetic = 0.0_dp
        do i = 1,ndim
            calc_kinetic = calc_kinetic + cos((2*pi*k(i))/Length(i))
        end do

        calc_kinetic = -2*calc_kinetic

    end function calc_kinetic

end module kpoints

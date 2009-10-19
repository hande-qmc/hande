module kpoints

use const
use system

implicit none

type kpoint
    integer :: k(3)
    real(dp) :: eigv
end type kpoint

contains

    type(kpoint) function gen_kpoint(k)

        integer, intent(in)  :: k(ndim)

        gen_kpoint = kpoint(k,calc_eigv(k))

    end function gen_kpoint

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

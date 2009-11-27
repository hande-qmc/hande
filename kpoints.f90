module kpoints

use const
use system

implicit none

type kpoint
    integer, pointer :: k(:) => NULL()
    integer :: ms
    real(dp) :: kinetic
end type kpoint

contains

    pure subroutine init_kpoint(kp,k,ms)

        type(kpoint), intent(out) :: kp
        integer, intent(in), optional  :: k(ndim)
        integer, intent(in), optional  :: ms
        integer :: ierr

        if (.not.associated(kp%k)) then
            allocate(kp%k(ndim),stat=ierr)
        end if

        if (present(k)) then
            kp%k = k
            kp%kinetic = calc_kinetic(k)
        end if

        if (present(ms)) kp%ms = ms

    end subroutine init_kpoint

    pure function calc_kinetic(k) result(kinetic)

        real(dp) :: kinetic

        integer, intent(in) :: k(ndim)
        integer :: i
        real(dp) :: kc(ndim)

        forall (i=1:ndim) kc(i) = sum(k*rlattice(i,:))

        kinetic = -2*sum(cos(2*pi*kc))*hubt

    end function calc_kinetic

    subroutine write_kpoint(k)

        type(kpoint), intent(in) :: k
        integer :: i

        write (6,'(1X,"(")', advance='no')
        write (6,'(i2)',advance='no') k%k(1)
        do i = 2,ndim
            write (6,'(",",i2)',advance='no') k%k(i)
        end do
        write (6,'(")")', advance='no')
        write (6,'(5X,i2)', advance='no') k%ms
        write (6,'(4X,f12.8)') k%kinetic

    end subroutine write_kpoint

end module kpoints

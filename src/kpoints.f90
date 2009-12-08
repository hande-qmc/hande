module kpoints

! Module for handling wavevectors.

use const
use system

implicit none

! The kpoint type is used to specify a spin orbital in momentum space.
type kpoint
    ! Wavevector in terms of the reciprocal lattice vectors of the crystal cell.
    integer, pointer :: k(:) => NULL()
    ! Spin of the electron (1 or -1).
    integer :: ms
    ! Kinetic energy.
    real(dp) :: kinetic
end type kpoint

contains

    pure subroutine init_kpoint(kp,k,ms)

        ! Initialise a variable of type kpoint.
        ! In:
        !   k (optional): wavevector in units of the reciprocal lattice vectors
        !                 of the crystal cell.
        !   ms (optional): set spin of an electron occupying the basis function.
        ! Out:
        !   kp: initialsed kp.  The wavevector and kinetic energy components are
        !       set if the k arguments is given and the ms component is set if
        !       the ms argument is given.  If no optional arguments are
        !       specified then a completely blank variable is returned.
        !
        ! This should be called even if k and ms are not specified so that the
        ! k component can be correctly allocated.

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

        ! In:
        !   k: wavevector in terms of the reciprocal lattice vectors of the
        !      crystal cell.
        ! Returns:
        !   The kinetic energy associated with a given wavevector.

        real(dp) :: kinetic

        integer, intent(in) :: k(ndim)
        integer :: i
        real(dp) :: kc(ndim)

        forall (i=1:ndim) kc(i) = sum(k*rlattice(i,:))

        ! For a square lattice the kinetic energy of a wavevector is given by
        !    -2t \sum_i cos(k.x_i)
        ! where x_i is the i-th reciprocal lattice vector of the primitive unit
        ! cell.
        kinetic = -2*sum(cos(2*pi*kc))*hubt

    end function calc_kinetic

    subroutine write_kpoint(k, iunit, new_line)

        ! Print out information stored in k.
        ! In:
        !    k: kpoint variable.
        !    iunit (optional): io unit to which the output is written.
        !        Default: 6 (stdout).
        !    new_line (optional): if true, then a new line is written at
        !        the end of the list of occupied orbitals.  Default: no
        !        new line.

        type(kpoint), intent(in) :: k
        integer, intent(in), optional :: iunit
        logical, intent(in), optional :: new_line
        integer :: i, io

        if (present(iunit)) then
            io = iunit
        else
            io = 6
        end if

        write (io,'(1X,"(")', advance='no')
        write (io,'(i2)',advance='no') k%k(1)
        do i = 2,ndim
            write (io,'(",",i2)',advance='no') k%k(i)
        end do
        write (io,'(")")', advance='no')
        write (io,'(5X,i2)', advance='no') k%ms
        write (io,'(4X,f12.8)', advance='no') k%kinetic
        if (present(new_line)) then
            if (new_line) write (io,'()')
        end if

    end subroutine write_kpoint

end module kpoints

module kpoints

! Module for handling wavevectors.

use const
use system

implicit none

contains

    pure function calc_kinetic(k) result(kinetic)

        ! In:
        !    k: wavevector in terms of the reciprocal lattice vectors of the
        !       crystal cell.
        ! Returns:
        !    The kinetic energy associated with a given wavevector.

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

    pure function is_reciprocal_lattice_vector(k) result(t_rlv)

        ! In:
        !    k: wavevector in terms of the reciprocal lattice vectors of the
        !       crystal cell.
        ! Returns:
        !    True if k is a reciprocal lattice vector of the *primitive* cell.

        logical :: t_rlv
        integer, intent(in) :: k(ndim)
        real(dp) :: kc(ndim)
        integer :: i

        if (all(k == 0)) then
            ! Easy!
            t_rlv = .true.
        else
            ! Test to see if delta_k is 0 up to a reciprocal lattice vector
            ! of the primitive cell.
            forall (i=1:ndim) kc(i) = sum(k*rlattice(i,:))
            t_rlv = all(abs(kc - nint(kc)) < depsilon)
        end if

    end function is_reciprocal_lattice_vector

end module kpoints

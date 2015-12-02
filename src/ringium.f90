module ringium_system

! Routines for ringium

! Ringium is a system of n electrons confined to a ring (J Chem Phys 138 (2013) 164124).
! As it is 1D, the different spin polarisations are degenerate, so we only consider the
! spin polarised case.

use const

implicit none

interface
    pure function digamma(arg) bind(c, name='psi')
        use, intrinsic :: iso_c_binding, only: c_double
        real(c_double), value, intent(in) :: arg
        real(c_double) :: digamma
    end function digamma
end interface

contains

    subroutine init_symmetry_ringium(sys)

        ! Initialise the symmetry components of sys (no spatial symmetry is implemented)

        use system, only: sys_t
        type(sys_t), intent(inout) :: sys

        sys%nsym = 1
        sys%sym0 = 1
        sys%sym_max = 1
        sys%nsym_tot = 1
        sys%sym0_tot = 1
        sys%sym_max_tot = 1

    end subroutine
 
    pure function get_two_e_int_ringium(sys, i, j, a, b) result(intgrl)

        ! In:
        !    sys: system being studied.
        !    i: index of basis function.
        !    j: index of basis function.
        !    a: index of basis function.
        !    b: index of basis function.

        ! Returns: 
        !   The anti-symmetrized integral <ij||ab> 

        use system, only: sys_t

        real(p) :: intgrl
        type(sys_t), intent(in) :: sys
        integer, intent(in) :: i, j, a, b
        real(p) :: x1, x2

        ! According to J Chem Phys 138, 164124 (2103)
        ! <ij||ab> = delta_(i+j,a+b) 1/(pi*R) (digamma(a-j+0.5)-digamma(a-i+0.5))

        if ((sys%basis%basis_fns(i)%l(1) + sys%basis%basis_fns(j)%l(1)) /= &
                & (sys%basis%basis_fns(a)%l(1) + sys%basis%basis_fns(b)%l(1))) then
            intgrl = 0.0_p
        else
            ! value of basis_fns%l is 2*lz
            x1 = (sys%basis%basis_fns(a)%l(1) - sys%basis%basis_fns(j)%l(1) + 1.0_p) * 0.5_p
            x2 = (sys%basis%basis_fns(a)%l(1) - sys%basis%basis_fns(i)%l(1) + 1.0_p) * 0.5_p
            intgrl = (digamma(x1) - digamma(x2)) / (pi * sys%ringium%radius)
        end if

    end function get_two_e_int_ringium

    pure function get_two_e_int_ringium_nonzero(sys, i, j, a) result(intgrl)

        ! In:
        !    sys: system being studied.
        !    i: index of basis function.
        !    j: index of basis function.
        !    a: index of basis function.

        ! Returns:
        !   The anti-symmetrized integral <ij||ab>

        ! NOTE:
        !   This assumes <ij||ab> is known to be non-zero by symmetry, otherwise an
        !   incorrect value will be returned.  If it might be zero, get_two_e_int_ringium
        !   must be used instead.  The index b is assumed to be such that angular momentum
        !   is conserved.

        use system, only: sys_t

        real(p) :: intgrl
        type(sys_t), intent(in) :: sys
        integer, intent(in) :: i, j, a
        real(p) :: x1, x2

        ! According to J Chem Phys 138, 164124 (2103)
        ! <ij||ab> = delta_(i+j,a+b) 1/(pi*R) (digamma(a-j+0.5)-digamma(a-i+0.5))

        ! value of basis_fns%l is 2*lz
        x1 = (sys%basis%basis_fns(a)%l(1) - sys%basis%basis_fns(j)%l(1) + 1.0) * 0.5
        x2 = (sys%basis%basis_fns(a)%l(1) - sys%basis%basis_fns(i)%l(1) + 1.0) * 0.5
        intgrl = (digamma(x1) - digamma(x2)) / (pi * sys%ringium%radius)

    end function get_two_e_int_ringium_nonzero

end module ringium_system

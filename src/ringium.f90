module ringium_system

! [review] - JSS: again, some brief notes and references on ringium would help, given it's relatively unusual.

use const

implicit none

! Routines for ringium

interface
    pure function psi(arg) bind(c)
        ! [review] - JSS: as psi is often interpreted as a wavefunction,
        ! [review] - JSS: comments (and perhaps a more meaningful name) would help.  
        ! [review] - JSS: e.g. function digamma(arg) bind(c, name='psi').
        use, intrinsic :: iso_c_binding, only: c_double
        real(c_double), value, intent(in) :: arg
        real(c_double) :: psi
    end function psi
end interface

contains

    subroutine init_symmetry_ringium(sys)

        ! [review] - JSS: no spatial symmetry (angular momentum apparently is, at least in the excitation generators).
        ! Initialise the symmetry components of sys (no symmetry is implemented)

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

        ! [review] - JSS: ref to where this is derived?
        ! <ij||ab> = delta_(i+j,a+b) 1/(pi*R) (Psi(a-j+0.5)-Psi(a-i+0.5))
        ! where Psi is the digamma function 

        if ((sys%basis%basis_fns(i)%l(1) + sys%basis%basis_fns(j)%l(1)) /= &
                & (sys%basis%basis_fns(a)%l(1) + sys%basis%basis_fns(b)%l(1))) then
            intgrl = 0.0_p
        else
            ! value of basis_fns%l is 2*lz
            ! [review] - JSS: is it worth having a fast integral function which disposes of symmetry checks?
            ! [review] - JSS: 0.5 is represented exactly as it's a power of 2 (so having a default single precision is ok...)
            x1 = (sys%basis%basis_fns(a)%l(1) - sys%basis%basis_fns(j)%l(1) + 1.0) * 0.5
            x2 = (sys%basis%basis_fns(a)%l(1) - sys%basis%basis_fns(i)%l(1) + 1.0) * 0.5
            intgrl = (psi(x1) - psi(x2)) / (pi * sys%ringium%radius)
        end if

    end function get_two_e_int_ringium

end module ringium_system

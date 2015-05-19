module kpoints

! Module for handling wavevectors.

implicit none

contains

    pure function calc_kinetic(sys, k) result(kinetic)

        ! In:
        !    sys: system of interest.  Requires lattice information and system
        !       type to be set and (if relevant) the Hubbard t parameter.
        !    k: wavevector in terms of the reciprocal lattice vectors of the
        !       crystal cell.
        ! Returns:
        !    The kinetic energy associated with a given wavevector.

        use const, only: p, pi
        use system, only: sys_t, hub_k, ueg, ringium

        real(p) :: kinetic

        type(sys_t), intent(in) :: sys
        integer, intent(in) :: k(sys%lattice%ndim)
        integer :: i
        real(p) :: kc(sys%lattice%ndim)

        ! Note sys%lattice%rlattice is stored in units of 2pi.
        forall (i=1:sys%lattice%ndim) kc(i) = sum(k*sys%lattice%rlattice(i,:)) + sys%k_lattice%ktwist(i)

        select case(sys%system)
        case(hub_k)
            ! For a square lattice the kinetic energy of a wavevector is given by
            !    -2t \sum_i cos(k.x_i)
            ! where x_i is the i-th reciprocal lattice vector of the primitive unit
            ! cell.
            kinetic = -2*sum(cos(2*pi*kc))*sys%hubbard%t
        case(ueg)
            ! For the UEG the kinetic energy of a wavevector is given by
            !    k^2/2
            ! in atomic units (the most sensible choice!)
            kinetic = 2*pi**2*dot_product(kc,kc)
        case(ringium)
            ! Not really a wavevector, but the kinetic energy of the function with
            ! angular momentum lz is
            !   l_z^2/2R^2
            ! but value in k is 2*l_z as l_z can be a half-integer so have
            !   k^2/8R^2
            kinetic = k(1)**2*0.125_p/sys%ringium%radius**2
        end select

    end function calc_kinetic

    pure function is_reciprocal_lattice_vector(sys, k) result(t_rlv)

        ! In:
        !    sys: system of interest.  Requires lattice information and system
        !       type to be set.
        !    k: wavevector in terms of the reciprocal lattice vectors of the
        !       crystal cell.
        ! Returns:
        !    True if k is a reciprocal lattice vector of the *primitive* cell.

        use const, only: p, depsilon
        use system, only: sys_t, ueg

        logical :: t_rlv
        type(sys_t), intent(in) :: sys
        integer, intent(in) :: k(sys%lattice%ndim)
        real(p) :: kc(sys%lattice%ndim)
        integer :: i

        select case(sys%system)
        case(ueg)
            ! Reciprocal primitive cell is infinitesimally small.
            t_rlv = all(k == 0)
        case default
            if (all(k == 0)) then
                ! Easy!
                t_rlv = .true.
            else
                ! Test to see if delta_k is 0 up to a reciprocal lattice vector
                ! of the primitive cell.
                forall (i=1:sys%lattice%ndim) kc(i) = sum(k*sys%lattice%rlattice(i,:))
                t_rlv = all(abs(kc - nint(kc)) < depsilon)
            end if
        end select

    end function is_reciprocal_lattice_vector

end module kpoints

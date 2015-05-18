module check_input
! Parse input options and check input for validity.

implicit none

contains

    subroutine check_sys(sys)

        ! Check the input parameters for the system.

        ! I don't pretend this is the most comprehensive of tests, but at least
        ! make sure a few things are not completely insane.

        ! In:
        !    sys: system object, as set in read_input

        use system

        use errors, only: stop_all, warning

        type(sys_t), intent(in) :: sys

        integer :: ivec, jvec
        character(*), parameter :: this='check_sys'

        if (sys%system /= heisenberg) then
            if (sys%nel <= 0) call stop_all(this,'Number of electrons must be positive.')
        end if

        if (sys%system /= ueg .and. sys%system /= read_in) then
            if (.not.(allocated(sys%lattice%lattice))) call stop_all(this, 'Lattice vectors not provided')
            do ivec = 1, sys%lattice%ndim
                do jvec = ivec+1, sys%lattice%ndim
                    if (dot_product(sys%lattice%lattice(:,ivec), sys%lattice%lattice(:,jvec)) /= 0) then
                        call stop_all(this, 'Lattice vectors are not orthogonal.')
                    end if
                end do
            end do
        end if

        if (sys%system == heisenberg) then
            if (sys%heisenberg%staggered_magnetic_field /= 0.0_p .and. (.not.sys%lattice%bipartite_lattice)) &
                call stop_all(this, 'Cannot set a staggered field&
                                    & for this lattice because it is frustrated.')
            if (sys%heisenberg%staggered_magnetic_field /= 0.0_p .and. sys%heisenberg%magnetic_field /= 0.0_p) &
                call stop_all(this, 'Cannot set a uniform and a staggered field at the same time.')
                        else if (sys%system == hub_k .or. sys%system == hub_real) then
            if (sys%nel > 2*sys%lattice%nsites) call stop_all(this, 'More than two electrons per site.')
        end if

        if (sys%system /= read_in) then
            if (sys%lattice%ndim > 3) call stop_all(this, 'Limited to 1, 2 or 3 dimensions')
            if (sys%system == ueg .and. sys%lattice%ndim == 1) call stop_all(this, 'UEG only functional in 2D and 3D')
        end if

    end subroutine check_sys

end module check_input

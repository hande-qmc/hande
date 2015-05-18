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

        if (sys%momentum_space .or. sys%system == read_in) then
            if (sys%real_lattice%finite_cluster) call stop_all(this,'FINITE_CLUSTER keyword only valid for&
                                      & calculations in real-space: ignoring keyword')
        end if

    end subroutine check_sys

    subroutine check_fciqmc_opts(sys, fciqmc_in)

        ! Check the FCIQMC specific options.

        ! In:
        !   sys: system being studied
        !   fciqmc_in: FCIQMC input options

        use qmc_data, only: fciqmc_in_t, single_basis, no_guiding, neel_singlet_guiding, neel_singlet
        use system, only: sys_t, heisenberg
        use errors, only: stop_all

        type(sys_t), intent(in) :: sys
        type(fciqmc_in_t), intent(in) :: fciqmc_in

        character(*), parameter :: this = 'check_fciqmc_opts'

        if (sys%system /= heisenberg) then
            if (fciqmc_in%trial_function /= single_basis) call stop_all(this, 'Only a single determinant can be used as the&
                                                     & reference state for this system. Other trial functions are not available.')
            if (fciqmc_in%guiding_function /= no_guiding) call stop_all(this, 'Importance sampling is only avaliable for the&
                                                                              & Heisenberg model currently.')
         else
            if (fciqmc_in%guiding_function == neel_singlet_guiding .and. fciqmc_in%trial_function /= neel_singlet) &
                call stop_all(this, 'This guiding function is only available when using the Neel singlet state&
                                    & as an energy estimator.')
        end if

        if (fciqmc_in%init_spin_inv_D0 .and. sys%Ms /= 0) then
            call stop_all(this, 'Flipping the reference state will give&
                                            & a state which has a different value of Ms and so cannot be used here.')
        end if

    end subroutine check_fciqmc_opts

    subroutine check_qmc_opts(qmc_in, simple_fciqmc)

        ! Check options common to QMC methods.

        ! In:
        !   qmc_in: QMC input options
        !   simple_fciqmc: true if using the simple FCIQMC algorithm.

        use qmc_data, only: qmc_in_t
        use errors, only: stop_all

        type(qmc_in_t), intent(in) :: qmc_in
        logical, intent(in) :: simple_fciqmc

        character(*), parameter :: this = 'check_qmc_opts'

        if (qmc_in%tau <= 0) call stop_all(this,'Tau not positive.')
        if (qmc_in%shift_damping <= 0) call stop_all(this,'Shift damping not positive.')

        if (.not. simple_fciqmc) then
            if (qmc_in%walker_length == 0) call stop_all(this,'Walker length zero.')
            if (qmc_in%spawned_walker_length == 0) call stop_all(this,'Spawned walker length zero.')
        end if

    end subroutine check_qmc_opts

end module check_input

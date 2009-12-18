program hubbard_fciqmc

implicit none

call init_calc()

call run_calc()

call end_calc()

contains

    subroutine init_calc()

        ! Initialise the calculation.
        ! Print out information about the compiled executable, 
        ! read input options and initialse the system and basis functions
        ! to be used.

        use report, only: environment_report
        use parse_input, only: read_input, check_input, distribute_input
        use system, only: init_system, system_type, hub_real
        use hubbard, only: init_basis_fns
        use determinants, only: init_determinants
        use parallel, only: init_parallel, parallel_report, nprocs, parent
        use hubbard_real, only: init_real_space_hub

        call init_parallel

        if (parent) then
            write (6,'(/,a8,/)') 'Hubbard'
            call environment_report()
        end if

        if (nprocs > 1 .and. parent) call parallel_report()

        if (parent) call read_input()

        call distribute_input()

        call init_system()

        call check_input()

        call init_basis_fns()

        call init_determinants()

        if (system_type == hub_real) call init_real_space_hub()

    end subroutine init_calc

    subroutine run_calc()

        ! Run the calculation based upon the input options.

        use hamiltonian
        use lanczos
        use full_diagonalisation

        use system
        use errors
        use determinants
        integer :: i

        do i = -nel, nel, 2
            call enumerate_determinants(i)
        end do
        call stop_all('run_calc','test')

        if (nprocs == 1) then
            if (t_exact .or. (t_lanczos .and. .not.direct_lanczos) ) then
                call generate_hamil(distribute_off)
            end if
        end if

        ! Lanczos.
        if (t_lanczos) then
            ! Construct the Hamiltonian matrix distributed over the processors
            ! if running in parallel.
            if (nprocs > 1 .and. .not.direct_lanczos) call generate_hamil(distribute_cols)
            call lanczos_diagonalisation()
        end if

        ! Exact diagonalisation.
        ! Warning: this destroys the Hamiltonian matrix...
        if (t_exact) then
            ! Construct the Hamiltonian matrix distributed over the processors
            ! if running in parallel.
            if (nprocs > 1) call generate_hamil(distribute_blocks)
            call exact_diagonalisation()
        end if

    end subroutine run_calc

    subroutine end_calc()

        ! Clean up time!

        use system, only: end_system, system_type, hub_real
        use hubbard, only: end_basis_fns
        use determinants, only: end_determinants
        use diagonalisation, only: end_hamil
        use parallel, only: end_parallel
        use hubbard_real, only: end_real_space_hub

        call end_system()
        call end_basis_fns()
        call end_determinants()
        call end_hamil()

        if (system_type == hub_real) call end_real_space_hub()

        call end_parallel()

    end subroutine end_calc

end program hubbard_fciqmc

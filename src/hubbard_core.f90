program hubbard_fciqmc


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
        use parse_input, only: read_input, check_input
        use system, only: init_system
        use hubbard, only: init_basis_fns
        use determinants, only: init_determinants
        use parallel

        if (parent == proc) write (6,'(/,a8,/)') 'Hubbard'

        call environment_report()

        call read_input()

        call init_system()

        call check_input()

        call init_basis_fns()

        call init_determinants()

    end subroutine init_calc

    subroutine run_calc()

        ! Run the calculation based upon the input options.

        use determinants, only: find_all_determinants
        use hamiltonian, only: t_exact, t_lanczos, generate_hamil, exact_diagonalisation, lanczos_diagonalisation

        call find_all_determinants()

        call generate_hamil()

        if (t_lanczos) call lanczos_diagonalisation()

        ! Warning: this destroys the Hamiltonian matrix...
        if (t_exact) call exact_diagonalisation()

    end subroutine run_calc

    subroutine end_calc()

        ! Clean up time!

        use system, only: end_system
        use hubbard, only: end_basis_fns
        use determinants, only: end_determinants
        use hamiltonian, only: end_hamil

        call end_system()
        call end_basis_fns()
        call end_determinants()
        call end_hamil()

    end subroutine end_calc

end program hubbard_fciqmc

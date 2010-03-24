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
        use symmetry, only: init_symmetry
        use calc, only: t_fciqmc, tsimple, seed
        use fciqmc, only: init_fciqmc
        use simple_fciqmc, only: init_simple_fciqmc
        use dSFMT_interface, only: dSFMT_init

        call init_parallel()

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

        call init_symmetry()

        if (system_type == hub_real) call init_real_space_hub()

        if (t_fciqmc) then
            call dSFMT_init(seed)
            if (tsimple) then
                call init_simple_fciqmc()
            else
                call init_fciqmc()
            end if
        end if

    end subroutine init_calc

    subroutine run_calc()

        ! Run the calculation based upon the input options.

        use calc, only: t_exact, t_lanczos, t_fciqmc, tsimple
        use diagonalisation, only: diagonalise
        use fciqmc, only: fciqmc_main
        use simple_fciqmc, only: do_simple_fciqmc

        if (t_exact .or. t_lanczos) call diagonalise()

        if (t_fciqmc) then
            if (tsimple) then
                call do_simple_fciqmc()
            else
                call fciqmc_main()
            end if
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
        use symmetry, only: end_symmetry

        call end_system()
        call end_basis_fns()
        call end_symmetry()
        call end_determinants()
        call end_hamil()

        if (system_type == hub_real) call end_real_space_hub()

        call end_parallel()

    end subroutine end_calc

end program hubbard_fciqmc

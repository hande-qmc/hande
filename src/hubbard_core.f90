program hubbard_fciqmc

implicit none

real :: start_time

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
        use excitations, only: init_excitations
        use parallel, only: init_parallel, parallel_report, iproc, nprocs, parent
        use hubbard_real, only: init_real_space_hub
        use momentum_symmetry, only: init_momentum_symmetry
        use calc

        call init_parallel()

        call cpu_time(start_time)

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

        call init_excitations()

        call init_momentum_symmetry()

        if (system_type == hub_real) call init_real_space_hub()

    end subroutine init_calc

    subroutine run_calc()

        ! Run the calculation based upon the input options.

        use calc
        use diagonalisation, only: diagonalise
        use dSFMT_interface, only: dSFMT_init
        use fciqmc_common, only: init_fciqmc
        use fciqmc, only: fciqmc_main
        use hilbert_space, only: estimate_hilbert_space
        use parallel, only: iproc, parent
        use simple_fciqmc, only: do_simple_fciqmc, init_simple_fciqmc
        use utils, only: int_fmt

        if (doing_calc(exact_diag+lanczos_diag)) call diagonalise()

        if (doing_calc(mc_hilbert_space)) then
            if (parent) then
                write (6,'(1X,a3,/,1X,3("-"),/)') 'RNG'
                write (6,'(1X,a51,'//int_fmt(seed,1)//',a1,/)') 'Initialised random number generator with a seed of:', seed, '.'
            end if
            call dSFMT_init(seed + iproc)
            call estimate_hilbert_space()
        end if

        if (doing_calc(fciqmc_calc+initiator_fciqmc+hfs_fciqmc_calc+ct_fciqmc_calc)) then
            if (parent) then
                write (6,'(1X,a3,/,1X,3("-"),/)') 'RNG'
                write (6,'(1X,a51,'//int_fmt(seed,1)//',a1,/)') 'Initialised random number generator with a seed of:', seed, '.'
            end if
            call dSFMT_init(seed + iproc)
            if (doing_calc(simple_fciqmc_calc)) then
                call init_simple_fciqmc()
                call do_simple_fciqmc()
            else
                call init_fciqmc()
                call fciqmc_main()
            end if
        end if

    end subroutine run_calc

    subroutine end_calc()

        ! Clean up time!

        use calc
        use system, only: end_system, system_type, hub_real
        use hubbard, only: end_basis_fns
        use determinants, only: end_determinants
        use excitations, only: end_excitations
        use diagonalisation, only: end_hamil
        use fciqmc_data, only: end_fciqmc
        use parallel, only: parent, end_parallel
        use hubbard_real, only: end_real_space_hub
        use momentum_symmetry, only: end_momentum_symmetry
        use report, only: end_report

        real :: end_time

        call end_system()
        call end_basis_fns()
        call end_momentum_symmetry()
        call end_determinants()
        call end_excitations()
        call end_hamil()

        if (system_type == hub_real) call end_real_space_hub()

        if (doing_calc(fciqmc_calc+initiator_fciqmc+hfs_fciqmc_calc+ct_fciqmc_calc)) call end_fciqmc()

        call cpu_time(end_time)

        if (parent) call end_report(end_time-start_time)

        call end_parallel()

    end subroutine end_calc

end program hubbard_fciqmc

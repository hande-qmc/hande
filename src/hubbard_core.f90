program hubbard_fciqmc

use report, only: environment_report
use parse_input, only: read_input, check_input
use system, only: init_system
use hubbard, only: init_basis_fns
use determinants, only: init_determinants, find_all_determinants
use hamiltonian, only: generate_hamil, exact_diagonalisation

call init_calc()

call run_calc()

contains

    subroutine init_calc()

        ! Initialise the calculation.
        ! Print out information about the compiled executable, 
        ! read input options and initialse the system and basis functions
        ! to be used.

        write (6,'(/,a8,/)') 'Hubbard'

        call environment_report()

        call read_input()

        call init_system()

        call check_input()

        call init_basis_fns()

        call init_determinants()

    end subroutine init_calc

    subroutine run_calc()

        ! Run the calculation based upon the input options.

        call find_all_determinants()

        call generate_hamil()

        call exact_diagonalisation()

    end subroutine run_calc

end program hubbard_fciqmc

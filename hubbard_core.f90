program hubbard_fciqmc

use report, only: environment_report
use parse_input, only: read_input
use hubbard, only: init_basis_fns

call init_calc()

contains

    subroutine init_calc()

        write (6,'(/,a8,/)') 'Hubbard'

        call environment_report()

        call read_input()

        call init_basis_fns()

    end subroutine init_calc

end program hubbard_fciqmc

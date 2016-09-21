module logging

! Module to contain functions handling logging outputs.

contains

    subroutine init_logging(logging_in, logging)
        use qmc_data, only: logging_t, logging_in_t

        type(logging_t), intent(inout) :: logging
        type(logging_in_t), intent(in) :: logging_in

        if (logging_in%calculation > 0) call init_logging_calc(logging_in, logging)
        !if (logging_in%spawning > 0) call init_logging_calc(logging_in, logging)
        !if (logging_in%death > 0) call init_logging_calc(logging_in, logging)
        !if (logging_in%annihilation > 0) call init_logging_calc(logging_in, logging)

    end subroutine init_logging

    subroutine end_logging(logging)
        use qmc_data, only: logging_t

        type(logging_t), intent(in) :: logging

        call end_logging_calc(logging)

    end subroutine end_logging

    subroutine init_logging_calc(logging_in, logging)

        use qmc_data, only: logging_t, logging_in_t
        use utils, only: get_free_unit
        use calc, only: calc_type, fciqmc_calc, ccmc_calc

        type(logging_t), intent(inout) :: logging
        type(logging_in_t), intent(in) :: logging_in
        integer :: iunit

        logging%calc_unit = get_free_unit()

        print *, 'Opening file ', logging_in%calc_filename
        open(logging%calc_unit, file=logging_in%calc_filename, status='unknown')

        call write_logging_calc_header(logging)

    end subroutine init_logging_calc

    subroutine end_logging_calc(logging)
        use qmc_data, only: logging_t

        type(logging_t), intent(in) :: logging

        if (logging%calc_unit /= huge(1)) close(logging%calc_unit, status='keep')

    end subroutine end_logging_calc

    subroutine write_logging_calc_header(logging)
        use calc, only: calc_type, fciqmc_calc, ccmc_calc
        use qmc_data, only: logging_t
        type(logging_t), intent(in) :: logging

        select case (calc_type)

        case(fciqmc_calc)
            write (logging%calc_unit, '(1X,"HANDE QMC Calculation Log File")')
            write (logging%calc_unit, '(1X,"==============================")')
            write (logging%calc_unit, '(1X,"Calculation type: FCIQMC")')
            write (logging%calc_unit, '(1X,"Verbosity Settings:")')
        case(ccmc_calc)
            write (logging%calc_unit, *) 'dsfgsg'
        end select
    end subroutine write_logging_calc_header

end module logging

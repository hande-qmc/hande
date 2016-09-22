module logging

! Module to contain functions handling logging outputs.

contains

    subroutine init_logging(logging_in, logging, sys, qs)

        ! Subroutine to initialise logs within HANDE, calling more specific
        ! functions for each specific type of logging activated by the user.


        use qmc_data, only: logging_t, logging_in_t, qmc_state_t
        use system, only: sys_t

        type(logging_t), intent(inout) :: logging
        type(logging_in_t), intent(in) :: logging_in
        type(qmc_state_t), intent(in) :: qs
        type(sys_t), intent(in) :: sys

        if (logging_in%calc > 0) call init_logging_calc(logging_in, logging)
        if (logging_in%spawn > 0) call init_logging_spawn(logging_in, logging)
        if (logging_in%death > 0) call init_logging_death(logging_in, logging)
        !if (logging_in%annihilation > 0) call init_logging_calc(logging_in, logging)

    end subroutine init_logging

    subroutine prep_logging_mc_cycle(iter, logging_in, logging_info, ndeath_tot, qn, cmplx_wfn)

        ! Subroutine to perform updates to logs required each iteration.
        ! This currently includes:
        !   - checking if within the range of iterations where logging is required.
        !   - printing updates for iteration number in all active logs.

        use qmc_data, only: logging_t, logging_in_t
        use const, only: int_p

        integer, intent(in) :: iter
        type(logging_t), intent(inout) :: logging_info
        type(logging_in_t), intent(in) :: logging_in
        integer(int_p), intent(inout) :: ndeath_tot
        logical, intent(in) :: qn, cmplx_wfn

        ! Check if we're within the required range to print logging info.
        logging_info%write_logging = (logging_in%start_iter <= iter .and. &
                                    iter <= logging_in%end_iter)
        ! Zero total ndeath accumulation
        ndeath_tot = 0_int_p
        ndeath_tot_im = 0_int_p

        if (logging_info%write_logging) then
            !if (logging_info%calc_unit /= huge(1)) call write_iter_to_log(iter, logging_info%calc_unit)
            if (logging_info%spawn_unit /= huge(1)) then
                call write_iter_to_log(iter, logging_info%spawn_unit)
                call write_logging_spawn_header(logging_info, qn, cmplx_wfn)
            end if
            if (logging_info%death_unit /= huge(1)) call write_iter_to_log(iter, logging_info%death_unit)
        end if

    end subroutine prep_logging_mc_cycle

    subroutine write_iter_to_log(iter, iunit)

        ! Writes iteration number to given io unit for logging.

        integer :: iter, iunit

        write(iunit, '("#",1X,20("="))')
        write(iunit, '("#",1X,"Iteration",1X,i10)') iter
        write(iunit, '("#",1X,20("="))')

    end subroutine write_iter_to_log

    subroutine end_logging(logging_info)

        ! Closes all active log files within fortran.
        ! In future may also print message to demonstrate end of logs.

        use qmc_data, only: logging_t

        type(logging_t), intent(in) :: logging_info

        if (logging_info%calc_unit /= huge(1)) close(logging_info%calc_unit, status='keep')
        if (logging_info%spawn_unit /= huge(1)) close(logging_info%spawn_unit, status='keep')
        if (logging_info%death_unit /= huge(1)) close(logging_info%death_unit, status='keep')

    end subroutine end_logging

    subroutine init_logging_calc(logging_in, logging)

        ! Initialises logging of high-level calculation information.
        ! This includes:
        !   - opening filename given and obtaining unit identifier.
        !   - converting from verbosity level given into specific
        !       information required in logs.
        !   - writing preamble information in log.

        use qmc_data, only: logging_t, logging_in_t

        type(logging_t), intent(inout) :: logging
        type(logging_in_t), intent(in) :: logging_in
        integer :: iunit

        print *, 'Opening file ', logging_in%calc_filename
        open(newunit=logging%calc_unit, file=logging_in%calc_filename, status='unknown')

        if (logging_in%calc > 0) logging%write_highlevel_values = .true.
        if (logging_in%calc > 1) logging%write_highlevel_calculations = .true.

        call write_logging_calc_header(logging)

    end subroutine init_logging_calc

    subroutine init_logging_spawn(logging_in, logging_info)

        ! Initialises logging of spawning information.
        ! This includes:
        !   - opening filename given and obtaining unit identifier.
        !   - converting from verbosity level given into specific
        !       information required in logs.
        !   - writing preamble information in log.

        use qmc_data, only: logging_t, logging_in_t
        use calc, only: calc_type, fciqmc_calc, ccmc_calc

        type(logging_t), intent(inout) :: logging_info
        type(logging_in_t), intent(in) :: logging_in
        integer :: iunit

        print *, 'Opening file ', logging_in%spawn_filename
        open(newunit=logging_info%spawn_unit, file=logging_in%spawn_filename, status='unknown')

        if (logging_in%spawn > 0) logging_info%write_successful_spawn = .true.
        if (logging_in%spawn > 1) logging_info%write_spawn_dets = .true.
        if (logging_in%spawn > 2) logging_info%write_failed_spawn = .true.

        call write_logging_spawn_preamble(logging_info)

    end subroutine init_logging_spawn

    subroutine init_logging_death(logging_in, logging)

        ! Initialises logging of death information.
        ! This includes:
        !   - opening filename given and obtaining unit identifier.
        !   - converting from verbosity level given into specific
        !       information required in logs.
        !   - writing preamble information in log.

        use qmc_data, only: logging_t, logging_in_t
        use calc, only: calc_type, fciqmc_calc, ccmc_calc

        type(logging_t), intent(inout) :: logging
        type(logging_in_t), intent(in) :: logging_in
        integer :: iunit

        print *, 'Opening file ', logging_in%death_filename
        open(newunit=logging%death_unit, file=logging_in%death_filename, status='unknown')

        if (logging_in%death > 0) logging%write_successful_death = .true.
        if (logging_in%death > 1) logging%write_death_det = .true.
        if (logging_in%death > 2) logging%write_failed_death = .true.

        call write_logging_death_header(logging)

    end subroutine init_logging_death

    subroutine write_logging_calc_header(logging_info)
        ! Write header and preamble for calculation log file.

        ! In:
        !    logging_info: contains information on logging settings.

        use calc, only: calc_type, fciqmc_calc, ccmc_calc
        use qmc_data, only: logging_t
        use qmc_io, only: write_column_title
        use report, only: environment_report

        type(logging_t), intent(in) :: logging_info

        write (logging_info%calc_unit, '(1X,"HANDE QMC Calculation Log File")')
        write (logging_info%calc_unit,'()')

        call environment_report(logging_info%calc_unit)

        select case (calc_type)
        case(fciqmc_calc)
            write (logging_info%calc_unit, '(1X,"Calculation type: FCIQMC")')
        case(ccmc_calc)
            write (logging_info%calc_unit, '(1X,"Calculation type: CCMC")')
        end select

        write (logging_info%calc_unit, '(1X,"Verbosity Settings:")')
        write (logging_info%calc_unit, '(1X,10X,"Write Calculation Values:",2X,L)') &
                        logging_info%write_highlevel_values
        write (logging_info%calc_unit, '(1X,10X,"Write Calculation Accumulation:",2X,L)') &
                        logging_info%write_highlevel_calculations

        write (logging_info%calc_unit,'()')

        write (logging_info%calc_unit,'("#")', advance='no')
        select case (calc_type)
        case(fciqmc_calc)
            call write_column_title(logging_info%calc_unit, "iter", int_val=.true., justify=1)
            call write_column_title(logging_info%calc_unit, "# spawn events", int_val=.true., justify=1)
            call write_column_title(logging_info%calc_unit, "# death events", int_val=.true., justify=1)
            call write_column_title(logging_info%calc_unit, "# attempts", int_val=.true., justify=1)
            write (logging_info%calc_unit,'()')
        case(ccmc_calc)
            write (logging_info%calc_unit,'()')

        end select
    end subroutine write_logging_calc_header

    subroutine write_logging_spawn_preamble(logging_info)

        ! Write initial preamble for top of spawn logging file.

        use qmc_data, only: logging_t
        use report, only: environment_report
        use calc, only: calc_type, fciqmc_calc, ccmc_calc

        type(logging_t), intent(in) :: logging_info

        write (logging_info%spawn_unit, '(1X,"HANDE QMC Spawning Log File")')
        write (logging_info%spawn_unit,'()')

        call environment_report(logging_info%spawn_unit)

        select case (calc_type)
        case(fciqmc_calc)
            write (logging_info%spawn_unit, '(1X,"Calculation type: FCIQMC")')
        case(ccmc_calc)
            write (logging_info%spawn_unit, '(1X,"Calculation type: CCMC")')
        end select

        write (logging_info%spawn_unit, '(1X,"Verbosity Settings:")')
        write (logging_info%spawn_unit, '(1X,10X,"Write Successful Spawns:",2X,L)') &
                        logging_info%write_successful_spawn
        write (logging_info%spawn_unit, '(1X,10X,"Write Failed Spawns:",2X,L)') &
                        logging_info%write_failed_spawn
        write (logging_info%spawn_unit, '(1X,10X,"Write Spawn Determinants:",2X,L)') &
                        logging_info%write_spawn_dets

        write (logging_info%spawn_unit,'()')


    end subroutine write_logging_spawn_preamble

    subroutine write_logging_spawn_header(logging_info, qn, cmplx_wfn)

        ! Write column headers for spawn logging information output.

        use qmc_data, only: logging_t
        use qmc_io, only: write_column_title

        type(logging_t), intent(in) :: logging_info
        logical, intent(in) :: qn, cmplx_wfn


        write (logging_info%spawn_unit,'("#")', advance='no')
        if (cmplx_wfn) then
            call write_column_title(logging_info%spawn_unit, "Re{H_ij}", justify=-1)
            call write_column_title(logging_info%spawn_unit, "Im{H_ij}", justify=-1)
        else
            call write_column_title(logging_info%spawn_unit, "H_ij", justify=-1)
        end if

        call write_column_title(logging_info%spawn_unit, "pgen", justify=-1)
        call write_column_title(logging_info%spawn_unit, "qn weighting", justify=-1)
        call write_column_title(logging_info%spawn_unit, "parent_sign", int_val=.true., justify=1)

        if (cmplx_wfn) then
            call write_column_title(logging_info%spawn_unit, "# spawn", int_val=.true., justify=1)
            call write_column_title(logging_info%spawn_unit, "# spawn im", int_val=.true., justify=1)
        else
            call write_column_title(logging_info%spawn_unit, "# spawn", int_val=.true., justify=1)
        end if
        write (logging_info%spawn_unit,'()')

    end subroutine write_logging_spawn_header

    subroutine write_logging_death_header(logging)

        use calc, only: calc_type, fciqmc_calc, ccmc_calc
        use qmc_data, only: logging_t

        type(logging_t), intent(in) :: logging

        select case (calc_type)
        case(fciqmc_calc)
            write (logging%death_unit, '(1X,"HANDE QMC Death Log File")')
            write (logging%death_unit, '(1X,"========================")')
            write (logging%death_unit, '(1X,"Calculation type: FCIQMC")')
            write (logging%death_unit, '(1X,"Verbosity Settings:")')
        case(ccmc_calc)
            write (logging%death_unit, *) 'dsfgsg'
        end select

    end subroutine write_logging_death_header

    subroutine write_logging_calc_fciqmc(logging_info, iter, nspawn_events, ndeath_events, nattempts)

        use qmc_io, only: write_qmc_var
        use const, only: int_p, int_64
        use qmc_data, only: logging_t

        integer, intent(in) :: iter
        integer(int_p), intent(in) :: nspawn_events, ndeath_events
        integer(int_64), intent(in) :: nattempts
        type(logging_t), intent(in) :: logging_info

        if (logging_info%write_logging .and. logging_info%write_highlevel_values) then
            write (logging_info%calc_unit,'(1X)', advance='no')
            call write_qmc_var(logging_info%calc_unit, iter)
            call write_qmc_var(logging_info%calc_unit, nspawn_events)
            call write_qmc_var(logging_info%calc_unit, ndeath_events)
            call write_qmc_var(logging_info%calc_unit, nattempts)
            write (logging_info%calc_unit,'()')
        end if

    end subroutine write_logging_calc_fciqmc

    subroutine write_logging_spawn(logging_info, hmatel, pgen, qn_weighting, nspawned, parent_sign, cmplx_wfn)

        ! Write log entry for a single spawning event.

        use qmc_io, only: write_qmc_var
        use const, only: int_p, int_64, p
        use qmc_data, only: logging_t
        use hamiltonian_data, only: hmatel_t

        type(logging_t), intent(in) :: logging_info
        type(hmatel_t), intent(in) :: hmatel
        real(p), intent(in) :: pgen, qn_weighting
        integer(int_p), intent(in) :: nspawned(:), parent_sign
        logical, intent(in) ::  cmplx_wfn

        if (logging_info%write_logging) then
            if ((logging_info%write_successful_spawn .and. any(abs(nspawned) > 0)) .or. &
                 logging_info%write_failed_spawn) then

                write (logging_info%spawn_unit,'(1X)', advance='no')

                if (cmplx_wfn) then
                    call write_qmc_var(logging_info%spawn_unit, real(hmatel%c))
                    call write_qmc_var(logging_info%spawn_unit, aimag(hmatel%c))
                else
                    call write_qmc_var(logging_info%spawn_unit, hmatel%r)
                end if

                call write_qmc_var(logging_info%spawn_unit, pgen)
                call write_qmc_var(logging_info%spawn_unit, qn_weighting)
                call write_qmc_var(logging_info%spawn_unit, parent_sign)
                call write_qmc_var(logging_info%spawn_unit, nspawned(1))

                if (cmplx_wfn) call write_qmc_var(logging_info%spawn_unit, nspawned(2))

                write (logging_info%spawn_unit,'()')
            end if
        end if

    end subroutine write_logging_spawn

end module logging

module logging

! Module to contain functions handling logging outputs.

! All calls to functions within this module should be preceded by 'if (debug)'.
! This ensures that in an optimised build all overhead involved with logging
! be removed at compile-time. However, this does mean that using this functionality
! requires a debug build (without some changes to the source code).

! Logging output is controlled by verbosity levels for each area of logging. For
! details of currently implemented output please see the manual.

! Currently only a rudimentary implementation as an initial effort. A more elegant
! approach is probably possible, for instance using generalised functions to initialise+
! write reports and enabling configurations to be set simply by changing initialisation or
! switching between preset  combination of variables.
! However, actually implementing this cleanly is no mean feat.

! To add a new level of logging output to a calculation to an existing log file, you must:
! 1) Add a (preferably descriptive) logical flag within
!       logging_t (+update json).
! 2) Add setting of flag appropriately in init function +
!       update checks in check_logging_inputs.
! 3) Update log preamble, header and log writer subroutine
!       appropriately to give new output.
!
! Adding a new type of log requires writing+calling versions of various functions
! (init_logging_x, write_logging_x_preamble, write_logging_x_header, write_logging_x)
! as well as addition of new parameters to logging_in_t and logging_t, but should be
! able to follow a similar pattern to that already established here.

use const, only: int_32, int_64

implicit none

type logging_in_t
    ! High-level debugging flag (at level of calculation running).
    integer(int_32) :: calc = 0
    character(255) :: calc_filename = 'CALC'
    ! Spawning flag.
    integer(int_32) :: spawn = 0
    character(255) :: spawn_filename = 'SPAWN'
    ! Death flag.
    integer(int_32) :: death = 0
    character(255) :: death_filename = 'DEATH'
    ! Iteration to start outputting logs from.
    integer(int_64) :: start_iter = 0_int_64
    ! Iteration to stop outputting logs from.
    integer(int_64) :: end_iter = huge(0_int_64)
end type logging_in_t

! Derived type to contain debugging flags and avoid passing lots of different
! flags to the various procedures.

type logging_t
    ! High-level debugging flag (at level of calculation running).
    logical :: write_highlevel_values = .false.
    integer :: calc_unit = huge(1_int_32)

    ! Spawning flags.
    logical :: write_successful_spawn = .false.
    logical :: write_failed_spawn = .false.
    integer :: spawn_unit = huge(1_int_32)

    ! Death flag.
    logical :: write_successful_death = .false.
    logical :: write_failed_death = .false.
    integer :: death_unit = huge(1_int_32)

    ! Whether within iterations required to output logging info.
    logical :: write_logging = .false.
end type logging_t

contains

! --- General functions for initialising+ending logging and writing reports ---

    subroutine init_logging(logging_in, logging_info)

        ! Subroutine to initialise logs within HANDE, calling more specific
        ! functions for each specific type of logging activated by the user.

        ! In:
        !   logging_in: input options relating to logging.
        ! In/Out:
        !   logging_info: derived type to be used during calculation to
        !       control logging output

        type(logging_t), intent(inout) :: logging_info
        type(logging_in_t), intent(in) :: logging_in

        call check_logging_inputs(logging_in)

        if (logging_in%calc > 0) call init_logging_calc(logging_in, logging_info)
        if (logging_in%spawn > 0) call init_logging_spawn(logging_in, logging_info)
        if (logging_in%death > 0) call init_logging_death(logging_in, logging_info)

    end subroutine init_logging

    subroutine check_logging_inputs(logging_in)

        ! Subroutine checking if verbosity level requested has been implemented
        ! and calling a warning if not.
        ! As additional functionality is added should be updated appropriately.

        use calc, only: fciqmc_calc, ccmc_calc, calc_type

        type(logging_in_t), intent(in) :: logging_in

        select case(calc_type)
        case(fciqmc_calc, ccmc_calc)
            continue ! All available logging implemented.
        case default
            if (logging_in%calc > 0) call write_logging_warning(1, logging_in%calc)
            if (logging_in%spawn > 0) call write_logging_warning(2, logging_in%spawn)
            if (logging_in%death > 0) call write_logging_warning(3, logging_in%death)
        end select

    end subroutine check_logging_inputs

    subroutine write_logging_warning(log_type, verbosity_level)

        ! Subroutine to write warning when logging verbosity level requested
        ! is not yet implemented.
        ! Changing this function will change all responses to such errors.

        use errors, only: warning
        use calc, only: calc_type, get_calculation_string
        use parallel, only: parent

        integer, intent(in) :: log_type, verbosity_level
        character(255) :: message, calc_name, log_name

        if (parent) then

            calc_name = get_calculation_string(calc_type)

            select case(log_type)
            case(1)
                log_name = "CALC"
            case(2)
                log_name = "SPAWN"
            case(3)
                log_name = "DEATH"
            end select

            write(message, '("Verbosity level ",i0," not yet implemented for ",a," logging with ",a)') &
                            verbosity_level, trim(log_name), trim(calc_name)

            call warning('check_logging_inputs',message)
        end if

    end subroutine write_logging_warning

    subroutine prep_logging_mc_cycle(iter, logging_in, logging_info, cmplx_wfn)

        ! Subroutine to perform updates to logs required each iteration.
        ! This currently includes:
        !   - checking if within the range of iterations where logging is required.
        !   - printing updates for iteration number in all active logs.

        ! In:
        !   iter: current iteration number.
        !   logging_in: input options relating to logging.
        !   cmplx_wfn: logical. True if using a complex wavefunction, false if not.
        ! In/Out:
        !   logging_info: derived type to be used during calculation to
        !       control logging output

        use const, only: int_p

        integer, intent(in) :: iter
        type(logging_t), intent(inout) :: logging_info
        type(logging_in_t), intent(in) :: logging_in
        logical, intent(in) :: cmplx_wfn

        ! Check if we're within the required range to print logging info.
        logging_info%write_logging = (logging_in%start_iter <= iter .and. &
                                    iter <= logging_in%end_iter)
        if (logging_info%write_logging) then
            if (logging_info%spawn_unit /= huge(1)) then
                call write_iter_to_log(iter, logging_info%spawn_unit)
                call write_logging_spawn_header(logging_info, cmplx_wfn)
            end if
            if (logging_info%death_unit /= huge(1)) then
                call write_iter_to_log(iter, logging_info%death_unit)
                call write_logging_death_header(logging_info)
            end if
        end if

    end subroutine prep_logging_mc_cycle

    subroutine end_logging(logging_info)

        ! Closes all active log files within fortran.
        ! In future may also print message to demonstrate end of logs.

        use report, only: end_report

        integer :: date_values(8)
        type(logging_t), intent(in) :: logging_info

        call date_and_time(VALUES=date_values)

        call write_date_time_close(logging_info%calc_unit, date_values)
        call write_date_time_close(logging_info%spawn_unit, date_values)
        call write_date_time_close(logging_info%death_unit, date_values)

    end subroutine end_logging

! --- Log-specific initialisation functions ---

    subroutine init_logging_calc(logging_in, logging_info)

        ! Initialises logging of high-level calculation information.
        ! This includes:
        !   - opening filename given and obtaining unit identifier.
        !   - converting from ierbosity level given into specific
        !       information required in logs.
        !   - writing preamble information in log.

        type(logging_t), intent(inout) :: logging_info
        type(logging_in_t), intent(in) :: logging_in

        open(newunit=logging_info%calc_unit, file=get_log_filename(logging_in%calc_filename), &
                status='unknown')

        if (logging_in%calc > 0) logging_info%write_highlevel_values = .true.

        call write_logging_calc_header(logging_info)

    end subroutine init_logging_calc

    subroutine init_logging_spawn(logging_in, logging_info)

        ! Initialises logging of spawning information.
        ! This includes:
        !   - opening filename given and obtaining unit identifier.
        !   - converting from verbosity level given into specific
        !       information required in logs.
        !   - writing preamble information in log.

        type(logging_t), intent(inout) :: logging_info
        type(logging_in_t), intent(in) :: logging_in

        open(newunit=logging_info%spawn_unit, file=get_log_filename(logging_in%spawn_filename), &
                status='unknown')

        if (logging_in%spawn > 0) logging_info%write_successful_spawn = .true.
        if (logging_in%spawn > 1) logging_info%write_failed_spawn = .true.

        call write_logging_spawn_preamble(logging_info)

    end subroutine init_logging_spawn

    subroutine init_logging_death(logging_in, logging_info)

        ! Initialises logging of death information.
        ! This includes:
        !   - opening filename given and obtaining unit identifier.
        !   - converting from verbosity level given into specific
        !       information required in logs.
        !   - writing preamble information in log.

        type(logging_t), intent(inout) :: logging_info
        type(logging_in_t), intent(in) :: logging_in

        open(newunit=logging_info%death_unit, file=get_log_filename(logging_in%death_filename), &
                status='unknown')

        if (logging_in%death > 0) logging_info%write_successful_death = .true.
        if (logging_in%death > 1) logging_info%write_failed_death = .true.

        call write_logging_death_preamble(logging_info)

    end subroutine init_logging_death

! --- Log-specific functions to write log preambles and headers ---

    subroutine write_logging_calc_header(logging_info)
        ! Write header and preamble for calculation log file.

        ! In:
        !    logging_info: contains information on logging settings.

        use calc, only: calc_type, fciqmc_calc, ccmc_calc
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
        write (logging_info%calc_unit, '(1X,10X,"Write Calculation Values:",2X,L1)') &
                        logging_info%write_highlevel_values

        write (logging_info%calc_unit,'()')

        write (logging_info%calc_unit,'("#")', advance='no')
        select case (calc_type)
        case(fciqmc_calc)
            call write_column_title(logging_info%calc_unit, "iter", int_val=.true., justify=1, sep=',')
            call write_column_title(logging_info%calc_unit, "# spawn events", int_val=.true., justify=1, sep=',')
            call write_column_title(logging_info%calc_unit, "# death particles", int_val=.true., justify=1, sep=',')
            call write_column_title(logging_info%calc_unit, "# attempts", int_val=.true., justify=1, sep=',')
            write (logging_info%calc_unit,'()')
        case(ccmc_calc)
            call write_column_title(logging_info%calc_unit, "iter", int_val=.true., justify=1, sep=',')
            call write_column_title(logging_info%calc_unit, "# spawn events", int_val=.true., justify=1, sep=',')
            call write_column_title(logging_info%calc_unit, "# death particles", int_val=.true., justify=1, sep=',')
            call write_column_title(logging_info%calc_unit, "# attempts", int_val=.true., justify=1, sep=',')
            call write_column_title(logging_info%calc_unit, "# D0 select", int_val=.true., justify=1, sep=',')
            call write_column_title(logging_info%calc_unit, "# stochastic", int_val=.true., justify=1, sep=',')
            call write_column_title(logging_info%calc_unit, "# single excit", int_val=.true., justify=1, sep=',')
            write (logging_info%calc_unit,'()')
        end select

    end subroutine write_logging_calc_header

    subroutine write_logging_spawn_preamble(logging_info)

        ! Write initial preamble for top of spawn logging file.

        ! In:
        !    logging_info: contains information on logging settings.

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
        write (logging_info%spawn_unit, '(1X,10X,"Write Successful Spawns:",2X,L1)') &
                        logging_info%write_successful_spawn
        write (logging_info%spawn_unit, '(1X,10X,"Write Failed Spawns:",2X,L1)') &
                        logging_info%write_failed_spawn

        write (logging_info%spawn_unit,'()')

    end subroutine write_logging_spawn_preamble

    subroutine write_logging_spawn_header(logging_info, cmplx_wfn)

        ! Write column headers for spawn logging information output.

        ! In:
        !    logging_info: derived type containing information on logging settings.
        !    cmplx_wfn: logical. True if using complex-valued wavefunction, false if not.

        use qmc_io, only: write_column_title

        type(logging_t), intent(in) :: logging_info
        logical, intent(in) :: cmplx_wfn


        write (logging_info%spawn_unit,'("#")', advance='no')
        if (cmplx_wfn) then
            call write_column_title(logging_info%spawn_unit, "Re{H_ij}", justify=-1, sep=',')
            call write_column_title(logging_info%spawn_unit, "Im{H_ij}", justify=-1, sep=',')
        else
            call write_column_title(logging_info%spawn_unit, "H_ij", justify=-1, sep=',')
        end if

        call write_column_title(logging_info%spawn_unit, "pgen", justify=-1, sep=',')
        call write_column_title(logging_info%spawn_unit, "qn weighting", justify=-1, sep=',')
        call write_column_title(logging_info%spawn_unit, "parent amplitude", justify=-1, sep=',')

        if (cmplx_wfn) then
            call write_column_title(logging_info%spawn_unit, "# spawn", int_val=.true., justify=1, sep=',')
            call write_column_title(logging_info%spawn_unit, "# spawn im", int_val=.true., justify=1, sep=',')
        else
            call write_column_title(logging_info%spawn_unit, "# spawn", int_val=.true., justify=1, sep=',')
        end if
        write (logging_info%spawn_unit,'()')

    end subroutine write_logging_spawn_header

    subroutine write_logging_death_preamble(logging_info)

        ! Write initial preamble for top of death logging file.

        ! In:
        !    logging_info: contains information on logging settings.

        use calc, only: calc_type, fciqmc_calc, ccmc_calc
        use report, only: environment_report

        type(logging_t), intent(in) :: logging_info

        write (logging_info%death_unit, '(1X,"HANDE QMC Death Log File")')
        write (logging_info%death_unit, '()')

        call environment_report(logging_info%death_unit)

        select case (calc_type)
        case(fciqmc_calc)
            write (logging_info%death_unit, '(1X,"Calculation type: FCIQMC")')
        case(ccmc_calc)
            write (logging_info%death_unit, '(1X,"Calculation type: CCMC")')
        end select

        write (logging_info%death_unit, '(1X,"Verbosity Settings:")')
        write (logging_info%death_unit, '(1X,10X,"Write Successful Deaths:",2X,L1)') &
                        logging_info%write_successful_death
        write (logging_info%death_unit, '(1X,10X,"Write Failed Deaths:",2X,L1)') &
                        logging_info%write_failed_death

        write (logging_info%spawn_unit,'()')

    end subroutine write_logging_death_preamble

    subroutine write_logging_death_header(logging_info)

        ! Write column headers for spawn logging information output.

        ! In:
        !    logging_info: derived type containing information on logging settings.

        use qmc_io, only: write_column_title

        type(logging_t), intent(in) :: logging_info

        write (logging_info%death_unit,'("#")', advance='no')

        call write_column_title(logging_info%death_unit, "Kii", sep=',')
        call write_column_title(logging_info%death_unit, "proj_energy", sep=',')
        call write_column_title(logging_info%death_unit, "loc_shift", sep=',')
        call write_column_title(logging_info%death_unit, "qn_weight", sep=',')
        call write_column_title(logging_info%death_unit, "p_death", sep=',')

        call write_column_title(logging_info%death_unit, "nkill", int_val = .true., justify=1, sep=',')

        call write_column_title(logging_info%death_unit, "init pop", justify=-1, sep=',')
        call write_column_title(logging_info%death_unit, "fin pop", justify=-1, sep=',')

        write (logging_info%death_unit,'()')

    end subroutine write_logging_death_header

! --- Log-specific subroutines to write log entries ---

    subroutine write_logging_calc_fciqmc(logging_info, iter, nspawn_events, ndeath_tot, nattempts)

        ! Write a single log entry for FCIQMC calculation-level information.

        ! In:
        !   logging_info: derived type containing information on logging.
        !   nspawn_events: integer. Number of spawning events in previous iteration.
        !   ndeath_tot: total number of particles created or destroyed by death in
        !       the previous iteration.
        !   nattempts: total number of spawning attempts to be made on this processor.

        use qmc_io, only: write_qmc_var
        use const, only: int_p, int_64

        integer, intent(in) :: iter, nspawn_events
        integer(int_p), intent(in) :: ndeath_tot
        integer(int_64), intent(in) :: nattempts
        type(logging_t), intent(in) :: logging_info

        if (logging_info%write_logging .and. logging_info%write_highlevel_values) then
            write (logging_info%calc_unit,'(1X)', advance='no')
            call write_qmc_var(logging_info%calc_unit, iter, sep=',')
            call write_qmc_var(logging_info%calc_unit, nspawn_events, sep=',')
            call write_qmc_var(logging_info%calc_unit, ndeath_tot, sep=',')
            call write_qmc_var(logging_info%calc_unit, nattempts, sep=',')
            write (logging_info%calc_unit,'()')
        end if

    end subroutine write_logging_calc_fciqmc

    subroutine write_logging_calc_ccmc(logging_info, iter, nspawn_events, ndeath_tot, nD0_select, &
                                        nclusters, nstochastic_clusters, nsingle_excitors)

        ! Write a single log entry for CCMC calculation-level information.

        ! In:
        !   logging_info: derived type containing information on logging.
        !   nspawn_events: integer. Number of spawning events in previous iteration.
        !   ndeath_tot: total number of particles created or destroyed by death in
        !       the previous iteration.
        !   nD0_select: total number of selection of the reference made this iteration.
        !   nclusters: total number of selections made this iteration.
        !   nstochastic_clusters: total number of stochastic selections made this iteration.
        !   nsingle_excitors: total number of deterministic selections made this iteration.

        use qmc_io, only: write_qmc_var
        use const, only: int_p, int_64

        type(logging_t), intent(in) :: logging_info
        integer, intent(in) :: iter
        integer, intent(in) :: nspawn_events
        integer(int_p), intent(in) :: ndeath_tot
        integer(int_64), intent(in) :: nD0_select, nclusters, nstochastic_clusters, nsingle_excitors

        if (logging_info%write_logging .and. logging_info%write_highlevel_values) then
            write (logging_info%calc_unit, '(1X)', advance='no')
            call write_qmc_var(logging_info%calc_unit, iter, sep=',')
            call write_qmc_var(logging_info%calc_unit, nspawn_events, sep=',')
            call write_qmc_var(logging_info%calc_unit, ndeath_tot, sep=',')
            call write_qmc_var(logging_info%calc_unit, nclusters, sep=',')
            call write_qmc_var(logging_info%calc_unit, nD0_select, sep=',')
            call write_qmc_var(logging_info%calc_unit, nstochastic_clusters, sep=',')
            call write_qmc_var(logging_info%calc_unit, nsingle_excitors, sep=',')
            write (logging_info%calc_unit, '()')
        end if

    end subroutine write_logging_calc_ccmc

    subroutine write_logging_spawn(logging_info, hmatel, pgen, qn_weighting, nspawned, parent_sign, cmplx_wfn)

        ! Write log entry for a single spawning event.

        ! In:
        !   logging_info: Derived type containing information on logging.
        !   hmatel: Derived type containing hamiltonian coupling element spawned across.
        !   pgen: real. Normalised probability of generating given excitation.
        !   qn_weighting: real. Reweighting factor for quasi-newton solvers. Always printed
        !       as still multiplied by but should be 1.000...
        !   nspawned: integer. Total signed walkers spawned in this event.
        !   parent_sign: real. Total signed population on parent determinant.
        !   cmplx_wfn: logical. True if using complex wavefunction, false if not.

        use qmc_io, only: write_qmc_var
        use const, only: int_p, p
        use hamiltonian_data, only: hmatel_t

        type(logging_t), intent(in) :: logging_info
        type(hmatel_t), intent(in) :: hmatel
        real(p), intent(in) :: pgen, qn_weighting, parent_sign
        integer(int_p), intent(in) :: nspawned(:)
        logical, intent(in) ::  cmplx_wfn

        if (logging_info%write_logging) then
            if ((logging_info%write_successful_spawn .and. any(abs(nspawned) > 0)) .or. &
                 (logging_info%write_failed_spawn .and. all(nspawned == 0))) then

                write (logging_info%spawn_unit,'(1X)', advance='no')

                if (cmplx_wfn) then
                    call write_qmc_var(logging_info%spawn_unit, real(hmatel%c), sep=',')
                    call write_qmc_var(logging_info%spawn_unit, aimag(hmatel%c), sep=',')
                else
                    call write_qmc_var(logging_info%spawn_unit, hmatel%r, sep=',')
                end if

                call write_qmc_var(logging_info%spawn_unit, pgen, sep=',')
                call write_qmc_var(logging_info%spawn_unit, qn_weighting, sep=',')
                call write_qmc_var(logging_info%spawn_unit, parent_sign, sep=',')
                call write_qmc_var(logging_info%spawn_unit, nspawned(1), sep=',')

                if (cmplx_wfn) call write_qmc_var(logging_info%spawn_unit, nspawned(2), sep=',')

                write (logging_info%spawn_unit,'()')
            end if
        end if

    end subroutine write_logging_spawn

    subroutine write_logging_death(logging_info, Kii, proj_energy, loc_shift, qn_weight, &
                nkill, pd, init_pop, fin_pop)

        ! Write log entry for a single death event.

        ! In:
        !   logging_info: Derived type containing information on logging.
        !   Kii: real. Diagonal Hamiltonian element of occupied determinant.
        !   proj_energy: real. Instantaneous projected energy.
        !   loc_shift: real. Instantaneous shift.
        !   qn_weighting: real. Reweighting factor for quasi-newton solvers. Always printed
        !       as still multiplied by but should be 1.000...
        !   nkill: integer. Total particle change in event.
        !   pd: real. pdeath of a single particle on same determinant.
        !   init_pop: real. Scaled initial population on determinant (ie. non-integer).
        !   fin_pop: real. Scaled final population on determinant (ie. non-integer).

        use qmc_io, only: write_qmc_var
        use const, only: int_p, p

        type(logging_t), intent(in) :: logging_info
        real(p), intent(in) :: Kii, proj_energy, qn_weight, loc_shift, pd, init_pop, fin_pop
        integer(int_p), intent(in) :: nkill

        if (logging_info%write_logging) then
            if ((logging_info%write_successful_death .and. abs(nkill) > 0) .or. &
                    logging_info%write_failed_death .and. nkill == 0) then

                write (logging_info%death_unit,'(1X)', advance='no')

                call write_qmc_var(logging_info%death_unit, Kii, sep=',')
                call write_qmc_var(logging_info%death_unit, proj_energy, sep=',')
                call write_qmc_var(logging_info%death_unit, loc_shift, sep=',')
                call write_qmc_var(logging_info%death_unit, qn_weight, sep=',')
                call write_qmc_var(logging_info%death_unit, pd, sep=',')

                call write_qmc_var(logging_info%death_unit, nkill, sep=',')

                call write_qmc_var(logging_info%death_unit, init_pop, sep=',')
                call write_qmc_var(logging_info%death_unit, fin_pop, sep=',')

                write (logging_info%death_unit,'()')
            end if
        end if

    end subroutine write_logging_death

! --- Generic helper functions ---

    function get_log_filename(in_name) result(out_name)

        ! Helper function to generate filenames to use for a generic log file.
        ! In:
        !   in_name: base name of log to give (eg. CALC)
        ! Out:
        !   out_name: filename to be used. Will be of
        !       form in_name.Y.pX.log for process number X.

        use parallel, only: iproc
        use utils, only: get_unique_filename
        character(255), intent(in) :: in_name
        character(255) :: suffix
        character(255) :: out_name
        integer :: id

        write(suffix,'(".p",i0,".log")') iproc
        call get_unique_filename(trim(in_name), suffix, .true., 0, out_name, id, .true.)

    end function get_log_filename

    subroutine write_iter_to_log(iter, iunit)

        ! Writes iteration number to given io unit for logging.

        ! In:
        !   iter: iteration number.
        !   iunit: io unit.

        integer :: iter, iunit

        write(iunit, '("#",1X,20("="))')
        write(iunit, '("#",1X,"Iteration",1X,i10)') iter
        write(iunit, '("#",1X,20("="))')

    end subroutine write_iter_to_log

    subroutine write_date_time_close(iunit, date_values)

        ! Write final message to log and close

        integer, intent(in) :: iunit, date_values(8)

        if (iunit /= huge(1)) then
            write (iunit,'(1X,64("="))')
            write (iunit,'(1X,a19,1X,i2.2,"/",i2.2,"/",i4.4,1X,a2,1X,i2.2,2(":",i2.2))') &
                      "Finished running on", date_values(3:1:-1), "at", date_values(5:7)
            write (iunit,'(1X,64("="))')

            close(iunit, status='keep')
        end if

    end subroutine write_date_time_close

! --- Json functions for logging types ---

    subroutine logging_in_t_json(js, logging_in, terminal)

        ! Serialise a logging_in_t object in JSON format.

        ! In/Out:
        !   js: json_out_t controlling the output unit and handling JSON internal state.  Unchanged on output.
        ! In:
        !   logging_in: logging_in_t object containing logging input values (including any defaults set).
        !   terminal (optional): if true, this is the last entry in the enclosing JSON object.  Default: false.

        use json_out

        type(json_out_t), intent(inout) :: js
        type(logging_in_t), intent(in) :: logging_in
        logical, intent(in), optional :: terminal

        call json_object_init(js, 'logging_in')
        call json_write_key(js, 'calc', logging_in%calc)
        call json_write_key(js, 'calc_file', logging_in%calc_filename)
        call json_write_key(js, 'spawn', logging_in%spawn)
        call json_write_key(js, 'spawn_file', logging_in%spawn_filename)
        call json_write_key(js, 'death', logging_in%death)
        call json_write_key(js, 'death_file', logging_in%death_filename)

        call json_write_key(js, 'start_iter', logging_in%start_iter)
        call json_write_key(js, 'end_iter', logging_in%end_iter, .true.)

        call json_object_end(js, terminal)

    end subroutine logging_in_t_json

    subroutine logging_t_json(js, logging_info, terminal)

        ! Serialise a logging_t object in JSON format.

        ! In/Out:
        !   js: json_out_t controlling the output unit and handling JSON internal state.  Unchanged on output.
        ! In:
        !   logging_info: logging_t object containing logging values used in calculation (including any
        !               defaults set).
        !   terminal (optional): if true, this is the last entry in the enclosing JSON object.  Default: false.

        use json_out

        type(json_out_t), intent(inout) :: js
        type(logging_t), intent(in) :: logging_info
        logical, intent(in), optional :: terminal

        call json_object_init(js, 'logging')
        call json_write_key(js, 'write_highlevel_values', logging_info%write_highlevel_values)
        call json_write_key(js, 'calc_unit', logging_info%calc_unit)

        call json_write_key(js, 'write_successful_spawn', logging_info%write_successful_spawn)
        call json_write_key(js, 'write_failed_spawn', logging_info%write_failed_spawn)
        call json_write_key(js, 'spawn_unit', logging_info%spawn_unit)

        call json_write_key(js, 'write_successful_death', logging_info%write_successful_death)
        call json_write_key(js, 'write_failed_death', logging_info%write_failed_death)
        call json_write_key(js, 'death_unit', logging_info%death_unit, .true.)

        call json_object_end(js, terminal)

    end subroutine logging_t_json

end module logging

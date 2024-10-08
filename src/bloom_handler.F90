module bloom_handler

    ! Module for handling blooms (single spawning events which create multiple particles).

    use const, only: p, int_p, dp, int_64

    implicit none

    enum,bind(c)
        ! Two options for defining blooms:
        ! 1) a set number of particles are created in a single spawning event.
        enumerator :: bloom_mode_fixedn
        ! 2) a number of particles which exceeds a fraction of the total population are
        !    created in a single spawning event.
        enumerator :: bloom_mode_fractionn
    end enum

    ! Stats about QMC blooming events.
    ! In general in CCMC a bloom is defined to have occured if a population of particles larger than 5% of
    ! the total population is spawned by any one cluster.
    ! In FCIQMC, a fixed particle population is used. The optimal value to use is not trivial to determine.
    type bloom_stats_t
        ! Note the mode is for information only; the QMC routine is responsible for
        ! calling the bloom handler procedures when a bloom is detected (the definition of
        ! which is left to the QMC routine).
        integer :: mode

        ! The factor by which the true populations are multiplied to create their
        ! encoded representations. This is not equal to 1, for example, when using
        ! real walker amplitudes.
        integer(int_p) :: encoding_factor = 1_int_p

        ! The proportion of the total number of particles which has to be spawned in one
        ! event to define a bloom.   Used (by convention) if mode = bloom_mode_fractionn.
        real(p) :: prop = 0.05_p

        ! The absolute number of particles which has to be spawned in one event to define
        ! a bloom. Used regardless of mode, but must be updated every iteration if using
        ! mode = bloom_mode_fractionn.
        integer :: threshold = 3
        ! nparticles in its encoded form, nparticles*encoding_factor.
        integer(int_p) :: threshold_encoded

        ! The maximum number of verbose warnings to printed out by a processor.
        integer :: nverbose_warnings = 1

        ! -- Report loop quantities (summed using energy_estimator comms) --
        ! Data is for the current processor before and all processors after energy_estimator comms.
        ! The total number of particles spawned by all blooms.
        real(p) :: tot_bloom_curr = 0.0_p
        ! The number of blooms in the current iteration (should be zeroed at the start of
        ! each iteration).
        integer(int_64) :: nblooms_curr = 0

        ! -- Current processor quantities --
        ! The maximum number of particles spawned by a single bloom on the current processor.
        real(p) :: max_bloom_proc = 0.0_p
        ! The number of warnings printed out on the current processor.
        integer :: nwarnings = 0

        ! -- All processor quantities (updated in write_bloom_report) --
        ! The total number of blooms.
        integer(int_64) :: nblooms = 0
        ! The total number of particles spawned by all blooms on the current processor.
        real(p) :: tot_bloom = 0.0_p
        ! Largest bloom
        real(p) :: max_bloom = 0.0_p

    end type bloom_stats_t

    contains

        subroutine init_bloom_stats_t(bloom_stats, mode, prop, threshold, nverbose_warnings, encoding_factor)

            ! Initialise a bloom_stats_t object.

            ! In:
            !    mode, prop, threshold, nverbose_warnings, encoding_factor: see
            !        description of matching components in bloom_stats_t object.
            ! Out:
            !    bloom_stats: initialised object for QMC blooming stats.

            integer, intent(in) :: mode
            integer, intent(in), optional :: threshold, nverbose_warnings
            real(p), intent(in), optional :: prop
            integer(int_p), intent(in), optional :: encoding_factor
            type(bloom_stats_t), intent(out) :: bloom_stats

            type(bloom_stats_t) :: bs_tmp

            ! Initialise to default values.
            bloom_stats = bs_tmp

            bloom_stats%mode = mode
            if (present(prop)) bloom_stats%prop = prop
            if (present(threshold)) bloom_stats%threshold = threshold
            if (present(nverbose_warnings)) bloom_stats%nverbose_warnings = nverbose_warnings
            if (present(encoding_factor)) bloom_stats%encoding_factor = encoding_factor

            bloom_stats%threshold_encoded = int(bloom_stats%threshold, int_p)*bloom_stats%encoding_factor

        end subroutine init_bloom_stats_t

        subroutine accumulate_bloom_stats(bloom_stats, nspawn)

            ! Detect if a blooming event has occured and accumulate/print data
            ! about it if so.
            ! Note that it left to the caller to ensure bloom_stats%threshold_encoded
            ! is set to appropriately detect blooms.

            ! In/Out:
            !     bloom_stats: stats object to update with the blooming event.
            ! In:
            !     nspawn: number of particles created in the blooming event, in
            !         the 'encoded' representation whereby the true number of
            !         spawns has been multiplied by bloom_stats%encoding_factor.

            use errors, only: stop_all

            type(bloom_stats_t), intent(inout) :: bloom_stats
            integer(int_p), intent(in) :: nspawn
            real(p) :: true_nspawn

            ! Check if a bloom event has occurred using previously set bloom threshold.
            if (abs(nspawn) >= bloom_stats%threshold_encoded) then

                ! Not thread safe (unless each thread has its own bloom_stats_t object,
                ! but blooms should not be frequent!
                !$omp critical

                ! Remove the encoding factor to create the true population, to be
                ! printed to the user.
                true_nspawn = real(nspawn,p)/bloom_stats%encoding_factor

                ! If the number of warnings exceeds the size of an integer then the
                ! population may have exploded.
                if ( bloom_stats%nblooms+bloom_stats%nblooms_curr == huge(bloom_stats%nblooms) ) then
                    call stop_all('accumulate_bloom_stats', 'Number of bloom events exceeded &
                         &the size of an integer: the population may have exploded')
                end if

                bloom_stats%nblooms_curr = bloom_stats%nblooms_curr + 1
                bloom_stats%tot_bloom_curr = bloom_stats%tot_bloom_curr + abs(true_nspawn)
                if(abs(true_nspawn) > bloom_stats%max_bloom_proc) then
                    bloom_stats%max_bloom_proc = abs(true_nspawn)
                end if

                !$omp end critical
            end if

        end subroutine accumulate_bloom_stats

        subroutine bloom_stats_init_report_loop(bloom_stats)

            ! Zero report-loop quantities in bloom_stats_t object.

            ! In/Out:
            !    bloom_stats: blooming stats.

            type(bloom_stats_t), intent(inout) :: bloom_stats

            bloom_stats%tot_bloom_curr = 0.0_p
            bloom_stats%nblooms_curr = 0

        end subroutine bloom_stats_init_report_loop

        subroutine bloom_stats_warning(bloom_stats, io_unit)

            ! Write out a bloom warning if we haven't hit the warning limit already.

            ! In/Out:
            !    bloom_stats: blooming stats.  nwarnings is incremented if a warning is printed out.
            ! In (optional):
            !    io_unit: io unit to write warning to.

            use utils, only: int_fmt

            type(bloom_stats_t), intent(inout) :: bloom_stats
            integer :: iunit
            integer, intent(in), optional :: io_unit

            iunit = 6
            if (present(io_unit)) iunit = io_unit

            if ( bloom_stats%nwarnings < bloom_stats%nverbose_warnings ) then
                select case(bloom_stats%mode)
                case(bloom_mode_fixedn)
                    write (iunit,'(1X, "# WARNING: more than", '//int_fmt(bloom_stats%threshold,1)//',&
                       &" particles spawned in a single event", '//int_fmt(bloom_stats%nblooms_curr,1)//',&
                       &" times in the last report loop.")') bloom_stats%threshold, bloom_stats%nblooms_curr
                case(bloom_mode_fractionn)
                    write (iunit,'(1X, "# WARNING: more than", '//int_fmt(int(bloom_stats%prop)*100,1)//',&
                       &"% of the total population spawned in a single event", '//int_fmt(bloom_stats%nblooms_curr,1)//',&
                       &" times in the last report loop.")') int(bloom_stats%prop*100), bloom_stats%nblooms_curr
                end select
                write (iunit,'(1X, "# Mean number of particles created in blooms: ",f8.1)') &
                    bloom_stats%tot_bloom_curr/bloom_stats%nblooms_curr
                write (iunit,'(1X,"# This warning only prints",'//int_fmt(bloom_stats%nverbose_warnings)//',&
                    &" time(s). You may wish to reduce the time step.")') bloom_stats%nverbose_warnings
                bloom_stats%nwarnings = bloom_stats%nwarnings + 1
            end if

        end subroutine bloom_stats_warning

        subroutine write_bloom_report(bloom_stats, io_unit)

            ! Prints a nice report on any blooming events which have occured.

            ! In:
            !     bloom_stats: stats about blooming to print
            ! In (optional):
            !    io_unit: io unit to write warning to.

            use parallel
            use utils, only: int_fmt

            type(bloom_stats_t), intent(inout) :: bloom_stats

            integer :: iunit
            integer, intent(in), optional :: io_unit
#ifdef PARALLEL
            integer :: ierr

            call mpi_reduce(bloom_stats%max_bloom_proc , bloom_stats%max_bloom, 1, &
                            MPI_preal, MPI_MAX, root, MPI_COMM_WORLD, ierr)
#else
            bloom_stats%max_bloom = bloom_stats%max_bloom_proc
#endif

            iunit = 6
            if (present(io_unit)) iunit = io_unit

            if (bloom_stats%nblooms > 0 .and. parent) then
                write (iunit, '(1X, "Blooming events occured: a more efficent calulation may be possible &
                    &with a smaller timestep.")')
                write (iunit, '(1X, "Total number of blooming events:", '//int_fmt(bloom_stats%nblooms,1)//')') &
                    bloom_stats%nblooms
                write (iunit, '(1X, "Maximum number of particles spawned in a blooming event:",f11.2)') &
                    bloom_stats%max_bloom
                write (iunit, '(1X, "Mean number of particles spawned in a blooming event:", 2X, f11.2, /)') &
                    bloom_stats%tot_bloom/bloom_stats%nblooms
            end if

        end subroutine write_bloom_report

        subroutine update_bloom_threshold_prop(bs, nparticles)

            ! Updates threshold and threshold_encoded in bloom_stats_t object
            ! according to previous population and preset proportion of the
            ! population that constitutes a bloom.

            ! In:
            !   nparticles: total population on previous iteration.
            ! In/Out:
            !   bs: bloom_stats_t object to be updated.

            type(bloom_stats_t), intent(inout) :: bs
            real(dp), intent(in) :: nparticles

            bs%threshold = int(real(nparticles,p)*bs%prop)
            bs%threshold_encoded = int(real(nparticles,p)*bs%prop*bs%encoding_factor, int_p)

        end subroutine update_bloom_threshold_prop

end module bloom_handler

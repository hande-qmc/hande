module bloom_handler

    ! Module for handling blooms (single spawning events which create multiple particles).

    use const, only: p, int_p

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
    ! A bloom is defined to have occured if a population of excips larger than 5% of
    ! the total population is spawned by any one cluster.
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
        ! event to define a bloom.   Used (by convention) if mode = bloom_mode_fixedn.
        real(p) :: prop = 0.05_p

        ! The absolute number of particles which has to be spawned in one event to define
        ! a bloom.   Used (by convention) if mode = bloom_mode_fractionn.
        integer :: n_bloom = 3
        ! n_bloom in its encoded form, n_bloom*encoding_factor.
        integer(int_p) :: n_bloom_encoded

        ! The number of blooms.
        integer :: nwarnings = 0
        ! The number of verbose warnings to print out (per processor, to avoid MPI comms).
        integer :: nverbose_warnings = 1
        ! The number of blooms in the current iteration (should be zeroed at the start of
        ! each iteration).
        integer :: nwarnings_curr = 0

        ! The maximum number of excips spawned by a single bloom.
        real(p) :: max_bloom = 0.0_p
        ! The total number of excips spawned by all blooms.
        real(p) :: tot_bloom = 0.0_p

    end type bloom_stats_t 

    contains

        subroutine init_bloom_stats_t(bloom_stats, mode, prop, n_bloom, nverbose_warnings, encoding_factor)

            ! Initialise a bloom_stats_t object.

            ! In:
            !    mode, prop, n_bloom, nverbose_warnings, encoding_factor: see
            !        description of matching components in bloom_stats_t object.
            ! Out:
            !    bloom_stats: initialised object for QMC blooming stats.

            integer, intent(in) :: mode
            integer, intent(in), optional :: n_bloom, nverbose_warnings
            real(p), intent(in), optional :: prop
            integer(int_p), intent(in), optional :: encoding_factor
            type(bloom_stats_t), intent(out) :: bloom_stats

            type(bloom_stats_t) :: bs_tmp

            ! Initialise to default values.
            bloom_stats = bs_tmp

            bloom_stats%mode = mode
            if (present(prop)) bloom_stats%prop = prop
            if (present(n_bloom)) bloom_stats%n_bloom = n_bloom
            if (present(nverbose_warnings)) bloom_stats%nverbose_warnings = nverbose_warnings
            if (present(encoding_factor)) bloom_stats%encoding_factor = encoding_factor

            bloom_stats%n_bloom_encoded = int(bloom_stats%n_bloom, int_p)*bloom_stats%encoding_factor

        end subroutine init_bloom_stats_t

        subroutine accumulate_bloom_stats(bloom_stats, nspawn)

            ! Accumulate/print data about a blooming event.  Note that it left to the
            ! caller to determine if a bloom has occured.

            ! In/Out:
            !     bloom_stats: stats object to update with the blooming event.
            ! In:
            !     nspawn: number of particles created in the blooming event, in
            !         the 'encoded' representation whereby the true number of
            !         spawns has been multiplied by bloom_stats%encoding_factor.

            use utils, only: int_fmt
            use errors, only: stop_all

            type(bloom_stats_t), intent(inout) :: bloom_stats
            integer(int_p), intent(in) :: nspawn
            real(p) :: true_nspawn

            ! Not thread safe (unless each thread has its own bloom_stats_t object, 
            ! but blooms should not be frequent!
            !$omp critical

            ! Remove the encoding factor to create the true population, to be
            ! printed to the user.
            true_nspawn = real(nspawn,p)/bloom_stats%encoding_factor

            ! If the number of warnings exceeds the size of an integer then the
            ! population may have exploded.
            if ( bloom_stats%nwarnings == huge(bloom_stats%nwarnings) ) then
                call stop_all('accumulate_bloom_stats', 'Number of bloom events exceeded &
                     &the size of an integer: the population may have exploded')
            end if

            ! Print out a warning message the first nverbose warning times only.
            if ( bloom_stats%nwarnings < bloom_stats%nverbose_warnings ) then
                select case(bloom_stats%mode)
                case(bloom_mode_fixedn)
                    write (6,'(1X, "# WARNING more than", '//int_fmt(bloom_stats%n_bloom,1)//',&
                       &" particles spawned in a single event.  Spawned: ", f8.1)') &
                        bloom_stats%n_bloom, true_nspawn
                case(bloom_mode_fractionn)
                    write (6,'(1X, "# WARNING more than", '//int_fmt(int(bloom_stats%prop*100))//',&
                       &" % of the total number of particles spawned in a single event.  # spawned: ", f8.1)') &
                        int(bloom_stats%prop*100), true_nspawn
                end select
                write (6,'(1X,"# This warning only prints",'//int_fmt(bloom_stats%nverbose_warnings)//',& 
                    &" time(s). You may wish to reduce the time step.")') bloom_stats%nverbose_warnings
            end if

            bloom_stats%nwarnings = bloom_stats%nwarnings + 1
            bloom_stats%nwarnings_curr = bloom_stats%nwarnings_curr + 1
            if(abs(true_nspawn) > bloom_stats%max_bloom) then
                bloom_stats%max_bloom = abs(true_nspawn)
            end if
            bloom_stats%tot_bloom = bloom_stats%tot_bloom + abs(true_nspawn)

            !$omp end critical

        end subroutine accumulate_bloom_stats
     
        subroutine write_bloom_report(bloom_stats)

            ! Prints a nice report on any blooming events which have occured.

            ! In:
            !     bloom_stats: stats about blooming to print

            use parallel
            use utils, only: int_fmt

            type(bloom_stats_t), intent(in) :: bloom_stats

            type(bloom_stats_t) :: bloom_stats_sum

#ifdef PARALLEL
            integer :: ierr

            bloom_stats_sum = bloom_stats
            call mpi_reduce(bloom_stats%nwarnings, bloom_stats_sum%nwarnings, 1, &
                            MPI_INTEGER, MPI_SUM, root, MPI_COMM_WORLD, ierr)
            call mpi_reduce(bloom_stats%tot_bloom, bloom_stats_sum%tot_bloom, 1, &
                            mpi_preal, MPI_SUM, root, MPI_COMM_WORLD, ierr)
            call mpi_reduce(bloom_stats%max_bloom, bloom_stats_sum%max_bloom, 1, &
                            MPI_preal, MPI_MAX, root, MPI_COMM_WORLD, ierr)
#else
            bloom_stats_sum = bloom_stats
#endif

            if(bloom_stats%nwarnings > 0 .and. parent) then
                write (6, '(1X, "Blooming events occured: a more efficent calulation may be possible &
                    &with a smaller timestep.")')
                write (6, '(1X, "Total number of blooming events:", '//int_fmt(bloom_stats%nwarnings,1)//')') &
                    bloom_stats%nwarnings
                write (6, '(1X, "Maxium number of excips spawned in a blooming event:",f11.2)') &
                    bloom_stats%max_bloom
                write (6, '(1X, "Mean number of excips spawned in a blooming event:", 2X, f11.2, /)') &
                    bloom_stats%tot_bloom/bloom_stats%nwarnings
            end if

        end subroutine write_bloom_report

end module bloom_handler

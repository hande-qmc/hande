module bloom_handler

    ! Module for handling blooms (single spawning events which create multiple particles).

    use const, only: p

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
        integer :: mode = bloom_mode_fixedn

        ! The proportion of the total number of particles which has to be spawned in one
        ! event to define a bloom.   Used (by convention) if mode = bloom_mode_fixedn.
        real(p) :: prop = 0.05_p

        ! The absolute number of particles which has to be spawned in one event to define
        ! a bloom.   Used (by convention) if mode = bloom_mode_fractionn.
        integer :: n_bloom = 3

        ! The number of blooms.
        integer :: nwarnings = 0
        ! The number of verbose warnings to print out (per processor, to avoid MPI comms).
        integer :: nverbose_warnings = 1
        ! The number of blooms in the current iteration (should be zeroed at the start of
        ! each iteration).
        integer :: nwarnings_curr = 0

        ! The maximum number of excips spawned by a single bloom.
        integer :: max_bloom = 0
        ! The total number of excips spawned by all blooms.
        real(p) :: tot_bloom = 0

    end type bloom_stats_t 

    contains

        subroutine accumulate_bloom_stats(bloom_stats, nspawn)

            ! Accumulate/print data about a blooming event.  Note that it left to the
            ! caller to determine if a bloom has occured.

            ! In/Out:
            !     bloom_stats: stats object to update with the blooming event.
            ! In:
            !     nspawn: number of particles created in the blooming event.

            use utils, only: int_fmt
            use errors, only: stop_all

            type(bloom_stats_t), intent(inout) :: bloom_stats
            integer, intent(in) :: nspawn

            ! Not thread safe (unless each thread has its own bloom_stats_t object, 
            ! but blooms should not be frequent!
            !$omp critical

            ! If the number of warnings exceeds the size of an integer then the
            ! population may have exploded.
            if ( bloom_stats%nwarnings == huge(bloom_stats%nwarnings) ) then
                call stop_all('accumulate_bloom_stats', 'Number of bloom events exceeded &
                     the size of an integer: the population may have exploded')
            end if

            ! Print out a warning message the first nverbose warning times only.
            if ( bloom_stats%nwarnings < bloom_stats%nverbose_warnings ) then
                select case(bloom_stats%mode)
                case(bloom_mode_fixedn)
                    write (6,'(1X, "# WARNING more than", '//int_fmt(bloom_stats%n_bloom,1)//',&
                        " particles spawned in a single event.  Spawned: ", '//int_fmt(nspawn,1)//')') &
                        bloom_stats%n_bloom, nspawn
                case(bloom_mode_fractionn)
                    write (6,'(1X, "# WARNING more than", '//int_fmt(int(bloom_stats%prop*100))//',&
                        " % of the total number of particles spawned in a single event.  # spawned: ", &
                        '//int_fmt(nspawn,1)//')') &
                        int(bloom_stats%prop*100), nspawn
                end select
                write (6,'(1X,"# This warning only prints",'//int_fmt(bloom_stats%nverbose_warnings)//',& 
                    " time(s). You may wish to reduce the time step.")'), bloom_stats%nverbose_warnings
            end if

            bloom_stats%nwarnings = bloom_stats%nwarnings + 1
            bloom_stats%nwarnings_curr = bloom_stats%nwarnings_curr + 1
            if(abs(nspawn) > bloom_stats%max_bloom) then
                bloom_stats%max_bloom = abs(nspawn)
            end if
            bloom_stats%tot_bloom = bloom_stats%tot_bloom + abs(nspawn)

            !$omp end critical

        end subroutine accumulate_bloom_stats
     
        subroutine write_bloom_report(bloom_stats)

            ! Prints a nice report on any blooming events which have occured.

            ! In:
            !     bloom_stats: stats about blooming to print

            use utils, only: int_fmt
            type(bloom_stats_t), intent(in) :: bloom_stats

            if(bloom_stats%nwarnings > 0) then
                write (6,'()')
                write (6, '(1X, "Blooming events occured: a more efficent calulation may be possible &
                    with a smaller timestep.")')
                write (6, '(1X, "Total number of blooming events:", '//int_fmt(bloom_stats%nwarnings,1)//')'), &
                    bloom_stats%nwarnings
                write (6, '(1X, "Maxium number of excips spawned in a blooming event:",&
                    '//int_fmt(bloom_stats%max_bloom,1)//')'), bloom_stats%max_bloom
                write (6, '(1X, "Mean number of excips spawned in a blooming event:", f11.2)'),&
                    bloom_stats%tot_bloom/bloom_stats%nwarnings
            end if

        end subroutine write_bloom_report

end module bloom_handler

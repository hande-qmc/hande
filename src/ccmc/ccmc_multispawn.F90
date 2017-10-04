module ccmc_multispawn

use const, only: i0, p

implicit none

type multispawn_stats_t
    integer :: nevents = 0
    integer :: nspawnings = 0
    integer :: nspawnings_max = 0
end type multispawn_stats_t

contains

    subroutine ms_stats_update(nspawnings, ms_stats)

        ! Update a multispawn_stats_t object.

        ! In:
        !    nspawnings: number of spawning attempts to be made from a single cluster selection.
        ! In/Out:
        !    ms_stats: On output, nevents is incremented and the nspawnings counts are increased accordingly.

        integer, intent(in) :: nspawnings
        type(multispawn_stats_t), intent(inout) :: ms_stats

        if (nspawnings > 1) then
            ms_stats%nevents = ms_stats%nevents + 1
            ms_stats%nspawnings = ms_stats%nspawnings + nspawnings
            ms_stats%nspawnings_max = max(ms_stats%nspawnings_max,nspawnings)
        end if

    end subroutine ms_stats_update

    pure function ms_stats_reduction(ms_stats_arr) result(ms_stats)

        ! In:
        !    ms_stats_arr: array of multispawn_stats_t objects.
        ! Returns:
        !    The multispawn_stats_t overall object (ie totals and max values from each component).

        type(multispawn_stats_t) :: ms_stats
        type(multispawn_stats_t), intent(in) :: ms_stats_arr(:)

        ms_stats%nevents = sum(ms_stats_arr%nevents)
        ms_stats%nspawnings = sum(ms_stats_arr%nspawnings)
        ms_stats%nspawnings_max = maxval(ms_stats_arr%nspawnings_max)

    end function ms_stats_reduction

    subroutine multispawn_stats_report(ms_stats)

        ! Print out a report of the cluster multispawn events.

        ! In:
        !    ms_stats: array of multispawn_stats_t objects.  The reduced version is printed out.

        use parallel
        use utils, only: int_fmt

        type(multispawn_stats_t), intent(in) :: ms_stats(:)
        type(multispawn_stats_t) :: ms_stats_total
#ifdef PARALLEL
        type(multispawn_stats_t) :: ms_stats_local
        integer :: ierr

        ms_stats_local = ms_stats_reduction(ms_stats)
        call mpi_reduce(ms_stats_local%nevents, ms_stats_total%nevents, 1, MPI_INTEGER, &
                        MPI_SUM, root, MPI_COMM_WORLD, ierr)
        call mpi_reduce(ms_stats_local%nspawnings, ms_stats_total%nspawnings, 1, MPI_INTEGER, &
                        MPI_SUM, root, MPI_COMM_WORLD, ierr)
        call mpi_reduce(ms_stats_local%nspawnings_max, ms_stats_total%nspawnings_max, 1, MPI_INTEGER, &
                        MPI_MAX, root, MPI_COMM_WORLD, ierr)
#else
        ms_stats_total = ms_stats_reduction(ms_stats)
#endif
        if (ms_stats_total%nevents > 0 .and. parent) then
            write (6,'(1X,"Multiple spawning events occurred.")')
            write (6,'(1X,"Number of multiple spawning events:",'//int_fmt(ms_stats_total%nevents,1)//')') &
                ms_stats_total%nevents
            write (6,'(1X,"Mean number of multiple spawning attempts per event:",2X,f11.2)') &
                real(ms_stats_total%nspawnings)/ms_stats_total%nevents
            write (6,'(1X,"Largest multiple spawning in a single event:",'//int_fmt(ms_stats_total%nspawnings_max,1)//',/)') &
                ms_stats_total%nspawnings_max
        end if

    end subroutine multispawn_stats_report

end module ccmc_multispawn

module ccmc_data

use const, only: i0, p, dp, int_64, depsilon
use determinants, only: det_info_t
use base_types, only: alloc_int2d, alloc_rdp1d

implicit none

! Work around Fortran's lack of arrays of pointers...
type bit_string_ptr
    integer(i0), pointer :: f(:)
end type bit_string_ptr

! Information about a cluster of excitors in addition to that stored in
! a det_info_t variable.
type cluster_t
    ! Pointers to the determinants formed by applying individual excitors to the
    ! reference determinant.
    type(bit_string_ptr), allocatable :: excitors(:) ! max: ex_level+2
    ! Number of excitors in cluster
    integer :: nexcitors
    ! Excitation level relative to the reference determinant of the determinant
    ! formed by applying the cluster of excitors to the reference determinant.
    integer :: excitation_level
    ! Overall complex component of amplitude of the cluster.
    complex(p) :: amplitude
    ! < D_i | a_i D_0 >, where D_i is the determinant formed by applying the
    ! cluster of excitors to the reference determinant.  Equal to +1 or -1.
    integer :: cluster_to_det_sign
    ! probability of selecting this cluster at random
    real(p) :: pselect
end type cluster_t

type multispawn_stats_t
    integer :: nevents = 0
    integer :: nspawnings = 0
    integer :: nspawnings_max = 0
end type multispawn_stats_t

! type to collect temporary information for qmc_state%excit_gen_data%p_single_double
! during the OpenMP loop in CCMC.
type p_single_double_stats_t
    real(p) :: h_pgen_singles_sum = 0.0_p ! hmatel/pgen sum for singles
    real(p) :: h_pgen_doubles_sum = 0.0_p ! hamtel/pgen sum for doubles
    real(p) :: excit_gen_singles = 0.0_p ! number of valid singles excitations created
    real(p) :: excit_gen_doubles = 0.0_p ! number of valid doubles excitations created
    logical :: overflow_loc = .false. ! true if precision has become too low for changes to be recognised.
end type p_single_double_stats_t

! Derived type to store all required information about the contribution selected
! (cluster, the determinant resulting from collapse and the partition in linked
! CCMC).
type wfn_contrib_t
    ! All cluster types we could possibly use in calculation.
    ! left_cluster and right_cluster will only be initialised if doing
    ! linked propogation.
    type(cluster_t) :: cluster
    type(cluster_t) :: left_cluster
    type(cluster_t) :: right_cluster
    ! All determinant types we could possibly use in calculation.
    ! ldet and rdet will only be initialised if doing
    ! linked propogation.
    type(det_info_t) :: cdet
    type(det_info_t) :: ldet
    type(det_info_t) :: rdet
end type wfn_contrib_t

! Information about possible clusters for each cluster size.
type selection_data_t
    ! Number of selections of each type within a calculation.
    integer(int_64) :: nD0_select
    integer(int_64) :: nclusters
    integer(int_64) :: nstochastic_clusters
    integer(int_64) :: nsingle_excitors
    ! Contains information on allowable cluster selections for each size.
    ! cluster_sizes_info(i) contains info for clusters of size i.
    ! cluster_sizes_info(i)%v(j,0) contains proportion of selections of size
    ! i that should use combination j.
    ! cluster_sizes_info(i)%v(j,k) contains number of excitors of excitation
    ! level k in allowed combination j.
    type(alloc_int2d), allocatable :: cluster_sizes_info(:) ! (ex_level+2)
    ! Info on the proportion of selections at a particular excitation level
    ! that should use a particular combination of excitor sizes.
    type(alloc_rdp1d), allocatable :: cluster_sizes_proportion(:) ! (ex_level+2)
    ! Proportion of selections at each cluster size.
    real(dp), allocatable :: size_weighting(:) ! (0:ex_level+2)
    ! Cumulative probability of selecting cluster of a given size.
    real(dp), allocatable :: cumulative_size_weighting(:) ! (0:ex_level+2)
    ! Weightings between stochastically selected clusters.
    real(dp), allocatable :: stoch_size_weighting(:) ! (2:ex_level+2)
    real(dp), allocatable :: cumulative_stoch_size_weighting(:) ! (2:ex_level+2)
    ! Average amplitude of \prod_i |N_I| / p_{select} for clusters of a given
    ! size.
    real(dp), allocatable :: average_amplitude(:) ! (0:ex_level+2)
    real(dp), allocatable :: variance_amplitude(:) ! (0:ex_level+2)
    ! Total number of successful selections of a given cluster size.
    integer(int_64), allocatable :: nsuccessful(:) ! (0:ex_level+2)
end type selection_data_t

type ex_lvl_dist_t
    ! Number of occupied excitors at each excitation level.
    integer, allocatable :: nstates_ex_lvl(:)
    ! Cumulative number of occupied excitors at each excitation level.
    integer, allocatable :: cumulative_nstates_ex_lvl(:)
    ! Population of excips on occupied excitors at each excitation level.
    real(p), allocatable :: pop_ex_lvl(:)
    ! Cumulative population of excips on occupied excitors at each excitation level.
    real(p), allocatable :: cumulative_pop_ex_lvl(:)
end type ex_lvl_dist_t

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

    subroutine multispawn_stats_report(ms_stats, io_unit)

        ! Print out a report of the cluster multispawn events.

        ! In:
        !    ms_stats: array of multispawn_stats_t objects.  The reduced version is printed out.
        ! In (optional):
        !    io_unit: io unit to write report to.

        use parallel
        use utils, only: int_fmt

        type(multispawn_stats_t), intent(in) :: ms_stats(:)
        type(multispawn_stats_t) :: ms_stats_total
        integer, intent(in), optional :: io_unit
        integer :: iunit
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
        iunit = 6
        if (present(io_unit)) iunit = io_unit

        if (ms_stats_total%nevents > 0 .and. parent) then
            write (iunit,'(1X,"Multiple spawning events occurred.")')
            write (iunit,'(1X,"Number of multiple spawning events:",'//int_fmt(ms_stats_total%nevents,1)//')') &
                ms_stats_total%nevents
            write (iunit,'(1X,"Mean number of multiple spawning attempts per event:",2X,f11.2)') &
                real(ms_stats_total%nspawnings)/ms_stats_total%nevents
            write (iunit,'(1X,"Largest multiple spawning in a single event:",'//int_fmt(ms_stats_total%nspawnings_max,1)//',/)') &
                ms_stats_total%nspawnings_max
        end if

    end subroutine multispawn_stats_report

    subroutine zero_ps_stats(ps_stats)

        ! Zero the ps_stats(:) components.

        ! In/Out:
        !    ps_stats: array of p_single_double_stats_t objects.

        type(p_single_double_stats_t), intent(inout) :: ps_stats(:)

        ps_stats%h_pgen_singles_sum = 0.0_p
        ps_stats%excit_gen_singles = 0.0_p
        ps_stats%h_pgen_doubles_sum = 0.0_p
        ps_stats%excit_gen_doubles = 0.0_p

    end subroutine zero_ps_stats
    
    subroutine ps_stats_reduction_update(ps, ps_stats)

        ! Reduce the ps_stats(:) components after the OpenMP loop and update ps.

        ! In:
        !    ps_stats: array of p_single_double_stats_t objects.
        ! In/Out:
        !    ps: the qmc_state%excit_gen_data%p_single_double_t object to be updated.

        use excit_gens, only: p_single_double_t

        type(p_single_double_stats_t), intent(in) :: ps_stats(:)
        type(p_single_double_t), intent(inout) :: ps

        real(p) :: excit_gen_singles_old, excit_gen_doubles_old

        excit_gen_singles_old = ps%excit_gen_singles
        excit_gen_doubles_old = ps%excit_gen_doubles

        ps%h_pgen_singles_sum = ps%h_pgen_singles_sum + sum(ps_stats%h_pgen_singles_sum)
        ps%excit_gen_singles = ps%excit_gen_singles + sum(ps_stats%excit_gen_singles)
        ps%h_pgen_doubles_sum = ps%h_pgen_doubles_sum + sum(ps_stats%h_pgen_doubles_sum)
        ps%excit_gen_doubles = ps%excit_gen_doubles + sum(ps_stats%excit_gen_doubles)
        
        ! check for overflow
        if (((abs(ps%excit_gen_singles-excit_gen_singles_old) < depsilon) .and. (sum(ps_stats%excit_gen_singles) > 0.0_p)) .or. &
            ((abs(ps%excit_gen_doubles-excit_gen_doubles_old) < depsilon) .and. (sum(ps_stats%excit_gen_doubles) > 0.0_p)) .or. &
             (any(ps_stats%overflow_loc))) then
            ps%overflow_loc = .true.
        end if

    end subroutine ps_stats_reduction_update
    
    subroutine end_selection_data(sd)

        ! Fully deallocate a passed in selection_data_t object.

        ! In:
        !   sd: selection_data_t object. Deallocated on output.


        use checking, only: check_deallocate
        use base_types, only: dealloc_int2d, dealloc_rdp1d

        type(selection_data_t), intent(inout) :: sd
        integer :: ierr

        if (allocated(sd%cluster_sizes_info)) then
            call dealloc_int2d(sd%cluster_sizes_info)
        end if
        if (allocated(sd%cluster_sizes_proportion)) then
            call dealloc_rdp1d(sd%cluster_sizes_proportion)
        end if
        if (allocated(sd%size_weighting)) then
            deallocate(sd%size_weighting, stat=ierr)
            call check_deallocate('sd%size_weighting', ierr)
        end if
        if (allocated(sd%cumulative_size_weighting)) then
            deallocate(sd%cumulative_size_weighting, stat=ierr)
            call check_deallocate('sd%cumulative_size_weighting', ierr)
        end if
        if (allocated(sd%stoch_size_weighting)) then
            deallocate(sd%stoch_size_weighting, stat=ierr)
            call check_deallocate('sd%stoch_size_weighting', ierr)
        end if
        if (allocated(sd%cumulative_stoch_size_weighting)) then
            deallocate(sd%cumulative_stoch_size_weighting, stat=ierr)
            call check_deallocate('sd%cumulative_stoch_size_weighting', ierr)
        end if
        if (allocated(sd%average_amplitude)) then
            deallocate(sd%average_amplitude, stat=ierr)
            call check_deallocate('sd%average_amplitude', ierr)
        end if
        if (allocated(sd%variance_amplitude)) then
            deallocate(sd%variance_amplitude, stat=ierr)
            call check_deallocate('sd%variance_amplitude', ierr)
        end if
        if (allocated(sd%nsuccessful)) then
            deallocate(sd%nsuccessful, stat=ierr)
            call check_deallocate('sd%nsuccessful', ierr)
        end if

    end subroutine end_selection_data

end module ccmc_data

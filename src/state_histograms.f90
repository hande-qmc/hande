module state_histograms

! Routines and functions to perform analysis on the
! excitation distribution of walkers in FCIQMC and DMQMC.
! [TODO] - WZV
!   1) Need to make a seperate readin/lua block for these
!       parameters probably to make the code more distinct?
!   2) clean up the various pieces of code to be more readable
!   3) Work the output into the stdout?
!   4) Rebase with main.
!   5) Does derived type location work here, does the name
!       need to be more clear? I.E. histogram_excit_dist.f90?
!

use const
implicit none

type state_histograms_data_t

    ! File names which will be used to write the data to.
    character(1024) :: states_per_bin_file
    character(1024) :: states_per_bin_per_exlevel_file

    ! Store the current seed, as in DMQMC this changes per
    ! beta loop. See 'init_dmqmc_beta_loop' for more information.
    integer :: current_seed

    ! Stores the number of mc cycles between the update
    ! and reporting of histogram data.
    integer :: histogram_frequency = -1

    ! The number of bins between each decade of psips
    integer :: bins_per_decade = 5

    ! Writing the psips (Nw) population as the general form:
    !   Nw = a \times 10^{b}
    ! we want to bin the population between decades based on these
    ! powers to encode the psips population. Reasonable
    ! values are 1 - 9 for a as this covers an entire decade,
    ! and 0 - 30 for the power b as this covers nearly all
    ! the possible machine precision (for 32 bit atleast)
    integer :: bmin = 0
    integer :: bmax = 30
    integer :: amin = 1
    integer :: amax = 9

    ! Some variables for tracking the largest 'a' and 'b' used
    ! in the calculation based on the input target population and
    ! the bins_per_decade.
    integer :: max_calc_a = 0
    integer :: max_calc_b = 0

    ! We want to keep track of the excitation levels psips can
    ! be contained on. So we store the minimum and maximum values.
    ! To start we just set them to values which should throw
    ! cause an error if they are never set for the calculation.
    integer :: max_ex_level = -1
    integer :: min_ex_level = 0

    ! Stores for the psips population, indexed by the
    ! a, b which define the bin the determinant's population belongs to.
    ! See bin, bmax, amin, amax above for more detail.
    real(p), allocatable :: particle_bins(:,:)

    ! Similar to particle_bins, used as the recieving array
    ! from MPI communication of the values on each proc.
    real(p), allocatable :: comm_bins(:,:)

    ! This stores the number of psips in a given excitation
    ! and particle number indexed again by "a" and "b" and the excitation
    ! level to the reference, and in DMQMC as well between the outer-product
    ! bitstring labels.
    real(p), allocatable :: excit_bins(:,:,:,:)

    ! Similar to excit_bins, used as the recieving array
    ! from MPI communication of the values on each proc.
    real(p), allocatable :: comm_excit_bins(:,:,:,:)

end type state_histograms_data_t

contains

    subroutine init_histogram_t(qmc_in, reference_in, hist, iunit)

        ! [TODO] WZV - Come back and add useful comments so the code
        ! can be digested by others and yourself in the future :).
        ! Don't forget the docstring!

        use parallel
        use qmc_data, only: qmc_in_t
        use reference_determinant, only: reference_t
        use checking, only: check_allocate
        use errors, only: stop_all

        type(qmc_in_t), intent(in) :: qmc_in
        type(reference_t), intent(in) :: reference_in

        type(state_histograms_data_t), intent(inout) :: hist

        integer, intent(in) :: iunit

        integer :: detpop_bytes_est, exlevelpop_bytes_est
        integer :: ierr

        ierr = 0

        hist%histogram_frequency = qmc_in%nreport
        if (qmc_in%state_histograms_freq /= -1) hist%histogram_frequency = qmc_in%state_histograms_freq

        hist%bins_per_decade = qmc_in%state_histograms_bpd

        hist%max_calc_a = hist%bins_per_decade
        hist%max_calc_b = floor(log10(qmc_in%target_particles)) + 3

        if (reference_in%max_ex_level == -1) then
            hist%max_ex_level = reference_in%ex_level
        else
            hist%max_ex_level = reference_in%max_ex_level
        end if

        ! Do a simple estimate of the memory requirment to write
        ! the bin and histogram data. Report to the user the size estimate
        ! and if its large exit. Should be an overestimate on the size.
        ! We skip this check if they provided the skip option.
        if (.not. qmc_in%state_histograms_mchk .and. parent) then
            ! Memory estimates in bytes for each file.
            ! bytes = columns * lines * bytes/(column*line)
            detpop_bytes_est = (hist%max_calc_a*(hist%max_calc_b + 2))*(qmc_in%nreport + 1)*18
            exlevelpop_bytes_est = ((hist%max_ex_level + 1)**2)*(hist%max_calc_b + 2)&
                                    *(mod(qmc_in%nreport, hist%histogram_frequency) + 1)*18
            write(iunit, '()')
            write(iunit, '(1X,"Memory estimate for state histograms DETPOPS files (MB):", f9.2)') 0.000001_p*detpop_bytes_est
            write(iunit, '(1X,"Memory estimate for state histograms EXLEVELPOPS files (MB):", f9.2)') 0.000001_p*exlevelpop_bytes_est
            write(iunit, '()')
            if (0.000001_p * (detpop_bytes_est + exlevelpop_bytes_est) > 1000.0_p) then
                call stop_all('init_histogram_t', 'The memory estimate for state&
                                histogram files is over 1000 (MB), if you acknowledge this warning and wish to&
                                proceed use state_histograms_mchk = true!')
            end if
        end if

        hist%current_seed = qmc_in%seed
        write(hist%states_per_bin_file, '("DETPOPS_RNG",I0)') qmc_in%seed
        write(hist%states_per_bin_per_exlevel_file, '("EXLEVELPOPS_RNG",I0,"_IREPORT",I0)') qmc_in%seed, 0

        allocate(hist%particle_bins(hist%amin:hist%amax, hist%bmin:hist%bmax), stat=ierr)
        call check_allocate('hist%particle_bins', (hist%amax-hist%amin)*(hist%bmax-hist%bmin), ierr)
        hist%particle_bins = 0.0_p

        allocate(hist%comm_bins(hist%amin:hist%amax, hist%bmin:hist%bmax), stat=ierr)
        call check_allocate('hist%comm_bins', (hist%amax-hist%amin)*(hist%bmax-hist%bmin), ierr)
        hist%comm_bins = 0.0_p

        allocate(hist%excit_bins(hist%amin:hist%amax, hist%bmin:hist%bmax, &
                    hist%min_ex_level:hist%max_ex_level, hist%min_ex_level:hist%max_ex_level), stat=ierr)
        call check_allocate('hist%excit_bins', (hist%amax-hist%amin)*(hist%bmax-hist%bmin)&
                                *((hist%max_ex_level-hist%min_ex_level)**2), ierr)
        hist%excit_bins = 0.0_p

        allocate(hist%comm_excit_bins(hist%amin:hist%amax, hist%bmin:hist%bmax, &
                    hist%min_ex_level:hist%max_ex_level, hist%min_ex_level:hist%max_ex_level), stat=ierr)
        call check_allocate('hist%comm_excit_bins', (hist%amax-hist%amin)*(hist%bmax-hist%bmin)&
                                *((hist%max_ex_level-hist%min_ex_level)**2), ierr)
        hist%comm_excit_bins = 0.0_p

    end subroutine init_histogram_t

    subroutine update_histogram_excitation_distribution(qs, f1, f2, real_pops, hist)

        ! [TODO] WZV - Come back and add useful comments so the code
        ! can be digested by others and yourself in the future :).
        ! Don't forget the docstring!

        use parallel
        use qmc_data, only: qmc_state_t
        use excitations, only: get_excitation_level

        type(qmc_state_t), intent(in) :: qs

        type(state_histograms_data_t), intent(inout) :: hist

        real(p), intent(in) :: real_pops

        integer(i0), intent(in) :: f1(:), f2(:)

        integer :: exlevel_1, exlevel_2, afac, bpow

        exlevel_1 = get_excitation_level(qs%ref%f0, f2)
        exlevel_2 = get_excitation_level(f1, f2)

        ! [TODO] WZV - Matches the original code but is really hard
        ! to read, see if there is a way to clean this up.
        bpow = floor(log10(abs(real_pops)))
        afac = floor((log10(abs(real_pops))-bpow)*hist%bins_per_decade) + 1

        associate(exbins => hist%excit_bins, bins => hist%particle_bins)
            exbins(afac, bpow, exlevel_1, exlevel_2) = exbins(afac, bpow, exlevel_1, exlevel_2) + 1.0_p
            bins(afac, bpow) = bins(afac, bpow) + 1.0_p
        end associate

    end subroutine update_histogram_excitation_distribution

    subroutine comm_and_report_histogram_excitation_distribution(hist, qmc_seed, ireport, start_calc, write_histogram, lazy_shift)

        ! [TODO] WZV - Come back and add useful comments so the code
        ! can be digested by others and yourself in the future :).
        ! Don't forget the docstring!

        use parallel

        type(state_histograms_data_t), intent(inout) :: hist

        logical, intent(in) :: start_calc, write_histogram

        integer, intent(in) :: qmc_seed, ireport
        integer, intent(in), optional :: lazy_shift

        integer :: ierr, afac, bpow, ex1, ex2, lazy_trunc
        real(p) :: detpop

        ierr = 0

        lazy_trunc = 0
        if (present(lazy_shift)) lazy_trunc = lazy_shift

        if (ireport > 0) then
            write(hist%states_per_bin_per_exlevel_file, '("EXLEVELPOPS_RNG",I0,"_IREPORT",I0)') hist%current_seed, ireport
        else if (hist%current_seed /= qmc_seed .and. ireport == 0) then
            write(hist%states_per_bin_file, '("DETPOPS_RNG",I0)') hist%current_seed
            write(hist%states_per_bin_per_exlevel_file, '("EXLEVELPOPS_RNG",I0,"_IREPORT",I0)') hist%current_seed, 0
        end if

        if (parent) then
            if (start_calc) then
                open(unit=117, file=trim(hist%states_per_bin_file), action='write')
            else
                open(unit=117, file=trim(hist%states_per_bin_file), access='append')
            end if

            if (write_histogram) then
                open(unit=343, file=trim(hist%states_per_bin_per_exlevel_file), action='write')
            end if
        end if

        associate(exbins => hist%excit_bins, bins => hist%particle_bins, &
                  comm_exbins => hist%comm_excit_bins, comm_bins => hist%comm_bins, &
                  bpd => hist%bins_per_decade)

#ifdef PARALLEL
            call mpi_reduce(bins, comm_bins, size(comm_bins), &
                            mpi_preal, mpi_sum, root, MPI_COMM_WORLD, ierr)
            call mpi_reduce(exbins, comm_exbins, size(comm_exbins), &
                            mpi_preal, mpi_sum, root, MPI_COMM_WORLD, ierr)
#else
            comm_bins = bins
            comm_exbins = exbins
#endif

            if (parent) then
                if (start_calc) then
                    do bpow = 0, hist%max_calc_b
                        do afac = 1, hist%max_calc_a

                            detpop = 10.0_p**(bpow + (1.0_p/bpd)*(afac - 1.0_p))
                            write(117, '(es18.10)', advance='no') detpop

                        end do
                    end do

                    write(117, '()')

                end if

                do bpow = 0, hist%max_calc_b
                    do afac = 1, hist%max_calc_a

                        write(117, '(es18.10)', advance='no') comm_bins(afac, bpow)

                    end do
                end do

                write(117, '()')

            end if

            ! Only want to print the full excitation level spectrums at the
            ! end of the calculation. 
            if (parent .and. write_histogram) then
                do ex1 = 0, hist%max_ex_level - lazy_trunc
                    do ex2 = 0, hist%max_ex_level

                        ! [TODO] WZV - Ask HP/James if add the walker column?
                        write(343, '(6X,A6,I3,I3)', advance='no') 'Ex.Lvl', ex1, ex2

                    end do
                end do

                write(343, '()')

                do bpow = 0, hist%max_calc_b
                    do afac = 1, hist%max_calc_a

                        ! [TODO] WZV - Ask HP/James if add the walker column?
                        write(343, '()')

                        do ex1 = 0, hist%max_ex_level - lazy_trunc
                            do ex2 = 0, hist%max_ex_level

                                write(343, '(es18.10)', advance='no') comm_exbins(afac, bpow, ex1, ex2)

                            end do
                        end do
                    end do
                end do

            end if

            bins = 0.0_p
            comm_bins = 0.0_p
            exbins = 0.0_p
            comm_exbins = 0.0_p

        end associate

        if (parent) then
            close(117)
            if (write_histogram) then
                close(343)
            end if
        end if

    end subroutine comm_and_report_histogram_excitation_distribution

    subroutine deallocate_histogram_t(hist)

        ! [TODO] WZV - Come back and add useful comments so the code
        ! can be digested by others and yourself in the future :).
        ! Don't forget the docstring!

        use parallel
        use checking, only: check_deallocate

        type(state_histograms_data_t), intent(inout) :: hist

        integer :: ierr

        ierr = 0

        if (allocated(hist%particle_bins)) then
            deallocate(hist%particle_bins, stat=ierr)
            call check_deallocate('hist%particle_bins', ierr)

            deallocate(hist%comm_bins, stat=ierr)
            call check_deallocate('hist%comm_bins', ierr)

            deallocate(hist%excit_bins, stat=ierr)
            call check_deallocate('hist%excit_bins', ierr)

            deallocate(hist%comm_excit_bins, stat=ierr)
            call check_deallocate('hist%comm_excit_bins', ierr)
        end if

    end subroutine deallocate_histogram_t

end module state_histograms

module state_histograms

! In this code we consider in general the excitation level
! between the orthoganol slater determinants used to label states on
! the density matrix (DMQMC) or wavefunction (FCIQMC) of a given system.
! This is a generalization of the excitation distribution analysis originally
! proposed by Cleland et al. See: https://doi.org/10.1021/ct300504f
!
! A brief explaination is given using H4/STO-3G in its ground state symmetry 
! block as an example.  The Hartree--Fock slater determinant bitstring would be:
!     | D_0 > = | 0 0 0 0 1 1 1 1 >
! For the density matrix, sites labeled using the slater determinants are:
!     \rho_ij = < D_i | \hat{\rho} | D_j >
! We can then define a matrix of excitations between those slater
! determinants < D_i | and | D_j > which would give us the following matrix:
!
!   0  1  1  2  1  1  2  2  2  2  2  2  2  2  3  3  2  3  3  4
!   1  0  2  2  2  2  1  3  1  1  3  2  2  3  2  2  2  4  2  3
!   1  2  0  2  2  2  3  1  1  3  1  2  2  3  2  2  2  2  4  3
!   2  2  2  0  2  2  2  2  2  2  2  2  2  2  2  2  4  2  2  2
!   1  2  2  2  0  2  1  3  3  3  1  2  2  1  2  4  2  2  2  3
!   1  2  2  2  2  0  3  1  3  1  3  2  2  1  4  2  2  2  2  3
!   2  1  3  2  1  3  0  4  2  2  2  2  2  2  1  3  2  3  1  2
!   2  3  1  2  3  1  4  0  2  2  2  2  2  2  3  1  2  1  3  2
!   2  1  1  2  3  3  2  2  0  2  2  2  2  4  1  1  2  3  3  2
!   2  1  3  2  3  1  2  2  2  0  4  2  2  2  3  1  2  3  1  2
!   2  3  1  2  1  3  2  2  2  4  0  2  2  2  1  3  2  1  3  2
!   2  2  2  2  2  2  2  2  2  2  2  0  4  2  2  2  2  2  2  2
!   2  2  2  2  2  2  2  2  2  2  2  4  0  2  2  2  2  2  2  2
!   2  3  3  2  1  1  2  2  4  2  2  2  2  0  3  3  2  1  1  2
!   3  2  2  2  2  4  1  3  1  3  1  2  2  3  0  2  2  2  2  1
!   3  2  2  2  4  2  3  1  1  1  3  2  2  3  2  0  2  2  2  1
!   2  2  2  4  2  2  2  2  2  2  2  2  2  2  2  2  0  2  2  2
!   3  4  2  2  2  2  3  1  3  3  1  2  2  1  2  2  2  0  2  1
!   3  2  4  2  2  2  1  3  3  1  3  2  2  1  2  2  2  2  0  1
!   4  3  3  2  3  3  2  2  2  2  2  2  2  2  1  1  2  1  1  0
!
! We refer to these excitation labels as exlevel_2 in the code.
!
! It is also relevent to consider the excitation between the Hartree--Fock
! and our site ij. This is useful as it allows us to preserve information
! about how far we are from the ground state row of the density matrix.
! So we write a matrix of excitations between < D_0 | and < D_i |, giving:
!
!   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
!   1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1
!   1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1
!   2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2
!   1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1
!   1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1
!   2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2
!   2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2
!   2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2
!   2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2
!   2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2
!   2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2
!   2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2
!   2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2
!   3  3  3  3  3  3  3  3  3  3  3  3  3  3  3  3  3  3  3  3
!   3  3  3  3  3  3  3  3  3  3  3  3  3  3  3  3  3  3  3  3
!   2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2
!   3  3  3  3  3  3  3  3  3  3  3  3  3  3  3  3  3  3  3  3
!   3  3  3  3  3  3  3  3  3  3  3  3  3  3  3  3  3  3  3  3
!   4  4  4  4  4  4  4  4  4  4  4  4  4  4  4  4  4  4  4  4
!
! We refer to these excitation labels as exlevel_1 in the code.
!
! Then concatenating matrix two with matrix one, we find
! a matrix of excitation labels which label each site on our
! density matrix. For the wavefunction (FCIQMC), the same analysis follows
! but we only consider the ground state "row" of the various matricies.
!
!   0,0  0,1  0,1  0,2  0,1  0,1  0,2  0,2  0,2  0,2  0,2  0,2  0,2  0,2  0,3  0,3  0,2  0,3  0,3  0,4
!   1,1  1,0  1,2  1,2  1,2  1,2  1,1  1,3  1,1  1,1  1,3  1,2  1,2  1,3  1,2  1,2  1,2  1,4  1,2  1,3
!   1,1  1,2  1,0  1,2  1,2  1,2  1,3  1,1  1,1  1,3  1,1  1,2  1,2  1,3  1,2  1,2  1,2  1,2  1,4  1,3
!   2,2  2,2  2,2  2,0  2,2  2,2  2,2  2,2  2,2  2,2  2,2  2,2  2,2  2,2  2,2  2,2  2,4  2,2  2,2  2,2
!   1,1  1,2  1,2  1,2  1,0  1,2  1,1  1,3  1,3  1,3  1,1  1,2  1,2  1,1  1,2  1,4  1,2  1,2  1,2  1,3
!   1,1  1,2  1,2  1,2  1,2  1,0  1,3  1,1  1,3  1,1  1,3  1,2  1,2  1,1  1,4  1,2  1,2  1,2  1,2  1,3
!   2,2  2,1  2,3  2,2  2,1  2,3  2,0  2,4  2,2  2,2  2,2  2,2  2,2  2,2  2,1  2,3  2,2  2,3  2,1  2,2
!   2,2  2,3  2,1  2,2  2,3  2,1  2,4  2,0  2,2  2,2  2,2  2,2  2,2  2,2  2,3  2,1  2,2  2,1  2,3  2,2
!   2,2  2,1  2,1  2,2  2,3  2,3  2,2  2,2  2,0  2,2  2,2  2,2  2,2  2,4  2,1  2,1  2,2  2,3  2,3  2,2
!   2,2  2,1  2,3  2,2  2,3  2,1  2,2  2,2  2,2  2,0  2,4  2,2  2,2  2,2  2,3  2,1  2,2  2,3  2,1  2,2
!   2,2  2,3  2,1  2,2  2,1  2,3  2,2  2,2  2,2  2,4  2,0  2,2  2,2  2,2  2,1  2,3  2,2  2,1  2,3  2,2
!   2,2  2,2  2,2  2,2  2,2  2,2  2,2  2,2  2,2  2,2  2,2  2,0  2,4  2,2  2,2  2,2  2,2  2,2  2,2  2,2
!   2,2  2,2  2,2  2,2  2,2  2,2  2,2  2,2  2,2  2,2  2,2  2,4  2,0  2,2  2,2  2,2  2,2  2,2  2,2  2,2
!   2,2  2,3  2,3  2,2  2,1  2,1  2,2  2,2  2,4  2,2  2,2  2,2  2,2  2,0  2,3  2,3  2,2  2,1  2,1  2,2
!   3,3  3,2  3,2  3,2  3,2  3,4  3,1  3,3  3,1  3,3  3,1  3,2  3,2  3,3  3,0  3,2  3,2  3,2  3,2  3,1
!   3,3  3,2  3,2  3,2  3,4  3,2  3,3  3,1  3,1  3,1  3,3  3,2  3,2  3,3  3,2  3,0  3,2  3,2  3,2  3,1
!   2,2  2,2  2,2  2,4  2,2  2,2  2,2  2,2  2,2  2,2  2,2  2,2  2,2  2,2  2,2  2,2  2,0  2,2  2,2  2,2
!   3,3  3,4  3,2  3,2  3,2  3,2  3,3  3,1  3,3  3,3  3,1  3,2  3,2  3,1  3,2  3,2  3,2  3,0  3,2  3,1
!   3,3  3,2  3,4  3,2  3,2  3,2  3,1  3,3  3,3  3,1  3,3  3,2  3,2  3,1  3,2  3,2  3,2  3,2  3,0  3,1
!   4,4  4,3  4,3  4,2  4,3  4,3  4,2  4,2  4,2  4,2  4,2  4,2  4,2  4,2  4,1  4,1  4,2  4,1  4,1  4,0
!
! When running this code then, we accumulate histograms for each unique site label
! generated from the concatenation of ex_level1 and ex_level2.
!     I.e., {0,0}, {0,1}, {0,2}, ... , {2,0}, ..., {4,2}, {4,3}, {4,4}
! The histograms store the count of determinants in a given population range
! indexed by the concatenated labels. We then report histograms for site labels
! in set intervals which can be controlled with the input parameters.

! [TODO] - WZV
!   1) Need to make a seperate readin/lua block for these
!       parameters probably to make the code more distinct?
!   2) clean up the various pieces of code to be more readable
!

use const
implicit none

type state_histograms_data_t

    ! File names which will be used to write the data to.
    character(1024) :: state_histograms_output_file

    ! Store the current seed, as in DMQMC this changes per
    ! beta loop. See 'init_dmqmc_beta_loop' for more information.
    integer :: current_seed

    ! Stores the number of mc cycles between the update
    ! and reporting of histogram data.
    integer :: histogram_frequency = -1

    ! Writing the psips (Nw) population as the general form:
    !   Nw = a \times 10^{b}
    ! Then we have a range of 'a' and 'b' values for a given calculation.
    ! The minimum 'a' is 1, and the minimum 'b' is 0 as a simulation
    ! will store at least 1 walker on a site in the simulation.
    ! Some variables for tracking the largest 'a' and 'b' used
    ! in the calculation are based on the input target population and
    ! the max_calc_a. These are set in 'init_histogram_t'.
    integer :: max_calc_a = 0
    integer :: max_calc_b = 0

    ! We want to keep track of the excitation levels psips can
    ! be contained on. So we store the maximum excitation level
    ! and set it to a value which will cause error's if it is never set.
    integer :: max_ex_level = -1

    ! This stores the number of psips in a given excitation
    ! and particle number indexed again by "a" and "b" and the excitation
    ! level to the reference, and in DMQMC the excitation level between
    ! the finite basis bitstring labels.
    ! shape = (max_calc_a, max_calc_b, max_ex_level + 1, max_ex_level + 1)
    real(p), allocatable :: excit_bins(:,:,:,:)

    ! Similar to excit_bins, used as the recieving array
    ! from MPI communication of the excit_bins from each proc.
    ! shape = (max_calc_a, max_calc_b, max_ex_level + 1, max_ex_level + 1)
    real(p), allocatable :: comm_excit_bins(:,:,:,:)

end type state_histograms_data_t

contains

    subroutine init_histogram_t(iunit, qmc_in, reference_in, hist, dmqmc_in)

        ! [TODO] WZV - Come back and add useful comments so the code
        ! can be digested by others and yourself in the future :).
        ! Don't forget the docstring!

        use parallel, only: parent
        use qmc_data, only: qmc_in_t
        use dmqmc_data, only: dmqmc_in_t
        use reference_determinant, only: reference_t
        use checking, only: check_allocate
        use errors, only: stop_all

        type(qmc_in_t), intent(in) :: qmc_in
        type(reference_t), intent(in) :: reference_in
        type(dmqmc_in_t), intent(in), optional :: dmqmc_in

        type(state_histograms_data_t), intent(inout) :: hist

        integer, intent(in) :: iunit

        integer :: bytes_est, nfiles
        integer :: ierr

        ierr = 0

        if (qmc_in%state_histograms_freq /= -1) then
            ! The input file has a report frequency
            hist%histogram_frequency = qmc_in%state_histograms_freq
        else
            ! If not set in the input, default to the end of calculation.
            hist%histogram_frequency = qmc_in%nreport
        end if

        hist%max_calc_a = qmc_in%state_histograms_bpd
        hist%max_calc_b = floor(log10(qmc_in%target_particles)) + 3

        if (reference_in%max_ex_level == -1) then
            hist%max_ex_level = reference_in%ex_level
        else
            hist%max_ex_level = reference_in%max_ex_level
        end if

        associate(amax => hist%max_calc_a, bmax => hist%max_calc_b, max_nex => hist%max_ex_level)

            ! Do a simple (and technically incorrect) estimate of the memory
            ! requirment to write all the histogram data. Report to the user
            ! the size estimate and if its large exit. This should be an
            ! overestimate on the size. One can skip this check if the user
            ! supplies state_histograms_mchk = false in the lua.
            if (parent) then
                ! Memory estimates in bytes for the file(s).
                ! bytes = columns * lines * bytes/(column*line)
                nfiles = (qmc_in%nreport / hist%histogram_frequency) + 1
                bytes_est = amax*(bmax+1)*(max_nex+1)*(max_nex+1)*nfiles*18
                if (present(dmqmc_in)) bytes_est = bytes_est*dmqmc_in%beta_loops

                write(iunit, '()')
                write(iunit, '(1X, a60, f9.2)') &
                    "Memory estimate for state histograms EXLEVELPOPS files (MB):", &
                    0.000001_p*bytes_est
                write(iunit, '()')

                if (qmc_in%state_histograms_mchk .and. 0.000001_p * bytes_est > 1000.0_p) then
                    call stop_all('init_histogram_t', 'The memory estimate for the state &
                                   histogram files is over 1000 (MB), if you acknowledge this &
                                   warning and wish to proceed add "state_histograms_mchk = false," &
                                   to the qmc lua block.')
                end if
            end if

            allocate(hist%excit_bins(1:amax, 0:bmax, 0:max_nex, 0:max_nex), stat=ierr)
            call check_allocate('hist%excit_bins', amax*(bmax+1)*(max_nex+1)*(max_nex+1), ierr)
            hist%excit_bins = 0.0_p

            allocate(hist%comm_excit_bins(1:amax, 0:bmax, 0:max_nex, 0:max_nex), stat=ierr)
            call check_allocate('hist%comm_excit_bins', amax*(bmax+1)*(max_nex+1)*(max_nex+1), ierr)
            hist%comm_excit_bins = 0.0_p

        end associate

    end subroutine init_histogram_t

    subroutine update_histogram_excitation_distribution(qs, f1, f2, real_pops, hist)

        ! [TODO] WZV - Come back and add useful comments so the code
        ! can be digested by others and yourself in the future :).
        ! Don't forget the docstring!

        use qmc_data, only: qmc_state_t
        use excitations, only: get_excitation_level

        type(qmc_state_t), intent(in) :: qs

        type(state_histograms_data_t), intent(inout) :: hist

        real(p), intent(in) :: real_pops

        integer(i0), intent(in) :: f1(:), f2(:)

        integer :: exlevel_1, exlevel_2, afac, bpow

        exlevel_1 = get_excitation_level(qs%ref%f0, f2)
        exlevel_2 = get_excitation_level(f1, f2)

        ! We want to write the current walker population in the form:
        !     Nw = a \times 10^{b}
        ! Then we can use a and b to find the histogram bin where this
        ! determinant belongs.
        ! The case of b is simple in that we take the floor of the base 10
        ! logarithm applied to the walker population.
        ! To get a, we find where a falls in the range (0, 1), and scale by
        ! the number of bins per decade (max_calc_a) then apply the floor.
        ! This results in the range [0, max_calc_a - 1], however we want to
        ! index the range [1, max_calc_a] so we add 1.
        bpow = floor(log10(abs(real_pops)))
        afac = floor((log10(abs(real_pops))-bpow)*hist%max_calc_a) + 1

        hist%excit_bins(afac, bpow, exlevel_1, exlevel_2) = hist%excit_bins(afac, bpow, exlevel_1, exlevel_2) + 1.0_p

    end subroutine update_histogram_excitation_distribution

    subroutine comm_and_report_histogram_excitation_distribution(hist, qmc_seed, ireport, lazy_shift)

        ! [TODO] WZV - Come back and add useful comments so the code
        ! can be digested by others and yourself in the future :).
        ! Don't forget the docstring!

        use parallel

        type(state_histograms_data_t), intent(inout) :: hist

        integer, intent(in) :: qmc_seed, ireport
        integer, intent(in), optional :: lazy_shift

        integer :: ierr, afac, bpow, ex1, ex2, lazy_trunc
        real(p) :: detpop

        ierr = 0

        lazy_trunc = 0
        if (present(lazy_shift)) lazy_trunc = lazy_shift

        write(hist%state_histograms_output_file, '("EXLEVELPOPS_RNG",I0,"_IREPORT",I0)') hist%current_seed, ireport

        if (parent) then
            open(unit=343, file=trim(hist%state_histograms_output_file), action='write')
        end if

        associate(exbins => hist%excit_bins, comm_exbins => hist%comm_excit_bins, bpd => hist%max_calc_a)

#ifdef PARALLEL
            call mpi_reduce(exbins, comm_exbins, size(exbins), &
                            mpi_preal, mpi_sum, root, MPI_COMM_WORLD, ierr)
#else
            comm_exbins = exbins
#endif

            if (parent) then

                write(343, '(9X, a9)', advance='no') 'bin_edges'

                do ex1 = 0, hist%max_ex_level - lazy_trunc
                    do ex2 = 0, hist%max_ex_level

                        write(343, '(4X,A6,I4,I4)', advance='no') 'Ex.Lvl', ex1, ex2

                    end do
                end do

                write(343, '()')

                do bpow = 0, hist%max_calc_b
                    do afac = 1, hist%max_calc_a

                        detpop = 10.0_p**(bpow + (1.0_p/bpd)*(afac - 1.0_p))
                        write(343, '(es18.10)', advance='no') detpop

                        do ex1 = 0, hist%max_ex_level - lazy_trunc
                            do ex2 = 0, hist%max_ex_level

                                write(343, '(es18.10)', advance='no') comm_exbins(afac, bpow, ex1, ex2)

                            end do
                        end do

                        write(343, '()')

                    end do
                end do

            end if

            exbins = 0.0_p
            comm_exbins = 0.0_p

        end associate

        if (parent) then
            close(343)
        end if

    end subroutine comm_and_report_histogram_excitation_distribution

    subroutine deallocate_histogram_t(hist)

        ! [TODO] WZV - Come back and add useful comments so the code
        ! can be digested by others and yourself in the future :).
        ! Don't forget the docstring!

        use checking, only: check_deallocate

        type(state_histograms_data_t), intent(inout) :: hist

        integer :: ierr

        ierr = 0

        deallocate(hist%excit_bins, stat=ierr)
        call check_deallocate('hist%excit_bins', ierr)

        deallocate(hist%comm_excit_bins, stat=ierr)
        call check_deallocate('hist%comm_excit_bins', ierr)

    end subroutine deallocate_histogram_t

end module state_histograms

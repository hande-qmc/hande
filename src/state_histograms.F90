module state_histograms

! In this code we consider in general the excitation level
! between the orthogonal Slater determinants used to label states on
! the density matrix (DMQMC) or wavefunction (FCIQMC) of a given system.
! This is a generalization of the excitation distribution analysis originally
! proposed by Cleland et al. See: https://doi.org/10.1021/ct300504f
!
! A brief explanation is given using H4/STO-3G in its ground state symmetry
! block as an example.  The Hartree--Fock Slater determinant bitstring would be:
!     | D_0 > = | 0 0 0 0 1 1 1 1 >
! For the density matrix, sites labeled using the Slater determinants are:
!     \rho_ij = < D_i | \hat{\rho} | D_j >
! We can then define a matrix of excitations between those Slater
! determinants < D_i | and | D_j > which would give us the following matrix:
! (*Note I believe that in the case of DMQMC, this is precisely what
!   the excitation distribution calculates when for example doing
!   Importance Sampling)
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
! It is also relevant to consider the excitation between the Hartree--Fock
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
! but we only consider the ground state "row" of the various matrices.
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

use const, only: p, int_64, i0
implicit none

type state_histograms_data_t

    ! File names which will be used to write the data to.
    character(1024) :: state_histograms_output_file

    ! Pad the report cycle in the file name for nice directory sorting
    character(1024) :: report_formatter

    ! Store the current seed, as in DMQMC this changes per
    ! beta loop. See 'init_dmqmc_beta_loop' for more information.
    ! We use the seed to label output files.
    integer :: current_seed

    ! Stores the number of reports between the reporting of histogram data.
    integer :: histogram_frequency = -1

    ! Writing the psips (Nw) population as the general form:
    !   Nw = a \times 10^{b}
    ! Then we have a range of 'a' and 'b' values for a given calculation.
    ! Some variables for tracking the largest 'a' and 'b' used
    ! in the calculation are based on the input target population and
    ! the number of bins we are using for the calculation.
    ! See 'init_state_histogram_t' for the full explanation.
    integer :: max_calc_a
    integer :: max_calc_b

    ! We want to keep track of the excitation levels psips can
    ! be contained on. So we store the maximum excitation level
    integer :: max_ex_level

    ! This stores the number of psips in a given excitation
    ! and particle number range indexed again by 'a' and 'b' and the excitation
    ! level to the reference, and in DMQMC the excitation level between
    ! the finite basis bitstring labels.
    ! shape = (max_calc_a, max_calc_b, max_ex_level + 1, max_ex_level + 1)
    integer(int_64), allocatable :: excit_bins(:,:,:,:)

    ! Similar to excit_bins, used as the receiving array
    ! from MPI communication of the excit_bins from each proc.
    ! shape = (max_calc_a, max_calc_b, max_ex_level + 1, max_ex_level + 1)
    integer(int_64), allocatable :: comm_excit_bins(:,:,:,:)

end type state_histograms_data_t

contains

    subroutine init_state_histogram_t(iunit, qmc_in, reference_in, hist, dmqmc_in)

        ! Set up the parameters and arrays for supporting the
        ! state histogram data collection and reporting.

        ! In:
        !    iunit: io unit to print information to.
        !    qmc_in: A type containing input options relating to QMC methods.
        !    reference_in: A type for the reference determinant, contains
        !       information about the systems maximum excitation level
        !    dmqmc_in (optional): if present, beta_loops is used to scale
        !       the memory estimate of the output files
        ! In/Out:
        !    hist: type containing all the state histograms information,
        !       this will have the system and calculation information as well
        !       as the arrays set-up upon exiting this routine.

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

        integer :: bytes_est, nfiles, ierr

        ierr = 0

        hist%current_seed = qmc_in%seed
        write(hist%report_formatter, '("(I0.",I0,")")') floor(log10(real(qmc_in%nreport))) + 3

        ! Set the report frequency, i.e. the number of report cycles
        ! between dumping the state histograms.
        if (qmc_in%state_histograms_nreport /= -1) then
            ! The user set report frequency
            hist%histogram_frequency = qmc_in%state_histograms_nreport
        else
            ! If not set in the input, default to the end of calculation.
            hist%histogram_frequency = qmc_in%nreport
        end if

        ! We want to write the walker population on states in the form:
        !     Nw = a \times 10^{b}
        ! It turns out to be convenient to write 'a' as the remainder of 'b',
        ! then we can bin the states based on the values of 'a' and 'b'.
        !   a: exists in the range [0, 1)
        !   b: {1, 2, 3, 4, ..., \inf}
        ! The right edge of a must be non-inclusive as 'a' and 'b' could simply
        ! be rewritten as b + 1 and a = 0.
        !
        ! Then we sub-divide the range of 'a' into bins, hence the histogram bit.
        ! The number of bins can be controlled via the input, but typically
        ! 5 bins is reasonable. This results in the following bins for 'a':
        !
        !   [0 -> 0.2), [0.2 -> 0.4), [0.4 -> 0.6], [0.6 -> 0.8), [0.8 -> 1)
        !
        ! Then to convert these bin ranges to an integer index, simply
        ! multiply by 5 and add 1, resulting in:
        !
        !   bin index : bin range
        ! -------------------------
        !           1 : [0 -> 0.2)
        !           2 : [0.2 -> 0.4)
        !           3 : [0.4 -> 0.6)
        !           4 : [0.6 -> 0.8)
        !           5 : [0.8 -> 1)
        !
        ! Then we can define the values 'a' can take based on the bins per decade
        ! or value of 'b'. Then the range of possible 'b' values falls within
        ! the walkers used in the simulation. Typically the maximum number of
        ! walkers plus three decades is a reasonable estimate.
        !
        ! Then along with these values of 'a' and 'b', all that remains to
        ! define the information for the state histograms are the possible
        ! excitation levels the exlevel_1 and exlevel_2 can take on.

        hist%max_calc_a = qmc_in%state_histograms_nbins
        hist%max_calc_b = floor(log10(qmc_in%target_particles)) + 3

        if (reference_in%max_ex_level == -1) then
            hist%max_ex_level = reference_in%ex_level
        else
            hist%max_ex_level = reference_in%max_ex_level
        end if

        associate(amax => hist%max_calc_a, bmax => hist%max_calc_b, max_nex => hist%max_ex_level)

            ! Do a simple (and technically incorrect) calculation of the memory
            ! requirement to write all the histogram data. Report to the user
            ! the size estimate and if its large exit. This should be an
            ! overestimate on the size. One can skip this check if the user
            ! supplies skip_histograms_mem_chk = false in the lua.
            if (parent) then
                ! Memory estimates in bytes for the file(s).
                nfiles = (qmc_in%nreport / hist%histogram_frequency) + 1
                bytes_est = amax*(bmax+1)*(max_nex+1)*(max_nex+1)*nfiles*18
                if (present(dmqmc_in)) bytes_est = bytes_est*dmqmc_in%beta_loops

                write(iunit, '()')
                write(iunit, '(1X, a60, f9.2)') &
                    "Memory estimate for state histograms EXLEVELPOPS files (MB):", &
                    0.000001_p*bytes_est
                write(iunit, '()')

                if (.not. qmc_in%skip_histograms_mem_chk .and. 0.000001_p * bytes_est > 1000.0_p) then
                    call stop_all('init_state_histogram_t', 'The memory estimate for the state &
                                   histogram files is over 1000 (MB), if you acknowledge this &
                                   warning and wish to proceed add "skip_histograms_mem_chk = false," &
                                   to the qmc lua block.')
                end if
            end if

            allocate(hist%excit_bins(1:amax, 0:bmax, 0:max_nex, 0:max_nex), stat=ierr)
            call check_allocate('hist%excit_bins', amax*(bmax+1)*(max_nex+1)*(max_nex+1), ierr)
            hist%excit_bins = 0_int_64

            allocate(hist%comm_excit_bins(1:amax, 0:bmax, 0:max_nex, 0:max_nex), stat=ierr)
            call check_allocate('hist%comm_excit_bins', amax*(bmax+1)*(max_nex+1)*(max_nex+1), ierr)
            hist%comm_excit_bins = 0_int_64

        end associate

    end subroutine init_state_histogram_t

    subroutine update_state_histogram(qs, f1, f2, real_pops, hist, icycle, ireport, final_report)

        ! Given the bitstrings and population for our determinant,
        ! use them to find find the population bin ('a' and 'b' index)
        ! as well as the excitation level indexes that a determinant belongs to.

        ! In:
        !   qs: qmc_state_t derived type with information on
        !       current calculation. Used for the reference bitstring.
        !    f1: The first bitstring label of \rho_ij (DMQMC)
        !       or the bitstring of C_i (FCIQMC)
        !    f2: The second bitstring label of \rho_ij (DMQMC)
        !       or the reference bitstring (FCIQMC)
        !    real_pops: The non-encoded walker population of \rho_ij/C_i
        !    icycle: The current mc cycle of the QMC method
        !    ireport: The current report cycle of the QMC method
        !    final_report (optional): A boolean used to ensure an update is
        !       performed regardless of the current ireport.
        ! In/Out:
        !    hist: type containing all the state histograms information,
        !       upon exiting, the excitation bins should be updated with
        !       the location of the current determinant.

        use qmc_data, only: qmc_state_t
        use excitations, only: get_excitation_level
        use errors, only: stop_all

        type(qmc_state_t), intent(in) :: qs
        type(state_histograms_data_t), intent(inout) :: hist

        real(p), intent(in) :: real_pops
        integer(i0), intent(in) :: f1(:), f2(:)
        integer, intent(in) :: icycle, ireport
        logical, intent(in), optional :: final_report

        logical :: skip_update_check
        integer :: exlevel_1, exlevel_2, iafac, ibpow

        skip_update_check = .false.
        if (present(final_report)) skip_update_check = final_report

        ! We only update the state histograms in certain intervals.
        if ((icycle /= 1 .or. mod(ireport - 1, hist%histogram_frequency) /= 0) .and. .not. skip_update_check) return

        ! exlevel_1 is the distance from the ground state row in DMQMC,
        ! while in FCIQMC we set f2 to qs%ref%f0, thus always results in 0.
        exlevel_1 = get_excitation_level(qs%ref%f0, f2)
        ! exlevel_2 is the excitation between the two labels of \rho_ij in DMQMC.
        ! while in FCIQMC it is the excitation from the ground state for C_i.
        exlevel_2 = get_excitation_level(f1, f2)

        ! See `init_state_histogram_t` for a more in depth explanation on why this works.
        ibpow = floor( log10( abs(real_pops)))
        iafac = floor( (log10( abs(real_pops)) - ibpow) * hist%max_calc_a) + 1

        if (ibpow > hist%max_calc_b) then
            call stop_all('update_state_histogram', 'The walker population has grown several orders &
                          of magnitude beyond the target population, this will cause a segfault!')
        end if

        hist%excit_bins(iafac, ibpow, exlevel_1, exlevel_2) = hist%excit_bins(iafac, ibpow, exlevel_1, exlevel_2) + 1_int_64

    end subroutine update_state_histogram

    subroutine comm_and_report_state_histogram(hist, ireport, iunit, index_shift, final_report)

        ! Share the histogram data with the root, and write out the information
        ! then zero the arrays so we can report again.

        ! In:
        !    ireport: The report cycle we are dumping the histogram data from.
        !    iunit: io unit to write all output to.
        !    index_shift (optional): Used in FCIQMC, equal to the maximum excitation
        !       of the simulation. A lazy way to truncate the histogram printing
        !       to only span those values relevent for FCIQMC.
        !    final_report (optional): A boolean used to ensure a report is
        !       performed regardless of the current ireport.
        ! In/Out:
        !    hist: type containing all the state histograms information,
        !       upon exiting, the excitation bins should be zerod.

        use parallel

        type(state_histograms_data_t), intent(inout) :: hist

        integer, intent(in) :: ireport, iunit
        integer, intent(in), optional :: index_shift
        logical, intent(in), optional :: final_report

        logical :: skip_report_check, file_exists
        integer :: ierr, iafac, ibpow, exlevel_1, exlevel_2, lazy_trunc, ifile
        real(p) :: detpop
        character(1024) :: ireport_string, tmp_output_file

        ierr = 0

        lazy_trunc = 0
        if (present(index_shift)) lazy_trunc = index_shift

        skip_report_check = .false.
        if (present(final_report)) skip_report_check = final_report

        ! We only report out the state histograms in certain intervals
        if (mod(ireport - 1, hist%histogram_frequency) /= 0 .and. .not. skip_report_check) return

        write(ireport_string, fmt=hist%report_formatter) ireport - 1
        write(hist%state_histograms_output_file, '("EXLEVELPOPS-RNG",I0,"-IREPORT",A)') hist%current_seed, trim(ireport_string)

        ifile = 0
        file_exists = .true.
        tmp_output_file = trim(hist%state_histograms_output_file)
        do while (file_exists)
            inquire(file=tmp_output_file, exist=file_exists)
            if (file_exists) then
                ifile = ifile + 1
                write(tmp_output_file, '(A, "-", I0)') trim(hist%state_histograms_output_file), ifile
            else
                hist%state_histograms_output_file = trim(tmp_output_file)
            end if
        end do

        if (parent) then
            write(iunit, '(1X,"# Writing state histogram to: ", A)') trim(hist%state_histograms_output_file)
            open(unit=343, file=trim(hist%state_histograms_output_file), action='write')
        end if

        associate(exbins => hist%excit_bins, comm_exbins => hist%comm_excit_bins, bpd => hist%max_calc_a)

#ifdef PARALLEL
            call mpi_reduce(exbins, comm_exbins, size(exbins), &
                            mpi_pop_integer, mpi_sum, root, MPI_COMM_WORLD, ierr)
#else
            comm_exbins = exbins
#endif

            if (parent) then

                write(343, '(9X, a9)', advance='no') 'bin_edges'

                do exlevel_1 = 0, hist%max_ex_level - lazy_trunc
                    do exlevel_2 = 0, hist%max_ex_level

                        write(343, '(4X,A6,I4,I4)', advance='no') 'Ex.Lvl', exlevel_1, exlevel_2

                    end do
                end do

                write(343, '()')

                do ibpow = 0, hist%max_calc_b
                    do iafac = 1, hist%max_calc_a

                        detpop = 10.0_p**(ibpow + (1.0_p/bpd)*(iafac - 1.0_p))
                        write(343, '(es18.10)', advance='no') detpop

                        do exlevel_1 = 0, hist%max_ex_level - lazy_trunc
                            do exlevel_2 = 0, hist%max_ex_level

                                write(343, '(es18.10)', advance='no') &
                                    real(comm_exbins(iafac, ibpow, exlevel_1, exlevel_2), p)

                            end do
                        end do

                        write(343, '()')

                    end do
                end do

            end if

            exbins = 0_int_64
            comm_exbins = 0_int_64

        end associate

        if (parent) then
            close(343)
        end if

    end subroutine comm_and_report_state_histogram

    subroutine fciqmc_state_histogram_final_report(qs, hist, ireport, iunit)

        ! For FCIQMC, we do one final update and report at the end of a calculation.
        ! This makes FCIQMC and DMQMC report state histograms in a similar fashion.
        ! Normally they would not due to DMQMC reporting data for the previous cycle.

        ! In:
        !    qs: qmc_state_t derived type with information on
        !        current calculation. Used for the reference bitstring.
        !    ireport: The current report cycle of the QMC method
        !    iunit: io unit to write all output to.
        ! In/Out:
        !    hist: type containing all the state histograms information,
        !        upon exiting, the excitation bins should be updated with
        !        the location of the current determinant.

        use qmc_data, only: qmc_state_t

        type(qmc_state_t), intent(in) :: qs
        type(state_histograms_data_t), intent(inout) :: hist

        integer, intent(in) :: ireport, iunit

        integer :: idet
        real(p) :: real_pops

        do idet = 1, qs%psip_list%nstates
            real_pops = real(qs%psip_list%pops(1,idet),p)/qs%psip_list%pop_real_factor
            call update_state_histogram(qs, qs%psip_list%states(:,idet), qs%ref%f0, &
                                       real_pops, hist, 1, ireport, final_report=.true.)
        end do

        call comm_and_report_state_histogram(hist, ireport + 1, iunit, &
                                             index_shift = hist%max_ex_level, final_report=.true.)

    end subroutine fciqmc_state_histogram_final_report

    subroutine deallocate_histogram_t(hist)

        ! Free the memory allocated for the arrays used to store histogram data.

        ! In/Out:
        !    hist: type containing all the state histograms information,
        !       upon exiting the excit_bins and comm_excit_bins should be
        !       freed from the memory.

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

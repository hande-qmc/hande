module fciqmc_data

! Data for fciqmc calculations and procedures which manipulate fciqmc and only
! fciqmc data.

use const
implicit none

!--- Input data: FCIQMC ---

! number of monte carlo cycles/report loop
integer :: ncycles
! number of report cycles
! the shift is updated and the calculation information printed out
! at the end of each report cycle.
integer :: nreport

! For DMQMC, beta_loops specifies the number of times
! the program will loop over each value of beta in the main loop.
integer :: beta_loops = 100

! timestep
real(p) :: tau

! Array sizes
! If these are < 0, then the values represent the number of MB to be used to
! store the main walker and spawned walker data respectively.
integer :: walker_length
integer :: spawned_walker_length

! Number of particles before which varyshift mode is turned on.
integer(lint) :: target_particles = 10000

! Don't bother renormalising generation probabilities; instead allow forbidden
! excitations to be generated and then rejected.
logical :: no_renorm = .false.

! probability of attempting single or double excitations...
! set to be nonsense value so can easily detect if it's given as an input option
real(p) :: pattempt_single = -1, pattempt_double = -1

!--- Input data: initiator-FCIQMC ---

! Complete active space within which a determinant is an initiator.
! (0,0) corresponds to the reference determinant only.
integer :: initiator_cas(2) = (/ 0,0 /)

! Population above which a determinant is an initiator.
integer :: initiator_population = 3

!--- Energy data ---

! shift: the shift is held constant at the initial value (from input) unless
! vary_shift is true.
! vary_shift_from_proje: if true, then the when variable shift mode is entered
! the shift is set to be the current projected energy.
! vary_shift_from: if vary_shift_from_proje is false, then the shift is set to
! this value when variable shift mode is entered.
! warning: if both initial_shift and vary_shift_from are set, then we expect the
! user to have been sensible.
real(p) :: shift = 0.0_p, vary_shift_from = 0.0_p
logical :: vary_shift_from_proje = .false.

! Initial shift, needed in DMQMC to reset the shift at the start of each
! beta loop.
real(p) :: initial_shift = 0.0_p

! Factor by which the changes in the population are damped when updating the
! shift.
real(p) :: shift_damping = 0.050_dp

! projected energy
! This stores during an FCIQMC report loop
!   \sum_{i/=0} <D_0|H|D_i> N_i
! where D_0 is the reference determinants and N_i is the walker population on
! determinant D_i.
! The projected energy is given as
!   <D_0|H|D_0> + \sum_{i/=0} <D_0|H|D_i> N_i/N_0
! and so proj_energy must be 'normalised' and averaged over the report loops
! accordingly.
! For convenience (especially in CCMC), where the numerator and denominator in
! the above equation are sampled and hence have to be weighted, we accumulate
! over each MC cycle and then over each report loop.
real(p) :: proj_energy, proj_energy_cycle

!--- Walker data ---

! Current number of walkers stored in the main list (processor dependent).
! This is updated during annihilation and merging of the spawned walkers into
! the main list.
integer :: tot_walkers

! Total number of particles on all walkers/determinants (processor dependent)
! Updated during death and annihilation and merging.
! The first element is the number of normal (Hamiltonian) particles.
! Subsequent elements are the number of Hellmann--Feynamnn particles.
integer(lint), allocatable :: nparticles(:) ! (sampling_size)
! Total number of particles across *all* processors, i.e. \sum_{proc} nparticles_{proc}
integer(lint), allocatable :: tot_nparticles(:) ! (sampling_size)

! Walker information: main list.
! sampling_size is one for each quantity sampled (i.e. 1 for standard
! FCIQMC/initiator-FCIQMC, 2 for FCIQMC+Hellmann--Feynman sampling).
integer :: sampling_size
! number of additional elements stored for each determinant in walker_data for
! (e.g.) importance sampling.
integer :: info_size
! a) determinants
integer(i0), allocatable, target :: walker_dets(:,:) ! (basis_length, walker_length)
! b) walker population
integer, allocatable, target :: walker_population(:,:) ! (sampling_size,walker_length)
! c) Walker information.  This contains:
! * Diagonal matrix elements, K_ii.  Storing them avoids recalculation.
!   K_ii = < D_i | H | D_i > - E_0, where E_0 = <D_0 | H | D_0> and |D_0> is the
!   reference determinant.  Always the first element.
! * Diagonal matrix elements for Hellmann--Feynmann sampling in 2:sampling_size
!   elements.
! * Further data in sampling_size+1:sampling_size:info_size.  For example, when
!   calculating the projected energy with various trial wavefunctions, it is
!   useful to store quantites which are expensive to calculate and which are
!   instead of recalculating them. For the Neel singlet state, the first component
!   gives the total number of spins up on the first sublattice. The second
!   component gives the number of 0-1 bonds where the 1 is on the first
!   sublattice.
real(p), allocatable, target :: walker_data(:,:) ! (sampling_size+info_size,walker_length)

! Walker information: spawned list.
! By combining the info in with the determinant, we can reduce the number of MPI
! communication calls during annihilation.
! a) array size.
! The size of each element in the spawned_walkers arrays depend upon what
! calculation is being done.  Each element has at least basis_length elements.
! * FCIQMC requires an additional element to store the population of the spawned
! walker.
! * initiator-FCIQMC requires a further additional element for information
! about the parent of the spawned walker.
! * Hellmann--Feynman sampling requires a further additional element for the
! population of the spawned Hellmann--Feynman walkers.

! spawned_walkers*(:basis_length,i) gives the determinant of the spawned walker.
! spawned_walkers*(spawned_pop,i) gives the population of the spawned walker.
! spawned_walkers*(spawned_hf_pop,i) gives the population of the spawned walker
! (Hellmann--Feynman sampling only).
! spawned_walkers*(spawned_parent,i) gives information about the parent of the
! spawned walker (initiator-FCIQMC only).
! spawned_hf_pop (if it exists) will always be equal to spawned_pop+1.

! In simple_fciqmc we only need to store the walker populations, so spawned_size
! is 1.
integer :: spawned_size
integer :: spawned_pop, spawned_parent, spawned_hf_pop
! b) determinants and the spawn times of the progeny (only used for ct_fciqmc)
integer(i0), allocatable, target :: spawned_walkers1(:,:) ! (spawned_size, spawned_walker_length)
integer(i0), allocatable, target :: spawned_walkers2(:,:) ! (spawned_size, spawned_walker_length)
real(p), allocatable :: spawn_times(:) ! (spawned_walker_length)
! c) pointers.
! In serial we only use spawned_walker_*1.  In parallel it is useful to have two
! arrays (one for receiving data and one for sending data when we need to
! communicate).  To avoid copying, we use pointers.
! spawned_walkers points at the current data,
! spawned_walkers_recvd is only used in data communication (see
! distribute_walkers in the annihilation module).
integer(i0), pointer :: spawned_walkers(:,:), spawned_walkers_recvd(:,:)
! d) current (filled) slot in the spawning arrays.
! In parallel we divide the spawning lists into blocks (one for each processor).
! spawning_head(j,i) gives the current filled slot in the spawning arrays for the
! block associated with the j-th thread on the i-th processor.
! After compress_spawned is called, spawning_head(0,i) gives the current filled
! slot on the i-th processor (all other elements are not meaningful).
! After distribute_walkers is called in the annihilation algorithm,
! spawning_head(0,0) is the number of spawned_walkers on the *current* processor
! and all other elements are not meaningful.
! It is convenient if the minimum size of spawning_head and spawning_block_start
! are both 0:1 along the processor dimension.
integer, allocatable :: spawning_head(:,:) ! (0:nthreads, max_0:(max(1,nprocs-1))
! spawning_block_start(i) contains the first position to be used in the spawning
! lists for storing a walker which is to be sent to the i-th processor.
integer, allocatable :: spawning_block_start(:,:) ! (0:nthreads,0:max(1,nprocs-1))

! Rate of spawning.  This is a running total over MC cycles on each processor
! until it is summed over processors and averaged over cycles in
! update_energy_estimators.
real(p) :: rspawn

!--- Reference determinant ---

! Bit string of reference determinant.
integer(i0), allocatable :: f0(:)

! List of occupied orbitals in reference determinant.
integer, allocatable :: occ_list0(:)

! Bit string of reference determinant used to generate the Hilbert space.
! This is usually identical to f0, but not necessarily (e.g. if we're doing
! a spin-flip calculation).
integer(i0), allocatable :: hs_f0(:)

! hs_f0:hs_occ_list0 as f0:occ_list0.
integer, allocatable :: hs_occ_list0(:)

! Population of walkers on reference determinant.
! The initial value can be overridden by a restart file or input option.
! For DMQMC, this variable stores the initial number of psips to be
! randomly distributed along the diagonal elements of the density matrix.
real(p) :: D0_population = 10.0_p

! Store value population of reference over a Monte Carlo cycle rather than
! accumulating it as we go.  This makes it easier to have multiple spawning
! events from the same determinant (as required in CCMC).
! D0_population is accumulated when updating the energy estimators.
real(p) :: D0_population_cycle

! Also start with D0_population on i_s|D_0>, where i_s is the spin-version
! operator.  This is only done if no restart file is used *and* |D_0> is not
! a closed shell determinant.
logical :: init_spin_inv_D0 = .false.

! When performing dmqmc calculations, dmqmc_factor = 2.0. This factor is
! required because in DMQMC calculations, instead of spawning from one end with
! the full probability, we spawn from two different ends with half probability each.
! Hence, tau is set to tau/2 in DMQMC calculations, so that an extra factor is not
! required in every spawning routine. In the death step however, we use
! walker_energies(1,idet), which has a factor of 1/2 included for convenience
! already, for conveniece elsewhere. Hence we have to multiply by an extra factor
! of 2 to account for the extra 1/2 in tau. dmqmc_factor is set to 1.0 when not
! performing a DMQMC calculation, and so can be ignored in these cases.
real(p) :: dmqmc_factor = 1.0_p

! This variable stores the number of estimators which are to be
! calculated and printed out in a DMQMC calculation.
integer :: number_dmqmc_estimators = 0

! In DMQMC the trace of the density matrix is an important quantity
! used in calculating all thermal estimators. This quantity stores
! the this value, Tr(\rho), where rho is the density matrix which
! the DMQMC algorithm calculates stochastically.
integer(i0) :: trace

! estimator_numerators stores all the numerators for the estimators in DMQMC
! which the user has asked to be calculated. These are, for a general
! operator O which we wish to find the thermal average of:
! \sum_{i,j} \rho_{ij} * O_{ji}
! This variabe will store this value from the first iteration of each
! report loop. At the end of a report loop, the values from each
! processor are combined and stored in estimator_numerators on the parent
! processor. This is then output, and the values of estimator_numerators
! are reset on each processor to start the next report loop.
real(p), allocatable :: estimator_numerators(:) !(number_dmqmc_estimators)
! The integers below store the index of the element in the array
! estimator_numerators, in which the corresponding thermal quantity is
! stored. In general, a different number and combination of estimators
! may be calculated, and hence we need a way of knowing which element
! refers to which operator. These indices give labels to operators.
! They will be set to 0 if no used, or else their positions in the array,
! from 1-number_dmqmc_estimators.
integer :: energy_index = 0
integer :: energy_squared_index = 0
integer :: correlation_index = 0
integer :: staggered_mag_index = 0

! If this logical is true then the program runs the DMQMC algorithm with
! importance sampling.
! dmqmc_sampling_prob stores the factors by which the probabilities of
! spawning to a larger excitation are reduced by. So, when spawning from
! a diagonal element to a element with one excitation, the probability
! of spawning is reduced by a factor dmqmc_sampling_probs(1).
! dmqmc_accumulated_probs(i) stores the multiplication of all the elements
! of dmqmc_sampling_probs up to the ith element. This quantity is often
! needed, so it is stored.
logical :: dmqmc_weighted_sampling
real(p), allocatable :: dmqmc_sampling_probs(:) ! (min(nel, nsites-nel))
real(p), allocatable :: dmqmc_accumulated_probs(:) ! (min(nel, nsites-nel) + 1)
! If dmqmc_vary_weights is true, then instead of using the final sampling
! weights for all the iterations, the weights will be gradually increased
! until finish_varying_weights, at which point they will be held constant.
! weight_altering_factors stores the factors by which each weight is
! multiplied at each step. 
logical :: dmqmc_vary_weights = .false.
integer :: finish_varying_weights = 0
real(dp), allocatable :: weight_altering_factors(:)
! If this logical is true then the program will calculate the ratios
! of the numbers of the psips on neighbouring excitation levels. These
! are output so that they can be used when doing importance sampling
! for DMQMC, so that each level will have roughly equal numbers of psips.
! The resulting new weights are used in the next beta loop.
logical :: dmqmc_find_weights

! If half_density_matrix is true then half the density matrix will be 
! calculated by reflecting spawning onto the lower triangle into the
! upper triangle. This is allowed because the density matrix is 
! symmetric.
logical :: half_density_matrix = .false.

! Calculate replicas (ie evolve two wavefunctions/density matrices at once)?
! Currently only implemented for DMQMC.
logical :: replica_tricks = .false.

! For DMQMC: If this locial is true then the fraction of psips at each
! excitation level will be output at each report loop. These fractions
! will be stored in the array below.
logical :: calculate_excit_distribution = .false.
real(p), allocatable :: excit_distribution(:) ! (min(nel, nsites-nel) + 1)

! If true, then the reduced density matrix will be calulated
! for the subsystem A specified by the user.
logical :: doing_reduced_dm = .false.

! If true then calculate the concurrence for reduced density matrix of two sites
logical :: doing_concurrence = .false.

! If true then calculate the Von-Neumann entanglement entropy for specified subsystem
logical :: doing_von_neumann_entropy = .false.

! If true then, if doing an exact diagonalisation, calculate and output the
! eigenvalues of the reduced density matrix requested.
logical :: doing_exact_rdm_eigv

! This stores the reduces matrix, which is slowly accumulated over time
! (on each processor).
real(p), allocatable :: reduced_density_matrix(:,:)

! If true then the reduced density matrix is output to a file, 'reduced_dm'
! each beta loop.
logical :: output_rdm
! The unit of the file reduced_dm.
integer :: rdm_unit

! This will store the 4x4 flip spin matrix \sigma_y \otimes \sigma_y if
! concurrence is to be calculated
real(p), allocatable :: flip_spin_matrix(:,:)

! When calculating certain DMQMC properties, we only want to start
! averaging once the ground state is reached. The below integer is input
! by the user, and gives the iteration at which data should start being
! accumulated for the quantity. This is currently only used for the
! reduced density matrix and calculating importance sampling weights.
integer :: start_averaging = 0

! In DMQMC, the user may want want the shift as a function of beta to be
! the same for each beta loop. If average_shift_until is non-zero then
! shift_profile is allocated, and for the first average_shift_until
! beta loops, the shift is stored at each beta value and then averaged
! over all these beta loops afterwards. These shift values are stored
! in shift profile, and then used in future beta loops as the shift
! profile for each one.
! When average_shift_until is set equal to -1, all the averaging has
! finished, and the algorithm uses the averaged values stored in
! shift_profile.
integer :: average_shift_until = 0
real(p), allocatable :: shift_profile(:) ! (nreport)

! correlation_mask is a bit string with a 1 at positions i and j which
! are considered when finding the spin correlation function, C(r_{i,j}).
! All other bits are set to 0. i and j are chosen by the user initially.
integer(i0), allocatable :: correlation_mask(:)
! correlation_sites stores the site positions specified by the users
! initially (as orbital labels).
integer, allocatable :: correlation_sites(:)

! When using the Neel singlet trial wavefunction, it is convenient
! to store all possible amplitudes in the wavefunction, since
! there are relativley few of them and they are expensive to calculate
real(dp), allocatable :: neel_singlet_amp(:) ! (nsites/2) + 1

! Energy of reference determinant.
real(p) :: H00

! Processor on which the reference determinant is kept.
integer :: D0_proc

! How often do we change the reference determinant to the determinant with
! greatest population?
! Default: we don't.
integer :: select_ref_det_every_nreports = huge(1)
! Factor by which the population on a determinant must exceed the reference
! determinant's population in order to be accepted as the new reference
! determinant.
real(p) :: ref_det_factor = 1.50_p

!--- Simple FCIQMC ---

! Data used *only* in the simple_fciqmc algorithm.
! Not set in the optimised algorithm.

! Location of reference determinant in dets_list.
integer :: ref_det

!--- Calculation modes ---

! The shift is updated at the end of each report loop when vary_shift is true.
logical :: vary_shift = .false.

!--- Restart data ---

! Restart calculation from file.
logical :: restart = .false.

! Print out restart file.
logical :: dump_restart_file = .false.

! Restart data.
integer :: mc_cycles_done = 0

!--- Folded spectrum data ---

! The line about which you are folding i.e. eps in (H-eps)^2 - E_0
real(p) :: fold_line = 0

! The generation probabilities of a dual excitation type
real(p) :: P__=0.05, Po_=0.475, P_o=0.475

! The split generation normalisations
real(p) :: X__=0, Xo_=0, X_o=0

contains

    !--- Statistics. ---

    function spawning_rate(ndeath, nattempts) result(rate)

        use parallel, only: nprocs

        ! Calculate the rate of spawning on the current processor.
        ! In:
        !    ndeath: number of particles that were killed/cloned during the MC
        !    cycle.
        !    nattempts: The number of attempts to spawn made in order to
        !    generate the current population of walkers in the spawned arrays.
        !    Note that this is *not* the same as nparticles as nparticles is
        !    updated during the Monte Carlo cycle as particles die.
        !    It is, however, identical to the number of particles on the
        !    processor at the beginning of the Monte Carlo cycle (multiplied by
        !    2 for the timestep algorithm).

        real(p) :: rate
        integer, intent(in) :: ndeath
        integer(lint), intent(in) :: nattempts
        integer :: nspawn

        nspawn = sum(spawning_head(0,:nprocs-1) - spawning_block_start(0,:nprocs-1))
        ! The total spawning rate is
        !   (nspawn + ndeath) / nattempts
        ! In the timestep algorithm each particle has 2 attempts (one to spawn on a different
        ! determinant and one to clone/die).
        rate = real(nspawn+ndeath,p)/nattempts

    end function spawning_rate

    !--- Operations on the spawned lists. ---

    subroutine sort_spawned_lists()

        ! Sort spawned_walkers according to the determinant list using
        ! quicksort.

        ! Uses the sample code in Numerical Recipies as a base.

        use basis, only: basis_length
        use determinants

        ! Threshold.  When a sublist gets to this length, switch to using
        ! insertion sort to sort the sublist.
        integer, parameter :: switch_threshold = 7

        ! sort needs auxiliary storage of length 2*log_2(n).
        integer, parameter :: stack_max = 50

        integer :: pivot, lo, hi, i, j
        integer(i0) :: tmp_spawned(spawned_size)

        ! Stack.  This is the auxilliary memory required by quicksort.
        integer, save :: stack(2,stack_max), nstack

        nstack = 0
        lo = 1
        hi = spawning_head(0,0)
        do
            ! If the section/partition we are looking at is smaller than
            ! switch_threshold then perform an insertion sort.
            if (hi - lo < switch_threshold) then
                do j = lo + 1, hi
                    tmp_spawned = spawned_walkers(:,j)
                    do i = j - 1, 1, -1
                        if (tmp_spawned(1:basis_length) .detgt. spawned_walkers(1:basis_length,i)) exit
                        spawned_walkers(:,i+1) = spawned_walkers(:,i)
                    end do
                    spawned_walkers(:,i+1) = tmp_spawned
                end do

                if (nstack == 0) exit
                hi = stack(2,nstack)
                lo = stack(1,nstack)
                nstack = nstack - 1

            else
                ! Otherwise start partitioning with quicksort.

                ! Pick the pivot element to be the median of spawned_walkers(:,lo), spawned_walkers(:,hi)
                ! and spawned_walkers(:,(lo+hi)/2).
                ! This largely overcomes a major problem with quicksort, where it
                ! degrades if the pivot is always the smallest element.
                pivot = (lo + hi)/2
                call swap_spawned(spawned_walkers(:,pivot), spawned_walkers(:,lo + 1))
                if (spawned_walkers(1:basis_length,lo) .detgt. spawned_walkers(1:basis_length,hi)) then
                    call swap_spawned(spawned_walkers(:,lo), spawned_walkers(:,hi))
                end if
                if (spawned_walkers(1:basis_length,lo+1) .detgt. spawned_walkers(1:basis_length,hi)) then
                    call swap_spawned(spawned_walkers(:,lo+1), spawned_walkers(:,hi))
                end if
                if (spawned_walkers(1:basis_length,lo) .detgt. spawned_walkers(1:basis_length,lo+1)) then
                    call swap_spawned(spawned_walkers(:,lo), spawned_walkers(:,lo+1))
                end if

                i = lo + 1
                j = hi
                tmp_spawned = spawned_walkers(:,lo + 1) ! a is the pivot value
                do while (.true.)
                    ! Scan down list to find element > a.
                    i = i + 1
                    do while (tmp_spawned(1:basis_length) .detgt. spawned_walkers(1:basis_length,i))
                        i = i + 1
                    end do

                    ! Scan down list to find element < a.
                    j = j - 1
                    do while (spawned_walkers(1:basis_length,j) .detgt.  tmp_spawned(1:basis_length))
                        j = j - 1
                    end do

                    ! When the pointers crossed, partitioning is complete.
                    if (j < i) exit

                    ! Swap the elements, so that all elements < a end up
                    ! in lower indexed variables.
                    call swap_spawned(spawned_walkers(:,i), spawned_walkers(:,j))
                end do

                ! Insert partitioning element
                spawned_walkers(:,lo + 1) = spawned_walkers(:,j)
                spawned_walkers(:,j) = tmp_spawned

                ! Push the larger of the partitioned sections onto the stack
                ! of sections to look at later.
                ! --> need fewest stack elements.
                nstack = nstack + 1

                ! With a stack_max of 50, we can sort arrays of length
                ! 1125899906842624.  It is safe to say this will never be
                ! exceeded, and so this test can be skipped.
!                if (nstack > stack_max) call stop_all('sort_spawned_lists', "parameter stack_max too small")

                if (hi - i + 1 >= j - lo) then
                    stack(2,nstack) = hi
                    stack(1,nstack) = i
                    hi = j - 1
                else
                    stack(2,nstack) = j - 1
                    stack(1,nstack) = lo
                    lo = i
                end if

            end if
        end do

    contains

        subroutine swap_spawned(s1,s2)

            integer(i0), intent(inout) :: s1(spawned_size), s2(spawned_size)
            integer(i0) :: tmp(spawned_size)

            tmp = s1
            s1 = s2
            s2 = tmp

        end subroutine swap_spawned

    end subroutine sort_spawned_lists

    !--- Output procedures ---

    subroutine write_fciqmc_report_header()

        use calc, only: doing_calc, hfs_fciqmc_calc, dmqmc_calc, doing_dmqmc_calc, dmqmc_correlation
        use calc, only: dmqmc_energy, dmqmc_energy_squared, dmqmc_staggered_magnetisation

        integer :: i

        if (doing_calc(dmqmc_calc)) then
           write (6,'(1X,a12,3X,a13,8X,a5)', advance = 'no') &
           '# iterations','Instant shift','Trace'

            if (doing_dmqmc_calc(dmqmc_energy)) then
                write (6, '(2X,a19)', advance = 'no') '\sum\rho_{ij}H_{ji}'
            end if
            if (doing_dmqmc_calc(dmqmc_energy_squared)) then
                write (6, '(2X,a19)', advance = 'no') '\sum\rho_{ij}H2{ji}'
            end if
            if (doing_dmqmc_calc(dmqmc_correlation)) then
                write (6, '(2X,a19)', advance = 'no') '\sum\rho_{ij}C_{ji}'
            end if
            if (doing_dmqmc_calc(dmqmc_staggered_magnetisation)) then
                write (6, '(2X,a19)', advance = 'no') '\sum\rho_{ij}M2{ji}'
            end if
            if (calculate_excit_distribution) then
                do i = 0, ubound(excit_distribution,1)
                    write (6, '(4X,a13,1X,i3)', advance = 'no') 'Excit. level:', i
                end do
            end if

            write (6, '(2X,a11,2X,a7,2X,a4)') '# particles', 'R_spawn', 'time'

        else
            write (6,'(1X,a13,3(2X,a17))', advance='no') &
                     "# iterations ", "Shift            ", "\sum H_0j N_j    ", "N_0              "
            if (doing_calc(hfs_fciqmc_calc)) then
                write (6,'(4(2X,a17),3X,a11,2X,a10)', advance='no') &
                    "H.F. Shift       ","\sum O_0j N_j    ","\sum H_0j N'_j   ","N'_0              ", &
                    "# H psips","# HF psips"
            else
                write (6,'(3X,a11)', advance='no') "# H psips"
            end if
            write (6,'(2X,a8,2X,a4)') "R_spawn ",  "time"
        end if

    end subroutine write_fciqmc_report_header

    subroutine write_fciqmc_report(ireport, ntot_particles, elapsed_time, comment)

        ! Write the report line at the end of a report loop.
        ! In:
        !    ireport: index of the report loop.
        !    ntot_particles: total number of particles in main walker list.
        !    elapsed_time: time taken for the report loop.
        !    comment: if true, then prefix the line with a #.

        use calc, only: doing_calc, dmqmc_calc, hfs_fciqmc_calc
        use hfs_data, only: proj_hf_O_hpsip, proj_hf_H_hfpsip, D0_hf_population, hf_shift

        integer, intent(in) :: ireport
        integer(lint), intent(in) :: ntot_particles(:)
        real, intent(in) :: elapsed_time
        logical :: comment
        integer :: mc_cycles, i

        mc_cycles = ireport*ncycles

        if (comment) then
            write (6,'(1X,"#",3X)', advance='no')
        else
            write (6,'(5X)', advance='no')
        end if

        ! See also the format used in inital_fciqmc_status if this is changed.
        if (doing_calc(dmqmc_calc)) then
            write (6,'(i8,2X,es17.10,i10)',advance = 'no') &
                                             (mc_cycles_done+mc_cycles-ncycles), shift, trace
            ! Perform a loop which outputs the numerators for each of the different
            ! estimators, as stored in total_estimator_numerators.
            do i = 1, number_dmqmc_estimators
                write (6, '(4X,es17.10)', advance = 'no') estimator_numerators(i)
            end do
            if (calculate_excit_distribution) then
                excit_distribution = excit_distribution/ntot_particles
                do i = 0, ubound(excit_distribution,1)
                    write (6, '(4X,es17.10)', advance = 'no') excit_distribution(i)
                end do
            end if
            write (6, '(2X, i11,3X,f6.4,2X,f4.2)') ntot_particles, rspawn, elapsed_time/ncycles
        else if (doing_calc(hfs_fciqmc_calc)) then
            write (6,'(i8,2X,6(es17.10,2X),es17.10,4X,i11,X,i11)', advance = 'no') &
                                             mc_cycles_done+mc_cycles, shift,   &
                                             proj_energy, D0_population, &
                                             hf_shift, proj_hf_O_hpsip, proj_hf_H_hfpsip, &
                                             D0_hf_population, &
                                             ntot_particles
        else
            write (6,'(i8,2X,2(es17.10,2X),es17.10,4X,i11)', advance='no') &
                                             mc_cycles_done+mc_cycles, shift,   &
                                             proj_energy, D0_population, &
                                             ntot_particles
        end if
        write (6,'(2X,f7.4,2X,f6.3)') rspawn, elapsed_time/ncycles

    end subroutine write_fciqmc_report

    subroutine write_fciqmc_final(ireport)

        ! Write out the energies (shift and projected energy) at the end of an
        ! FCIQMC calculation.
        ! In:
        !    ireport: index of the report loop after the report loop has been
        !    exited.

        integer, intent(in) :: ireport
        integer :: report_cycles_done

        if (ireport /= nreport+1) then
            ! exited calculation early via softexit.
            ! number of report loops done that actually were done is ireport.
            report_cycles_done = ireport
        else
            ! terminated the report loop cycle after reaching the last index.
            report_cycles_done = nreport
        end if

        write (6,'(/,1X,a13,10X,f22.12)') 'final shift =', shift
        write (6,'(1X,a20,3X,f22.12)') 'final proj. energy =', proj_energy/D0_population
        write (6,'(1X,a12,11X,f22.12)') 'E0 + shift =', shift+H00
        write (6,'(1X,a19,4X,f22.12)') 'E0 + proj. energy =', proj_energy/D0_population+H00

    end subroutine write_fciqmc_final

    subroutine end_fciqmc()

        ! Deallocate fciqmc data arrays.

        use checking, only: check_deallocate

        integer :: ierr

        if (allocated(occ_list0)) then
            deallocate(occ_list0, stat=ierr)
            call check_deallocate('occ_list0',ierr)
        end if
        if (allocated(nparticles)) then
            deallocate(nparticles, stat=ierr)
            call check_deallocate('nparticles',ierr)
        end if
        if (allocated(walker_dets)) then
            deallocate(walker_dets, stat=ierr)
            call check_deallocate('walker_dets',ierr)
        end if
        if (allocated(walker_population)) then
            deallocate(walker_population, stat=ierr)
            call check_deallocate('walker_population',ierr)
        end if
        if (allocated(walker_data)) then
            deallocate(walker_data, stat=ierr)
            call check_deallocate('walker_data',ierr)
        end if
        if (allocated(spawned_walkers1)) then
            deallocate(spawned_walkers1, stat=ierr)
            call check_deallocate('spawned_walkers1',ierr)
        end if
        if (allocated(spawned_walkers2)) then
            deallocate(spawned_walkers2, stat=ierr)
            call check_deallocate('spawned_walkers2',ierr)
        end if
        if (allocated(spawning_head)) then
            deallocate(spawning_head, stat=ierr)
            call check_deallocate('spawning_head',ierr)
        end if
        if (allocated(spawning_block_start)) then
            deallocate(spawning_block_start, stat=ierr)
            call check_deallocate('spawning_block_start',ierr)
        end if
        if (allocated(f0)) then
            deallocate(f0, stat=ierr)
            call check_deallocate('f0',ierr)
        end if
        if (allocated(neel_singlet_amp)) then
            deallocate(neel_singlet_amp, stat=ierr)
            call check_deallocate('neel_singlet_amp',ierr)
        end if
        if (allocated(estimator_numerators)) then
            deallocate(estimator_numerators, stat=ierr)
            call check_deallocate('estimator_numerators', ierr)
        end if

    end subroutine end_fciqmc

end module fciqmc_data

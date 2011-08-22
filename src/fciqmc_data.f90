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

! timestep
real(p) :: tau

! Array sizes
! If these are < 0, then the values represent the number of MB to be used to
! store the main walker and spawned walker data respectively.
integer :: walker_length
integer :: spawned_walker_length

! Number of particles before which varyshift mode is turned on.
integer :: target_particles = 10000

!--- Input data: initiator-FCIQMC ---

integer :: CAS(2) = (/ 0,0 /)

integer :: initiator_population = 3

!--- Energy data ---

! shift
real(p) :: shift = 0.0_p

! shift averaged over the calculation, once varyshift mode has been entered.
! This is really a running total: the average is only taken at output time (in
! write_fciqmc_report).
real(p) :: av_shift = 0.0_p

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
real(p) :: proj_energy

! This is \sum_{i/=0} <D_i|M_z^2|D_i> N_i^2 which is numerator of the expectation
! value of the staggered magnetisation squared (this is normalised by dividing
! through by population_squared).
real(p) :: average_magnetisation = 0.0_p

! projected energy averaged over the calculation.
! This is really a running total of \sum_{i/=0} <D_0|H|D_i> N_i: the average is
! only taken at output time (in write_fciqmc_report) by considering the ratio
!   av_proj_energy/av_D0_population.
real(p) :: av_proj_energy = 0.0_p
! Similarly a running total for the population on the reference determinant.
real(p) :: av_D0_population = 0.0_p

! Report loop at which the averages are set to 0.
integer :: start_averaging_from = 0

!--- Walker data ---

! Current number of walkers stored in the main list (processor dependent).
! This is updated during annihilation and merging of the spawned walkers into
! the main list.
integer :: tot_walkers

! Total number of particles on all walkers/determinants (processor dependent)
! Updated during death and annihilation and merging.
! The first element is the number of normal (Hamiltonian) particles.
! Subsequent elements are the number of Hellmann--Feynamnn particles.
integer, allocatable :: nparticles(:) ! (sampling_size)

! Walker information: main list.
! sampling_size is one for each quantity sampled (i.e. 1 for standard
! FCIQMC/initiator-FCIQMC, 2 for FCIQMC+Hellmann--Feynman sampling).
integer :: sampling_size
! a) determinants
integer(i0), allocatable :: walker_dets(:,:) ! (basis_length, walker_length)
! b) walker population
integer, allocatable :: walker_population(:,:) ! (sampling_size,walker_length)
! c) Diagonal matrix elements, K_ii.  Storing them avoids recalculation.
! K_ii = < D_i | H | D_i > - E_0, where E_0 = <D_0 | H | D_0> and |D_0> is the
! reference determinant.
real(p), allocatable :: walker_energies(:,:) ! (sampling_size,walker_length)
! When calculating the projected energy with various trial wavefunctions, it
! is useful to store quantites which are expensive to calculate and which are
! instead of recalculating them. For the Neel singlet state, the first component
! gives stores the total number of spins up on the first sublattice. The second
! component gives the number of 0-1 bonds where the 1 is on the first sublattice.
integer, allocatable :: walker_reference_data(:,:) ! (2,walker_length)

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
! spawning_head(i) gives the current filled slot in the spawning arrays for the
! block associated with the i-th processor.
! After distribute_walkers is called in the annihilation algorithm,
! spawning_head(0) is the number of spawned_walkers on the *current* processor
! and all other elements are not meaningful.
! It is convenient if the minimum size of spawning_head and spawning_block_start
! are both 0:1.
integer, allocatable :: spawning_head(:) ! (0:(max(1,nprocs-1))
! spawning_block_start(i) contains the first position to be used in the spawning
! lists for storing a walker which is to be sent to the i-th processor.
integer, allocatable :: spawning_block_start(:) ! (0:max(1,nprocs-1))

! Rate of spawning.  This is a running total over MC cycles on each processor
! until it is summed over processors and averaged over cycles in
! update_energy_estimators.
real(p) :: rspawn

!--- Reference determinant ---

! Bit string of reference determinant.
integer(i0), allocatable :: f0(:)

! List of occupied orbitals in reference determinant.
integer, allocatable :: occ_list0(:)

! Population of walkers on reference determinant.
! The initial value can be overridden by a restart file or input option.
real(p) :: D0_population = 10.0_p

! The modulus squared of the wavefunction which the psips represent
! This is used in calculating the expectation value of the
! staggered magnetisation.
real(p) :: population_squared = 0.0_p

! When using the Neel singlet trial wavefunction, it is convenient
! to store all possible amplitudes in the wavefunction, since
! there are relativley few of them and they are expensive to calculate
real(dp), allocatable :: neel_singlet_amp(:) ! (nsites/2) + 1

! Energy of reference determinant.
real(p) :: H00

! Processor on which the reference determinant is kept.
integer :: D0_proc

!--- Simple FCIQMC ---

! Data used *only* in the simple_fciqmc algorithm.
! Not set in the optimised algorithm.

! Location of reference determinant in dets_list.
integer :: ref_det

!--- Calculation modes ---

! The shift is updated at the end of each report loop when vary_shift is true.
logical :: vary_shift = .false.
! The number of report loops after which vary_shift mode was entered.
integer :: start_vary_shift
! True if the staggered magnetisation is to be calculated in the Heisenberg model
logical :: calculate_magnetisation = .false.

!--- Restart data ---

! Restart calculation from file.
logical :: restart = .false.

! Print out restart file.
logical :: dump_restart_file = .false.

! Restart data.
integer :: mc_cycles_done = 0, nparticles_old_restart = 0

contains

    !--- Initialisation. ---

    subroutine set_reference_det()

        ! Set the list of occupied orbitals in the reference determinant to be
        ! the spin-orbitals with the lowest kinetic energy which satisfy the
        ! spin polarisation.

        ! Note: this is for testing only!  The symmetry input is currently
        ! ignored.

        ! This should be used as a last resort if the user doesn't specify
        ! a reference determinant.

        use checking, only: check_allocate
        use errors, only: stop_all
        use system, only: nalpha, nbeta, nel, system_type, hub_k, hub_real, nsites, &
                          heisenberg, J_coupling
        use basis, only: bit_lookup
        use hubbard_real, only: connected_orbs
        
        integer :: i, j, ierr, spins_set, connections
        integer :: bit_element, bit_pos

        ! Leave the reference determinant unchanged if it's already been
        ! allocated (and presumably set).
        
        if (allocated(occ_list0)) then
            if (size(occ_list0) /= nel) then
                select case(system_type)
                case(heisenberg)
                    call stop_all('set_reference_det', &
                        'Reference determinant supplied does not contain the &
                        &specified number of up electrons.')
                case default
                    call stop_all('set_reference_det', &
                        'Reference determinant supplied does not contain the &
                        &specified number of electrons.')
                end select
            end if
        else
            allocate(occ_list0(nel), stat=ierr)
            call check_allocate('occ_list0',nel,ierr)
            select case(system_type)
            case(hub_k)
                ! Occupy the Fermi sphere.
                forall (i=1:nalpha) occ_list0(i) = 2*i-1
                forall (i=1:nbeta) occ_list0(i+nalpha) = 2*i
            case(hub_real)
                ! Attempt to keep electrons on different sites where possible.
                ! Sites 1, 3, 5, ... (occupy every other alpha orbital first, ie
                ! place a max of nsites/2 electrons.  (nsites+1)/2 accounts for
                ! the possibility that we have an odd number of sites.)
                forall (i=1:min(nalpha,(nsites+1)/2)) occ_list0(i) = 4*i-3
                ! now occupy the alternate alpha orbitals
                forall (i=1:nalpha-min(nalpha,(nsites+1)/2)) &
                    occ_list0(i+min(nalpha,(nsites+1)/2)) = 4*i-1
                ! Similarly for beta, but now occupying orbitals sites 2, 4,
                ! ..., preferentially.
                forall (i=1:min(nbeta,nsites/2)) occ_list0(i+nalpha) = 4*i
                forall (i=1:nbeta-min(nbeta,nsites/2)) &
                    occ_list0(i+nalpha+min(nbeta,nsites/2)) = 4*i-2
            case(heisenberg)
                if (J_coupling >= 0) then
                    forall (i=1:nel) occ_list0(i) = i
                ! For the antiferromagnetic case, below. This is messy but should 
                ! give a reasonable reference determinant for general cases, even
                ! for bizarre lattices. For bipartite lattices (eg 4x4, 6x6...)
                ! it will give the best possible reference determinant.
                else if (J_coupling < 0) then
                    ! Always set the first spin up
                    occ_list0(1) = 1
                    spins_set = 1
                    ! Loop over other sites to find orbitals which are not connected to
                    ! the other sites previously chosen.
                    do i=2,nsites
                        bit_pos = bit_lookup(1,i)
                        bit_element = bit_lookup(2,i)
                        connections = 0
                        ! Loop over all chosen sites to see if they neighbour this site.
                        do j=1,spins_set
                            if (btest(connected_orbs(bit_element, occ_list0(j)), bit_pos)) then
                                  connections = connections + 1
                            end if
                        end do
                        ! If this site has no neighbours which have been previously added
                        ! to the reference determinant, then we include it.
                        if (connections == 0) then
                            spins_set = spins_set + 1
                            occ_list0(spins_set) = i
                        end if
                    end do
                    ! If, after finding all the sites which are not connected, we still haven't
                    ! chosen enough sites, we accept that we must have some neigbouring sites
                    ! included in the reference determinant and start choosing the remaining sites.
                    if (spins_set /= nel) then
                        ! Loop over all sites looking for extra spins to include in the
                        ! reference detereminant.
                        fill_sites: do i=2,nsites
                            connections = 0
                            ! Check if this site is already included.
                            do j=1,spins_set
                                if (occ_list0(j) == i) connections = connections + 1
                            end do
                            ! If connection = 0, this site is not currently included in the
                            ! reference determinant, so add it.
                            if (connections == 0) then
                                spins_set = spins_set + 1
                                occ_list0(spins_set) = i
                            end if
                            ! When the correct number of spins have been chosen to be up,
                            ! we are finished.
                            if (spins_set == nel) exit fill_sites
                        end do fill_sites
                    end if
                end if
            end select
        end if

    end subroutine set_reference_det

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
        !    processor at the beginning of the Monte Carlo cycle (miltiplied by
        !    2 for the timestep algorithm).

        real(p) :: rate
        integer, intent(in) :: ndeath, nattempts
        integer :: nspawn

        nspawn = sum(spawning_head(:nprocs-1) - spawning_block_start(:nprocs-1))
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
        hi = spawning_head(0)
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

    pure subroutine search_walker_list(f, istart, iend, hit, pos)

        ! Find where a determinant belongs in the main walker list.
        ! Only elements between istart and iend are examined (use the 
        ! array boundaries in the worst case).
        !
        ! In:
        !    f: bit string representation of the Slater determinant.
        !    istart: first position to examine in the walker list.
        !    iend: last position to examine in the walker list.
        ! Out:
        !    hit: true if found f in the main walker list.
        !    pos : the corresponding position in the main walker list
        !        where the determinant belongs.  If hit is true, then
        !        the determinant in this position is the same as f, else
        !        this is where f should go to keep the main walker list sorted.

        use basis, only: basis_length
        use determinants, only: det_compare

        integer(i0), intent(in) :: f(basis_length)
        integer, intent(in) :: istart, iend
        logical, intent(out) :: hit
        integer, intent(out) :: pos

        integer :: hi, lo, compare

        if (istart > iend) then

            ! Already know the element has to be appended to the list.
            ! This should only occur if istart = iend + 1.
            pos = istart
            hit = .false.
            
        else

            ! Search range.
            lo = istart
            hi = iend

            ! Assume f doesn't exist in the main walkers list initially.
            hit = .false.

            do while (hi /= lo)
                ! Narrow the search range down in steps.

                ! Mid-point.
                ! We shift one of the search limits to be the mid-point.
                ! The successive dividing the search range by 2 gives a O[log N]
                ! search algorithm.
                pos = (hi+lo)/2

                compare = det_compare(walker_dets(:,pos), f)
                select case(compare)
                case (0)
                    ! hit!
                    hit = .true.
                    exit
                case(1)
                    ! walker_dets(:,pos) is "smaller" than f.
                    ! The lowest position f can take is hence pos + 1 (i.e. if
                    ! f is greater than pos by smaller than pos + 1).
                    lo = pos + 1
                case(-1)
                    ! walker_dets(:,pos) is "greater" than f.
                    ! The highest position f can take is hence pos (i.e. if f is
                    ! smaller than pos but greater than pos - 1).  This is why
                    ! we differ slightly from a standard binary search (where lo
                    ! is set to be pos+1 and hi to be pos-1 accordingly), as
                    ! a standard binary search assumes that the element you are
                    ! searching for actually appears in the array being
                    ! searched...
                    hi = pos
                end select

            end do

            ! If hi == lo, then we have narrowed the search down to one position but
            ! not checked if that position is the item we're hunting for.
            ! Because walker_dets can expand (i.e. we might be searching for an
            ! element which doesn't exist yet) the binary search can find either
            ! the element before or after where f should be placed.
            if (hi == lo) then
                compare = det_compare(walker_dets(:,hi), f)
                select case(compare)
                case (0)
                    ! hit!
                    hit = .true.
                    pos = hi
                case(1)
                    ! walker_dets(:,pos) is "smaller" than f.
                    ! f should be placed in the next slot.
                    pos = hi + 1
                case(-1)
                    ! walker_dets(:,pos) is "greater" than f.
                    ! f should ber placed here.
                    pos = hi
                end select
            end if

        end if

    end subroutine search_walker_list

    !--- Output procedures ---

    subroutine write_fciqmc_report_header()
        
        if (calculate_magnetisation) then
            write (6,'(1X,a12,6X,a13,6X,a9,9X,a12,7X,a11,12X,a4,3X,a11,3X,a16,11X,a9,3X,a7,3X,a4)') &
           '# iterations','Instant shift','Av. shift','\sum H_0j Nj', &
           'Av. Proj. E','# D0','# particles','\sum M_ii^2 Ni^2','\sum Ni^2','R_spawn','time'
       else
           write (6,'(1X,a12,3X,a13,6X,a9,10X,a12,7X,a11,8X,a4,16X,a11,2X,a7,2X,a4)') &
           '# iterations','Instant shift','Av. shift','\sum H_0j Nj',    &
           'Av. Proj. E','# D0','# particles','R_spawn','time'
       end if

    end subroutine write_fciqmc_report_header

    subroutine write_fciqmc_report(ireport, ntot_particles, elapsed_time)

        ! Write the report line at the end of a report loop.
        ! In:
        !    ireport: index of the report loop.
        !    ntot_particles: total number of particles in main walker list.
        !    elapsed_time: time taken for the report loop.

        integer, intent(in) :: ireport, ntot_particles
        real, intent(in) :: elapsed_time
        integer :: mc_cycles, vary_shift_reports

        mc_cycles = ireport*ncycles

        vary_shift_reports = ireport - start_averaging_from - start_vary_shift

        ! See also the format used in inital_fciqmc_status if this is changed.
        if (calculate_magnetisation) then
            write (6,'(5X,i8,4(2X,es17.10),2X,f11.4,5X,i9,2X,es17.10,3X,es17.10,4X,f6.4,3X,f4.2)') &
                                             mc_cycles_done+mc_cycles, shift, &
                                             av_shift/vary_shift_reports, proj_energy, &
                                             av_proj_energy/av_D0_population, D0_population, &
                                             ntot_particles,average_magnetisation, &
                                             population_squared, rspawn, elapsed_time/ncycles
        else if (.not.calculate_magnetisation) then                                    
            write (6,'(5X,i8,2X,4(es17.10,2X),es17.10,4X,i11,3X,f6.4,2X,f4.2)') &
                                             mc_cycles_done+mc_cycles, shift,   &
                                             av_shift/vary_shift_reports, proj_energy,       &
                                             av_proj_energy/av_D0_population, D0_population, & 
                                             ntot_particles, rspawn, elapsed_time/ncycles
        end if

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

        av_shift = av_shift/(report_cycles_done - start_vary_shift - start_averaging_from)

        write (6,'(/,1X,a13,10X,f22.12)') 'final shift =', shift
        write (6,'(1X,a20,3X,f22.12)') 'final proj. energy =', proj_energy/D0_population
        write (6,'(1X,a11,12X,f22.12)') 'av. shift =', av_shift
        write (6,'(1X,a18,5X,f22.12)') 'av. proj. energy =', av_proj_energy/av_D0_population
        write (6,'(1X,a12,11X,f22.12)') 'E0 + shift =', shift+H00
        write (6,'(1X,a19,4X,f22.12)') 'E0 + proj. energy =', proj_energy/D0_population+H00
        write (6,'(1X,a16,7X,f22.12)') 'E0 + av. shift =', av_shift+H00
        write (6,'(1X,a23,f22.12)') 'E0 + av. proj. energy =', av_proj_energy/av_D0_population+H00

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
        if (allocated(walker_energies)) then
            deallocate(walker_energies, stat=ierr)
            call check_deallocate('walker_energies',ierr)
        end if
        if (allocated(walker_reference_data)) then
            deallocate(walker_reference_data, stat=ierr)
            call check_deallocate('walker_reference_data',ierr)
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

    end subroutine end_fciqmc

end module fciqmc_data

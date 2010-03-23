module fciqmc_data

! Data for fciqmc calculations and procedures which manipulate fciqmc and only
! fciqmc data.

use const
implicit none

!--- Input data ---

! number of monte carlo cycles/report loop
integer :: ncycles
! number of report cycles
! the shift is updated and the calculation information printed out
! at the end of each report cycle.
integer :: nreport

! timestep
real(p) :: tau

! Array sizes
integer :: walker_length
integer :: spawned_walker_length

! Current number of walkers stored in the main list.
integer :: tot_walkers

! Number of particles before which varyshift mode is turned on.
integer :: target_particles = 10000

!--- Energy data ---

! shift
real(p) :: shift = 0.0_p

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

!--- Walker data ---

! Walker information: main list.
! a) determinants
integer(i0), allocatable :: walker_dets(:,:) ! (basis_length, walker_length)
! b) walker population
integer, allocatable :: walker_population(:) ! (walker_length)
! c) Diagonal matrix elements, K_ii.  Storing them avoids recalculation.
! K_ii = < D_i | H | D_i > - E_0, where E_0 = <D_0 | H | D_0> and |D_0> is the
! reference determinant.
real(p), allocatable :: walker_energies(:)

! Walker information: spawned list.
! a) determinants.
integer(i0), allocatable :: spawned_walker_dets(:,:) ! (basis_length, spawned_walker_length)
! b) walker population.
integer, allocatable :: spawned_walker_population(:) ! (spawned_walker_length)
! c) current (filled) slot in the spawning arrays.
!    if 0, then no elements are in the spawning arrays.
integer :: spawning_head

!--- Reference determinant ---

! Bit string of reference determinant.
integer(i0), allocatable :: f0(:)

! List of occupied orbitals in reference determinant.
integer, allocatable :: occ_list0(:)

! Population of walkers on reference determinant.
integer :: D0_population

! Energy of reference determinant.
real(p) :: H00

!--- Simple FCIQMC ---

! Data used *only* in the simple_fciqmc algorithm.
! Not set in the optimised algorithm.

! Location of reference determinant in dets_list.
integer :: ref_det

!--- Calculation modes ---

! The shift is updated at the end of each report loop when vary_shift is true.
logical :: vary_shift = .false.

contains

    subroutine sort_spawned_lists()

        ! Sort spawned_walker_dets and spawned_walker_populations according to
        ! the determinant list using quicksort.

        ! Uses the sample code in Numerical Recipies as a base.

        use basis, only: basis_length
        use determinants

        ! Threshold.  When a sublist gets to this length, switch to using
        ! insertion sort to sort the sublist.
        integer, parameter :: switch_threshold = 7

        ! sort needs auxiliary storage of length 2*log_2(n).
        integer, parameter :: stack_max = 50

        integer :: pivot, lo, hi, i, j, tmp_pop
        integer :: tmp_det(basis_length)

        ! Stack.  This is the auxilliary memory required by quicksort.
        integer :: stack(2,stack_max), nstack

        write (6,*) 'bsort',spawning_head

        nstack = 0
        lo = 1
        hi = spawning_head
        do
            ! If the section/partition we are looking at is smaller than
            ! switch_threshold then perform an insertion sort.
            if (hi - lo < switch_threshold) then
                do j = lo + 1, hi
                    tmp_det = spawned_walker_dets(:,j)
                    tmp_pop = spawned_walker_population(j)
                    do i = j - 1, 1, -1
                        if (tmp_det .detgt. spawned_walker_dets(:,i)) exit
                        spawned_walker_dets(:,i+1) = spawned_walker_dets(:,i)
                        spawned_walker_population(i+1) = spawned_walker_population(i)
                    end do
                    spawned_walker_dets(:,i+1) = tmp_det
                    spawned_walker_population(i+1) = tmp_pop
                end do

                if (nstack == 0) exit
                hi = stack(2,nstack)
                lo = stack(1,nstack)
                nstack = nstack - 1

            else
                ! Otherwise start partitioning with quicksort.

                ! Pick the pivot element to be the median of spawned_walker_dets(:,lo), spawned_walker_dets(:,hi)
                ! and spawned_walker_dets(:,(lo+hi)/2).
                ! This largely overcomes a major problem with quicksort, where it
                ! degrades if the pivot is always the smallest element.
                pivot = (lo + hi)/2
                call swap_dets(spawned_walker_dets(:,pivot), spawned_walker_dets(:,lo + 1))
                call swap_int(spawned_walker_population(pivot), spawned_walker_population(lo + 1))
                if (spawned_walker_dets(:,lo) .detgt. spawned_walker_dets(:,hi)) then
                    call swap_dets(spawned_walker_dets(:,lo), spawned_walker_dets(:,hi))
                    call swap_int(spawned_walker_population(lo), spawned_walker_population(hi))
                end if
                if (spawned_walker_dets(:,lo+1) .detgt. spawned_walker_dets(:,hi)) then
                    call swap_dets(spawned_walker_dets(:,lo+1), spawned_walker_dets(:,hi))
                    call swap_int(spawned_walker_population(lo+1), spawned_walker_population(hi))
                end if
                if (spawned_walker_dets(:,lo) .detgt. spawned_walker_dets(:,lo+1)) then
                    call swap_dets(spawned_walker_dets(:,lo), spawned_walker_dets(:,lo+1))
                    call swap_int(spawned_walker_population(lo), spawned_walker_population(lo+1))
                end if

                i = lo + 1
                j = hi
                tmp_det = spawned_walker_dets(:,lo + 1) ! a is the pivot value
                do while (.true.)
                    ! Scan down list to find element > a.
                    i = i + 1
                    do while (tmp_det .detgt. spawned_walker_dets(:,i))
                        i = i + 1
                    end do

                    ! Scan down list to find element < a.
                    j = j - 1
                    do while (spawned_walker_dets(:,j) .detgt. tmp_det)
                        j = j - 1
                    end do

                    ! When the pointers crossed, partitioning is complete.
                    if (j < i) exit

                    ! Swap the elements, so that all elements < a end up
                    ! in lower indexed variables.
                    call swap_dets(spawned_walker_dets(:,i), spawned_walker_dets(:,j))
                    call swap_int(spawned_walker_population(i), spawned_walker_population(j))
                end do

                ! Insert partitioning element
                spawned_walker_dets(:,lo + 1) = spawned_walker_dets(:,j)
                spawned_walker_dets(:,j) = tmp_det
                call swap_int(spawned_walker_population(j), spawned_walker_population(lo + 1))

                ! Push the larger of the partitioned sections onto the stack
                ! of sections to look at later.
                ! --> need fewest stack elements.
                nstack = nstack + 1

                ! With a stack_max of 50, we can sort arrays of length 
                ! 1125899906842624.  It is safe to say this will never be
                ! exceeded, and so this test can be skipped.
!                if (nstack > stack_max) call stop_all('sort_spawned_lists', "parameter stack_max too small")

                if (hi - i + 1 >= j - 1) then
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

        ! verify
        tmp_det = spawned_walker_dets(:,1)
        do i = 2, spawning_head
            if (tmp_det .detgt. spawned_walker_dets(:,i)) then
                write (6,*) 'error sorting'
                stop
            end if
            tmp_det = spawned_walker_dets(:,i)
        end do

        write (6,*) 'asort',spawning_head

    contains

        subroutine swap_int(i1, i2)

            integer, intent(inout) :: i1, i2
            integer :: tmp

            tmp = i1
            i1 = i2
            i2 = tmp

        end subroutine swap_int

        subroutine swap_dets(f1,f2)

            integer(i0), intent(inout) :: f1(basis_length), f2(basis_length)
            integer(i0) :: tmp(basis_length)

            tmp = f1
            f1 = f2
            f2 = tmp

        end subroutine swap_dets

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
                    lo = pos + 1
                case(-1)
                    ! walker_dets(:,pos) is "greater" than f.
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

        write (6,'(1X,a12,7X,a13,10X,a12,2X,a11)') '# iterations','Instant shift','Proj. Energy','# particles'

    end subroutine write_fciqmc_report_header

    subroutine write_fciqmc_report(ireport, nparticles)

        ! Write the report line at the end of a report loop.
        ! In:
        !    ireport: index of the report loop.
        !    nparticles: total number of particles in main walker list.

        integer, intent(in) :: ireport, nparticles

        write (6,'(5X,i8,2(f20.10,2X),i11)') ireport*ncycles, shift, proj_energy, nparticles

    end subroutine write_fciqmc_report

    subroutine write_fciqmc_final()

        write (6,'(/,1X,a13,7X,f22.12)') 'final_shift =', shift
        write (6,'(1X,a20,f22.12)') 'final proj. energy =', proj_energy
        write (6,'(1X,a12,8X,f22.12)') 'E0 + shift =', shift+H00
        write (6,'(1X,a19,1X,f22.12)') 'E0 + proj. energy =', proj_energy+H00

    end subroutine write_fciqmc_final

    subroutine end_fciqmc()

        ! Deallocate fciqmc data arrays.

        integer :: ierr

        if (allocated(walker_dets)) deallocate(walker_dets, stat=ierr)
        if (allocated(walker_population)) deallocate(walker_population, stat=ierr)
        if (allocated(walker_energies)) deallocate(walker_energies, stat=ierr)
        if (allocated(spawned_walker_dets)) deallocate(spawned_walker_dets, stat=ierr)
        if (allocated(spawned_walker_population)) deallocate(spawned_walker_population, stat=ierr)
        if (allocated(f0)) deallocate(f0, stat=ierr)

    end subroutine end_fciqmc

end module fciqmc_data

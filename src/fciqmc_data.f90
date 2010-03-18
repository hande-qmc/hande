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
real(dp) :: tau

! Array sizes
integer :: walker_length
integer :: spawned_walker_length

! Current number of walkers stored in the main list.
integer :: tot_walkers

! Number of particles before which varyshift mode is turned on.
integer :: target_particles = 10000

!--- Energy data ---

! shift
real(dp) :: shift = 0.0_dp

! projected energy
! This stores during a FCIQMC cycle
!   \sum_{i/=0} <D_0|H|D_i> N_i
! where D_0 is the reference determinants and N_i is the walker population on
! determinant D_i.
! The projected energy is given as
!   <D_0|H|D_0> + \sum_{i/=0} <D_0|H|D_i> N_i/N_0
! and so proj_energy must be converted accordingly.
real(dp) :: proj_energy

!--- Walker data ---

! Walker information: main list.
! a) determinants
integer(i0), allocatable :: walker_dets(:,:) ! (basis_length, walker_length)
! b) walker population
integer, allocatable :: walker_population(:) ! (walker_length)
! c) Diagonal matrix elements, K_ii.  Storing them avoids recalculation.
! K_ii = < D_i | H | D_i > - E_0, where E_0 = <D_0 | H | D_0> and |D_0> is the
! reference determinant.
real(dp), allocatable :: walker_energies(:)

! Walker information: spawned list.
! a) determinants.
integer(i0), allocatable :: spawned_walker_dets(:,:) ! (basis_length, spawned_walker_length)
! b) walker population.
integer, allocatable :: spawned_walker_population(:) ! (spawned_walker_length)
! c) next empty slot in the spawning array.
integer :: spawning_head

!--- Reference determinant ---

! Bit string of reference determinant.
integer(i0), allocatable :: f0(:)

! List of occupied orbitals in reference determinant.
integer, allocatable :: occ_list0(:)

! Population of walkers on reference determinant.
integer :: D0_population

! Energy of reference determinant.
real(dp) :: H00

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
        ! the determinant list.

        ! This is a simple insertion sort so should be modified to quicksort
        ! for larger lists.

        use basis, only: basis_length
        use determinants, only: det_gt

        integer :: i, j
        integer(i0) :: tmp_det(basis_length)
        integer :: tmp_pop

        do i = 2, spawning_head
            j = i - 1
            tmp_det = spawned_walker_dets(:,i)
            tmp_pop = spawned_walker_population(i)
            do while (j >= 1 .and. det_gt(spawned_walker_dets(:,j),tmp_det))
                spawned_walker_dets(:,j+1) = spawned_walker_dets(:,j)
                spawned_walker_population(j+1) = spawned_walker_population(j)
                j = j - 1
            end do
            spawned_walker_dets(:,j+1) = tmp_det
            spawned_walker_population(j+1) = tmp_pop
        end do

    end subroutine sort_spawned_lists

    subroutine search_walker_list(f, istart, iend, hit, pos)

        ! Find where a determinant belongs in the main walker list.
        ! Only elements between istart and iend are examined (use the 
        ! array boundaries in the worst case).
        !
        ! This currently uses a linear search and should be changed to a binary
        ! search algorithm for efficiency.
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
        use determinants, only: det_gt

        integer(i0), intent(in) :: f(basis_length)
        integer, intent(in) :: istart, iend
        logical, intent(out) :: hit
        integer, intent(out) :: pos

        integer :: i

        ! Assume it should go at the end to start with...
        pos = iend + 1
        hit = .false.
        do i = istart, iend
            if (all(walker_dets(:,i) == f)) then
                hit = .true.
                pos = i
                exit
            else if (det_gt(walker_dets(:,i),f)) then
                hit = .false.
                pos = i
                exit
            end if
        end do

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

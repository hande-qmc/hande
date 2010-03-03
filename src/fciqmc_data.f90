module fciqmc_data

! Data for fciqmc calculations and procedures which manipulate fciqmc and only
! fciqmc data.

use const
implicit none

!--- Input data ---

! number of monte carlo cycles/report loop
integer :: ncycles
! number of report cycles
integer :: nreport

! timestep
real(dp) :: tau = 0.10_dp

! shift
real(dp) :: shift = 0.0_dp

! Array sizes
integer :: walker_length
integer :: spawned_walker_length

! Current number of walkers stored in the main list.
integer :: tot_walkers

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

! Energy of reference determinant.
real(dp) :: H00

contains

    subroutine sort_spawning_list()

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

    end subroutine sort_spawning_list

    subroutine search_walker_list(f, hit, pos)

        ! Find where a determinant belongs in the main walker list.
        !
        ! This currently uses a linear search and should be changed to a binary
        ! search algorithm for efficiency.
        !
        ! In:
        !    f: bit string representation of the Slater determinant.
        ! Out:
        !    hit: true if found f in the main walker list.
        !    pos : the corresponding position in the main walker list
        !        where the determinant belongs.  If hit is true, then
        !        the determinant in this position is the same as f, else
        !        this is where f should go to keep the main walker list sorted.

        use basis, only: basis_length
        use determinants, only: det_gt

        integer(i0), intent(in) :: f(basis_length)
        logical, intent(out) :: hit
        integer, intent(out) :: pos

        integer :: i

        do i = 1, tot_walkers
            if (all(f == walker_dets(:,i))) then
                hit = .true.
                pos = i
                exit
            else if (det_gt(f, walker_dets(:,i))) then
                hit = .false.
                pos = i
                exit
            end if
        end do

    end subroutine search_walker_list

end module fciqmc_data

module ranking

! Module of ranking procedures.

! "Half my kingdom for a sane implementation of generic programming in Fortran!"
!   Richard III, shortly after falling into a time vortex at the Battle
!   of Agincourt.

use const

implicit none

contains

    pure subroutine insertion_rank_rp(arr, rank)

        ! Rank a real(p) array in increasing order using the insertion sort
        ! algorithm.
        !
        ! Resultant ranking is *stable* and insertion sort is really pretty
        ! decent, especially when the array is small.  Naturally one should
        ! investigate quicksort and the like for large arrays, but insertion
        ! sort is a good compromise and easy to code.
        !
        ! Based upon the F90 insertion sort implementation in Rosetta Stone
        ! (TODO reference).
        !
        ! ***WARNING***
        ! We assume that the arrays are 1-indexed due to a feature of how array
        ! bounds are handled in procedure arguments in Fortran.
        !
        ! In:
        !   arr: array of real values.
        ! In/Out:
        !    rank: on output rank(i) contains the ranked index (in increasing
        !    order) of the value in arr(i), that is arr(rank(i)) returns the
        !    i-th element in the sorted list of values of arr.
        !    NOTE: rank must have at least the dimensions of arr on input.

        real(p), intent(in) :: arr(:)
        integer, intent(inout) :: rank(:) ! inout to avoid automatic deallocation 
                                          ! of an allocatable array on entry

        integer :: i, j, tmp

        forall (i=1:size(arr)) rank(i) = i

        do i = 2, size(arr)
            j = i - 1
            tmp = rank(i)
            do while ( j >= 1 .and. arr(rank(j)) > arr(tmp) )
                rank(j+1) = rank(j)
                j = j - 1
            end do
            rank(j+1) = tmp
        end do

    end subroutine insertion_rank_rp

end module ranking

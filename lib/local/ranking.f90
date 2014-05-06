module ranking

! Module of ranking procedures.

! "Half my kingdom for a sane implementation of generic programming in Fortran!"
!   Richard III, shortly after falling into a time vortex at the Battle
!   of Agincourt.

use const

implicit none

interface insertion_rank
    module procedure insertion_rank_rp
    module procedure insertion_rank_int
end interface

contains

    pure subroutine insertion_rank_rp(arr, rank, tolerance)

        ! Rank a real(p) array in increasing order using the insertion sort
        ! algorithm.
        !
        ! Resultant ranking is *stable* and insertion sort is really pretty
        ! decent, especially when the array is small.  Naturally one should
        ! investigate quicksort and the like for large arrays, but insertion
        ! sort is a good compromise and easy to code.
        !
        ! Based upon the F90 insertion sort implementation in Rosetta Stone
        ! http://rosettacode.org/wiki/Sorting_algorithms/Insertion_sort#Fortran.
        !
        ! ***WARNING***
        ! We assume that the arrays are 1-indexed due to a feature of how array
        ! bounds are handled in procedure arguments in Fortran.
        !
        ! In:
        !   arr: array of real values.
        !   tolerance (optional, default 0.0): tolerance for comparing values.
        !   If present, then two values which differ by less than tolerance are
        !   treated as equal.  This allows one to avoid changing the order of
        !   items that have identical values to within a desired precision.
        !   Naturally, one should have this as small as possible.
        ! In/Out:
        !    rank: on output rank(i) contains the ranked index (in increasing
        !    order) of the value in arr(i), that is arr(rank(i)) returns the
        !    i-th element in the sorted list of values of arr.
        !    NOTE: rank must have at least the dimensions of arr on input.

        real(p), intent(in) :: arr(:)
        real(p), intent(in), optional :: tolerance
        integer, intent(inout) :: rank(:) ! inout to avoid automatic deallocation
                                          ! of an allocatable array on entry

        integer :: i, j, tmp
        real(p) :: tol

        if (present(tolerance)) then
            tol = tolerance
        else
            tol = 0.0_p
        end if

        forall (i=1:size(arr)) rank(i) = i

        do i = 2, size(arr)
            j = i - 1
            tmp = rank(i)
            do while ( j >= 1 )
                if (arr(rank(j)) - arr(tmp) < tol) exit
                rank(j+1) = rank(j)
                j = j - 1
            end do
            rank(j+1) = tmp
        end do

    end subroutine insertion_rank_rp

    ! [review] - JSS: as in load balancing branch, tolerance is not needed here.
    pure subroutine insertion_rank_int(arr, rank, tolerance)

        ! Rank an int array in increasing order using the insertion sort
        ! algorithm.
        !
        ! Resultant ranking is *stable* and insertion sort is really pretty
        ! decent, especially when the array is small.  Naturally one should
        ! investigate quicksort and the like for large arrays, but insertion
        ! sort is a good compromise and easy to code.
        !
        ! Based upon the F90 insertion sort implementation in Rosetta Stone
        ! http://rosettacode.org/wiki/Sorting_algorithms/Insertion_sort#Fortran.
        !
        ! ***WARNING***
        ! We assume that the arrays are 1-indexed due to a feature of how array
        ! bounds are handled in procedure arguments in Fortran.
        !
        ! In:
        !   arr: array of real values.
        !   tolerance (optional, default 0.0): tolerance for comparing values.
        !   If present, then two values which differ by less than tolerance are
        !   treated as equal.  This allows one to avoid changing the order of
        !   items that have identical values to within a desired precision.
        !   Naturally, one should have this as small as possible.
        ! In/Out:
        !    rank: on output rank(i) contains the ranked index (in increasing
        !    order) of the value in arr(i), that is arr(rank(i)) returns the
        !    i-th element in the sorted list of values of arr.
        !    NOTE: rank must have at least the dimensions of arr on input.

        integer, intent(in) :: arr(:)
        real(p), intent(in), optional :: tolerance
        integer, intent(inout) :: rank(:) ! inout to avoid automatic deallocation
                                          ! of an allocatable array on entry

        integer :: i, j, tmp
        real(p) :: tol

        if (present(tolerance)) then
            tol = tolerance
        else
            tol = 0.0_p
        end if

        forall (i=1:size(arr)) rank(i) = i

        do i = 2, size(arr)
            j = i - 1
            tmp = rank(i)
            do while ( j >= 1 )
                if (arr(rank(j)) - arr(tmp) < tol) exit
                rank(j+1) = rank(j)
                j = j - 1
            end do
            rank(j+1) = tmp
        end do

    end subroutine insertion_rank_int

end module ranking

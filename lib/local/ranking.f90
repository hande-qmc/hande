module ranking

! Module of ranking procedures.

! "Half my kingdom for a sane implementation of generic programming in Fortran!"
!   Richard III, shortly after falling into a time vortex at the Battle
!   of Agincourt.

use const

implicit none

interface insertion_rank
    module procedure insertion_rank_sp
    module procedure insertion_rank_dp
    module procedure insertion_rank_int_32
    module procedure insertion_rank_int_64
end interface

contains

    pure subroutine insertion_rank_sp(arr, rank, tolerance)

        ! Rank a real(sp) array in increasing order using the insertion sort
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

        real(sp), intent(in) :: arr(:)
        real(sp), intent(in), optional :: tolerance
        integer, intent(inout) :: rank(:) ! inout to avoid automatic deallocation
                                          ! of an allocatable array on entry

        integer :: i, j, tmp
        real(sp) :: tol

        if (present(tolerance)) then
            tol = tolerance
        else
            tol = 0.0_sp
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

    end subroutine insertion_rank_sp

    pure subroutine insertion_rank_dp(arr, rank, tolerance)

        ! Rank a real(dp) array in increasing order using the insertion sort
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

        real(dp), intent(in) :: arr(:)
        real(dp), intent(in), optional :: tolerance
        integer, intent(inout) :: rank(:) ! inout to avoid automatic deallocation
                                          ! of an allocatable array on entry

        integer :: i, j, tmp
        real(dp) :: tol

        if (present(tolerance)) then
            tol = tolerance
        else
            tol = 0.0_dp
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

    end subroutine insertion_rank_dp

    pure subroutine insertion_rank_int_32(arr, rank)

        ! Rank an array of integers in increasing order using the insertion sort
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
        ! In/Out:
        !    rank: on output rank(i) contains the ranked index (in increasing
        !    order) of the value in arr(i), that is arr(rank(i)) returns the
        !    i-th element in the sorted list of values of arr.
        !    NOTE: rank must have at least the dimensions of arr on input.

        integer, intent(in) :: arr(:)
        integer, intent(inout) :: rank(:) ! inout to avoid automatic deallocation
                                          ! of an allocatable array on entry
        integer :: i, j, tmp

        forall (i=1:size(arr)) rank(i) = i

        do i = 2, size(arr)
            j = i - 1
            tmp = rank(i)
            do while ( j >= 1 )
                if (arr(rank(j)) - arr(tmp) < 0) exit
                rank(j+1) = rank(j)
                j = j - 1
            end do
            rank(j+1) = tmp
        end do


    end subroutine insertion_rank_int_32

    pure subroutine insertion_rank_int_64(arr, rank)

        ! Rank an array of long integers in increasing order using the insertion sort
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
        ! In/Out:
        !    rank: on output rank(i) contains the ranked index (in increasing
        !    order) of the value in arr(i), that is arr(rank(i)) returns the
        !    i-th element in the sorted list of values of arr.
        !    NOTE: rank must have at least the dimensions of arr on input.

        integer(int_64), intent(in) :: arr(:)
        integer, intent(inout) :: rank(:) ! inout to avoid automatic deallocation
                                          ! of an allocatable array on entry
        integer :: i, j, tmp

        forall (i=1:size(arr)) rank(i) = i

        do i = 2, size(arr)
            j = i - 1
            tmp = rank(i)
            do while ( j >= 1 )
                if (arr(rank(j)) - arr(tmp) < 0) exit
                rank(j+1) = rank(j)
                j = j - 1
            end do
            rank(j+1) = tmp
        end do

    end subroutine insertion_rank_int_64

end module ranking

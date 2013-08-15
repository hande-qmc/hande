module sort

use const, only: i0, p

private
public :: insert_sort

! Insertion sort...
interface insert_sort
    module procedure insert_sort_int
end interface insert_sort

contains

    pure subroutine insert_sort_int(list)

        ! Sort in-place an integer list by the insertion sort algorithm.

        ! In/Out:
        !    list: list of integers.  List is sorted on output.

        integer, intent(inout) :: list(:)
        real :: tmp
        integer :: i, j

        ! Based on http://rosettacode.org/wiki/Sorting_algorithms/Insertion_sort#Fortran.
        ! Essentially a direct translation of the pseudo-code for the algorithm.

        do i = 2, ubound(list,dim=1)
            j = i - 1
            tmp = list(i)
            do while (j>=1 .and. list(j)>tmp)
                list(j+1) = list(j)
                j = j - 1
            end do
            list(j+1) = tmp
        end do

    end subroutine insert_sort_int

end module sort

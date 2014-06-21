module sort

use const, only: int_4, int_8, p

private
public :: qsort
public :: insert_sort

! Quicksort...
interface qsort
    module procedure qsort_int_32_list
    module procedure qsort_int_64_list
end interface qsort

! Insertion sort...
interface insert_sort
    module procedure insert_sort_int_32
    module procedure insert_sort_int_64
end interface insert_sort

contains

!--- In-place quicksort.  Uses the sample code in Numerical Recipies as a base.  ---

    pure subroutine qsort_int_32_list(list, head, nsort)

        ! Sort a 2D array of int_4 integers.

        ! list(:,i) is regarded as greater than list(:,j) if the first
        ! non-identical element between list(:,i) and list(:,j) is greater in
        ! list(:,i).

        ! In/Out:
        !    list: 2D array of int_s integers.  Sorted on output.
        ! In:
        !    head (optional): sort list up to and including list(:,:head) and
        !        leave the rest of the array untouched.  Default: sort the
        !        entire array.
        !    nsort (optional): sort list only using the first nsort elements in
        !        each 1D slice to compare entries (ie compare list(:nsort,i) and
        !        list(:nsort,j), so list is sorted according to list(:nsort,:)).
        !        Default: use entire slice.

        use bit_utils, only: operator(.bitstrgt.)

        integer(int_4), intent(inout) :: list(:,:)
        integer, intent(in), optional :: head, nsort

        ! Threshold.  When a sublist gets to this length, switch to using
        ! insertion sort to sort the sublist.
        integer, parameter :: switch_threshold = 7

        ! sort needs auxiliary storage of length 2*log_2(n).
        integer, parameter :: stack_max = 50

        integer :: pivot, lo, hi, i, j, ns
        integer(int_4) :: tmp(ubound(list,dim=1))

        ! Stack.  This is the auxilliary memory required by quicksort.
        integer :: stack(2,stack_max), nstack

        if (present(nsort)) then
            ns = nsort
        else
            ns = ubound(list, dim=1)
        end if

        nstack = 0
        lo = 1
        if (present(head)) then
            hi = head
        else
            hi = ubound(list, dim=2)
        end if
        do
            ! If the section/partition we are looking at is smaller than
            ! switch_threshold then perform an insertion sort.
            if (hi - lo < switch_threshold) then
                do j = lo + 1, hi
                    tmp = list(:,j)
                    do i = j - 1, 1, -1
                        if (tmp(1:ns) .bitstrgt. list(1:ns,i)) exit
                        list(:,i+1) = list(:,i)
                    end do
                    list(:,i+1) = tmp
                end do

                if (nstack == 0) exit
                hi = stack(2,nstack)
                lo = stack(1,nstack)
                nstack = nstack - 1

            else
                ! Otherwise start partitioning with quicksort.

                ! Pick the pivot element to be the median of list(:,lo), list(:,hi)
                ! and list(:,(lo+hi)/2).
                ! This largely overcomes a major problem with quicksort, where it
                ! degrades if the pivot is always the smallest element.
                pivot = (lo + hi)/2
                call swap_sublist(list(:,pivot), list(:,lo + 1))
                if (list(1:ns,lo) .bitstrgt. list(1:ns,hi)) then
                    call swap_sublist(list(:,lo), list(:,hi))
                end if
                if (list(1:ns,lo+1) .bitstrgt. list(1:ns,hi)) then
                    call swap_sublist(list(:,lo+1), list(:,hi))
                end if
                if (list(1:ns,lo) .bitstrgt. list(1:ns,lo+1)) then
                    call swap_sublist(list(:,lo), list(:,lo+1))
                end if

                i = lo + 1
                j = hi
                tmp = list(:,lo + 1) ! a is the pivot value
                do while (.true.)
                    ! Scan down list to find element > a.
                    i = i + 1
                    do while (tmp(1:ns) .bitstrgt. list(1:ns,i))
                        i = i + 1
                    end do

                    ! Scan down list to find element < a.
                    j = j - 1
                    do while (list(1:ns,j) .bitstrgt. tmp(1:ns))
                        j = j - 1
                    end do

                    ! When the pointers crossed, partitioning is complete.
                    if (j < i) exit

                    ! Swap the elements, so that all elements < a end up
                    ! in lower indexed variables.
                    call swap_sublist(list(:,i), list(:,j))
                end do

                ! Insert partitioning element
                list(:,lo + 1) = list(:,j)
                list(:,j) = tmp

                ! Push the larger of the partitioned sections onto the stack
                ! of sections to look at later.
                ! --> need fewest stack elements.
                nstack = nstack + 1

                ! With a stack_max of 50, we can sort arrays of length
                ! 1125899906842624.  It is safe to say this will never be
                ! exceeded, and so this test can be skipped.
!                if (nstack > stack_max) call stop_all('qsort_int_4_list', "parameter stack_max too small")

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

        pure subroutine swap_sublist(s1,s2)

            integer(int_4), intent(inout) :: s1(:), s2(:)
            integer(int_4) :: tmp(ubound(s1,dim=1))

            tmp = s1
            s1 = s2
            s2 = tmp

        end subroutine swap_sublist

    end subroutine qsort_int_32_list

    pure subroutine qsort_int_64_list(list, head, nsort)

        ! Sort a 2D array of int_8 integers.

        ! list(:,i) is regarded as greater than list(:,j) if the first
        ! non-identical element between list(:,i) and list(:,j) is greater in
        ! list(:,i).

        ! In/Out:
        !    list: 2D array of int_s integers.  Sorted on output.
        ! In:
        !    head (optional): sort list up to and including list(:,:head) and
        !        leave the rest of the array untouched.  Default: sort the
        !        entire array.
        !    nsort (optional): sort list only using the first nsort elements in
        !        each 1D slice to compare entries (ie compare list(:nsort,i) and
        !        list(:nsort,j), so list is sorted according to list(:nsort,:)).
        !        Default: use entire slice.

        use bit_utils, only: operator(.bitstrgt.)

        integer(int_8), intent(inout) :: list(:,:)
        integer, intent(in), optional :: head, nsort

        ! Threshold.  When a sublist gets to this length, switch to using
        ! insertion sort to sort the sublist.
        integer, parameter :: switch_threshold = 7

        ! sort needs auxiliary storage of length 2*log_2(n).
        integer, parameter :: stack_max = 50

        integer :: pivot, lo, hi, i, j, ns
        integer(int_8) :: tmp(ubound(list,dim=1))

        ! Stack.  This is the auxilliary memory required by quicksort.
        integer :: stack(2,stack_max), nstack

        if (present(nsort)) then
            ns = nsort
        else
            ns = ubound(list, dim=1)
        end if

        nstack = 0
        lo = 1
        if (present(head)) then
            hi = head
        else
            hi = ubound(list, dim=2)
        end if
        do
            ! If the section/partition we are looking at is smaller than
            ! switch_threshold then perform an insertion sort.
            if (hi - lo < switch_threshold) then
                do j = lo + 1, hi
                    tmp = list(:,j)
                    do i = j - 1, 1, -1
                        if (tmp(1:ns) .bitstrgt. list(1:ns,i)) exit
                        list(:,i+1) = list(:,i)
                    end do
                    list(:,i+1) = tmp
                end do

                if (nstack == 0) exit
                hi = stack(2,nstack)
                lo = stack(1,nstack)
                nstack = nstack - 1

            else
                ! Otherwise start partitioning with quicksort.

                ! Pick the pivot element to be the median of list(:,lo), list(:,hi)
                ! and list(:,(lo+hi)/2).
                ! This largely overcomes a major problem with quicksort, where it
                ! degrades if the pivot is always the smallest element.
                pivot = (lo + hi)/2
                call swap_sublist(list(:,pivot), list(:,lo + 1))
                if (list(1:ns,lo) .bitstrgt. list(1:ns,hi)) then
                    call swap_sublist(list(:,lo), list(:,hi))
                end if
                if (list(1:ns,lo+1) .bitstrgt. list(1:ns,hi)) then
                    call swap_sublist(list(:,lo+1), list(:,hi))
                end if
                if (list(1:ns,lo) .bitstrgt. list(1:ns,lo+1)) then
                    call swap_sublist(list(:,lo), list(:,lo+1))
                end if

                i = lo + 1
                j = hi
                tmp = list(:,lo + 1) ! a is the pivot value
                do while (.true.)
                    ! Scan down list to find element > a.
                    i = i + 1
                    do while (tmp(1:ns) .bitstrgt. list(1:ns,i))
                        i = i + 1
                    end do

                    ! Scan down list to find element < a.
                    j = j - 1
                    do while (list(1:ns,j) .bitstrgt. tmp(1:ns))
                        j = j - 1
                    end do

                    ! When the pointers crossed, partitioning is complete.
                    if (j < i) exit

                    ! Swap the elements, so that all elements < a end up
                    ! in lower indexed variables.
                    call swap_sublist(list(:,i), list(:,j))
                end do

                ! Insert partitioning element
                list(:,lo + 1) = list(:,j)
                list(:,j) = tmp

                ! Push the larger of the partitioned sections onto the stack
                ! of sections to look at later.
                ! --> need fewest stack elements.
                nstack = nstack + 1

                ! With a stack_max of 50, we can sort arrays of length
                ! 1125899906842624.  It is safe to say this will never be
                ! exceeded, and so this test can be skipped.
!                if (nstack > stack_max) call stop_all('qsort_int_8_list', "parameter stack_max too small")

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

        pure subroutine swap_sublist(s1,s2)

            integer(int_8), intent(inout) :: s1(:), s2(:)
            integer(int_8) :: tmp(ubound(s1,dim=1))

            tmp = s1
            s1 = s2
            s2 = tmp

        end subroutine swap_sublist

    end subroutine qsort_int_64_list

!--- In-place insertion sort.  ---

    pure subroutine insert_sort_int_32(list)

        ! Sort in-place an integer list by the insertion sort algorithm.

        ! In/Out:
        !    list: list of integers.  List is sorted on output.

        integer(int_4), intent(inout) :: list(:)
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

    end subroutine insert_sort_int_32

    pure subroutine insert_sort_int_64(list)

        ! Sort in-place an integer list by the insertion sort algorithm.

        ! In/Out:
        !    list: list of integers.  List is sorted on output.

        integer(int_8), intent(inout) :: list(:)
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

    end subroutine insert_sort_int_64

end module sort

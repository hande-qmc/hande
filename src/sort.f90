module sort

use const, only: int_32, int_64, p

implicit none

private
public :: qsort
public :: insert_sort

! Quicksort...
interface qsort
    module procedure qsort_int_32_list
    module procedure qsort_int_64_list
    module procedure qsort_psip_info
end interface qsort

! Insertion sort...
interface insert_sort
    module procedure insert_sort_int_32
    module procedure insert_sort_int_64
end interface insert_sort

contains

!--- In-place quicksort.  Uses the sample code in Numerical Recipies as a base.  ---

    pure subroutine qsort_int_32_list(list, head, nsort)

        ! Sort a 2D array of int_32 integers.

        ! list(:,i) is regarded as greater than list(:,j) if the first
        ! non-identical element between list(:,i) and list(:,j) is greater in
        ! list(:,i).

        ! In/Out:
        !    list: 2D array of int_32 integers.  Sorted on output.
        ! In:
        !    head (optional): sort list up to and including list(:,:head) and
        !        leave the rest of the array untouched.  Default: sort the
        !        entire array.
        !    nsort (optional): sort list only using the first nsort elements in
        !        each 1D slice to compare entries (ie compare list(:nsort,i) and
        !        list(:nsort,j), so list is sorted according to list(:nsort,:)).
        !        Default: use entire slice.

        use bit_utils, only: operator(.bitstrge.), operator(.bitstrgt.)

        integer(int_32), intent(inout) :: list(:,:)
        integer, intent(in), optional :: head, nsort

        ! Threshold.  When a sublist gets to this length, switch to using
        ! insertion sort to sort the sublist.
        integer, parameter :: switch_threshold = 7

        ! sort needs auxiliary storage of length 2*log_2(n).
        integer, parameter :: stack_max = 50

        integer :: pivot, lo, hi, i, j, ns
        integer(int_32) :: tmp(ubound(list,dim=1))

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
                        if (tmp(1:ns) .bitstrge. list(1:ns,i)) exit
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
!                if (nstack > stack_max) call stop_all('qsort_int_32_list', "parameter stack_max too small")

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

            integer(int_32), intent(inout) :: s1(:), s2(:)
            integer(int_32) :: tmp(ubound(s1,dim=1))

            tmp = s1
            s1 = s2
            s2 = tmp

        end subroutine swap_sublist

    end subroutine qsort_int_32_list

    pure subroutine qsort_int_64_list(list, head, nsort)

        ! Sort a 2D array of int_64 integers.

        ! list(:,i) is regarded as greater than list(:,j) if the first
        ! non-identical element between list(:,i) and list(:,j) is greater in
        ! list(:,i).

        ! In/Out:
        !    list: 2D array of int_64 integers.  Sorted on output.
        ! In:
        !    head (optional): sort list up to and including list(:,:head) and
        !        leave the rest of the array untouched.  Default: sort the
        !        entire array.
        !    nsort (optional): sort list only using the first nsort elements in
        !        each 1D slice to compare entries (ie compare list(:nsort,i) and
        !        list(:nsort,j), so list is sorted according to list(:nsort,:)).
        !        Default: use entire slice.

        use bit_utils, only: operator(.bitstrge.), operator(.bitstrgt.)

        integer(int_64), intent(inout) :: list(:,:)
        integer, intent(in), optional :: head, nsort

        ! Threshold.  When a sublist gets to this length, switch to using
        ! insertion sort to sort the sublist.
        integer, parameter :: switch_threshold = 7

        ! sort needs auxiliary storage of length 2*log_2(n).
        integer, parameter :: stack_max = 50

        integer :: pivot, lo, hi, i, j, ns
        integer(int_64) :: tmp(ubound(list,dim=1))

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
                        if (tmp(1:ns) .bitstrge. list(1:ns,i)) exit
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
!                if (nstack > stack_max) call stop_all('qsort_int_64_list', "parameter stack_max too small")

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

            integer(int_64), intent(inout) :: s1(:), s2(:)
            integer(int_64) :: tmp(ubound(s1,dim=1))

            tmp = s1
            s1 = s2
            s2 = tmp

        end subroutine swap_sublist

    end subroutine qsort_int_64_list

    pure subroutine qsort_psip_info(nstates, states, pops, dat)

        ! Sort a set of psip information (states, populations and data) in order according
        ! to the state labels.

        ! states(:,i) is regarded as greater than states(:,j) if the first
        ! non-identical element between states(:,i) and states(:,j) is greater in
        ! states(:,i).

        ! In:
        !    nstates: number of occupied states.
        ! In/Out:
        !    states: 2D array of i0 integers containing the state label for each occupied state.
        !        Sorted on output.
        !    pops, dat: population and data arrays for each state.  Sorted by states on output.

        ! Note: the size of the first dimension of states is immaterial, as we do a comparison
        ! based on the entire slice.  The second dimensions of pops, dat and states must be >= nstates.

        use const, only: int_p, i0, p
        use bit_utils, only: operator(.bitstrge.), operator(.bitstrgt.)

        integer, intent(in) :: nstates
        integer(i0), intent(inout) :: states(:,:)
        integer(int_p), intent(inout) :: pops(:,:)
        real(p), intent(inout) :: dat(:,:)

        ! Threshold.  When a substates gets to this length, switch to using
        ! insertion sort to sort the substates.
        integer, parameter :: switch_threshold = 7

        ! sort needs auxiliary storage of length 2*log_2(n).
        integer, parameter :: stack_max = 50

        integer :: pivot, lo, hi, i, j
        integer(i0) :: tmp_state(ubound(states,dim=1))
        integer(int_p) :: tmp_pop(ubound(pops,dim=1))
        real(p) :: tmp_dat(ubound(dat,dim=1))

        ! Stack.  This is the auxilliary memory required by quicksort.
        integer :: stack(2,stack_max), nstack

        hi = nstates

        nstack = 0
        lo = 1

        do
            ! If the section/partition we are looking at is smaller than
            ! switch_threshold then perform an insertion sort.
            if (hi - lo < switch_threshold) then
                do j = lo + 1, hi
                    tmp_state = states(:,j)
                    tmp_pop = pops(:,j)
                    tmp_dat = dat(:,j)
                    do i = j - 1, 1, -1
                        if (tmp_state .bitstrge. states(:,i)) exit
                        states(:,i+1) = states(:,i)
                        pops(:,i+1) = pops(:,i)
                        dat(:,i+1) = dat(:,i)
                    end do
                    states(:,i+1) = tmp_state
                    pops(:,i+1) = tmp_pop
                    dat(:,i+1) = tmp_dat
                end do

                if (nstack == 0) exit
                hi = stack(2,nstack)
                lo = stack(1,nstack)
                nstack = nstack - 1

            else
                ! Otherwise start partitioning with quicksort.

                ! Pick the pivot element to be the median of states(:,lo), states(:,hi)
                ! and states(:,(lo+hi)/2).
                ! This largely overcomes a major problem with quicksort, where it
                ! degrades if the pivot is always the smallest element.
                pivot = (lo + hi)/2
                call swap_states(states(:,pivot), pops(:,pivot), dat(:,pivot), states(:,lo+1), pops(:,lo+1), dat(:,lo+1))
                if (states(:,lo) .bitstrgt. states(:,hi)) then
                    call swap_states(states(:,lo), pops(:,lo), dat(:,lo), states(:,hi), pops(:,hi), dat(:,hi))
                end if
                if (states(:,lo+1) .bitstrgt. states(:,hi)) then
                    call swap_states(states(:,lo+1), pops(:,lo+1), dat(:,lo+1), states(:,hi), pops(:,hi), dat(:,hi))
                end if
                if (states(:,lo) .bitstrgt. states(:,lo+1)) then
                    call swap_states(states(:,lo), pops(:,lo), dat(:,lo), states(:,lo+1), pops(:,lo+1), dat(:,lo+1))
                end if

                i = lo + 1
                j = hi
                tmp_state = states(:,lo+1) ! a is the pivot value
                tmp_pop = pops(:,lo+1)
                tmp_dat = dat(:,lo+1)
                do while (.true.)
                    ! Scan down states to find element > a.
                    i = i + 1
                    do while (tmp_state .bitstrgt. states(:,i))
                        i = i + 1
                    end do

                    ! Scan down states to find element < a.
                    j = j - 1
                    do while (states(:,j) .bitstrgt. tmp_state)
                        j = j - 1
                    end do

                    ! When the pointers crossed, partitioning is complete.
                    if (j < i) exit

                    ! Swap the elements, so that all elements < a end up
                    ! in lower indexed variables.
                    call swap_states(states(:,i), pops(:,i), dat(:,i), states(:,j), pops(:,j), dat(:,j))
                end do

                ! Insert partitioning element
                states(:,lo + 1) = states(:,j)
                pops(:,lo + 1) = pops(:,j)
                dat(:,lo + 1) = dat(:,j)
                states(:,j) = tmp_state
                pops(:,j) = tmp_pop
                dat(:,j) = tmp_dat

                ! Push the larger of the partitioned sections onto the stack
                ! of sections to look at later.
                ! --> need fewest stack elements.
                nstack = nstack + 1

                ! With a stack_max of 50, we can sort arrays of length
                ! 1125899906842624.  It is safe to say this will never be
                ! exceeded, and so this test can be skipped.
!                if (nstack > stack_max) call stop_all('qsort_int_64_states', "parameter stack_max too small")

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

        pure subroutine swap_states(s1,p1,d1,s2,p2,d2)

            integer(i0), intent(inout) :: s1(:), s2(:)
            integer(int_p), intent(inout) :: p1(:), p2(:)
            real(p), intent(inout) :: d1(:), d2(:)
            integer(i0) :: tmp_state(ubound(s1,dim=1))
            integer(int_p) :: tmp_pop(ubound(p1,dim=1))
            real(p) :: tmp_dat(ubound(d1,dim=1))

            tmp_state = s1
            s1 = s2
            s2 = tmp_state

            tmp_pop = p1
            p1 = p2
            p2 = tmp_pop

            tmp_dat = d1
            d1 = d2
            d2 = tmp_dat

        end subroutine swap_states

    end subroutine qsort_psip_info

!--- In-place insertion sort.  ---

    pure subroutine insert_sort_int_32(list)

        ! Sort in-place an integer list by the insertion sort algorithm.

        ! In/Out:
        !    list: list of integers.  List is sorted on output.

        integer(int_32), intent(inout) :: list(:)
        integer(int_32) :: tmp
        integer :: i, j

        ! Based on http://rosettacode.org/wiki/Sorting_algorithms/Insertion_sort#Fortran.
        ! Essentially a direct translation of the pseudo-code for the algorithm.

        do i = 2, ubound(list,dim=1)
            j = i - 1
            tmp = list(i)
            do while (j>=1)
                if (list(j)<=tmp) exit ! Can't combine with while conditional as Fortran can evaluate boolean statements in any order.
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

        integer(int_64), intent(inout) :: list(:)
        integer(int_64) :: tmp
        integer :: i, j

        ! Based on http://rosettacode.org/wiki/Sorting_algorithms/Insertion_sort#Fortran.
        ! Essentially a direct translation of the pseudo-code for the algorithm.

        do i = 2, ubound(list,dim=1)
            j = i - 1
            tmp = list(i)
            do while (j>=1)
                if (list(j)<=tmp) exit ! Can't combine with while conditional as Fortran can evaluate boolean statements in any order.
                list(j+1) = list(j)
                j = j - 1
            end do
            list(j+1) = tmp
        end do

    end subroutine insert_sort_int_64

end module sort

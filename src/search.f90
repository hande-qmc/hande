module search

implicit none

private
public :: binary_search

interface binary_search
    ! All binary search procedures work in essentially the same way with the
    ! same interface (bar list to search & item to search for).
    module procedure binary_search_i0_list
    module procedure binary_search_integer
end interface binary_search

contains

    pure subroutine binary_search_i0_list(list, item, istart, iend, hit, pos)

        ! Find where an item resides in a list of such items.
        ! Only elements between istart and iend are examined (use the
        ! array boundaries in the worst case).
        !
        ! In:
        !    list: a sorted i0 integer 2D list/array; the first dimension
        !        corresponds to 1D arrays to compare to item.
        !    item: an i0 integer 1D list/array.
        !    istart: first position to examine in the list.
        !    iend: last position to examine in the list.
        ! Out:
        !    hit: true if found item in list.
        !    pos: the position corresponding to item in list.
        !        If hit is true, then the element in this position is the same
        !        as item, else this is where item should go to keep the list
        !        sorted.

        use const, only: i0

        integer(i0), intent(in) :: list(:,:), item(:)
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

            ! Assume item doesn't exist in the list initially.
            hit = .false.

            do while (hi /= lo)
                ! Narrow the search range down in steps.

                ! Mid-point.
                ! We shift one of the search limits to be the mid-point.
                ! The successive dividing the search range by 2 gives a O[log N]
                ! search algorithm.
                pos = (hi+lo)/2

                compare = list_compare(list(:,pos), item)
                select case(compare)
                case (0)
                    ! hit!
                    hit = .true.
                    exit
                case(1)
                    ! list(:,pos) is "smaller" than item.
                    ! The lowest position item can take is hence pos + 1 (i.e. if
                    ! item is greater than pos by smaller than pos + 1).
                    lo = pos + 1
                case(-1)
                    ! list(:,pos) is "greater" than item.
                    ! The highest position item can take is hence pos (i.e. if item is
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
            ! Because list can expand (i.e. we might be searching for an
            ! element which doesn't exist yet) the binary search can find either
            ! the element before or after where item should be placed.
            if (hi == lo) then
                compare = list_compare(list(:,hi), item)
                select case(compare)
                case (0)
                    ! hit!
                    hit = .true.
                    pos = hi
                case(1)
                    ! list(:,pos) is "smaller" than item.
                    ! item should be placed in the next slot.
                    pos = hi + 1
                case(-1)
                    ! list(:,pos) is "greater" than item.
                    ! item should ber placed here.
                    pos = hi
                end select
            end if

        end if

        contains

            pure function list_compare(item1, item2) result(compare)

                ! In:
                !    f1(total_basis_length): bit string representation of the Slater
                !        determinant.
                !    f2(total_basis_length): bit string representation of the Slater
                !        determinant.
                !    (For DMQMC this bitstring contains information for both determinants)
                ! Returns:
                !    0 if f1 and f2 are identical;
                !    1 if the first non-identical element in f1 is smaller than the
                !    corresponding element in f2;
                !    -1 if the first non-identical element in f1 is greater than the
                !    corresponding element in f2;

                integer :: compare
                integer(i0), intent(in) :: item1(:), item2(:)

                integer :: i

                compare = 0
                do i = 1, ubound(item1, dim=1)
                    if (item1(i) < item2(i)) then
                        compare = 1
                        exit
                    else if (item1(i) > item2(i)) then
                        compare = -1
                        exit
                    end if
                end do

            end function list_compare

    end subroutine binary_search_i0_list

    pure subroutine binary_search_integer(list, item, istart, iend, hit, pos)

        ! Find where an item resides in a list of such items.
        ! Only elements between istart and iend are examined (use the
        ! array boundaries in the worst case).
        !
        ! In:
        !    list: a sorted integer list/array.
        !    item: an integer to find in list.
        !    istart: first position to examine in the list.
        !    iend: last position to examine in the list.
        ! Out:
        !    hit: true if found item in list.
        !    pos: the position corresponding to item in list.
        !        If hit is true, then the element in this position is the same
        !        as item, else this is where item should go to keep the list
        !        sorted.

        integer, intent(in) :: list(:), item
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

            ! Assume item doesn't exist in the list initially.
            hit = .false.

            do while (hi /= lo)
                ! Narrow the search range down in steps.

                ! Mid-point.
                ! We shift one of the search limits to be the mid-point.
                ! The successive dividing the search range by 2 gives a O[log N]
                ! search algorithm.
                pos = (hi+lo)/2

                compare = item - list(pos)
                if (compare == 0) then
                    ! hit!
                    hit = .true.
                    exit
                else if (compare > 0) then
                    ! list(pos) is "smaller" than item.
                    ! The lowest position item can take is hence pos + 1 (i.e. if
                    ! item is greater than pos by smaller than pos + 1).
                    lo = pos + 1
                else
                    ! list(pos) is "greater" than item.
                    ! The highest position item can take is hence pos (i.e. if item is
                    ! smaller than pos but greater than pos - 1).  This is why
                    ! we differ slightly from a standard binary search (where lo
                    ! is set to be pos+1 and hi to be pos-1 accordingly), as
                    ! a standard binary search assumes that the element you are
                    ! searching for actually appears in the array being
                    ! searched...
                    hi = pos
                end if

            end do

            ! If hi == lo, then we have narrowed the search down to one position but
            ! not checked if that position is the item we're hunting for.
            ! Because list can expand (i.e. we might be searching for an
            ! element which doesn't exist yet) the binary search can find either
            ! the element before or after where item should be placed.
            if (hi == lo) then
                compare = item - list(hi)
                if (compare == 0) then
                    ! hit!
                    hit = .true.
                    pos = hi
                else if (compare > 0) then
                    ! list(pos) is smaller than item.
                    ! item should be placed in the next slot.
                    pos = hi + 1
                else
                    ! list(pos) is greater than item.
                    ! item should be placed here.
                    pos = hi
                end if
            end if

        end if

    end subroutine binary_search_integer

end module search

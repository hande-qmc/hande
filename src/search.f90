module search

implicit none

private
public :: binary_search

interface binary_search
    ! All binary search procedures work in essentially the same way with the
    ! same interface (bar list to search & item to search for).
    module procedure binary_search_i0_list
    module procedure binary_search_int_32
    module procedure binary_search_int_64
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
        use bit_utils, only: bit_str_cmp

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

                compare = bit_str_cmp(list(:,pos), item)
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
                compare = bit_str_cmp(list(:,hi), item)
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

    end subroutine binary_search_i0_list

    pure subroutine binary_search_int_32(list, item, istart, iend, hit, pos)

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

        use const, only: int_32

        integer(int_32), intent(in) :: list(:), item
        integer, intent(in) :: istart, iend
        logical, intent(out) :: hit
        integer, intent(out) :: pos

        integer :: hi, lo
        integer(int_32) :: compare

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
                if (compare == 0_int_32) then
                    ! hit!
                    hit = .true.
                    exit
                else if (compare > 0_int_32) then
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
                if (compare == 0_int_32) then
                    ! hit!
                    hit = .true.
                    pos = hi
                else if (compare > 0_int_32) then
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

    end subroutine binary_search_int_32

    pure subroutine binary_search_int_64(list, item, istart, iend, hit, pos)

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

        use const, only: int_64

        integer(int_64), intent(in) :: list(:), item
        integer, intent(in) :: istart, iend
        logical, intent(out) :: hit
        integer, intent(out) :: pos

        integer :: hi, lo
        integer(int_64) :: compare

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
                if (compare == 0_int_64) then
                    ! hit!
                    hit = .true.
                    exit
                else if (compare > 0_int_64) then
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
                if (compare == 0_int_64) then
                    ! hit!
                    hit = .true.
                    pos = hi
                else if (compare > 0_int_64) then
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

    end subroutine binary_search_int_64

end module search

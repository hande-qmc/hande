module search

use const
use excitations, only: get_excitation_level

implicit none

private
public :: binary_search, tree_t, node_t, tree_add, tree_search

type node_t
    ! For use in BK-tree search of secondary references
    ! See also comments below under tree_add and tree_search
    ! A node contains a bitstring, and an array of edges, which go from 1 to the maximum excitation level
    ! possible in the current Fock space
    ! Each of the edges point to another node
    integer(i0), allocatable :: bstring(:)
    type(node_t), pointer :: edges(:) => null()
end type node_t

type tree_t
    ! For use in BK-tree search of secondary references
    ! A tree stores a pointer to the root node, and some of the constants specific to the system
    type(node_t), pointer :: root => null()
    integer :: n_secondary_ref, ex_lvl, max_excit
end type tree_t

interface binary_search
    ! All binary search procedures work in essentially the same way with the
    ! same interface (bar list to search & item to search for).
    module procedure binary_search_i0_list
    module procedure binary_search_int_32
    module procedure binary_search_int_64
    module procedure binary_search_real_p
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

    pure subroutine binary_search_real_p(list, item, istart, iend, hit, pos)

        ! Find where an item resides in a list of such items.
        ! Only elements between istart and iend are examined (use the
        ! array boundaries in the worst case).
        !
        ! NB item is tested to be within depsilon of list values, rather than equality.
        !
        ! In:
        !    list: a sorted real list/array.
        !    item: a real to find in list.
        !    istart: first position to examine in the list.
        !    iend: last position to examine in the list.
        ! Out:
        !    hit: true if found item in list.
        !    pos: the position corresponding to item in list.
        !        If hit is true, then the element in this position is the same
        !        as item, else this is where item should go to keep the list
        !        sorted.

        use const, only: p, depsilon

        real(p), intent(in) :: list(:), item
        integer, intent(in) :: istart, iend
        logical, intent(out) :: hit
        integer, intent(out) :: pos

        integer :: hi, lo
        real(p) :: compare

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
                if (abs(compare) < depsilon) then
                    ! hit!
                    hit = .true.
                    exit
                else if (compare > 0.0_p) then
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
                if (abs(compare) < depsilon) then
                    ! hit!
                    hit = .true.
                    pos = hi
                else if (compare > 0.0_p) then
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

    end subroutine binary_search_real_p

    subroutine tree_add(this, next_bstring)

        ! This builds a Hamming-distance BK-tree for k-nearest neighbour search
        ! A very good explanation can be found here: 
        ! https://daniel-j-h.github.io/post/nearest-neighbors-in-metric-spaces/

        type(tree_t), intent(inout) :: this
        integer(i0), intent(in) :: next_bstring(:)
        type(node_t), pointer :: curr_node => null()
        integer(i0), pointer :: parent(:) => null()
        type(node_t), pointer :: child(:) => null()
        integer :: excit_lvl

        ! Initialises the root node
        if (.not. associated(this%root)) then
            allocate(this%root)
        end if
        curr_node => this%root
        if (.not. allocated(curr_node%bstring)) then
            ! If the node was initialised but empty, store the bitstring here
            curr_node%bstring = next_bstring
        else if (.not. associated(curr_node%edges)) then
            ! If the node has a bitstring, compare and allocate edges, and store the bitstring at the end 
            ! of the correct edge
            excit_lvl = get_excitation_level(curr_node%bstring, next_bstring)
            allocate(curr_node%edges(this%max_excit))
            curr_node%edges(excit_lvl)%bstring = next_bstring
        else
            do
                ! If the node has both a bitstring and allocated edges, that means
                ! you have a descend another level

                ! This looks redundant but it's not, since the do loop goes on forever
                ! you'll get to the end of the 'branch' and see a node with unallocated edges
                if (.not. associated(curr_node%edges)) then
                  allocate(curr_node%edges(this%max_excit))
                end if
                parent => curr_node%bstring
                child => curr_node%edges
                excit_lvl = get_excitation_level(parent, next_bstring)
                curr_node => child(excit_lvl)
                if (.not. allocated(curr_node%bstring)) then
                    curr_node%bstring = next_bstring
                    exit
                end if
            end do
        end if

    end subroutine tree_add

    recursive function tree_search(this, new_bstring, curr_node, offset) result(hit)
        ! This is a recursive, depth-first traversal of a Hamming-distance BK-tree.
        ! This implementation was in part inspired by the C++ implementation here
        ! https://www.geeksforgeeks.org/bk-tree-introduction-implementation/ 
        ! and the Python implementation here
        ! https://github.com/ahupp/bktree
        !
        ! [TODO]: A non-recursive search may improve performance, and an example in Python 
        ! using deque can be found at https://github.com/benhoyt/pybktree


        type(tree_t), intent(in) :: this
        integer(i0), intent(in) :: new_bstring(:)
        integer :: excit_lvl, endlvl, startlvl, lvl
        type(node_t), intent(in) :: curr_node
        integer, intent(in), optional :: offset
        logical :: hit
        integer :: ex_lvl

        excit_lvl = get_excitation_level(curr_node%bstring, new_bstring)
        
        if (.not. present(offset)) then
            ex_lvl = this%ex_lvl
        else
            ex_lvl = offset
        end if

        if (excit_lvl <= ex_lvl) then
            ! the current node is within ex_lvl of the new bitstring
            hit = .true.
            return
        end if
        
        ! If not, descend recursively into a subspace of the tree
        startlvl = excit_lvl - ex_lvl
        endlvl = excit_lvl + ex_lvl
        if (startlvl <= ex_lvl) then
            ! We already establishedd that this bitstring is not within ex_lvl of the current node
            startlvl = ex_lvl + 1
        end if
        if (endlvl > this%max_excit) then
            endlvl = this%max_excit
        end if

        do lvl = startlvl, endlvl
            if (.not. associated(curr_node%edges)) then
                ! reached the bottom, nothing found
                hit = .false.
                return
            end if
            if (.not. allocated(curr_node%edges(lvl)%bstring)) then
                ! Empty node, look into the neighbour now
                continue
            else
                ! If this node has information in it, dig into it
                ! This is the recursive bit, and it is doing the descending / traversal into the tree
                if (.not. present(offset)) then
                    hit = tree_search(this, new_bstring, curr_node%edges(lvl))
                else
                    hit = tree_search(this, new_bstring, curr_node%edges(lvl), offset)
                end if
                if (hit) then
                    ! Again might look redundant but this check makes sure once hit = .true.
                    ! the result gets propagated to the top level and breaks recursion
                    hit = .true.
                    return
                end if
            end if
        end do

        ! nothing found after recursion, return false
        hit = .false. 
        return

    end function tree_search
end module search

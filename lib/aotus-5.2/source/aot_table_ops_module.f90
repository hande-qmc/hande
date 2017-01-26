! Copyright (C) 2011-2013 German Research School for Simulation Sciences GmbH,
!                         Aachen and others.
!               2013-2014 University of Siegen.
! Please see the COPYRIGHT file in this directory for details.

!> This module provides general operations on Lua tables.
!!
!! These operations are a common set of actions, that are used by the various
!! type specific implementations.
module aot_table_ops_module
  use flu_binding
  use aot_kinds_module, only: double_k, single_k, long_k
  use aot_top_module, only: aot_top_get_val

  implicit none

  private

  public :: aot_table_open, aot_table_close
  public :: aot_table_top, aot_table_length, aot_table_first, aot_table_push

contains

  !> Return the position at the top of the stack as a
  !! table handle.
  !!
  !! If it actually exists and is a table, this handle can be used
  !! for further operations on that table.
  !! Otherwise a 0 will be returned.
  function aot_table_top(L) result(thandle)
    type(flu_state) :: L !< Handle for the Lua script.

    !> A handle for the table on the top of the stack to access it.
    integer :: thandle

    if (flu_isNoneOrNil(L, -1) .or. (.not. flu_isTable(L, -1))) then
      thandle = 0
      call flu_pop(L)
    else
      thandle = flu_gettop(L)
    end if
  end function aot_table_top


  !> This subroutine tries to open a table, and returns a handle for it.
  !!
  !! If parent is present, the table is tried to open within that table.
  !! Return its position in the stack as a handle for this
  !! table. If it does not exist or the table entry is not
  !! a table itself, the handle will be set to 0.
  !! The table can be looked up either by position or name.
  !!
  !! If a key is present but no parent, a global table is opened.
  !! If neither key nor parent is present, a new table is created.
  !!
  !! After the table is opened, the returned handle can be used to access its
  !! components.
  subroutine aot_table_open(L, parent, thandle, key, pos)
    type(flu_state) :: L !< Handle for the Lua script.

    !> Handle of the table containing the requested table.
    integer, intent(in), optional :: parent

    !> A handle for the table to access it, 0 if no table available.
    integer, intent(out) :: thandle

    !> Name of the entry in the parent table to access.
    !!
    !! The key takes precedence over the position, if both are provided.
    !! In this case the positional address is only tried, if the access to the
    !! key failed.
    character(len=*), intent(in), optional :: key

    !> Position of the entry in the parent table to access.
    integer, intent(in), optional :: pos

    if (present(parent)) then
      call aot_table_push(L, parent, key, pos)
      thandle = aot_table_top(L)
    else
      if (present(key)) then
        call flu_getglobal(L, key)
        thandle = aot_table_top(L)
      else
        call flu_createtable(L, 0, 0)
        thandle = flu_gettop(L)
      end if
    end if

  end subroutine aot_table_open


  !> Close a table again.
  !!
  !! This is done by popping all values above and itself from the stack.
  subroutine aot_table_close(L, thandle)
    type(flu_state) :: L !< Handle for the Lua script.

    !> Handle of the table to close.
    integer, intent(in) :: thandle

    if (thandle > 0) call flu_settop(L, thandle-1)
  end subroutine aot_table_close


  !> This subroutine tries to push the value of the entry given by key or pos
  !! within the table thandle on the lua stack.
  !!
  !! If no corresponding value is found, a nil value is pushed to the stack.
  !! Key and pos are both optional, but one of them has to be supplied. If one
  !! is supplied, the key is checked first and only if this fails the entry at
  !! pos will be looked up.
  subroutine aot_table_push(L, thandle, key, pos)
    type(flu_state) :: L !< Handle for the Lua script.

    !> Handle to the table to look in.
    integer :: thandle

    !> Name of the entry to push to the stack.
    character(len=*), intent(in), optional :: key

    !> Position of the entry to push to the stack.
    integer, intent(in), optional :: pos

    if (thandle /= 0) then
      ! Only proceed if thandle is actually a table
      ! (Should be received with aot_table_global or aot_table_top)

      if (present(key)) then
        ! Try to look up the given key first
        call flu_getfield(L, thandle, key)
        if (flu_isNoneOrNil(L, -1)) then
          ! If this is not found, try to retrieve
          ! the value at the given position
          if (present(pos)) then
            call flu_pop(L)
            call flu_pushInteger(L, pos)
            call flu_getTable(L, thandle)
          end if
        end if
      else
        ! No key to look up, just check the given position
        if (present(pos)) then
          call flu_pushInteger(L, pos)
          call flu_getTable(L, thandle)
        else
          ! Neither key nor pos present, nothing to look up
          ! Just push a NIL onto the stack as a result
          call flu_pushnil(L)
        end if
      end if

    else

      call flu_pushnil(L)

    end if

  end subroutine aot_table_push


  !> Load the first key-value pair of table thandle on the
  !! stack.
  !!
  !! This serves as an entry point, further traversal
  !! can be done by flu_next(L, thandle).
  !! If there are no entries in the table the function
  !! returns false, otherwise the result will be true.
  function aot_table_first(L, thandle) result(exists)
    type(flu_state) :: L !< Handle for the Lua script.

    !> Handle to the table to get the first entry of.
    integer, intent(in) :: thandle

    !> The return value signals, if there actually is such a first entry.
    logical :: exists

    if (thandle /= 0) then
      call flu_pushnil(L)
      exists = flu_next(L, thandle)
    else
      exists = .false.
    end if
  end function aot_table_first


  !> Count the entries in a lua table.
  function aot_table_length(L, thandle) result(length)
    type(flu_state) :: L !< Handle for the Lua script.

    !> Handle of the table to count the enries in.
    integer, intent(in) :: thandle

    !> Returns the number of entries in the table.
    integer :: length

    length = 0
    if (aot_table_first(L, thandle)) then
      do
        length = length + 1
        call flu_pop(L)
        if (.not. flu_next(L, thandle)) exit
      end do
    end if
  end function aot_table_length


end module aot_table_ops_module

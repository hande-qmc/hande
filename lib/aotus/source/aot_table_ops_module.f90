! Copyright (C) 2011-2013 German Research School for Simulation Sciences GmbH,
!                         Aachen and others.
!               2013-2016 University of Siegen.
! Please see the COPYRIGHT file in this directory for details.

!> This module provides general operations on Lua tables.
!!
!! These operations are a common set of actions, that are used by the various
!! type specific implementations.
module aot_table_ops_module
  use flu_binding
  use flu_kinds_module, only: double_k, single_k, long_k
  use aot_top_module, only: aot_top_get_val

  implicit none

  private

  public :: aot_table_open, aot_table_close
  public :: aot_table_top, aot_table_length, aot_table_first, aot_table_push
  public :: aot_push
  public :: aot_type_of

  interface aot_push
    module procedure aot_table_push
  end interface aot_push


contains


  !> Return the position at the top of the stack as a
  !! table handle.
  !!
  !! If it actually exists and is a table, this handle can be used
  !! for further operations on that table.
  !! Otherwise a 0 will be returned.
  function aot_table_top(L) result(thandle)
    type(flu_state) :: L !! Handle for the Lua script.

    !> A handle for the table on the top of the stack to access it.
    integer :: thandle

    if (.not. flu_isTable(L, -1)) then
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
  !! Only passing pos, without a thandle is erroneous and always
  !! results in a thandle = 0.
  !!
  !! After the table is opened, the returned handle can be used to access its
  !! components.
  subroutine aot_table_open(L, parent, thandle, key, pos)
    type(flu_state) :: L !! Handle for the Lua script.

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

    integer :: luatype

    thandle = 0

    if (present(parent)) then
      call aot_table_push(L, parent, key, pos)
      thandle = aot_table_top(L)
    else
      if (present(key)) then
        luatype = flu_getglobal(L, key)
        thandle = aot_table_top(L)
      else if (.not. present(pos)) then
        call flu_createtable(L, 0, 0)
        thandle = flu_gettop(L)
      end if
    end if

  end subroutine aot_table_open


  !> Close a table again.
  !!
  !! This is done by popping all values above and itself from the stack.
  subroutine aot_table_close(L, thandle)
    type(flu_state) :: L !! Handle for the Lua script.

    !> Handle of the table to close.
    integer, intent(in) :: thandle

    if (thandle > 0) call flu_settop(L, thandle-1)
  end subroutine aot_table_close


  !> This subroutine tries to push the value of the entry given by key or pos
  !! within the table thandle onto the Lua stack.
  !!
  !! If no corresponding value is found, a nil value is pushed to the stack.
  !! Key, pos and thandle are all optional.
  !! If no thandle is provided, the key will be obtained as a global variable.
  !! When none of thandle, key and pos are provided, the subroutine does
  !! nothing and the resulting type returned in toptype is the type of the
  !! current top entry in the Lua stack.
  !! Passing only pos without thandle is illegal and will result in a NIL
  !! value on the top of the stack.
  subroutine aot_table_push(L, thandle, key, pos, toptype)
    type(flu_state) :: L !! Handle for the Lua script.

    !> Handle to the table to look in.
    integer, intent(in), optional :: thandle

    !> Name of the entry to push to the stack.
    character(len=*), intent(in), optional :: key

    !> Position of the entry to push to the stack.
    integer, intent(in), optional :: pos

    integer, intent(out), optional :: toptype

    integer :: loctype

    loctype = FLU_TNIL

    istable: if (present(thandle)) then

      if (thandle /= 0) then
        ! Only proceed if thandle is actually a table
        ! (Should be received with aot_table_global or aot_table_top)

        if (present(key)) then
          ! Try to look up the given key first
          loctype = flu_getfield(L, thandle, key)
          if ((loctype == FLU_TNONE) .or. (loctype == FLU_TNIL)) then
            ! If this is not found, try to retrieve
            ! the value at the given position
            if (present(pos)) then
              call flu_pop(L)
              call flu_pushInteger(L, pos)
              loctype = flu_getTable(L, thandle)
            end if
          end if
        else
          ! No key to look up, just check the given position
          if (present(pos)) then
            call flu_pushInteger(L, pos)
            loctype = flu_getTable(L, thandle)
          else
            ! Neither key nor pos present, nothing to look up
            ! Just push a NIL onto the stack as a result
            call flu_pushnil(L)
          end if
        end if

      else

        call flu_pushnil(L)

      end if

    else istable

      if (present(key)) then
        ! Try to look up the given key as a global variable
        loctype = flu_getglobal(L, key)
      else
        ! No key, no thandle, treat this as a no-op if also no pos is provided
        ! and return the type of the current top of the Lua stack.
        if (present(pos)) then
          ! Passing pos without thandle is illegal, and we always push a NIL
          ! in this case.
          call flu_pushnil(L)
        else
          loctype = flu_type(L, -1)
        end if
      end if

    end if istable

    if (present(toptype)) then
      toptype = loctype
    end if

  end subroutine aot_table_push


  !> Get the Lua object in table thandle under the given key or pos on the
  !! top of the stack and return the Lua type of the gotten entry.
  !!
  !! This might be used to get a Lua entry to the top of the stack without
  !! knowing its type beforehand, and then deciding what to load, based on
  !! the type.
  !! Lua types are encoded as integer values and available in the
  !! [[flu_binding]] module.
  !!
  !! - FLU_TNONE    : not existing
  !! - FLU_TNIL     : not available
  !! - FLU_TBOOLEAN : logical value
  !! - FLU_TNUMBER  : a number
  !! - FLU_TSTRING  : a string
  !! - FLU_TTABLE   : a table
  !! - FLU_TFUNCTION: a function
  !!
  !! If none of key, pos or thandle are provided, the type of the current
  !! top of the stack will be returned. Just passing pos without a thandle
  !! is invalid and always returns FLU_TNONE.
  !!
  function aot_type_of(L, thandle, key, pos) result(luatype)
    type(flu_State) :: L !! Handle to the Lua script.

    !> Handle of the table to get the value from
    integer, intent(in), optional :: thandle

    !> Key of the value to find the type for.
    character(len=*), intent(in), optional :: key

    !> Position of the value to find the type for.
    integer, intent(in), optional :: pos

    !> Type of the Lua object found in L, thandle, key and pos
    integer :: luatype

    luatype = FLU_TNONE

    if (present(thandle)) then
      call aot_table_push( L       = L,       &
        &                  thandle = thandle, &
        &                  key     = key,     &
        &                  pos     = pos,     &
        &                  toptype = luatype  )
    else
      if (present(key)) then
        luatype = flu_getglobal(L, key)
      else if (.not. present(pos)) then
        luatype = flu_type(L, -1)
      end if
    end if

  end function aot_type_of


  !> Load the first key-value pair of table thandle on the
  !! stack.
  !!
  !! This serves as an entry point, further traversal
  !! can be done by flu_next(L, thandle).
  !! If there are no entries in the table the function
  !! returns false, otherwise the result will be true.
  function aot_table_first(L, thandle) result(exists)
    type(flu_state) :: L !! Handle for the Lua script.

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
    type(flu_state) :: L !! Handle for the Lua script.

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

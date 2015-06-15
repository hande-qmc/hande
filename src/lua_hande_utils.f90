module lua_hande_utils

! Utility procedures for working with the Lua API.

implicit none

contains

    subroutine warn_unused_args(lua_state, valid_keys, pos)

        ! Print a warning message for the keys not recognised in a table.

        ! In/Out:
        !    lua_state: lua state containing the table.
        ! In:
        !    keys: list of keys recognised in the table.
        !    pos (optional, default 1): position (ie handle) of the table in the
        !        lua stack.

        use flu_binding, only: flu_State, flu_pushnil, flu_next, flu_tolstring, flu_pop, flu_tolstring

        use, intrinsic :: iso_c_binding,  only: c_null_char

        use errors, only: warning
        use parallel, only: parent

        type(flu_State), intent(inout) :: lua_state
        character(*), intent(in) :: valid_keys(:)
        integer, intent(in), optional :: pos
        integer :: pos_loc, len, j
        character(:), allocatable :: key, key_list
        character, pointer :: str(:)

        pos_loc = 1
        if (present(pos)) pos_loc = pos

        ! Iterate through the table and print key, value pairs.
        ! See example code from the lua api manual: http://pgl.yoyo.org/luai/i/lua_next.
        call flu_pushnil(lua_state)
        do while (flu_next(lua_state, pos_loc))
            ! key is at index -2 and value at index -1
            str => flu_tolstring(lua_state, -2, len)
            allocate(character(size(str)) :: key)
            do j = 1, size(str)
                key(j:j) = str(j)
            end do
            if (all(valid_keys /= key)) then
                if (allocated(key_list)) then
                    key_list = key_list//', '//key
                else
                    key_list = key
                end if
            end if
            ! remove value, keep key for next iteration.
            call flu_pop(lua_state, 1)
            deallocate(key)
        end do

        if (allocated(key_list) .and. parent) then
            call warning('warn_unused_args', 'The following keywords are not recognised and have been ignored: '//key_list//'.')
        end if

    end subroutine warn_unused_args

    subroutine get_flag_and_id(lua_state, thandle, flag, id, key)

        ! Parse a combined boolean flag and integer field from lua, where a key
        ! is set to an boolean (which implies using the default integer value)
        ! or an integer (which implies the flag should be set to true)..

        ! In:
        !    thandle: handle of the table containing the key.
        !    key: key in the lua state to examine.
        ! In/Out:
        !    lua_state: flu/Lua state to which the HANDE API is added.
        !    id: unaltered if the value is boolean or the key doesn't exist and
        !       set to the value of the key if the value is an integer.
        !    flag: set to the value if the value is boolean, unaltered if the
        !       key does not exist in the table and set to true if the value is
        !       an integer.

        use flu_binding, only: flu_State
        use aot_table_module, only: aot_exists, aot_get_val

        type(flu_State), intent(inout) :: lua_state
        integer, intent(in) :: thandle
        logical, intent(out) :: flag
        integer, intent(inout) :: id
        character(*), intent(in) :: key

        integer :: err

        if (aot_exists(lua_state, thandle, key)) then
            call aot_get_val(flag, err, lua_state, thandle, key)
            if (err /= 0) then
                ! Passed an id instead.
                flag = .true.
                call aot_get_val(id, err, lua_state, thandle, key)
            end if
        end if

    end subroutine get_flag_and_id

end module lua_hande_utils

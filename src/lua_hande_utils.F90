module lua_hande_utils

! Utility procedures for working with the Lua API.

implicit none

integer, parameter :: lua_registryindex = -1001000

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
        integer :: pos_loc, strlen, j
        character(:), allocatable :: key, key_list
        character, pointer :: str(:)

        pos_loc = 1
        if (present(pos)) pos_loc = pos

        ! Iterate through the table and print key, value pairs.
        ! See example code from the lua api manual: http://pgl.yoyo.org/luai/i/lua_next.
        call flu_pushnil(lua_state)
        do while (flu_next(lua_state, pos_loc))
            ! key is at index -2 and value at index -1
            str => flu_tolstring(lua_state, -2, strlen)
#if ! defined(__GNUC__) || __GNUC__ > 4 || (__GNUC__ == 4 && (__GNUC_MINOR__ > 7))
            allocate(character(strlen) :: key)
#else
            allocate(character(size(str)) :: key)
#endif
            do j = 1, strlen
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
        ! or an integer (which implies the flag should be set to true).

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

    subroutine get_rng_seed(lua_state, table, rng_seed)

        ! Get the rng_seed value.  If set by user, read from table, otherwise generate one based upon the UUID and time.

        ! In/Out:
        !    lua_state: flu/Lua state to which the HANDE API is added.
        ! In:
        !    table: handle for the table (possibly) containing a 'rng_seed' value.
        ! Out:
        !    rng_seed: set to table['rng_seed'] if it exists and to a randomly-ish generated value otherwise.

        use flu_binding, only: flu_State
        use aot_table_module, only: aot_exists, aot_get_val
        use calc, only: gen_seed, GLOBAL_META

        type(flu_State), intent(inout) :: lua_state
        integer, intent(in) :: table
        integer, intent(out) :: rng_seed
        integer :: err

        if (aot_exists(lua_state, table, 'rng_seed')) then
            call aot_get_val(rng_seed, err, lua_state, table, 'rng_seed')
        else
            rng_seed = gen_seed(GLOBAL_META%uuid)
        end if

    end subroutine get_rng_seed

    subroutine get_userdata(lua_state, table, name, ptr)

        ! Check a table is the correct type (has the correct metatable) and return its userdata.

        ! In/Out:
        !    lua_state: flu/Lua state to which the HANDE API is added.
        ! In:
        !    table: handle for the table containing the userdata.
        !    name: expected name of metatable and key for userdata.
        ! Out:
        !    ptr: C pointer to userdata.

        use, intrinsic :: iso_c_binding, only: c_ptr
        use flu_binding, only: flu_State, flu_getmetatable, flu_pop
        use aot_table_module, only: aot_get_val

        use errors, only: stop_all

        type(flu_state), intent(inout) :: lua_state
        integer, intent(in) :: table
        character(*), intent(in) :: name
        type(c_ptr), intent(out) :: ptr

        integer :: ierr
        character(100) :: udata_name

        ! check metatable of object
        ierr = flu_getmetatable(lua_state, table)
        if (ierr == 0) call stop_all("get_userdata", "Problem receiving "//trim(name)//" object: no metatable.")
        call aot_get_val(udata_name, ierr, lua_state, -1, key="__metatable")
        if (trim(name) /= trim(udata_name)) call stop_all("get_userdata", trim(name)//" expected but " &
                                                          //trim(udata_name)//" recieved.")
        ! Pop metatable off stack
        call flu_pop(lua_state)

        ! If successful, return userdata.
        call aot_get_val(ptr, ierr, lua_state, table, key=trim(name))

    end subroutine get_userdata

    subroutine register_timing(lua_state, tag, time)

        ! Add information about time taken by a function to the timer object in the registry

        ! In/Out:
        !   lua_state: flu/lua state to which the HANDE API is added.
        ! In:
        !   tag: string naming the function being timed
        !   time: time taken by the function

        use flu_binding, only: flu_pushinteger, flu_settable, flu_State
        use aot_table_module, only: aot_table_length, aot_table_open, aot_table_set_val, aot_table_close

        type(flu_State), intent(inout) :: lua_state
        character(*), intent(in) :: tag
        real, intent(in) :: time

        integer :: timer, entries, handle, ierr

        call aot_table_open(lua_state, lua_registryindex, timer, "timer")

        entries = aot_table_length(lua_state, timer)
        call flu_pushinteger(lua_state, entries+1)
        call aot_table_open(lua_state, thandle=handle)
        call aot_table_set_val(tag, lua_state, handle, key="tag")
        call aot_table_set_val(time, lua_state, handle, key="time")
        call flu_settable(lua_state, timer)

        call aot_table_close(lua_state, timer)

    end subroutine register_timing

    subroutine timing_summary(lua_state)

        ! Print a breakdown of the time used by each system setup or calculation.

        ! In/Out:
        !   lua_state: flu/Lua state to which the HANDE API is added.

        use flu_binding!, only: flu_State, flu_pushnil, flu_next, flu_pop
        use aot_table_module, only: aot_table_open, aot_table_close, aot_get_val, aot_table_top

        type(flu_State), intent(inout) :: lua_state

        integer :: timer, ierr, strlen
        real :: time
!        character, pointer :: key(:)
        character(100) :: key

        call aot_table_open(lua_state, lua_registryindex, timer, "timer")
        
        ! Traversal of timer table
        call flu_pushnil(lua_state)
        do while (flu_next(lua_state, timer))
            call aot_get_val(time, ierr, lua_state, aot_table_top(lua_state), "time")
            call aot_get_val(key, ierr, lua_state, aot_table_top(lua_state), "tag")
            print *, trim(key), time
            call flu_pop(lua_state, 1)
        end do

        call aot_table_close(lua_state, timer)

    end subroutine timing_summary

end module lua_hande_utils

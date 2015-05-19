module lua_hande_fns

! Wrappers around other useful functions we wish to expose (ie those that don't 
! perform calculations nor create a system object).

implicit none

contains

    function lua_redistribute_restart(L) result(nresult) bind(c)

        ! Redistribute restart files over a different number of processors.

        ! In/Out:
        !    L: lua state (bare C pointer).

        ! Lua:
        !    redistribute {
        !       nprocs = N, -- required
        !       read = id,
        !       write = id,
        !    }

        use, intrinsic :: iso_c_binding, only: c_ptr, c_int
        use flu_binding, only: flu_State, flu_copyptr

        use aot_table_module, only: aot_table_top, aot_get_val, aot_exists

        use errors, only: stop_all
        use restart_hdf5, only: restart_info_t, init_restart_info_t, redistribute_restart_hdf5

        integer(c_int) :: nresult
        type(c_ptr), value :: L

        type(flu_State) :: lua_state

        integer :: opts, nprocs_target, read_id, write_id, err
        logical :: read_exists, write_exists
        type(restart_info_t) :: ri

        lua_state = flu_copyptr(l)
        opts = aot_table_top(lua_state)

        call aot_get_val(nprocs_target, err, lua_state, opts, 'nprocs')
        if (err /= 0) call stop_all('lua_redistribute_restart', 'nprocs: target number of processors not supplied.')

        read_exists = aot_exists(lua_state, opts, 'read')
        write_exists = aot_exists(lua_state, opts, 'write')

        if (read_exists .and. write_exists) then
            call aot_get_val(read_id, err, lua_state, opts, 'read')
            call aot_get_val(write_id, err, lua_state, opts, 'write')
            call init_restart_info_t(ri, write_id, read_id)
        else if (read_exists) then
            call aot_get_val(read_id, err, lua_state, opts, 'read')
            call init_restart_info_t(ri, read_id=read_id)
        else if (write_exists) then
            call aot_get_val(write_id, err, lua_state, opts, 'write')
            call init_restart_info_t(ri, write_id=write_id)
        else
            call init_restart_info_t(ri)
        end if

        call redistribute_restart_hdf5(ri, nprocs_target)

        nresult = 0

    end function lua_redistribute_restart

end module lua_hande_fns

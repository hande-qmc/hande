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
        !       nprocs = N,
        !       read = id,
        !       write = id,
        !       sys = system,
        !    }

        use, intrinsic :: iso_c_binding, only: c_ptr, c_int
        use flu_binding, only: flu_State, flu_copyptr

        use aot_table_module, only: aot_table_top, aot_get_val, aot_exists, aot_table_close

        use errors, only: stop_all
        use parallel, only: nprocs
        use restart_hdf5, only: restart_info_t, init_restart_info_t, redistribute_restart_hdf5
        use system, only: sys_t
        use lua_hande_system, only: get_sys_t
        use lua_hande_utils, only: warn_unused_args

        integer(c_int) :: nresult
        type(c_ptr), value :: L

        type(flu_State) :: lua_state

        integer :: opts, nprocs_target, read_id, write_id, err
        logical :: read_exists, write_exists
        type(restart_info_t) :: ri
        type(sys_t), pointer :: sys
        character(6), parameter :: keys(4) = [character(6) :: 'sys', 'read', 'write', 'nprocs']

        lua_state = flu_copyptr(l)
        opts = aot_table_top(lua_state)

        call aot_get_val(nprocs_target, err, lua_state, opts, 'nprocs', default=nprocs)

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

        call warn_unused_args(lua_state, keys, opts)
        call aot_table_close(lua_state, opts)

        if (aot_exists(lua_state, opts, 'sys')) then
            call get_sys_t(lua_state, sys)
            call redistribute_restart_hdf5(ri, nprocs_target, sys)
        else
            call redistribute_restart_hdf5(ri, nprocs_target)
        end if

        nresult = 0

    end function lua_redistribute_restart

    function lua_dump_hdf5_generic_system(L) result(nresult) bind(c)

        ! Write a read_in system to an HDF5 file, which can subsequently be used instead of an ASCII FCIDUMP file.

        ! In/Out:
        !    L: lua state (bare C pointer).

        ! Lua:
        !    dump_hdf5_system {
        !       sys = system,       -- required
        !       filename = filename,
        !    }
        ! Returns:
        !    name of HDF5 file created.  Only set on the root process and set
        !    to an empty string for all other processors.

        use, intrinsic :: iso_c_binding, only: c_ptr, c_int
        use flu_binding, only: flu_State, flu_copyptr, flu_pushstring

        use aot_table_module, only: aot_table_top, aot_get_val, aot_table_close
        use lua_hande_system, only: get_sys_t
        use lua_hande_utils, only: warn_unused_args

        use errors, only: warning
        use parallel, only: parent
        use system, only: sys_t, read_in

        use hdf5_system, only: dump_system_hdf5

        type(c_ptr), value :: L
        integer(c_int) :: nresult

        type(flu_State) :: lua_state
        type(sys_t), pointer :: sys

        integer :: opts, err
        character(255) :: filename

        character(8), parameter :: keys(2) = [character(8) :: 'sys', 'filename']

        lua_state = flu_copyptr(L)
        call get_sys_t(lua_state, sys)
        opts = aot_table_top(lua_state)
        filename = ''
        call aot_get_val(filename, err, lua_state, opts, 'filename')
        call warn_unused_args(lua_state, keys, opts)
        call aot_table_close(lua_state, opts)

        if (sys%system == read_in) then
            if (parent) call dump_system_hdf5(sys, filename)
        else
            call warning('lua_dump_hdf5_system', 'Cannot write systems other than read_in to an HDF5 file.')
        end if

        nresult = 1
        call flu_pushstring(lua_state, filename)

    end function lua_dump_hdf5_generic_system

end module lua_hande_fns

module lua_hande

! A lua interface to HANDE using the AOTUS library.

implicit none

contains

    subroutine run_lua_hande(ierr)

        ! Generic entry point from which a lua script is run.

        ! Out:
        !    ierr: error code.  Non-zero if there was any problem reading or running
        !          the lua script.

        use aotus_module, only: open_config_file
        use flu_binding, only: flu_State, fluL_newstate, flu_close
        use errors, only: stop_all, warning

        integer, intent(out) :: ierr

        character(255) :: inp_file, err_string
        type(flu_State) :: lua_state
        logical :: t_exists

        ! [todo] - MPI communication.  It is probably easiest if we just read the lua
        ! [todo] - input file into a string and broadcast that than parsing only on the
        ! [todo] - parent process and then having to figure out which function is being
        ! [todo] - called from the lua script.

        if (command_argument_count() > 0) then
            ! Input file specified on the command line.
            call get_command_argument(1, inp_file)
            inquire(file=inp_file, exist=t_exists)
            if (.not.t_exists) then
                call stop_all('read_input','File does not exist:'//trim(inp_file))
            end if

            lua_state = fluL_newstate()
            call register_lua_hande_api(lua_state)

            ! Attempt to run script.  If it fails (ie ierr is non-zero) then try
            ! parsing it as traditional input file for now.
            call open_config_file(lua_state, inp_file, ierr, err_string)
            ! [todo] - Change to stop once traditional input mode has been removed.
            if (ierr == 0) then
                call flu_close(lua_state)
            else
                call warning('run_lua_hande', trim(err_string)//'.  Assuming traditional input file...', 2)
            end if

        else
            ! Maybe read via STDIN?  Only in traditional mode for now...
            ! [todo] - STDIN functionality with lua?
            ierr = 1
        end if

    end subroutine run_lua_hande

    subroutine register_lua_hande_api(lua_state)

        ! Register HANDE functions which are exposed to Lua scripts with a Lua state
        ! so that they can be called from a Lua script loaded by the same Lua state.

        ! In/Out:
        !    lua_state: flu/Lua state to which the HANDE API is added.

        use flu_binding, only: flu_State, flu_register
        use tests, only: test_lua_api

        type(flu_State), intent(inout) :: lua_state

        call flu_register(lua_state, 'test_lua_api', test_lua_api)

    end subroutine register_lua_hande_api

end module lua_hande

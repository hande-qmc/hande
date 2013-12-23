module lua_hande

! A lua interface to HANDE using the AOTUS library.

! A very quick tutorial/summary/notes on interacting with lua:

! 1. Please read the Lua manual on the C API: http://www.lua.org/pil/24.html.
! 2. aotus has many useful wrappers around the Lua C API for interacting with Lua,
!    including transferring information from Fortran to Lua and vice versa.  See
!    lib/aotus/source for the wrapper functions.
! 3. Everything goes via the Lua stack, which is a C pointer (Fortran)/Lua_State* (C).
!    aotus wraps this into a flu_State object, which we should use wherever possible (ie
!    always).
! 4. Note that each function called by lua gets its own private stack.
! 5. The Lua stack is LIFO (last in, first out), so the first argument passed is at the
!    bottom of the stack and the last argument at the top.  Similarly the first argument
!    returned is at the bottom of the stack and the last argument at the top.
! 6. A function can be called from lua if it is 'registered' with the Lua stack and has
!    a very specific interface, e.g.:
!
!        function test_lua_api(L) result(nreturn) bind(c)
!
!            ! Example function callable from a Lua script.
!            ! Must be bind(C) so it can be called from C/Lua.
!            ! Must take a C pointer (which is a Lua stack) as the sole argument.
!            ! Must return an integer which is the number of variables the
!            ! function returns to Lua by pushing to the stack.
!
!            ! In/Out:
!            !    L: lua state (bare C pointer).
!
!            use flu_binding, only: flu_State, flu_copyptr
!            use, intrinsic :: iso_c_binding, only: c_ptr, c_int
!
!            integer(c_int) :: nreturn
!            type(c_ptr), value :: L
!
!            type(flu_State) :: lua_state
!
!            ! Create flu_state from lua state so we can use the nice bindings
!            ! provided by AOTUS.
!            lua_state = flu_copyptr(L)
!
!            ! Number of variables returned on lua stack.
!            nreturn = 0
!
!            ! Get arguments passed to us by lua (if appropriate).
!
!            ! Now do our work...
!            write (6,'(1X,"Hello from fortran!",/)')
!
!        end function test_lua_api
!
! 7. The function can the be called in lua using function_name(arg1, arg2, ...).  In the
!    example above, there are no arguments so it's simply test_lua_api(), assuming
!    test_lua_api is the name provided to lua in the flu_register call.
! 8. Please read the lua documentation and see the examples in the lua_hande* files for
!    more details!

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

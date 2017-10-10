module tests

! Handy boilerplate module to use for testing and debugging code.

! Some rules:
! * Do commit test procedures which are generic or frequently useful.
! * Don't commit (at least not permanently) code that tests development work.
!   If you do commit such code, please tidy up and remove it once it stops being
!   useful.
! * Don't keep test calls in production code, which should be as lean and as
!   fast as possible.
! * Keep the modules on which this module depends to a minimum to avoid
!   potential circular dependency problems.
! * Don't assume old code in here is compatible with the existing code---caveat
!   emptor!

use const

implicit none

contains

    subroutine assert_statement(test)

        ! Exit if a test is not met.
        ! If this proves useful, then it might best belong in the utils module.

        ! In:
        !    test: logical statement which ought(!) to be true.

        use errors, only: stop_all

        logical, intent(in) :: test

        if (.not. test) call stop_all('assert_statement','assertion is false', print_backtrace=.true.)

    end subroutine assert_statement

    function test_lua_api(L) result(nreturn) bind(c)

        ! Example function callable from a Lua script.

        ! In/Out:
        !    L: lua state (bare C pointer).

        use flu_binding, only: flu_State, flu_copyptr
        use iso_c_binding, only: c_ptr, c_int
        use parallel, only: iproc

        integer(c_int) :: nreturn
        type(c_ptr), value :: L

        type(flu_State) :: lua_state

        ! Number of variables returned on Lua stack.
        nreturn = 0

        ! Create flu_state from lua state so we can use the nice bindings
        ! provided by AOTUS.
        lua_state = flu_copyptr(L)

        ! Get arguments passed to us by lua (if appropriate).

        ! Now do our work...
        write (6,'(1X,"Hello from fortran! Processor: ",i2,/)') iproc

    end function test_lua_api

end module tests

!> This module provides a Fortran interface to the Lua dump routine.
module dump_lua_fif_module
  use, intrinsic :: iso_c_binding

  implicit none

  interface
    function dump_lua_toBuf(L, length, ierr) &
      &        bind(c, name='dump_lua_toBuf')
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: L
      integer(kind=c_int) :: length
      integer(kind=c_int) :: ierr
      type(c_ptr) :: dump_lua_toBuf
    end function dump_lua_toBuf
  end interface

end module dump_lua_fif_module

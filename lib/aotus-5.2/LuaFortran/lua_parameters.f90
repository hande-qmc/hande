!> This module provides some parameters defined in the
!! Lua header file that are needed in the wrapper
!! functions for the Fortran interface.
!!
!! Lua_int and Lua_num are system dependent, and you
!! might need to adapt them on your system.
!! The type constants have to be consistent with the
!! Lua header definition, and thus should be checked
!! after version upgrades of the Lua library.
module lua_parameters
  use, intrinsic :: iso_c_binding

  implicit none

  ! System dependent, might need to be adapted:
  integer, parameter :: lua_int = c_long
  integer, parameter :: lua_num = c_double

  ! Lua constants
  integer(kind=c_int), parameter :: LUA_TNONE = -1
  integer(kind=c_int), parameter :: LUA_TNIL = 0
  integer(kind=c_int), parameter :: LUA_TBOOLEAN = 1
  integer(kind=c_int), parameter :: LUA_TLIGHTUSERDATA = 2
  integer(kind=c_int), parameter :: LUA_TTABLE = 5
  integer(kind=c_int), parameter :: LUA_TFUNCTION = 6

end module lua_parameters

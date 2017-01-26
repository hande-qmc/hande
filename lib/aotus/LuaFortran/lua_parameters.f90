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

  ! Lua config constants (see luaconf.h)
  ! Attention: might need to be adapted!
  integer(kind=c_int), parameter :: LUAI_MAXSTACK = 1000000

  ! Lua constants (see lua.h)
  integer(kind=c_int), parameter :: LUA_TNONE          = -1
  integer(kind=c_int), parameter :: LUA_TNIL           =  0
  integer(kind=c_int), parameter :: LUA_TBOOLEAN       =  1
  integer(kind=c_int), parameter :: LUA_TLIGHTUSERDATA =  2
  integer(kind=c_int), parameter :: LUA_TNUMBER        =  3
  integer(kind=c_int), parameter :: LUA_TSTRING        =  4
  integer(kind=c_int), parameter :: LUA_TTABLE         =  5
  integer(kind=c_int), parameter :: LUA_TFUNCTION      =  6
  integer(kind=c_int), parameter :: LUA_TUSERDATA      =  7
  integer(kind=c_int), parameter :: LUA_TTHREAD        =  8

  integer(kind=c_int), parameter :: LUA_REGISTRYINDEX  = -LUAI_MAXSTACK - 1000

end module lua_parameters

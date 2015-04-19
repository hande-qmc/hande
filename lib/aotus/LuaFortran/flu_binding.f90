!> This module provides the Fortran Lua interface.
!!
!! It defines a flu_state which encapsulates the
!! Lua state and is used to reference a Lua script.
!! The main content are then the wrapper implementations
!! which ease the usage of the Lua functions declared
!! in the lua_fif module.
module flu_binding
  use, intrinsic :: iso_c_binding
  use lua_fif
  use dump_lua_fif_module

  implicit none

  private

  type :: flu_State
    private
    type(c_ptr) :: state = c_null_ptr
    logical :: opened_libs = .false.
  end type flu_State

  type :: cbuf_type
    type(c_ptr) :: ptr = c_null_ptr
    character, pointer :: buffer(:) => NULL()
  end type cbuf_type

  public :: flu_State
  public :: cbuf_type

  public :: flu_close, flu_isopen
  public :: flu_createTable
  public :: flu_getField, flu_getGlobal, flu_getTable, flu_getTop
  public :: flu_setGlobal
  public :: flu_insert
  public :: flu_isFunction, flu_isNumber, flu_isTable
  public :: flu_isNone, flu_isNoneOrNil, flu_isNil
  public :: flu_isBoolean, flu_islightuserdata
  public :: flu_pcall
  public :: flu_next
  public :: flu_setTop
  public :: flu_setTable, flu_setField
  public :: flu_todouble
  public :: flu_tolstring, flu_tonumber, flu_toboolean, flu_touserdata
  public :: flu_pop
  public :: flu_pushinteger, flu_pushnil, flu_pushnumber, flu_pushboolean
  public :: flu_pushstring, flu_pushvalue, flu_pushlightuserdata

  public :: flu_copyptr
  public :: flu_register

  public :: flu_dump
  public :: flu_free_cbuf

  public :: fluL_loadfile, fluL_newstate, fluL_openlibs, fluL_loadstring
  public :: fluL_loadbuffer

  interface flu_pushnumber
    module procedure flu_pushreal
    module procedure flu_pushdouble
  end interface flu_pushnumber

  interface flu_dump
    module procedure flu_dump_toBuf
  end interface flu_dump


  ! Interoperable interface required for a function that is callable from Lua.
  abstract interface
    function lua_Function(s) result(val) bind(c)
      use, intrinsic :: iso_c_binding
      integer(c_int) :: val 
      type(c_ptr), value :: s
    end function lua_Function
  end interface

  interface
    subroutine c_free(ptr) bind(c, name="free")
      use, intrinsic :: iso_c_binding, only: c_ptr
      type(c_ptr), value :: ptr
    end subroutine c_free
  end interface


contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Wrapper routines for the lua API
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine flu_close(L)
    type(flu_State) :: L

    call lua_close(L%state)
    L%state = c_null_ptr
  end subroutine flu_close


  subroutine flu_createtable(L, narr, nrec)
    type(flu_State) :: L
    integer :: narr
    integer :: nrec

    integer(kind=c_int) :: c_narr
    integer(kind=c_int) :: c_nrec

    c_narr = narr
    c_nrec = nrec

    call lua_createtable(L%state, c_narr, c_nrec)
  end subroutine flu_createtable


  subroutine flu_getfield(L, index, k)
    type(flu_State)  :: L
    integer          :: index
    character(len=*) :: k

    integer(kind=c_int) :: c_index
    character(len=len_trim(k)+1) :: c_k

    c_k = trim(k) // c_null_char
    c_index = index
    call lua_getfield(L%state, c_index, c_k)
  end subroutine flu_getfield


  subroutine flu_getglobal(L, k)
    type(flu_State)  :: L
    character(len=*) :: k

    character(len=len_trim(k)+1) :: c_k

    c_k = trim(k) // c_null_char
    call lua_getglobal(L%state, c_k)
  end subroutine flu_getglobal


  subroutine flu_gettable(L, index)
    type(flu_State) :: L
    integer :: index

    integer(kind=c_int) :: c_index

    c_index = index
    call lua_gettable(L%state, c_index)
  end subroutine flu_gettable


  function flu_gettop(L) result(stacktop)
    type(flu_state) :: L
    integer :: stacktop

    stacktop = int(lua_gettop(L%state), kind=kind(stacktop))
  end function flu_gettop


  subroutine flu_insert(L, index)
    type(flu_state) :: L
    integer :: index

    integer(kind=c_int) :: c_index

    c_index = int(index, kind = c_int)
    call lua_insert(L%state, c_index)
  end subroutine flu_insert


  function flu_isBoolean(L, index) result(is_boolean)
    type(flu_State) :: L
    integer         :: index
    logical         :: is_boolean

    integer(kind=c_int) :: c_index

    c_index = int(index, kind = c_int)
    !! Only defined as a Macro, using lua_type:
    is_boolean = (lua_type(L%state, c_index) == LUA_TBOOLEAN)
  end function flu_isBoolean


  function flu_isFunction(L, index) result(is_function)
    type(flu_State) :: L
    integer         :: index
    logical         :: is_function

    integer(kind=c_int) :: c_index

    c_index = int(index, kind = c_int)
    is_function = (lua_type(L%state, c_index) == LUA_TFUNCTION)
  end function flu_isFunction


  function flu_isnumber(L, index) result(is_number)
    type(flu_State) :: L
    integer         :: index
    logical         :: is_number

    integer(kind=c_int) :: c_index

    c_index = int(index, kind = c_int)
    is_number = (lua_isnumber(L%state, c_index) .eq. 1)
  end function flu_isnumber


  function flu_isTable(L, index) result(is_Table)
    type(flu_State) :: L
    integer         :: index
    logical         :: is_Table

    integer(kind=c_int) :: c_index

    c_index = int(index, kind = c_int)
    !! Only defined as a Macro, using lua_type:
    is_Table = (lua_type(L%state, c_index) == LUA_TTABLE)
  end function flu_isTable


  function flu_isNoneOrNil(L, index) result(is_NoneOrNil)
    type(flu_State) :: L
    integer         :: index
    logical         :: is_NoneOrNil

    integer(kind=c_int) :: c_index

    c_index = int(index, kind = c_int)
    !! Only defined as a Macro, using lua_type:
    is_NoneOrNil = (lua_Type(L%state, c_index) <= 0)
  end function flu_isNoneOrNil


  function flu_isNil(L, index) result(is_Nil)
    type(flu_State) :: L
    integer         :: index
    logical         :: is_Nil

    integer(kind=c_int) :: c_index

    c_index = int(index, kind = c_int)
    !! Only defined as a Macro, using lua_type:
    is_Nil = (lua_Type(L%state, c_index) .eq. LUA_TNIL)
  end function flu_isNil


  function flu_isNone(L, index) result(is_None)
    type(flu_State) :: L
    integer         :: index
    logical         :: is_None

    integer(kind=c_int) :: c_index

    c_index = int(index, kind = c_int)
    !! Only defined as a Macro, using lua_type:
    is_None = (lua_Type(L%state, c_index) .eq. LUA_TNONE)
  end function flu_isNone


  function flu_islightuserdata(L, index) result(is_lightuserdata)
    type(flu_State) :: L
    integer         :: index
    logical         :: is_lightuserdata

    integer(kind=c_int) :: c_index

    c_index = int(index, kind = c_int)
    is_lightuserdata = (lua_Type(L%state, c_index) .eq. LUA_TLIGHTUSERDATA)
  end function flu_islightuserdata


  function flu_next(L, index) result(exists)
    type(flu_State) :: L
    integer, intent(in) :: index
    logical :: exists

    integer(kind=c_int) :: retCode

    retCode = lua_next(L%state, index)
    exists = (retCode /= 0)
  end function flu_next


  function flu_pcall(L, nargs, nresults, errfunc) result(errcode)
    type(flu_State) :: L
    integer :: nargs
    integer :: nresults
    integer :: errfunc
    integer :: errcode

    integer(kind=c_int) :: c_nargs
    integer(kind=c_int) :: c_nresults
    integer(kind=c_int) :: c_errfunc
    integer(kind=c_int) :: c_errcode

    c_nargs = nargs
    c_nresults = nresults
    c_errfunc = errfunc

    c_errcode = lua_pcallk(L%state, c_nargs, c_nresults, c_errfunc, &
      &                    0_c_int, C_NULL_PTR)
    errcode = c_errcode
  end function flu_pcall


  subroutine flu_pop(L, n)
    type(flu_State) :: L
    integer, optional, intent(in) :: n

    integer(kind=c_int) :: n_c

    n_c = -2
    if (present(n)) n_c = -n-1
    call lua_settop(L%state, n_c)
  end subroutine flu_pop


  subroutine flu_pushinteger(L, n)
    type(flu_State) :: L
    integer :: n

    integer(kind=lua_int) :: n_c

    n_c = int(n, lua_int)
    call lua_pushinteger(L%state, n_c)
  end subroutine flu_pushinteger

  subroutine flu_pushboolean(L, b)
    type(flu_State) :: L
    logical :: b

    integer(kind=lua_int) :: n_c

    if (b) then
      n_c = 1_lua_int
    else
      n_c = 0_lua_int
    end if
    call lua_pushboolean(L%state, n_c)
  end subroutine flu_pushboolean

  subroutine flu_pushstring(L, string)
    type(flu_State) :: L
    character(len=*), intent(in) :: string

    integer(kind=c_size_t) :: c_len
    type(c_ptr) :: ret

    c_len = len(string)
    ret = lua_pushlstring(L%state, string, c_len)
  end subroutine flu_pushstring

  subroutine flu_pushreal(L, n)
    type(flu_State) :: L
    real :: n

    real(kind=c_double) :: n_c

    n_c = real(n, c_double)
    call lua_pushnumber(L%state, n_c)
  end subroutine flu_pushreal

  subroutine flu_pushdouble(L, n)
    type(flu_State) :: L
    real(kind=c_double) :: n

    call lua_pushnumber(L%state, n)
  end subroutine flu_pushdouble

  subroutine flu_pushnil(L)
    type(flu_State) :: L

    call lua_pushnil(L%state)
  end subroutine flu_pushnil

  subroutine flu_pushvalue(L, index)
    type(flu_State) :: L
    integer :: index

    integer(kind=c_int) :: c_index

    c_index = index
    call lua_pushvalue(L%state, c_index)
  end subroutine flu_pushvalue

  subroutine flu_pushlightuserdata(L, ptr)
    type(flu_State) :: L
    type(c_ptr) :: ptr

    call lua_pushlightuserdata(L%state, ptr)

  end subroutine flu_pushlightuserdata


  subroutine flu_settable(L, n)
    type(flu_State) :: L
    integer, intent(in) :: n

    integer(kind=c_int) :: n_c

    n_c = n
    call lua_settable(L%state, n_c)
  end subroutine flu_settable


  subroutine flu_settop(L, n)
    type(flu_State) :: L
    integer, intent(in) :: n

    integer(kind=c_int) :: n_c

    n_c = n
    call lua_settop(L%state, n_c)
  end subroutine flu_settop


  subroutine flu_setfield(L, index, k)
    type(flu_State)  :: L
    integer          :: index
    character(len=*) :: k

    integer(kind=c_int) :: c_index
    character(len=len_trim(k)+1) :: c_k

    c_k = trim(k) // c_null_char
    c_index = index
    call lua_setfield(L%state, c_index, c_k)
  end subroutine flu_setfield


  subroutine flu_setglobal(L, k)
      type(flu_State) :: L
      character(len=*), intent(in) :: k

      character(len=len_trim(k)+1) :: c_k

      c_k = trim(k) // c_null_char

      call lua_setglobal(L%state, c_k)

  end subroutine flu_setglobal


  function flu_tolstring(L, index, len) result(string)
    type(flu_State) :: L
    integer :: index
    integer :: len
    character,pointer,dimension(:) :: string

    integer :: string_shape(1)
    integer(kind=c_int) :: c_index
    integer(kind=c_size_t) :: c_len
    type(c_ptr) :: c_string

    c_index = index
    c_string = lua_tolstring(L%state, c_index, c_len)
    len = int(c_len,kind=kind(len))
    string_shape(1) = len
    call c_f_pointer(c_string, string, string_shape)
  end function flu_tolstring


  function flu_todouble(L, index) result(number)
    type(flu_State) :: L
    integer :: index
    real(kind=c_double) :: number

    integer(kind=c_int) :: c_index
    integer(kind=c_int) :: isnum

    c_index = index
    number = lua_tonumberx(L%state, c_index, isnum)
  end function flu_todouble


  function flu_tonumber(L, index) result(number)
    type(flu_State) :: L
    integer :: index
    real :: number

    integer(kind=c_int) :: c_index
    integer(kind=c_int) :: isnum

    c_index = index
    number = real(lua_tonumberx(L%state, c_index, isnum), &
      &           kind=kind(number))
  end function flu_tonumber


  function flu_toBoolean(L, index) result(bool)
    type(flu_State) :: L
    integer :: index
    logical :: bool

    integer(kind=c_int) :: c_index

    c_index = index
    bool = (lua_toBoolean(L%state, c_index) == 1)
  end function flu_toBoolean


  function flu_touserdata(L, index) result(ptr)
    type(flu_State) :: L
    integer :: index
    type(c_ptr) :: ptr

    integer(kind=c_int) :: c_index

    c_index = index
    ptr = lua_touserdata(L%state, c_index)
  end function flu_touserdata

  subroutine flu_pushcclosure(L, fn, n)
    type(flu_State), value :: L 
    procedure(lua_Function) :: fn
    integer :: n 

    integer(c_int) :: c_n
    type(c_funptr) :: c_fn

    c_n = n 
    c_fn = c_funloc(fn)

    call lua_pushcclosure(L%state, c_fn, c_n)

  end subroutine flu_pushcclosure

  subroutine flu_register(L, fn_name, fn) 

    ! lua_register is defined as a macro in lua.h and isn't accessible from Fortran.
    ! Re-implement macro explictly.

    type(flu_State) :: L
    character(len=*), intent(in) :: fn_name
    procedure(lua_Function) :: fn

    call flu_pushcclosure(L, fn, 0)
    call flu_setglobal(L, fn_name)

  end subroutine flu_register

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Wrapper routines for the auxiliary library 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function fluL_loadfile(L, filename) result(errcode)
    type(flu_State) :: L
    character(len=*) :: filename
    integer :: errcode

    character(len=len_trim(filename)+1) :: c_filename
    character(len=3) :: c_mode
    integer(kind=c_int) :: c_errcode

    c_filename = trim(filename) // c_null_char
    c_mode = "bt" // c_null_char
    c_errcode = luaL_loadfilex(L%state, c_filename, c_mode)
    errcode = c_errcode
  end function fluL_loadfile


  function fluL_loadbuffer(L, buffer, bufName) result(errcode)
    type(flu_State) :: L
    character :: buffer(:)
    character(len=*), optional :: bufName
    integer :: errcode

    character(len=33) :: label
    character(len=3) :: c_mode
    integer(kind=c_int) :: c_errcode
    integer(kind=c_size_t) :: nChars

    if (present(bufName)) then
      label = trim(bufName) // c_null_char
    else
      label = 'ScriptBuffer' // c_null_char
    end if
    nChars = int(size(buffer),kind=kind(nChars))
    c_mode = "bt" // c_null_char
    c_errcode = luaL_loadbufferx(L%state, buffer, nChars, label, c_mode)
    errcode = c_errcode
  end function fluL_loadbuffer


  function fluL_loadstring(L, string) result(errcode)
    type(flu_State) :: L
    character(len=*) :: string
    integer :: errcode

    character(len=len_trim(string)+1) :: c_string
    integer(kind=c_int) :: c_errcode

    c_string = trim(string) // c_null_char
    c_errcode = luaL_loadstring(L%state, c_string)
    errcode = c_errcode

  end function fluL_loadstring


  function fluL_newstate() result(new_state)
    type(flu_State) :: new_state

    new_state%state = luaL_newstate()
  end function fluL_newstate


  subroutine fluL_openlibs(L)
    type(flu_State) :: L
    
    if (.not. L%opened_libs) then
      call luaL_openlibs(L%state)
      L%opened_libs = .true.
    end if
  end subroutine fluL_openlibs

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Routines for using existing Lua states with 
  ! flu_binding
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function flu_copyptr(lua_state) result(L)
      ! WARNING: this copies the pointer to an existing Lua state, not the Lua
      ! state itself.  Modifying L via the flu bindings will modify the same Lua
      ! state as pointed to by lua_state.
      type(flu_State) :: L
      type(c_ptr), intent(in) :: lua_state
      L%state = lua_state
  end function flu_copyptr

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Routines for probing the Lua state
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function flu_isopen(L) result(is_open)
      logical :: is_open
      type(flu_State), intent(in) :: L

      is_open = c_associated(L%state)
  end function flu_isopen


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! Wrapper implementation for lua_dump !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Dump to a buffer and return the pointer to the resulting string.
  subroutine flu_dump_toBuf(L, buf, length, iError)
    type(flu_State) :: L
    type(cbuf_type), intent(out) :: buf
    integer :: length
    integer :: iError

    type(c_ptr) :: string_c
    integer(kind=c_int) :: length_c
    integer(kind=c_int) :: iErr

    string_c = dump_lua_toBuf(L%state, length_c, iErr)
    iError = int(iErr)
    if (iError == 0) then
      length = int(length_c)
      buf%ptr = string_c
      call c_f_pointer(string_c, buf%buffer, [length])
    else
      length = 0
    end if
  end subroutine flu_dump_toBuf


  !> Free an allocated cbuf.
  !!
  !! This is a helping routine to deallocate memory that was allocated for
  !! the cbuf by C. (Cray compiler complained about its deallocation in Fortran)
  subroutine flu_free_cbuf(buf)
    type(cbuf_type) :: buf

    call c_free(buf%ptr)
    nullify(buf%buffer)
  end subroutine flu_free_cbuf

end module flu_binding

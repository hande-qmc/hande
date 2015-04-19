! Copyright (C) 2011-2013 German Research School for Simulation Sciences GmbH,
!                         Aachen and others.
!               2013-2014 University of Siegen.
! Please see the COPYRIGHT file in this directory for details.

!> This module provides high level Fortran interfaces to retrieve values from a
!! Lua script.
!!
!! The central interface of the module is aot_table_module#aot_get_val, which is
!! a generic interface that allows access to scalars and vectors in global Lua
!! variables as well as nested tables.
!!
!! In the \ref aot_overview "overview page" there are some more general
!! remarks and further pointers.
module aotus_module
  use flu_binding
  use aot_kinds_module, only: double_k, single_k, long_k
  use aot_top_module, only: aot_top_get_val, aot_err_handler, &
    &                       aoterr_Fatal, aoterr_NonExistent, aoterr_WrongType
  use aot_table_module, only: aot_get_val, aot_table_set_val, &
    &                         aot_table_open, aot_table_close
  use aot_vector_module, only: aot_top_get_val, aot_get_val

  implicit none

  private

  public :: aot_get_val
  public :: open_config_file, close_config
  public :: open_config_chunk, open_config_buffer
  public :: aot_require_buffer
  public :: aot_file_to_buffer

  ! Entities inherited from aot_top_module, published here to
  ! allow most functionality by "use aotus_module".
  public :: aoterr_Fatal, aoterr_NonExistent, aoterr_WrongType
  public :: aot_err_handler
  public :: aot_top_get_val

  ! Inherited from the flu_binding module, publish for convenience.
  public :: flu_State

contains

  !> Subroutine to load and execute a script from a file.
  !!
  !! If you are using MPI for parallelization, have a look at the
  !! tem_open_distconf routine in the
  !! [treelm library](https://bitbucket.org/apesteam/treelm) instead.
  subroutine open_config_file(L, filename, ErrCode, ErrString, buffer)
    type(flu_State) :: L !< Handle to the Lua script

    !> Name of file to load the Lua code from
    character(len=*), intent(in) :: filename

    !> Error code returned by Lua during loading or executing the file.
    !!
    !! This optional parameter might be used to react on errors in the calling
    !! side. If neither ErrCode nor ErrString are given, this subroutine will
    !! stop the program execution and print the error message from Lua to the
    !! stdout.
    integer, intent(out), optional :: ErrCode

    !> Obtained error description from the Lua stack.
    !!
    !! This optional argument holds the Lua error message in case somehting
    !! went wrong. It can be used to provide some feedback to the user in the
    !! calling routine. If neither ErrCode nor ErrString are provided,
    !! open_config() will print the error message and stop program execution.
    character(len=*), intent(out), optional :: ErrString

    !> Optional argument to return the compiled script after loading it to
    !! the caller.
    !!
    !! It might be handy to reuse the loaded script later on, this argument
    !! allows you to obtain the script in compiled form, before it is executed.
    !! The buffer will be allocated and filled with the Lua data.
    !! It contains the actual string in buffer%buffer which is a character
    !! pointer, and the original c_ptr to this
    type(cbuf_type), intent(out), optional :: buffer

    integer :: err
    integer :: length

    if (.not.flu_isopen(L)) L = fluL_newstate()

    err = fluL_loadfile(L, filename)
    call aot_err_handler(L, err, 'Cannot load configuration file:', ErrString, &
      &                  ErrCode)

    if (err == 0) then
      if (present(buffer)) then
        call flu_dump(L, buffer, length, err)
      end if
      call fluL_openlibs(L)

      err = flu_pcall(L, 0, 0, 0)
      call aot_err_handler(L, err, 'Cannot run configuration file:',  &
        &                  ErrString, ErrCode)
    end if

  end subroutine open_config_file


  !> Subroutine to load and execute a script given in a string.
  subroutine open_config_chunk(L, chunk, ErrCode, ErrString)
    type(flu_State) :: L !< Handle to the Lua script

    !> String with Lua code to load.
    character(len=*), intent(in) :: chunk

    !> Error code returned by Lua during loading or executing the file.
    !!
    !! This optional parameter might be used to react on errors in the calling
    !! side. If neither ErrCode nor ErrString are given, this subroutine will
    !! stop the program execution and print the error message from Lua to the
    !! stdout.
    integer, intent(out), optional :: ErrCode

    !> Obtained error description from the Lua stack.
    !!
    !! This optional argument holds the Lua error message in case somehting
    !! went wrong. It can be used to provide some feedback to the user in the
    !! calling routine. If neither ErrCode nor ErrString are provided,
    !! open_config() will print the error message and stop program execution.
    character(len=*), intent(out), optional :: ErrString

    integer :: err

    if (.not.flu_isopen(L)) L = fluL_newstate()

    err = fluL_loadstring(L, chunk)

    call aot_err_handler(L, err, 'Cannot load chunk:', ErrString, ErrCode)

    if (err == 0) then
      call fluL_openlibs(L)

      err = flu_pcall(L, 0, 0, 0)

      call aot_err_handler(L, err, 'Cannot run chunk:', ErrString, ErrCode)
    end if

  end subroutine open_config_chunk


  !> Subroutine to load and execute a script given in a buffer
  !! (bytecode).
  subroutine open_config_buffer(L, buffer, bufName, ErrCode, ErrString)
    type(flu_State) :: L !< Handle to the Lua script

    !> String with Lua code to load.
    character, intent(in) :: buffer(:)

    !> Name for the buffer to use in debug messages.
    character(len=*), intent(in), optional :: bufName

    !> Error code returned by Lua during loading or executing the file.
    !!
    !! This optional parameter might be used to react on errors in the calling
    !! side. If neither ErrCode nor ErrString are given, this subroutine will
    !! stop the program execution and print the error message from Lua to the
    !! stdout.
    integer, intent(out), optional :: ErrCode

    !> Obtained error description from the Lua stack.
    !!
    !! This optional argument holds the Lua error message in case somehting
    !! went wrong. It can be used to provide some feedback to the user in the
    !! calling routine. If neither ErrCode nor ErrString are provided,
    !! open_config() will print the error message and stop program execution.
    character(len=*), intent(out), optional :: ErrString

    integer :: err

    if (.not.flu_isopen(L)) L = fluL_newstate()

    err = fluL_loadbuffer(L, buffer, bufName)

    call aot_err_handler(L, err, 'Cannot load buffer:', ErrString, ErrCode)

    if (err == 0) then
      call fluL_openlibs(L)

      err = flu_pcall(L, 0, 0, 0)

      call aot_err_handler(L, err, 'Cannot run buffer:', ErrString, ErrCode)
    end if

  end subroutine open_config_buffer


  !> Close an opened Lua script again.
  subroutine close_config(L)
    type(flu_State) :: L !< Handle to the Lua script to close.

    call flu_close(L)

  end subroutine close_config


  !> Subroutine to load a script from a file and put it into a character buffer.
  !!
  !! This is useful to rerun a given code in a file without the need to touch
  !! the file itself again.
  subroutine aot_file_to_buffer(filename, buffer, ErrCode, ErrString)
    !> Name of file to load the Lua code from
    character(len=*), intent(in) :: filename

    !> Buffer to store the script in the given file in
    type(cbuf_type), intent(out) :: buffer

    !> Error code returned by Lua during loading or executing the file.
    !!
    !! This optional parameter might be used to react on errors in the calling
    !! side. If neither ErrCode nor ErrString are given, this subroutine will
    !! stop the program execution and print the error message from Lua to the
    !! stdout.
    integer, intent(out), optional :: ErrCode

    !> Obtained error description from the Lua stack.
    !!
    !! This optional argument holds the Lua error message in case somehting
    !! went wrong. It can be used to provide some feedback to the user in the
    !! calling routine. If neither ErrCode nor ErrString are provided,
    !! open_config() will print the error message and stop program execution.
    character(len=*), intent(out), optional :: ErrString

    type(flu_State) :: L
    integer :: err
    integer :: buflen

    L = fluL_newstate()

    err = fluL_loadfile(L, filename)
    call aot_err_handler(L, err, 'Cannot load configuration file:', ErrString, &
      &                  ErrCode)

    if (err == 0) then

      call flu_dump(L = L, buf = buffer, length = buflen, iError = err)

      if (err /= 0) then
        if (present(ErrCode)) then
           ErrCode = err
           if (present(ErrString)) then
             ErrString = 'Error while dumping the Lua script into a buffer!'
           end if
        else
          write(*,*) 'Error while dumping the Lua script into a buffer!'
          write(*,*) 'STOPPING'
          STOP
        end if
      end if

    end if

    call close_config(L)

  end subroutine aot_file_to_buffer


  !> Load and execute a given buffer and register it in the package table as
  !! the given module name.
  subroutine aot_require_buffer(L, buffer, modname)
    type(flu_State) :: L !< Lua State to set load the buffer into.
    character, intent(in) :: buffer(:) !< Buffer to load.
    character(len=*), intent(in) :: modname !< Module name to set.

    integer :: pac_handle
    integer :: ld_handle

    call open_config_buffer(L = L, buffer = buffer, bufName = trim(modname))
    call aot_table_open(L, thandle = pac_handle, key = "package")
    call aot_table_open(L, parent = pac_handle, &
      &                 thandle = ld_handle, key = "loaded")
    call aot_table_set_val(val = .true., L = L, thandle = ld_handle, &
      &                    key = trim(modname))
    call aot_table_close(L, ld_handle)
    call aot_table_close(L, pac_handle)
  end subroutine aot_require_buffer

end module aotus_module

!> \page aot_overview Overview for Aotus
!!
!! Aotus stands for *Advanced Options in Tables and Universal Scripting*.
!!
!! It is a Fortran wrapper for the [Lua](http://www.lua.org/) scripting
!! language.
!! The aim of this wrapper is to provide flexible configuration files to Fortran
!! applications with the full user experience provided by Lua.
!! Aotus is also known as the
!! [night monkey](http://en.wikipedia.org/wiki/Night_monkey), living in south
!! america.
!! Thus, we saw the name fit as it interacts with the moon (Lua, provided
!! by the Pontifical Catholic University of Rio de Janeiro in Brazil).
!!
!! The most prominent data structure in Lua are
!! [tables](http://www.lua.org/manual/5.2/manual.html#2), which provide the
!! possibility to store complex structures.
!! Configuration is typically done with global variables in the Lua script
!! or some tables to distinguish subsections in the configuration.
!!
!! Aotus provides several layers, encapsulating the bare
!! [Lua C-API](http://www.lua.org/manual/5.2/manual.html#4):
!! - \ref lua_fif this just provides the
!!   [ISO_C_Binding](http://www.fortran.bcs.org/2002/interop.htm)
!!   interface declarations.
!! - \ref flu_binding this the actural Fortran binding wrapped around lua_fif,
!!   to provide a more Fortran like interface.
!!   Especially the flu_binding::flu_state type is declared which maintains the
!!   handle for the
!!   [Lua state](http://www.lua.org/manual/5.2/manual.html#lua_state).
!! - \ref aot_table_module provides some convenience functions to work on Lua
!!   tables in Fortran.
!! - \ref aot_fun_module provides some convenience functions to work with Lua
!!   functions in Fortran.
!! - \ref aotus_module provides the high end level to easily retrieve data from
!!   a Lua script.
!! - On top of those there is an additional \ref aot_vector_module, which allows
!!   the concise reading of values into arrays of rank one.
!! - Finally there is an \ref aot_out_module that allows output of Fortran values
!!   into nested Lua tables.
!!
!! The library can be compiled by various modern Fortran compilers as described
!! in \ref compiler_support "Compiler Support".
!!
!! An example showing the usage of the library in a Fortran application is given
!! in sample/aotus_sample.f90 in the Aotus main directory.
!! The corresponding Lua script used as input is given in sample/config.lua.
!!
!! Note on usage in parallel environments: Aotus itself is not providing parallel
!! facilities. But it can be nicely used in parallel aswell. However, for
!! massively parallel systems, it is advisable to minimize the access to config
!! files. To avoid excessive filesystem meta accesses it is recommended to
!! load required files only on one process.
!! An implementation of this for MPI can be found in TreElMs
!! [distconf](https://geb.sts.nt.uni-siegen.de/doxy/treelm/classtem__aux__module.html#a1a6bd9f747c89e6f00791131e3d169de).
!!
!! *Sources are available at <https://bitbucket.org/apesteam/aotus>.*

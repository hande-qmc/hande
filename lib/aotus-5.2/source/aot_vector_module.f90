! Copyright (C) 2011-2013 German Research School for Simulation Sciences GmbH,
!                         Aachen and others.
!               2013-2014 University of Siegen.
! Please see the COPYRIGHT file in this directory for details.

!> This module provides some convenience functions to access complete vectors
!! from a lua table at once.
!!
!! It provides two generic interfaces, one for vectors inside tables, and one
!! for vectors defined as global variables (get_config_val).
!! Vectors might be accessed with a variable length, to be defined by the
!! Lua table and allocated in the get_ routines or with a fixed length.
!! For the variable length vectors, a maximal length has to be provided
!! up to which the vector might be allocated.
!! Otherwise the interfaces correspond to the scalar retrieval operations.
module aot_vector_module
  use flu_binding
  use aot_kinds_module, only: double_k, single_k, long_k
  use aot_table_ops_module, only: aot_table_close, aot_table_top, &
    &                             aot_table_length, aot_table_push, &
    &                             aot_table_first
  use aot_top_module, only: aot_top_get_val, aoterr_NonExistent, aoterr_Fatal

  ! The following module enables an interface for quadruple precision numbers,
  ! if the compiler supports them. However, you should be aware, that this is
  ! merely a convenience interface, as the values provided by Lua are only
  ! double precision.
  use aot_quadruple_vector_module

  ! Support for extended double precision.
  use aot_extdouble_vector_module

  implicit none

  public :: aot_table_get_val, aot_get_val, aot_top_get_val

  !> Use these routines to obtain a vector whose length is unknown.
  !!
  !! Arrays will be allocated as needed to read the data from the
  !! Lua script with these routines. A maximal length has to be
  !! specified to limit the allocated memory by these routines (and make the
  !! interfaces distinguishable).
  interface aot_get_val
    module procedure get_table_real_vvect
    module procedure get_table_double_vvect
    module procedure get_table_integer_vvect
    module procedure get_table_long_vvect
    module procedure get_table_logical_vvect
  end interface

  interface aot_table_get_val
    module procedure get_table_real_vvect
    module procedure get_table_double_vvect
    module procedure get_table_integer_vvect
    module procedure get_table_long_vvect
    module procedure get_table_logical_vvect
  end interface

  interface aot_top_get_val
    module procedure get_top_real_vvect
    module procedure get_top_double_vvect
    module procedure get_top_integer_vvect
    module procedure get_top_long_vvect
    module procedure get_top_logical_vvect
  end interface


  !> Use these routines to obtain a vector of known length.
  !!
  !! The given vector has to exist already and will be filled by
  !! values from the Lua table, as far as they exist.
  !! If the Lua table is longer than the available elements in the array
  !! only the first elements from the table will be stored in the array.
  interface aot_get_val
    module procedure get_table_real_v
    module procedure get_table_double_v
    module procedure get_table_integer_v
    module procedure get_table_long_v
    module procedure get_table_logical_v
  end interface

  interface aot_table_get_val
    module procedure get_table_real_v
    module procedure get_table_double_v
    module procedure get_table_integer_v
    module procedure get_table_long_v
    module procedure get_table_logical_v
  end interface

  interface aot_top_get_val
    module procedure get_top_real_v
    module procedure get_top_double_v
    module procedure get_top_integer_v
    module procedure get_top_long_v
    module procedure get_top_logical_v
  end interface


  private


contains

  !> This routine obtains a vectorial quantity with variable length from a Lua
  !! table as a whole.
  !!
  !! It is intented to ease the reading of vectors on the Fortran side by
  !! capsulating the parsing of the Lua table internally.
  !! For the dynamically sized array, which will be allocated, a upper limit
  !! to allocate has to be specified.
  subroutine get_table_real_vvect(val, ErrCode, maxlength, L, thandle, &
    &                             key, pos, default)
    type(flu_State) :: L !< Handle to the lua script
    integer, intent(in), optional :: thandle !< Handle of the parent table

    !> Vector read from the Lua table, will have the same length as the table
    !! but not exceed maxlength, if provided.
    real(kind=single_k), intent(out), allocatable :: val(:)

    !> Error code describing problems encountered in each of the components.
    !! Will be allocated with the same length as the returned vector.
    !! If the complete vector is not given in the Lua script, and no default
    !! is provided, an zerosized array will be returned.
    integer, intent(out), allocatable :: ErrCode(:)

    !> Maximal length to allocate for the vector.
    integer, intent(in) :: maxlength

    !> Name of the variable (vector) to read.
    character(len=*), intent(in), optional :: key

    !> Position of the (vector) to read.
    integer, intent(in), optional :: pos

    !> A default vector to use, if no proper definition is found.
    !! Components will be filled with the help of this default definition.
    real(kind=single_k), intent(in), optional :: default(:)

    logical :: valid_args

    valid_args = .true.
    if (present(thandle)) then
      ! Get the requested value from the provided table
      call aot_table_push(L=L, thandle=thandle, &
        &                   key=key, pos=pos)
    else
      if (present(key)) then
        ! Get the requeseted global variable
        call flu_getglobal(L, key)
      else
        valid_args = .false.
      end if
    end if
    if (valid_args) then
      call aot_top_get_val(val, ErrCode, maxlength, L, default)
    else
      ! In case of invalid arguments return 0-sized arrays.
      ! (Equivalent of not found Lua tables.)
      allocate(Val(0))
      allocate(ErrCode(0))
    end if

  end subroutine get_table_real_vvect


  !> This routine obtains a vectorial quantity with variable length from a Lua
  !! table as a whole.
  !!
  !! It is intented to ease the reading of vectors on the Fortran side by
  !! capsulating the parsing of the Lua table internally.
  !! For the dynamically sized array, which will be allocated, a upper limit
  !! to allocate has to be specified.
  subroutine get_table_double_vvect(val, ErrCode, maxlength, L, thandle, &
    &                               key, pos, default)
    type(flu_State) :: L !< Handle to the lua script
    integer, intent(in), optional :: thandle !< Handle of the parent table

    !> Vector read from the Lua table, will have the same length as the table
    !! but not exceed maxlength, if provided.
    real(kind=double_k), intent(out), allocatable :: val(:)

    !> Error code describing problems encountered in each of the components.
    !! Will be allocated with the same length as the returned vector.
    !! If the complete vector is not given in the Lua script, and no default
    !! is provided, an zerosized array will be returned.
    integer, intent(out), allocatable :: ErrCode(:)

    !> Maximal length to allocate for the vector.
    integer, intent(in) :: maxlength

    !> Name of the variable (vector) to read.
    character(len=*), intent(in), optional :: key

    !> Position of the (vector) to read.
    integer, intent(in), optional :: pos

    !> A default vector to use, if no proper definition is found.
    !! Components will be filled with the help of this default definition.
    real(kind=double_k), intent(in), optional :: default(:)

    logical :: valid_args

    valid_args = .true.
    if (present(thandle)) then
      ! Get the requested value from the provided table
      call aot_table_push(L=L, thandle=thandle, &
        &                   key=key, pos=pos)
    else
      if (present(key)) then
        ! Get the requeseted global variable
        call flu_getglobal(L, key)
      else
        valid_args = .false.
      end if
    end if
    if (valid_args) then
      call aot_top_get_val(val, ErrCode, maxlength, L, default)
    else
      ! In case of invalid arguments return 0-sized arrays.
      ! (Equivalent of not found Lua tables.)
      allocate(Val(0))
      allocate(ErrCode(0))
    end if

  end subroutine get_table_double_vvect


  !> This routine obtains a vectorial quantity with variable length from a Lua
  !! table as a whole.
  !!
  !! It is intented to ease the reading of vectors on the Fortran side by
  !! capsulating the parsing of the Lua table internally.
  !! For the dynamically sized array, which will be allocated, a upper limit
  !! to allocate has to be specified.
  subroutine get_table_integer_vvect(val, ErrCode, maxlength, L, thandle, &
    &                                key, pos, default)
    type(flu_State) :: L !< Handle to the lua script
    integer, intent(in), optional :: thandle !< Handle of the parent table

    !> Vector read from the Lua table, will have the same length as the table
    !! but not exceed maxlength, if provided.
    integer, intent(out), allocatable :: val(:)

    !> Error code describing problems encountered in each of the components.
    !! Will be allocated with the same length as the returned vector.
    !! If the complete vector is not given in the Lua script, and no default
    !! is provided, an zerosized array will be returned.
    integer, intent(out), allocatable :: ErrCode(:)

    !> Maximal length to allocate for the vector.
    integer, intent(in) :: maxlength

    !> Name of the variable (vector) to read.
    character(len=*), intent(in), optional :: key

    !> Position of the (vector) to read.
    integer, intent(in), optional :: pos

    !> A default vector to use, if no proper definition is found.
    !! Components will be filled with the help of this default definition.
    integer, intent(in), optional :: default(:)

    logical :: valid_args

    valid_args = .true.
    if (present(thandle)) then
      ! Get the requested value from the provided table
      call aot_table_push(L=L, thandle=thandle, &
        &                   key=key, pos=pos)
    else
      if (present(key)) then
        ! Get the requeseted global variable
        call flu_getglobal(L, key)
      else
        valid_args = .false.
      end if
    end if
    if (valid_args) then
      call aot_top_get_val(val, ErrCode, maxlength, L, default)
    else
      ! In case of invalid arguments return 0-sized arrays.
      ! (Equivalent of not found Lua tables.)
      allocate(Val(0))
      allocate(ErrCode(0))
    end if

  end subroutine get_table_integer_vvect


  !> This routine obtains a vectorial quantity with variable length from a Lua
  !! table as a whole.
  !!
  !! It is intented to ease the reading of vectors on the Fortran side by
  !! capsulating the parsing of the Lua table internally.
  !! For the dynamically sized array, which will be allocated, a upper limit
  !! to allocate has to be specified.
  subroutine get_table_long_vvect(val, ErrCode, maxlength, L, thandle, &
    &                             key, pos, default)
    type(flu_State) :: L !< Handle to the lua script
    integer, intent(in), optional :: thandle !< Handle of the parent table

    !> Vector read from the Lua table, will have the same length as the table
    !! but not exceed maxlength, if provided.
    integer(kind=long_k), intent(out), allocatable :: val(:)

    !> Error code describing problems encountered in each of the components.
    !! Will be allocated with the same length as the returned vector.
    !! If the complete vector is not given in the Lua script, and no default
    !! is provided, an zerosized array will be returned.
    integer, intent(out), allocatable :: ErrCode(:)

    !> Maximal length to allocate for the vector.
    integer, intent(in) :: maxlength

    !> Name of the variable (vector) to read.
    character(len=*), intent(in), optional :: key

    !> Position of the (vector) to read.
    integer, intent(in), optional :: pos

    !> A default vector to use, if no proper definition is found.
    !! Components will be filled with the help of this default definition.
    integer(kind=long_k), intent(in), optional :: default(:)

    logical :: valid_args

    valid_args = .true.
    if (present(thandle)) then
      ! Get the requested value from the provided table
      call aot_table_push(L=L, thandle=thandle, &
        &                   key=key, pos=pos)
    else
      if (present(key)) then
        ! Get the requeseted global variable
        call flu_getglobal(L, key)
      else
        valid_args = .false.
      end if
    end if
    if (valid_args) then
      call aot_top_get_val(val, ErrCode, maxlength, L, default)
    else
      ! In case of invalid arguments return 0-sized arrays.
      ! (Equivalent of not found Lua tables.)
      allocate(Val(0))
      allocate(ErrCode(0))
    end if

  end subroutine get_table_long_vvect


  !> This routine obtains a vectorial quantity with variable length from a Lua
  !! table as a whole.
  !!
  !! It is intented to ease the reading of vectors on the Fortran side by
  !! capsulating the parsing of the Lua table internally.
  !! For the dynamically sized array, which will be allocated, a upper limit
  !! to allocate has to be specified.
  subroutine get_table_logical_vvect(val, ErrCode, maxlength, L, thandle, &
    &                                key, pos, default)
    type(flu_State) :: L !< Handle to the lua script
    integer, intent(in), optional :: thandle !< Handle of the parent table

    !> Vector read from the Lua table, will have the same length as the table
    !! but not exceed maxlength, if provided.
    logical, intent(out), allocatable :: val(:)

    !> Error code describing problems encountered in each of the components.
    !! Will be allocated with the same length as the returned vector.
    !! If the complete vector is not given in the Lua script, and no default
    !! is provided, an zerosized array will be returned.
    integer, intent(out), allocatable :: ErrCode(:)

    !> Maximal length to allocate for the vector.
    integer, intent(in) :: maxlength

    !> Name of the variable (vector) to read.
    character(len=*), intent(in), optional :: key

    !> Position of the (vector) to read.
    integer, intent(in), optional :: pos

    !> A default vector to use, if no proper definition is found.
    !! Components will be filled with the help of this default definition.
    logical, intent(in), optional :: default(:)

    logical :: valid_args

    valid_args = .true.
    if (present(thandle)) then
      ! Get the requested value from the provided table
      call aot_table_push(L=L, thandle=thandle, &
        &                   key=key, pos=pos)
    else
      if (present(key)) then
        ! Get the requeseted global variable
        call flu_getglobal(L, key)
      else
        valid_args = .false.
      end if
    end if
    if (valid_args) then
      call aot_top_get_val(val, ErrCode, maxlength, L, default)
    else
      ! In case of invalid arguments return 0-sized arrays.
      ! (Equivalent of not found Lua tables.)
      allocate(Val(0))
      allocate(ErrCode(0))
    end if

  end subroutine get_table_logical_vvect



  !> This routine obtains a vectorial quantity with fixed length from a Lua
  !! table as a whole.
  !!
  !! It is intented to ease the reading of vectors on the Fortran side by
  !! capsulating the parsing of the Lua table internally.
  !! Components which are not found are filled with the data given in
  !! the default vector. For each component an error code will be returned
  !! to indicate the success when reading it.
  !! If the vector is not defined at all, all components will be indicated
  !! as non-existent.
  !! Components, which are neither defined in the Lua script, nor in the
  !! default will be marked with the aoterr_Fatal flag.
  subroutine get_table_real_v(val, ErrCode, L, thandle, key, &
    &                         pos, default)
    type(flu_State) :: L !< Handle to the lua script
    integer, intent(in), optional :: thandle !< Handle of the parent table

    !> Vector read from the Lua table.
    real(kind=single_k), intent(out) :: val(:)

    !> Error code describing problems encountered in each of the components.
    !! This array has to have the same length as val.
    integer, intent(out) :: ErrCode(:)

    !> Name of the variable (vector) to read.
    character(len=*), intent(in), optional :: key

    !> Position of the (vector) to read.
    integer, intent(in), optional :: pos

    !> A default vector to use, if no proper definition is found.
    !! Components will be filled with the help of this default definition.
    real(kind=single_k), intent(in), optional :: default(:)

    logical :: valid_args

    valid_args = .true.
    if (present(thandle)) then
      ! Get the requested value from the provided table
      call aot_table_push(L=L, thandle=thandle, &
        &                 key=key, pos=pos)
    else
      if (present(key)) then
        ! Get the requeseted global variable
        call flu_getglobal(L, key)
      else
        valid_args = .false.
      end if
    end if
    if (valid_args) then
      call aot_top_get_val(val, ErrCode, L, default)
    else
      ErrCode = ibSet(0, aoterr_NonExistent)
      ErrCode = ibSet(ErrCode, aoterr_Fatal)
    end if

  end subroutine get_table_real_v


  !> This routine obtains a vectorial quantity with fixed length from a Lua
  !! table as a whole.
  !!
  !! It is intented to ease the reading of vectors on the Fortran side by
  !! capsulating the parsing of the Lua table internally.
  !! Components which are not found are filled with the data given in
  !! the default vector. For each component an error code will be returned
  !! to indicate the success when reading it.
  !! If the vector is not defined at all, all components will be indicated
  !! as non-existent.
  !! Components, which are neither defined in the Lua script, nor in the
  !! default will be marked with the aoterr_Fatal flag.
  subroutine get_table_double_v(val, ErrCode, L, thandle, key, &
    &                         pos, default)
    type(flu_State) :: L !< Handle to the lua script
    integer, intent(in), optional :: thandle !< Handle of the parent table

    !> Vector read from the Lua table.
    real(kind=double_k), intent(out) :: val(:)

    !> Error code describing problems encountered in each of the components.
    !! This array has to have the same length as val.
    integer, intent(out) :: ErrCode(:)

    !> Name of the variable (vector) to read.
    character(len=*), intent(in), optional :: key

    !> Position of the (vector) to read.
    integer, intent(in), optional :: pos

    !> A default vector to use, if no proper definition is found.
    !! Components will be filled with the help of this default definition.
    real(kind=double_k), intent(in), optional :: default(:)

    logical :: valid_args

    valid_args = .true.
    if (present(thandle)) then
      ! Get the requested value from the provided table
      call aot_table_push(L=L, thandle=thandle, &
        &                 key=key, pos=pos)
    else
      if (present(key)) then
        ! Get the requeseted global variable
        call flu_getglobal(L, key)
      else
        valid_args = .false.
      end if
    end if
    if (valid_args) then
      call aot_top_get_val(val, ErrCode, L, default)
    else
      ErrCode = ibSet(0, aoterr_NonExistent)
      ErrCode = ibSet(ErrCode, aoterr_Fatal)
    end if

  end subroutine get_table_double_v


  !> This routine obtains a vectorial quantity with fixed length from a Lua
  !! table as a whole.
  !!
  !! It is intented to ease the reading of vectors on the Fortran side by
  !! capsulating the parsing of the Lua table internally.
  !! Components which are not found are filled with the data given in
  !! the default vector. For each component an error code will be returned
  !! to indicate the success when reading it.
  !! If the vector is not defined at all, all components will be indicated
  !! as non-existent.
  !! Components, which are neither defined in the Lua script, nor in the
  !! default will be marked with the aoterr_Fatal flag.
  subroutine get_table_integer_v(val, ErrCode, L, thandle, key, &
    &                         pos, default)
    type(flu_State) :: L !< Handle to the lua script
    integer, intent(in), optional :: thandle !< Handle of the parent table

    !> Vector read from the Lua table.
    integer, intent(out) :: val(:)

    !> Error code describing problems encountered in each of the components.
    !! This array has to have the same length as val.
    integer, intent(out) :: ErrCode(:)

    !> Name of the variable (vector) to read.
    character(len=*), intent(in), optional :: key

    !> Position of the (vector) to read.
    integer, intent(in), optional :: pos

    !> A default vector to use, if no proper definition is found.
    !! Components will be filled with the help of this default definition.
    integer, intent(in), optional :: default(:)

    logical :: valid_args

    valid_args = .true.
    if (present(thandle)) then
      ! Get the requested value from the provided table
      call aot_table_push(L=L, thandle=thandle, &
        &                 key=key, pos=pos)
    else
      if (present(key)) then
        ! Get the requeseted global variable
        call flu_getglobal(L, key)
      else
        valid_args = .false.
      end if
    end if
    if (valid_args) then
      call aot_top_get_val(val, ErrCode, L, default)
    else
      ErrCode = ibSet(0, aoterr_NonExistent)
      ErrCode = ibSet(ErrCode, aoterr_Fatal)
    end if

  end subroutine get_table_integer_v


  !> This routine obtains a vectorial quantity with fixed length from a Lua
  !! table as a whole.
  !!
  !! It is intented to ease the reading of vectors on the Fortran side by
  !! capsulating the parsing of the Lua table internally.
  !! Components which are not found are filled with the data given in
  !! the default vector. For each component an error code will be returned
  !! to indicate the success when reading it.
  !! If the vector is not defined at all, all components will be indicated
  !! as non-existent.
  !! Components, which are neither defined in the Lua script, nor in the
  !! default will be marked with the aoterr_Fatal flag.
  subroutine get_table_long_v(val, ErrCode, L, thandle, key, &
    &                         pos, default)
    type(flu_State) :: L !< Handle to the lua script
    integer, intent(in), optional :: thandle !< Handle of the parent table

    !> Vector read from the Lua table.
    integer(kind=long_k), intent(out) :: val(:)

    !> Error code describing problems encountered in each of the components.
    !! This array has to have the same length as val.
    integer, intent(out) :: ErrCode(:)

    !> Name of the variable (vector) to read.
    character(len=*), intent(in), optional :: key

    !> Position of the (vector) to read.
    integer, intent(in), optional :: pos

    !> A default vector to use, if no proper definition is found.
    !! Components will be filled with the help of this default definition.
    integer(kind=long_k), intent(in), optional :: default(:)

    logical :: valid_args

    valid_args = .true.
    if (present(thandle)) then
      ! Get the requested value from the provided table
      call aot_table_push(L=L, thandle=thandle, &
        &                 key=key, pos=pos)
    else
      if (present(key)) then
        ! Get the requeseted global variable
        call flu_getglobal(L, key)
      else
        valid_args = .false.
      end if
    end if
    if (valid_args) then
      call aot_top_get_val(val, ErrCode, L, default)
    else
      ErrCode = ibSet(0, aoterr_NonExistent)
      ErrCode = ibSet(ErrCode, aoterr_Fatal)
    end if

  end subroutine get_table_long_v


  !> This routine obtains a vectorial quantity with fixed length from a Lua
  !! table as a whole.
  !!
  !! It is intented to ease the reading of vectors on the Fortran side by
  !! capsulating the parsing of the Lua table internally.
  !! Components which are not found are filled with the data given in
  !! the default vector. For each component an error code will be returned
  !! to indicate the success when reading it.
  !! If the vector is not defined at all, all components will be indicated
  !! as non-existent.
  !! Components, which are neither defined in the Lua script, nor in the
  !! default will be marked with the aoterr_Fatal flag.
  subroutine get_table_logical_v(val, ErrCode, L, thandle, key, &
    &                         pos, default)
    type(flu_State) :: L !< Handle to the lua script
    integer, intent(in), optional :: thandle !< Handle of the parent table

    !> Vector read from the Lua table.
    logical, intent(out) :: val(:)

    !> Error code describing problems encountered in each of the components.
    !! This array has to have the same length as val.
    integer, intent(out) :: ErrCode(:)

    !> Name of the variable (vector) to read.
    character(len=*), intent(in), optional :: key

    !> Position of the (vector) to read.
    integer, intent(in), optional :: pos

    !> A default vector to use, if no proper definition is found.
    !! Components will be filled with the help of this default definition.
    logical, intent(in), optional :: default(:)

    logical :: valid_args

    valid_args = .true.
    if (present(thandle)) then
      ! Get the requested value from the provided table
      call aot_table_push(L=L, thandle=thandle, &
        &                 key=key, pos=pos)
    else
      if (present(key)) then
        ! Get the requeseted global variable
        call flu_getglobal(L, key)
      else
        valid_args = .false.
      end if
    end if
    if (valid_args) then
      call aot_top_get_val(val, ErrCode, L, default)
    else
      ErrCode = ibSet(0, aoterr_NonExistent)
      ErrCode = ibSet(ErrCode, aoterr_Fatal)
    end if

  end subroutine get_table_logical_v




  subroutine get_top_real_vvect(val, ErrCode, maxlength, L, default)
    type(flu_State) :: L !< Handle to the lua script

    !> Vector read from the Lua table, will have the same length as the table
    !! but not exceed maxlength, if provided.
    real(kind=single_k), intent(out), allocatable :: val(:)

    !> Error code describing problems encountered in each of the components.
    !! Will be allocated with the same length as the returned vector.
    !! If the complete vector is not given in the Lua script, and no default
    !! is provided, an zerosized array will be returned.
    integer, intent(out), allocatable :: ErrCode(:)

    !> Maximal length to allocate for the vector.
    integer, intent(in) :: maxlength

    !> A default vector to use, if no proper definition is found.
    !! Components will be filled with the help of this default definition.
    real(kind=single_k), intent(in), optional :: default(:)

    integer :: vect_handle
    integer :: table_len, vect_len, def_len
    integer :: vect_lb
    integer :: iComp

    ! Try to interpret the top entry on the stack as a table
    vect_handle = aot_table_top(L=L)
    table_len = aot_table_length(L=L, thandle=vect_handle)

    ! The size of the vector is limited by maxlength.
    vect_len = min(maxlength, table_len)

    ! Find the length of the default value, if it is not provided, its 0.
    def_len = 0
    if (present(default)) def_len = size(default)

    ! Now parse the table with the vector entries.
    if (aot_table_first(L, vect_handle)) then
      allocate(val(vect_len))
      allocate(errCode(vect_len))

      ErrCode = 0

      ! Only if the vector table actually exists, and has at least one entry,
      ! this parsing has to be done.
      if (present(default).and.(def_len > 0)) then
        call aot_top_get_val(val(1), ErrCode(1), L, default(1))
      else
        call aot_top_get_val(val(1), ErrCode(1), L)
      end if

      ! Up to the length of the default value, provide the default settings.
      do iComp=2,def_len
        if (.not. flu_next(L, vect_handle)) exit
        call aot_top_get_val(val(iComp), ErrCode(iComp), L, &
          &                  default(iComp))
      end do

      vect_lb = max(2, def_len+1)
      ! After def_len entries no default values for the components are
      ! available anymore, proceed without a default setting for the rest.
      do iComp=vect_lb,vect_len
        if (.not. flu_next(L, vect_handle)) exit
        call aot_top_get_val(val(iComp), ErrCode(iComp), L)
      end do
    else
      ! No vector definition found in the Lua script, use the default.
      if (present(default)) then
        allocate(val(def_len))
        allocate(errCode(vect_len))
        val = default
        ErrCode = ibSet(0, aoterr_NonExistent)
      else
        ! No vector definition in the Lua script and no default provided,
        ! return an empty array.
        allocate(val(0))
        allocate(errCode(0))
      end if
    end if
    call aot_table_close(L, vect_handle)

  end subroutine get_top_real_vvect



  subroutine get_top_double_vvect(val, ErrCode, maxlength, L, default)
    type(flu_State) :: L !< Handle to the lua script

    !> Vector read from the Lua table, will have the same length as the table
    !! but not exceed maxlength, if provided.
    real(kind=double_k), intent(out), allocatable :: val(:)

    !> Error code describing problems encountered in each of the components.
    !! Will be allocated with the same length as the returned vector.
    !! If the complete vector is not given in the Lua script, and no default
    !! is provided, an zerosized array will be returned.
    integer, intent(out), allocatable :: ErrCode(:)

    !> Maximal length to allocate for the vector.
    integer, intent(in) :: maxlength

    !> A default vector to use, if no proper definition is found.
    !! Components will be filled with the help of this default definition.
    real(kind=double_k), intent(in), optional :: default(:)

    integer :: vect_handle
    integer :: table_len, vect_len, def_len
    integer :: vect_lb
    integer :: iComp

    ! Try to interpret the top entry on the stack as a table
    vect_handle = aot_table_top(L=L)
    table_len = aot_table_length(L=L, thandle=vect_handle)

    ! The size of the vector is limited by maxlength.
    vect_len = min(maxlength, table_len)

    ! Find the length of the default value, if it is not provided, its 0.
    def_len = 0
    if (present(default)) def_len = size(default)

    ! Now parse the table with the vector entries.
    if (aot_table_first(L, vect_handle)) then
      allocate(val(vect_len))
      allocate(errCode(vect_len))

      ErrCode = 0

      ! Only if the vector table actually exists, and has at least one entry,
      ! this parsing has to be done.
      if (present(default).and.(def_len > 0)) then
        call aot_top_get_val(val(1), ErrCode(1), L, default(1))
      else
        call aot_top_get_val(val(1), ErrCode(1), L)
      end if

      ! Up to the length of the default value, provide the default settings.
      do iComp=2,def_len
        if (.not. flu_next(L, vect_handle)) exit
        call aot_top_get_val(val(iComp), ErrCode(iComp), L, &
          &                  default(iComp))
      end do

      vect_lb = max(2, def_len+1)
      ! After def_len entries no default values for the components are
      ! available anymore, proceed without a default setting for the rest.
      do iComp=vect_lb,vect_len
        if (.not. flu_next(L, vect_handle)) exit
        call aot_top_get_val(val(iComp), ErrCode(iComp), L)
      end do
    else
      ! No vector definition found in the Lua script, use the default.
      if (present(default)) then
        allocate(val(def_len))
        allocate(errCode(vect_len))
        val = default
        ErrCode = ibSet(0, aoterr_NonExistent)
      else
        ! No vector definition in the Lua script and no default provided,
        ! return an empty array.
        allocate(val(0))
        allocate(errCode(0))
      end if
    end if
    call aot_table_close(L, vect_handle)

  end subroutine get_top_double_vvect


  subroutine get_top_integer_vvect(val, ErrCode, maxlength, L, default)
    type(flu_State) :: L !< Handle to the lua script

    !> Vector read from the Lua table, will have the same length as the table
    !! but not exceed maxlength, if provided.
    integer, intent(out), allocatable :: val(:)

    !> Error code describing problems encountered in each of the components.
    !! Will be allocated with the same length as the returned vector.
    !! If the complete vector is not given in the Lua script, and no default
    !! is provided, an zerosized array will be returned.
    integer, intent(out), allocatable :: ErrCode(:)

    !> Maximal length to allocate for the vector.
    integer, intent(in) :: maxlength

    !> A default vector to use, if no proper definition is found.
    !! Components will be filled with the help of this default definition.
    integer, intent(in), optional :: default(:)

    integer :: vect_handle
    integer :: table_len, vect_len, def_len
    integer :: vect_lb
    integer :: iComp

    ! Try to interpret the top entry on the stack as a table
    vect_handle = aot_table_top(L=L)
    table_len = aot_table_length(L=L, thandle=vect_handle)

    ! The size of the vector is limited by maxlength.
    vect_len = min(maxlength, table_len)

    ! Find the length of the default value, if it is not provided, its 0.
    def_len = 0
    if (present(default)) def_len = size(default)

    ! Now parse the table with the vector entries.
    if (aot_table_first(L, vect_handle)) then
      allocate(val(vect_len))
      allocate(errCode(vect_len))

      ErrCode = 0

      ! Only if the vector table actually exists, and has at least one entry,
      ! this parsing has to be done.
      if (present(default).and.(def_len > 0)) then
        call aot_top_get_val(val(1), ErrCode(1), L, default(1))
      else
        call aot_top_get_val(val(1), ErrCode(1), L)
      end if

      ! Up to the length of the default value, provide the default settings.
      do iComp=2,def_len
        if (.not. flu_next(L, vect_handle)) exit
        call aot_top_get_val(val(iComp), ErrCode(iComp), L, &
          &                  default(iComp))
      end do

      vect_lb = max(2, def_len+1)
      ! After def_len entries no default values for the components are
      ! available anymore, proceed without a default setting for the rest.
      do iComp=vect_lb,vect_len
        if (.not. flu_next(L, vect_handle)) exit
        call aot_top_get_val(val(iComp), ErrCode(iComp), L)
      end do
    else
      ! No vector definition found in the Lua script, use the default.
      if (present(default)) then
        allocate(val(def_len))
        allocate(errCode(vect_len))
        val = default
        ErrCode = ibSet(0, aoterr_NonExistent)
      else
        ! No vector definition in the Lua script and no default provided,
        ! return an empty array.
        allocate(val(0))
        allocate(errCode(0))
      end if
    end if
    call aot_table_close(L, vect_handle)

  end subroutine get_top_integer_vvect


  subroutine get_top_long_vvect(val, ErrCode, maxlength, L, default)
    type(flu_State) :: L !< Handle to the lua script

    !> Vector read from the Lua table, will have the same length as the table
    !! but not exceed maxlength, if provided.
    integer(kind=long_k), intent(out), allocatable :: val(:)

    !> Error code describing problems encountered in each of the components.
    !! Will be allocated with the same length as the returned vector.
    !! If the complete vector is not given in the Lua script, and no default
    !! is provided, an zerosized array will be returned.
    integer, intent(out), allocatable :: ErrCode(:)

    !> Maximal length to allocate for the vector.
    integer, intent(in) :: maxlength

    !> A default vector to use, if no proper definition is found.
    !! Components will be filled with the help of this default definition.
    integer(kind=long_k), intent(in), optional :: default(:)

    integer :: vect_handle
    integer :: table_len, vect_len, def_len
    integer :: vect_lb
    integer :: iComp

    ! Try to interpret the top entry on the stack as a table
    vect_handle = aot_table_top(L=L)
    table_len = aot_table_length(L=L, thandle=vect_handle)

    ! The size of the vector is limited by maxlength.
    vect_len = min(maxlength, table_len)

    ! Find the length of the default value, if it is not provided, its 0.
    def_len = 0
    if (present(default)) def_len = size(default)

    ! Now parse the table with the vector entries.
    if (aot_table_first(L, vect_handle)) then
      allocate(val(vect_len))
      allocate(errCode(vect_len))

      ErrCode = 0

      ! Only if the vector table actually exists, and has at least one entry,
      ! this parsing has to be done.
      if (present(default).and.(def_len > 0)) then
        call aot_top_get_val(val(1), ErrCode(1), L, default(1))
      else
        call aot_top_get_val(val(1), ErrCode(1), L)
      end if

      ! Up to the length of the default value, provide the default settings.
      do iComp=2,def_len
        if (.not. flu_next(L, vect_handle)) exit
        call aot_top_get_val(val(iComp), ErrCode(iComp), L, &
          &                  default(iComp))
      end do

      vect_lb = max(2, def_len+1)
      ! After def_len entries no default values for the components are
      ! available anymore, proceed without a default setting for the rest.
      do iComp=vect_lb,vect_len
        if (.not. flu_next(L, vect_handle)) exit
        call aot_top_get_val(val(iComp), ErrCode(iComp), L)
      end do
    else
      ! No vector definition found in the Lua script, use the default.
      if (present(default)) then
        allocate(val(def_len))
        allocate(errCode(vect_len))
        val = default
        ErrCode = ibSet(0, aoterr_NonExistent)
      else
        ! No vector definition in the Lua script and no default provided,
        ! return an empty array.
        allocate(val(0))
        allocate(errCode(0))
      end if
    end if
    call aot_table_close(L, vect_handle)

  end subroutine get_top_long_vvect


  subroutine get_top_logical_vvect(val, ErrCode, maxlength, L, default)
    type(flu_State) :: L !< Handle to the lua script

    !> Vector read from the Lua table, will have the same length as the table
    !! but not exceed maxlength, if provided.
    logical, intent(out), allocatable :: val(:)

    !> Error code describing problems encountered in each of the components.
    !! Will be allocated with the same length as the returned vector.
    !! If the complete vector is not given in the Lua script, and no default
    !! is provided, an zerosized array will be returned.
    integer, intent(out), allocatable :: ErrCode(:)

    !> Maximal length to allocate for the vector.
    integer, intent(in) :: maxlength

    !> A default vector to use, if no proper definition is found.
    !! Components will be filled with the help of this default definition.
    logical, intent(in), optional :: default(:)

    integer :: vect_handle
    integer :: table_len, vect_len, def_len
    integer :: vect_lb
    integer :: iComp

    ! Try to interpret the top entry on the stack as a table
    vect_handle = aot_table_top(L=L)
    table_len = aot_table_length(L=L, thandle=vect_handle)

    ! The size of the vector is limited by maxlength.
    vect_len = min(maxlength, table_len)

    ! Find the length of the default value, if it is not provided, its 0.
    def_len = 0
    if (present(default)) def_len = size(default)

    ! Now parse the table with the vector entries.
    if (aot_table_first(L, vect_handle)) then
      allocate(val(vect_len))
      allocate(errCode(vect_len))

      ErrCode = 0

      ! Only if the vector table actually exists, and has at least one entry,
      ! this parsing has to be done.
      if (present(default).and.(def_len > 0)) then
        call aot_top_get_val(val(1), ErrCode(1), L, default(1))
      else
        call aot_top_get_val(val(1), ErrCode(1), L)
      end if

      ! Up to the length of the default value, provide the default settings.
      do iComp=2,def_len
        if (.not. flu_next(L, vect_handle)) exit
        call aot_top_get_val(val(iComp), ErrCode(iComp), L, &
          &                  default(iComp))
      end do

      vect_lb = max(2, def_len+1)
      ! After def_len entries no default values for the components are
      ! available anymore, proceed without a default setting for the rest.
      do iComp=vect_lb,vect_len
        if (.not. flu_next(L, vect_handle)) exit
        call aot_top_get_val(val(iComp), ErrCode(iComp), L)
      end do
    else
      ! No vector definition found in the Lua script, use the default.
      if (present(default)) then
        allocate(val(def_len))
        allocate(errCode(vect_len))
        val = default
        ErrCode = ibSet(0, aoterr_NonExistent)
      else
        ! No vector definition in the Lua script and no default provided,
        ! return an empty array.
        allocate(val(0))
        allocate(errCode(0))
      end if
    end if
    call aot_table_close(L, vect_handle)

  end subroutine get_top_logical_vvect



  subroutine get_top_real_v(val, ErrCode, L,  default)
    type(flu_State) :: L !< Handle to the lua script

    !> Vector read from the Lua table.
    real(kind=single_k), intent(out) :: val(:)

    !> Error code describing problems encountered in each of the components.
    !! This array has to have the same length as val.
    integer, intent(out) :: ErrCode(:)

    !> A default vector to use, if no proper definition is found.
    !! Components will be filled with the help of this default definition.
    real(kind=single_k), intent(in), optional :: default(:)

    integer :: vect_handle
    integer :: table_len, vect_len, def_len
    integer :: vect_lb
    integer :: iComp

    ErrCode = 0

    ! Try to interpret it as table.
    vect_handle = aot_table_top(L=L)
    table_len = aot_table_length(L=L, thandle=vect_handle)

    vect_len = min(table_len, size(val))

    ! Find the length of the default value, if it is not provided, its 0.
    def_len = 0
    if (present(default)) def_len = size(default)

    ! Now parse the table with the vector entries.
    if (aot_table_first(L, vect_handle).and.(vect_len > 0)) then

      ! Only if the vector table actually exists, and has at least one entry,
      ! this parsing has to be done.
      if (present(default).and.(def_len > 0)) then
        call aot_top_get_val(val(1), ErrCode(1), L, default(1))
      else
        call aot_top_get_val(val(1), ErrCode(1), L)
      end if

      ! Up to the length of the default value, provide the default settings.
      do iComp=2,def_len
        if (.not. flu_next(L, vect_handle)) exit
        call aot_top_get_val(val(iComp), ErrCode(iComp), L, &
          &                  default(iComp))
      end do

      vect_lb = max(2, def_len+1)
      ! After def_len entries no default values for the components are
      ! available anymore, proceed without a default setting for the rest.
      do iComp=vect_lb,vect_len
        if (.not. flu_next(L, vect_handle)) exit
        call aot_top_get_val(val(iComp), ErrCode(iComp), L)
      end do

      ! If the table in the Lua script is not long enough, fill the remaining
      ! components with the default components, as far as they are defined.
      do iComp=vect_len+1,def_len
        ErrCode(iComp) = ibSet(ErrCode(iComp), aoterr_NonExistent)
        val(iComp) = default(iComp)
      end do
      vect_lb = max(vect_len+1, def_len)
      do iComp=vect_lb,vect_len
        ErrCode(iComp) = ibSet(ErrCode(iComp), aoterr_Fatal)
      end do
    else
      ! No vector definition found in the Lua script, use the default.
      ErrCode = ibSet(ErrCode, aoterr_NonExistent)
      if (present(default)) then
        val(:def_len) = default(:def_len)
        if (def_len < vect_len) then
          ErrCode(def_len+1:) = ibSet(ErrCode(def_len+1:), aoterr_Fatal)
        end if
      else
        ! No vector definition in the Lua script and no default provided.
        ErrCode = ibSet(ErrCode, aoterr_Fatal)
      end if
    end if
    call aot_table_close(L, vect_handle)

  end subroutine get_top_real_v


  subroutine get_top_double_v(val, ErrCode, L,  default)
    type(flu_State) :: L !< Handle to the lua script

    !> Vector read from the Lua table.
    real(kind=double_k), intent(out) :: val(:)

    !> Error code describing problems encountered in each of the components.
    !! This array has to have the same length as val.
    integer, intent(out) :: ErrCode(:)

    !> A default vector to use, if no proper definition is found.
    !! Components will be filled with the help of this default definition.
    real(kind=double_k), intent(in), optional :: default(:)

    integer :: vect_handle
    integer :: table_len, vect_len, def_len
    integer :: vect_lb
    integer :: iComp

    ! Try to interpret it as table.
    vect_handle = aot_table_top(L=L)
    table_len = aot_table_length(L=L, thandle=vect_handle)

    ErrCode = 0

    vect_len = min(table_len, size(val))

    ! Find the length of the default value, if it is not provided, its 0.
    def_len = 0
    if (present(default)) def_len = size(default)

    ! Now parse the table with the vector entries.
    if (aot_table_first(L, vect_handle).and.(vect_len > 0)) then

      ! Only if the vector table actually exists, and has at least one entry,
      ! this parsing has to be done.
      if (present(default).and.(def_len > 0)) then
        call aot_top_get_val(val(1), ErrCode(1), L, default(1))
      else
        call aot_top_get_val(val(1), ErrCode(1), L)
      end if

      ! Up to the length of the default value, provide the default settings.
      do iComp=2,def_len
        if (.not. flu_next(L, vect_handle)) exit
        call aot_top_get_val(val(iComp), ErrCode(iComp), L, &
          &                  default(iComp))
      end do

      vect_lb = max(2, def_len+1)
      ! After def_len entries no default values for the components are
      ! available anymore, proceed without a default setting for the rest.
      do iComp=vect_lb,vect_len
        if (.not. flu_next(L, vect_handle)) exit
        call aot_top_get_val(val(iComp), ErrCode(iComp), L)
      end do

      ! If the table in the Lua script is not long enough, fill the remaining
      ! components with the default components, as far as they are defined.
      do iComp=vect_len+1,def_len
        ErrCode(iComp) = ibSet(ErrCode(iComp), aoterr_NonExistent)
        val(iComp) = default(iComp)
      end do
      vect_lb = max(vect_len+1, def_len)
      do iComp=vect_lb,vect_len
        ErrCode(iComp) = ibSet(ErrCode(iComp), aoterr_Fatal)
      end do
    else
      ! No vector definition found in the Lua script, use the default.
      ErrCode = ibSet(ErrCode, aoterr_NonExistent)
      if (present(default)) then
        val(:def_len) = default(:def_len)
        if (def_len < vect_len) then
          ErrCode(def_len+1:) = ibSet(ErrCode(def_len+1:), aoterr_Fatal)
        end if
      else
        ! No vector definition in the Lua script and no default provided.
        ErrCode = ibSet(ErrCode, aoterr_Fatal)
      end if
    end if
    call aot_table_close(L, vect_handle)

  end subroutine get_top_double_v


  subroutine get_top_integer_v(val, ErrCode, L,  default)
    type(flu_State) :: L !< Handle to the lua script

    !> Vector read from the Lua table.
    integer, intent(out) :: val(:)

    !> Error code describing problems encountered in each of the components.
    !! This array has to have the same length as val.
    integer, intent(out) :: ErrCode(:)

    !> A default vector to use, if no proper definition is found.
    !! Components will be filled with the help of this default definition.
    integer, intent(in), optional :: default(:)

    integer :: vect_handle
    integer :: table_len, vect_len, def_len
    integer :: vect_lb
    integer :: iComp

    ErrCode = 0

    ! Try to interpret it as table.
    vect_handle = aot_table_top(L=L)
    table_len = aot_table_length(L=L, thandle=vect_handle)

    vect_len = min(table_len, size(val))

    ! Find the length of the default value, if it is not provided, its 0.
    def_len = 0
    if (present(default)) def_len = size(default)

    ! Now parse the table with the vector entries.
    if (aot_table_first(L, vect_handle).and.(vect_len > 0)) then

      ! Only if the vector table actually exists, and has at least one entry,
      ! this parsing has to be done.
      if (present(default).and.(def_len > 0)) then
        call aot_top_get_val(val(1), ErrCode(1), L, default(1))
      else
        call aot_top_get_val(val(1), ErrCode(1), L)
      end if

      ! Up to the length of the default value, provide the default settings.
      do iComp=2,def_len
        if (.not. flu_next(L, vect_handle)) exit
        call aot_top_get_val(val(iComp), ErrCode(iComp), L, &
          &                  default(iComp))
      end do

      vect_lb = max(2, def_len+1)
      ! After def_len entries no default values for the components are
      ! available anymore, proceed without a default setting for the rest.
      do iComp=vect_lb,vect_len
        if (.not. flu_next(L, vect_handle)) exit
        call aot_top_get_val(val(iComp), ErrCode(iComp), L)
      end do

      ! If the table in the Lua script is not long enough, fill the remaining
      ! components with the default components, as far as they are defined.
      do iComp=vect_len+1,def_len
        ErrCode(iComp) = ibSet(ErrCode(iComp), aoterr_NonExistent)
        val(iComp) = default(iComp)
      end do
      vect_lb = max(vect_len+1, def_len)
      do iComp=vect_lb,vect_len
        ErrCode(iComp) = ibSet(ErrCode(iComp), aoterr_Fatal)
      end do
    else
      ! No vector definition found in the Lua script, use the default.
      ErrCode = ibSet(ErrCode, aoterr_NonExistent)
      if (present(default)) then
        def_len = def_len
        val(:def_len) = default(:def_len)
        if (def_len < vect_len) then
          ErrCode(def_len+1:) = ibSet(ErrCode(def_len+1:), aoterr_Fatal)
        end if
      else
        ! No vector definition in the Lua script and no default provided.
        ErrCode = ibSet(ErrCode, aoterr_Fatal)
      end if
    end if
    call aot_table_close(L, vect_handle)

  end subroutine get_top_integer_v


  subroutine get_top_long_v(val, ErrCode, L,  default)
    type(flu_State) :: L !< Handle to the lua script

    !> Vector read from the Lua table.
    integer(kind=long_k), intent(out) :: val(:)

    !> Error code describing problems encountered in each of the components.
    !! This array has to have the same length as val.
    integer, intent(out) :: ErrCode(:)

    !> A default vector to use, if no proper definition is found.
    !! Components will be filled with the help of this default definition.
    integer(kind=long_k), intent(in), optional :: default(:)

    integer :: vect_handle
    integer :: table_len, vect_len, def_len
    integer :: vect_lb
    integer :: iComp

    ErrCode = 0

    ! Try to interpret it as table.
    vect_handle = aot_table_top(L=L)
    table_len = aot_table_length(L=L, thandle=vect_handle)

    vect_len = min(table_len, size(val))

    ! Find the length of the default value, if it is not provided, its 0.
    def_len = 0
    if (present(default)) def_len = size(default)

    ! Now parse the table with the vector entries.
    if (aot_table_first(L, vect_handle).and.(vect_len > 0)) then

      ! Only if the vector table actually exists, and has at least one entry,
      ! this parsing has to be done.
      if (present(default).and.(def_len > 0)) then
        call aot_top_get_val(val(1), ErrCode(1), L, default(1))
      else
        call aot_top_get_val(val(1), ErrCode(1), L)
      end if

      ! Up to the length of the default value, provide the default settings.
      do iComp=2,def_len
        if (.not. flu_next(L, vect_handle)) exit
        call aot_top_get_val(val(iComp), ErrCode(iComp), L, &
          &                  default(iComp))
      end do

      vect_lb = max(2, def_len+1)
      ! After def_len entries no default values for the components are
      ! available anymore, proceed without a default setting for the rest.
      do iComp=vect_lb,vect_len
        if (.not. flu_next(L, vect_handle)) exit
        call aot_top_get_val(val(iComp), ErrCode(iComp), L)
      end do

      ! If the table in the Lua script is not long enough, fill the remaining
      ! components with the default components, as far as they are defined.
      do iComp=vect_len+1,def_len
        ErrCode(iComp) = ibSet(ErrCode(iComp), aoterr_NonExistent)
        val(iComp) = default(iComp)
      end do
      vect_lb = max(vect_len+1, def_len)
      do iComp=vect_lb,vect_len
        ErrCode(iComp) = ibSet(ErrCode(iComp), aoterr_Fatal)
      end do
    else
      ! No vector definition found in the Lua script, use the default.
      ErrCode = ibSet(ErrCode, aoterr_NonExistent)
      if (present(default)) then
        val(:def_len) = default(:def_len)
        if (def_len < vect_len) then
          ErrCode(def_len+1:) = ibSet(ErrCode(def_len+1:), aoterr_Fatal)
        end if
      else
        ! No vector definition in the Lua script and no default provided.
        ErrCode = ibSet(ErrCode, aoterr_Fatal)
      end if
    end if
    call aot_table_close(L, vect_handle)

  end subroutine get_top_long_v


  subroutine get_top_logical_v(val, ErrCode, L,  default)
    type(flu_State) :: L !< Handle to the lua script

    !> Vector read from the Lua table.
    logical, intent(out) :: val(:)

    !> Error code describing problems encountered in each of the components.
    !! This array has to have the same length as val.
    integer, intent(out) :: ErrCode(:)

    !> A default vector to use, if no proper definition is found.
    !! Components will be filled with the help of this default definition.
    logical, intent(in), optional :: default(:)

    integer :: vect_handle
    integer :: table_len, vect_len, def_len
    integer :: vect_lb
    integer :: iComp

    ErrCode = 0

    ! Try to interpret it as table.
    vect_handle = aot_table_top(L=L)
    table_len = aot_table_length(L=L, thandle=vect_handle)

    vect_len = min(table_len, size(val))

    ! Find the length of the default value, if it is not provided, its 0.
    def_len = 0
    if (present(default)) def_len = size(default)

    ! Now parse the table with the vector entries.
    if (aot_table_first(L, vect_handle).and.(vect_len > 0)) then

      ! Only if the vector table actually exists, and has at least one entry,
      ! this parsing has to be done.
      if (present(default).and.(def_len > 0)) then
        call aot_top_get_val(val(1), ErrCode(1), L, default(1))
      else
        call aot_top_get_val(val(1), ErrCode(1), L)
      end if

      ! Up to the length of the default value, provide the default settings.
      do iComp=2,def_len
        if (.not. flu_next(L, vect_handle)) exit
        call aot_top_get_val(val(iComp), ErrCode(iComp), L, &
          &                  default(iComp))
      end do

      vect_lb = max(2, def_len+1)
      ! After def_len entries no default values for the components are
      ! available anymore, proceed without a default setting for the rest.
      do iComp=vect_lb,vect_len
        if (.not. flu_next(L, vect_handle)) exit
        call aot_top_get_val(val(iComp), ErrCode(iComp), L)
      end do

      ! If the table in the Lua script is not long enough, fill the remaining
      ! components with the default components, as far as they are defined.
      do iComp=vect_len+1,def_len
        ErrCode(iComp) = ibSet(ErrCode(iComp), aoterr_NonExistent)
        val(iComp) = default(iComp)
      end do
      vect_lb = max(vect_len+1, def_len)
      do iComp=vect_lb,vect_len
        ErrCode(iComp) = ibSet(ErrCode(iComp), aoterr_Fatal)
      end do
    else
      ! No vector definition found in the Lua script, use the default.
      ErrCode = ibSet(ErrCode, aoterr_NonExistent)
      if (present(default)) then
        val(:def_len) = default(:def_len)
        if (def_len < vect_len) then
          ErrCode(def_len+1:) = ibSet(ErrCode(def_len+1:), aoterr_Fatal)
        end if
      else
        ! No vector definition in the Lua script and no default provided.
        ErrCode = ibSet(ErrCode, aoterr_Fatal)
      end if
    end if
    call aot_table_close(L, vect_handle)

  end subroutine get_top_logical_v

end module aot_vector_module

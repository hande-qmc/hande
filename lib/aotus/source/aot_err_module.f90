! Copyright (C) 2011-2013 German Research School for Simulation Sciences GmbH,
!                         Aachen and others.
!               2013-2014 University of Siegen
! Please see the COPYRIGHT file in this directory for details.

!> This module provides the handling of errors.
module aot_err_module
  use flu_binding

  implicit none

  private

  public :: aoterr_Fatal, aoterr_NonExistent, aoterr_WrongType
  public :: aot_err_handler

  !> Some parameters for the error handling.
  !!
  !! They indicate the bits to set in case of the corresponding error, to allow
  !! appropiate reactions of the calling application.
  !! As a bitmask is used to encode the error any combination of them might be
  !! returned.
  integer, parameter :: aoterr_Fatal = 0
  integer, parameter :: aoterr_NonExistent = 1
  integer, parameter :: aoterr_WrongType = 2


contains


  !> Error handler to capture Lua errors.
  !!
  !! This routine encapsulates the retrieval of error messages from the Lua
  !! stack upon a failing Lua operation.
  !! It should be be used after all flu functions, that return an err as result.
  !! Such as flu_binding::fluL_loadfile and flu_binding::flu_pcall.
  !! The ErrString and ErrCode parameters are both optional if none of them are
  !! provided, the execution will be stopped if an error had occured and err is
  !! not 0. The error message will be written to standard output in this case.
  !!
  !! If either of them are provide, the application will continue, and the
  !! calling side has to deal with the occured error.
  subroutine aot_err_handler(L, err, msg, ErrString, ErrCode)
    type(flu_State) :: L !< Handle to the Lua script

    !> Lua error code to evaluate
    integer, intent(in) :: err

    !> Some additional message that should be prepended to the Lua error
    !! message.
    character(len=*), intent(in) :: msg

    !> Resulting error string obtained by combination of msg and the error
    !! description on the Lua stack.
    character(len=*), intent(out), optional :: ErrString

    !> The Lua error code, just the same as err.
    integer, intent(out), optional :: ErrCode

    logical :: stop_on_error
    character, pointer, dimension(:) :: string
    integer :: str_len
    integer :: i

    stop_on_error = .not.(present(ErrString) .or. present(ErrCode))

    if (present(ErrCode)) then
      ErrCode = err
    end if

    if (err .ne. 0) then

      string => flu_tolstring(L, -1, str_len)
      if (present(ErrString)) then
        ErrString = ''
        do i=1,min(str_len, len(ErrString))
          ErrString(i:i) = string(i)
        end do
      end if

      if (stop_on_error) then
        write(*,*) msg, string
        STOP
      end if

    end if

  end subroutine aot_err_handler
end module aot_err_module

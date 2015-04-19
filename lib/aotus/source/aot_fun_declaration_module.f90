! Copyright (C) 2011-2013 German Research School for Simulation Sciences GmbH,
!                         Aachen and others.
!              2013-2014 University of Siegen
! Please see the COPYRIGHT file in this directory for details.

module aot_fun_declaration_module

  implicit none

  type aot_fun_type
    integer :: handle = 0
    integer :: arg_count = 0
  end type

  public :: aot_fun_type

  private

end module aot_fun_declaration_module

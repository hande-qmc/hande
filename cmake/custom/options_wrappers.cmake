# CMake code originally written for the Psi4 project <https://github.com/psi4/psi4>
# by Lori A. Burns and Ryan M. Richard.

# Macro for printing an option in a consistent manner
#
# Syntax: print_option(<option to print> <was specified>)
#
macro(print_option variable default)
  if(NOT DEFINED ${variable} OR "${${variable}}" STREQUAL "")
    message(STATUS "Setting (unspecified) option ${variable}: ${default}")
  else()
    message(STATUS "Setting option ${variable}: ${${variable}}")
  endif()
endmacro()

# Wraps an option with default ON/OFF. Adds nice messaging to option()
#
# Syntax: option_with_print(<option name> <description> <default value>)
#
macro(option_with_print variable msge default)
  print_option(${variable} ${default})
  option(${variable} ${msge} ${default})
endmacro(option_with_print)

# Wraps an option with a default other than ON/OFF and prints it
# NOTE: Can't combine with above b/c CMake handles ON/OFF options specially
# NOTE2: CMAKE_BUILD_TYPE (and other CMake variables) are always defined so need
#       to further check for if they are the NULL string.  This is also why we
#       need the force
#
# Syntax: option_with_default(<option name> <description> <default value>)
#
macro(option_with_default variable msge default)
  print_option(${variable} ${default})
  if(NOT DEFINED ${variable} OR "${${variable}}" STREQUAL "")
    set(${variable} ${default} CACHE STRING ${msge} FORCE)
  endif()
endmacro(option_with_default)

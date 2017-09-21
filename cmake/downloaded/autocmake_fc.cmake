#.rst:
#
# Adds Fortran support.
# Appends EXTRA_FCFLAGS to CMAKE_Fortran_FLAGS.
# If environment variable FCFLAGS is set, then the FCFLAGS are used
# and no other flags are used or appended.
#
# Variables used::
#
#   EXTRA_FCFLAGS
#
# Variables defined::
#
#   CMAKE_Fortran_MODULE_DIRECTORY
#
# Variables modified::
#
#   CMAKE_Fortran_FLAGS
#
# Environment variables used::
#
#   FCFLAGS
#
# autocmake.yml configuration::
#
#   docopt:
#     - "--fc=<FC> Fortran compiler [default: gfortran]."
#     - "--extra-fc-flags=<EXTRA_FCFLAGS> Extra Fortran compiler flags [default: '']."
#   export: "'FC={0}'.format(arguments['--fc'])"
#   define: "'-DEXTRA_FCFLAGS=\"{0}\"'.format(arguments['--extra-fc-flags'])"

set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)
include_directories(${PROJECT_BINARY_DIR}/modules)

if(NOT DEFINED CMAKE_Fortran_COMPILER_ID)
    message(FATAL_ERROR "CMAKE_Fortran_COMPILER_ID variable is not defined!")
endif()

if(NOT CMAKE_Fortran_COMPILER_WORKS)
    message(FATAL_ERROR "CMAKE_Fortran_COMPILER_WORKS is false!")
endif()

if(DEFINED EXTRA_FCFLAGS)
  if(NOT EXTRA_FCFLAGS STREQUAL "")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${EXTRA_FCFLAGS}")
  endif()
endif()

if(DEFINED ENV{FCFLAGS})
    message(STATUS "FCFLAGS is set to '$ENV{FCFLAGS}'.")
    set(CMAKE_Fortran_FLAGS "$ENV{FCFLAGS}")
endif()

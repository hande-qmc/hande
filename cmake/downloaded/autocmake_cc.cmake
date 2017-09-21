#.rst:
#
# Adds C support.
# Appends EXTRA_CFLAGS to CMAKE_C_FLAGS.
# If environment variable CFLAGS is set, then the CFLAGS are used
# and no other flags are used or appended.
#
# Variables used::
#
#   EXTRA_CFLAGS
#
# Variables modified::
#
#   CMAKE_C_FLAGS
#
# Environment variables used::
#
#   CFLAGS
#
# autocmake.yml configuration::
#
#   docopt:
#     - "--cc=<CC> C compiler [default: gcc]."
#     - "--extra-cc-flags=<EXTRA_CFLAGS> Extra C compiler flags [default: '']."
#   export: "'CC={0}'.format(arguments['--cc'])"
#   define: "'-DEXTRA_CFLAGS=\"{0}\"'.format(arguments['--extra-cc-flags'])"

if(NOT DEFINED CMAKE_C_COMPILER_ID)
    message(FATAL_ERROR "CMAKE_C_COMPILER_ID variable is not defined!")
endif()

if(NOT CMAKE_C_COMPILER_WORKS)
    message(FATAL_ERROR "CMAKE_C_COMPILER_WORKS is false!")
endif()

if(DEFINED EXTRA_CFLAGS)
  if(NOT EXTRA_CFLAGS STREQUAL "")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${EXTRA_CFLAGS}")
  endif()
endif()

if(DEFINED ENV{CFLAGS})
    message(STATUS "CFLAGS is set to '$ENV{CFLAGS}'.")
    set(CMAKE_C_FLAGS "$ENV{CFLAGS}")
endif()

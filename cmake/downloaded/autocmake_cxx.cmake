#.rst:
#
# Adds C++ support.
# Appends EXTRA_CXXFLAGS to CMAKE_CXX_FLAGS.
# If environment variable CXXFLAGS is set, then the CXXFLAGS are used
# and no other flags are used or appended.
#
# Variables used::
#
#   EXTRA_CXXFLAGS
#
# Variables modified::
#
#   CMAKE_CXX_FLAGS
#
# Environment variables used::
#
#   CXXFLAGS
#
# autocmake.yml configuration::
#
#   docopt:
#     - "--cxx=<CXX> C++ compiler [default: g++]."
#     - "--extra-cxx-flags=<EXTRA_CXXFLAGS> Extra C++ compiler flags [default: '']."
#   export: "'CXX={0}'.format(arguments['--cxx'])"
#   define: "'-DEXTRA_CXXFLAGS=\"{0}\"'.format(arguments['--extra-cxx-flags'])"

if(NOT DEFINED CMAKE_CXX_COMPILER_ID)
    message(FATAL_ERROR "CMAKE_CXX_COMPILER_ID variable is not defined!")
endif()

if(NOT CMAKE_CXX_COMPILER_WORKS)
    message(FATAL_ERROR "CMAKE_CXX_COMPILER_WORKS is false!")
endif()

if(DEFINED EXTRA_CXXFLAGS)
  if(NOT EXTRA_CXXFLAGS STREQUAL "")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${EXTRA_CXXFLAGS}")
  endif()
endif()

if(DEFINED ENV{CXXFLAGS})
    message(STATUS "CXXFLAGS is set to '$ENV{CXXFLAGS}'.")
    set(CMAKE_CXX_FLAGS "$ENV{CXXFLAGS}")
endif()

#.rst:
#
# Finds LibUUID and enables it
#
# Variables modified::
#
#   DISABLE_UUID
#
# Variables defined::
#
#   UUID_FOUND
#   UUID_LIBRARIES
#   UUID_LIBRARY_DIRS
#   UUID_INCLUDE_DIRS
#
# autocmake.yml configuration::
#
#   docopt:
#     - "--libuuid=<LIBUUID> Whether to activate UUID generation [default: OFF]."
#   define: "'-DLIBUUID=\"{0}\"'.format(arguments['--libuuid'])"

if(LIBUUID)
  find_package(PkgConfig REQUIRED QUIET)
  pkg_search_module(UUID REQUIRED uuid)
  if(UUID_FOUND)
    message(STATUS "Found libuuid")
    set(USE_LIBUUID ON)
  endif()
endif()

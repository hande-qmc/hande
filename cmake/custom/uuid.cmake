#.rst:
#
# Finds LibUUID and enables it
#
# Variables used::
#
#   ENABLE_UUID
#
# Variables defined::
#
#   USE_UUID
#   UUID_FOUND
#   UUID_LIBRARIES
#   UUID_LIBRARY_DIRS
#   UUID_INCLUDE_DIRS
#
# autocmake.yml configuration::
#
#   docopt: "--uuid=<ENABLE_UUID> Whether to activate UUID generation [default: True]."
#   define: "'-DENABLE_UUID=\"{0}\"'.format(arguments['--uuid'])"

option_with_print(ENABLE_UUID "Enable usage of UUID" ON)
set(USE_UUID OFF)
if(ENABLE_UUID)
  find_package(PkgConfig REQUIRED QUIET)
  pkg_search_module(UUID REQUIRED uuid IMPORTED_TARGET)
  if(UUID_FOUND)
    message(STATUS "Found libuuid")
    set(USE_UUID ON)
  else()
    message(FATAL_ERROR "UUID requested, but libuuid was not found")
  endif()
endif()

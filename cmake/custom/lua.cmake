#.rst:
#
# Finds Lua
#
# Variables used::
#
#   LUA_ROOT
#
# autocmake.yml configuration::
#
#   docopt: "--lua=<LUA_ROOT> Specify the path to the Lua installation to use [default: '']."
#   define: "'-DLUA_ROOT=\"{0}\"'.format(arguments['--lua'])"

option_with_default(LUA_ROOT "Specify the path to the Lua installation to use" "")
if(LUA_ROOT)
  list(APPEND CMAKE_PREFIX_PATH ${LUA_ROOT})
endif()
find_package(Lua 5.3 REQUIRED)

# Append the dl library to the LUA_LIBRARIES list for a statically
# linked Lua library.
# This contraption is needed for CMake older than 3.6.0
# This is adapted from:
# https://github.com/Kitware/CMake/blob/v3.8.2/Modules/FindLua.cmake
if(CMAKE_VERSION VERSION_LESS 3.8.0)
  if(UNIX AND NOT APPLE AND NOT BEOS)
    # include dl library for statically-linked Lua library
    get_filename_component(LUA_LIB_EXT ${LUA_LIBRARY} EXT)
    if(LUA_LIB_EXT STREQUAL CMAKE_STATIC_LIBRARY_SUFFIX)
      list(APPEND LUA_LIBRARIES ${CMAKE_DL_LIBS})
    endif()
  endif()
endif()

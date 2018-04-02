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
find_package(Lua 5.3 EXACT REQUIRED)

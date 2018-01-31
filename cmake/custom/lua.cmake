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

find_package(Lua 5.3 EXACT REQUIRED)

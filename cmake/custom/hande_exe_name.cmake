#.rst:
#
# Name for the HANDE executable.
#
# Variables defined::
#
#   HANDE_EXE_NAME
#
# autocmake.yml configuration::
#
#   docopt:
#     - "--exe-name=<HANDE_EXE_NAME> [default: \"hande.cmake.x\"]."
#   define:
#     - "'-DHANDE_EXE_NAME=\"{0}\"'.format(arguments['--exe-name'])"

option_with_default(HANDE_EXE_NAME "Name for the HANDE executable" "hande.cmake.x")

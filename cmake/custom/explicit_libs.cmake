#.rst:
#
# Sets explicit linker options for extra libraries to be linked in
#
# Variables modified::
#
#   EXPLICIT_LIBS
#
# autocmake.yml configuration::
#
#   docopt:
#     - "--explicit-libs=<EXPLICIT_LIBS> Explicit linker options for extra libraries to be linked in [default: '']."
#   define:
#     - "'-DEXPLICIT_LIBS=\"{0}\"'.format(arguments['--explicit-libs'])"

option_with_default(EXPLICIT_LIBS "Explicit linker options for extra libraries to be linked in" "")
if(EXPLICIT_LIBS)
  message(STATUS "Explicit libraries: ${EXPLICIT_LIBS}")
  set(EXPLICIT_LIBS "${EXPLICIT_LIBS}" CACHE STRING "Explicit linker options for extra libraries to be linked in")
endif()

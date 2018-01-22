#.rst:
#
# Enables Lanczos diagonalisation, provided the TRLan library was found.
#
# Variables used::
#
#   ENABLE_LANCZOS
#
# Variables defined::
#
#   USE_LANCZOS
#
# autocmake.yml configuration::
#
#   docopt:
#     - "--lanczos=<TRLan_LIBRARIES> Set TRLan libraries to be linked in [default: ]."
#   define:
#     - "'-DENABLE_LANCZOS=\"{0}\"'.format(True if arguments['--lanczos'] else False)"
#     - "'-DTRLan_LIBRARIES=\"{0}\"'.format(arguments['--lanczos'])"

option_with_print(ENABLE_LANCZOS "Use Lanczos diagonalisation, requires TRLan library" OFF)
set(USE_LANCZOS OFF)
if(ENABLE_LANCZOS)
  if(NOT TRLan_LIBRARIES)
    # List of TRLan_LIBRARIES not set, abort
    message(FATAL_ERROR "You have to explicitly provide the TRLan library to link in!")
  else()
    set(USE_LANCZOS ON)
    message(STATUS "TRLan libraries: ${TRLan_LIBRARIES}")
  endif()
endif()

#.rst:
#
# Enables Lanczos diagonalisation, provided the TRLan library was found.
#
# Variables used::
#
#   TRLan_LIBRARIES
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
#     - "'-DTRLan_LIBRARIES=\"{0}\"'.format(arguments['--lanczos'])"

option_with_default(TRLan_LIBRARIES "Link line for TRLan" "")
set(USE_LANCZOS OFF)
if(TRLan_LIBRARIES)
  set(USE_LANCZOS ON)
  message(STATUS "TRLan libraries: ${TRLan_LIBRARIES}")
  set(TRLan_LIBRARIES "${TRLan_LIBRARIES}" CACHE STRING "Link line for TRLan")
endif()

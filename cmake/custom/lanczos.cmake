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
#     - "--lanczos Toggle use of Lanczos diagonalisation [default: False]."
#   define: "'-DENABLE_LANCZOS=\"{0}\"'.format(arguments['--lanczos'])"

option_with_print(ENABLE_LANCZOS "Use Lanczos diagonalisation, requires TRLan library" OFF)
set(USE_LANCZOS OFF)
if(ENABLE_LANCZOS)
  # find_package(TRLan REQUIRED)
  # if(TRLan_FOUND)
  #   set(USE_LANCZOS ON)
  # else()
  #   message(FATAL_ERROR "Lanczos diagonalisation requested, but TRLan library not found!")
  # endif()
endif()

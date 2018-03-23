#.rst:
#
# Sets the integer length HANDE will use to store walker populations.
#
# Variables modified::
#
#   HANDE_POP_SIZE
#
# autocmake.yml configuration::
#
#   docopt:
#     - "--pop-size=<HANDE_POP_SIZE> An integer among 32 or 64 [default: 32]."
#   define: "'-DHANDE_POP_SIZE=\"{0}\"'.format(arguments['--pop-size'])"

option_with_default(HANDE_POP_SIZE "Integer length for storing walker populations" 32)
list(APPEND _valid_pop_size 32 64)
if(DEFINED HANDE_POP_SIZE AND NOT HANDE_POP_SIZE IN_LIST _valid_pop_size)
  message(STATUS "${HANDE_POP_SIZE} not a valid integer length for storing walker populations, resetting to default value 32")
  set(HANDE_POP_SIZE 32 CACHE STRING "Integer length for storing walker populations" FORCE)
endif()

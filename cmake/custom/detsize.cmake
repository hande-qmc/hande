#.rst:
#
# Sets the integer length HANDE will use to represent Slater determinants as bit strings.
#
# Variables modified::
#
#   HANDE_DET_SIZE
#
# autocmake.yml configuration::
#
#   docopt:
#     - "--det-size=<HANDE_DET_SIZE> An integer among 32 or 64 [default: 32]."
#   define: "'-DHANDE_DET_SIZE=\"{0}\"'.format(arguments['--det-size'])"

option_with_default(HANDE_DET_SIZE "Integer length for representing Slater determinants as bit strings" 32)
list(APPEND _valid_det_size 32 64)
if(DEFINED HANDE_DET_SIZE AND NOT HANDE_DET_SIZE IN_LIST _valid_det_size)
  message(STATUS "${HANDE_DET_SIZE} not a valid integer length for representing Slater determinants, resetting to default value 32")
  set(HANDE_DET_SIZE 32 CACHE STRING "Integer length for representing Slater determinants as bit strings" FORCE)
endif()

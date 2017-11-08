#.rst:
#
# Enable single precision arithmetics, where appropriate.
#
# Variables used::
#
#   ENABLE_SINGLE_PRECISION
#
# Variables defined::
#
#   USE_SINGLE_PRECISION
#
# autocmake.yml configuration::
#
#   docopt:
#     - "--single Enable usage of single precision, where appropriate [default: False]."
#   define: "'-DENABLE_SINGLE_PRECISION=\"{0}\"'.format(arguments['--single'])"

option_with_print(ENABLE_SINGLE_PRECISION "Enable usage of single precision arithmetics, where appropriate" OFF)
set(USE_SINGLE_PRECISION OFF)
if(ENABLE_SINGLE_PRECISION)
  set(USE_SINGLE_PRECISION ON)
endif()

#.rst:
#
# Enable backtrace functionality.
#
# Variables used::
#
#   ENABLE_BACKTRACE
#
# Variables defined::
#
#   USE_BACKTRACE
#
# autocmake.yml configuration::
#
#   docopt:
#     - "--backtrace Enable backtrace functionality [default: False]."
#   define: "'-DENABLE_BACKTRACE=\"{0}\"'.format(arguments['--backtrace'])"

option_with_print(ENABLE_BACKTRACE "Enable use of backtrace functionality" OFF)
set(USE_BACKTRACE OFF)
if(ENABLE_BACKTRACE)
  find_package(Backtrace REQUIRED)
  if(Backtrace_FOUND)
    set(USE_BACKTRACE ON)
  else()
    message(FATAL_ERROR "Backtrace functionality requested, but unsupported or not found")
  endif()
endif()

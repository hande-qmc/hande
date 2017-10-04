#.rst:
#
# Finds TRLan and enables it.
#
# Variables modified::
#
#   DISABLE_LANCZOS
#
# autocmake.yml configuration::
#
#   docopt:
#     - "--trlan=<TRLAN> Toggle use of TRLan [default: OFF]."
#   define: "'-DTRLAN=\"{0}\"'.format(arguments['--trlan'])"

if(TRLAN)
  # Find TRLan and enable it
else()
  add_definitions(-DDISABLE_LANCZOS)
endif()

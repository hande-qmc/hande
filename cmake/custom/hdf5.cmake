#.rst:
#
# Enables and finds HDF5.
#
# Variables modified::
#
#   DISABLE_HDF5
#
# autocmake.yml configuration::
#
#   docopt:
#     - "--hdf5=<HDF5> C++ compiler [default: OFF]."
#   define: "'-DHDF5=\"{0}\"'.format(arguments['--hdf5'])"

if(HDF5)
  find_package(HDF5 COMPONENTS Fortran REQUIRED)
  set(USE_HDF5 ON)
else()
  add_definitions(-DDISABLE_HDF5)
endif()

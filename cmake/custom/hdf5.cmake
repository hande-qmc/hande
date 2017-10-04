#.rst:
#
# Finds HDF5 and enables it, provided that the Fortran 2003 interface was
# compiled.
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
  find_package(HDF5 1.8.15 COMPONENTS Fortran REQUIRED)
  set(USE_HDF5 ON)
  message(STATUS "Parallel HDF5? ${HDF5_IS_PARALLEL}")
  # Was the Fortran 2003 interface to HDF5 enabled?
  # Compile an example from the HDF5 website:
  # https://support.hdfgroup.org/HDF5/examples/f-src.html
  set(scratch_directory
    ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/hdf5)
  try_compile(has_Fortran2003 ${scratch_directory} ${PROJECT_SOURCE_DIR}/cmake/custom/compound_complex_fortran2003.f90)
  if(NOT has_Fortran2003)
    message(FATAL_ERROR "HDF5 was not compiled with --enable-fortran2003")
  endif()
else()
  add_definitions(-DDISABLE_HDF5)
endif()

#.rst:
#
# Finds HDF5 and enables it, provided that the Fortran 2003 interface was
# compiled.
#
# Variables used::
#
#   ENABLE_HDF5
#
# Variables defined::
#
#   USE_HDF5
#   HDF5_Fortran_LIBRARIES
#   HDF5_INCLUDE_DIRS
#
# autocmake.yml configuration::
#
#   docopt: "--hdf5=<HDF5> Enable HDF5 [default: True]."
#   define: "'-DENABLE_HDF5=\"{0}\"'.format(arguments['--hdf5'])"

option_with_print(ENABLE_HDF5 "Enable usage of HDF5 (requires Fortran 2003 bindings)" ON)
set(USE_HDF5 OFF)
if(ENABLE_HDF5)
  find_package(HDF5 1.8.15 COMPONENTS Fortran REQUIRED)
  if(HDF5_IS_PARALLEL)
    message(STATUS "Parallel HDF5 FOUND")
  else()
    message(STATUS "Parallel HDF5 NOT FOUND")
  endif()
  # Was the Fortran 2003 interface to HDF5 enabled?
  # Compile an example from the HDF5 website:
  # https://support.hdfgroup.org/HDF5/examples/f-src.html
  set(scratch_directory ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/hdf5-f03)
  try_compile(HDF5_Fortran2003
    ${scratch_directory}
    SOURCES
      ${CMAKE_CURRENT_LIST_DIR}/compound_complex_fortran2003.f90
    CMAKE_FLAGS
      -DINCLUDE_DIRECTORIES=${HDF5_INCLUDE_DIRS}
    LINK_LIBRARIES
      ${HDF5_Fortran_LIBRARIES}
    )
  if(NOT HDF5_Fortran2003)
    message(FATAL_ERROR "HDF5 requested, but library was not compiled with --enable-fortran2003")
  else()
    set(USE_HDF5 ON)
  endif()
endif()

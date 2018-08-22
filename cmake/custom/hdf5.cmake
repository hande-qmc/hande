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
  # Was the Fortran 2003 interface to HDF5 enabled?
  # Compile an example from the HDF5 website:
  # https://support.hdfgroup.org/HDF5/examples/f-src.html
  set(scratch_directory ${CMAKE_CURRENT_BINARY_DIR}/feature_tests/HDF5_HAS_Fortran2003-test)
  try_compile(HDF5_HAS_Fortran2003
    ${scratch_directory}
    SOURCES
      ${CMAKE_CURRENT_LIST_DIR}/test_hdf5_has_fortran2003.f90
    CMAKE_FLAGS
      -DINCLUDE_DIRECTORIES=${HDF5_INCLUDE_DIRS}
    LINK_LIBRARIES
      ${HDF5_Fortran_LIBRARIES}
    OUTPUT_VARIABLE
      HDF5_HAS_Fortran2003-test-output
    )
  if(NOT HDF5_HAS_Fortran2003)
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/HDF5_HAS_Fortran2003-test.log ${HDF5_HAS_Fortran2003-test-output})
    message(FATAL_ERROR "HDF5 requested, but library was not compiled with --enable-fortran2003 \
Compiling a simple test executable failed, consult the log file: ${CMAKE_CURRENT_BINARY_DIR}/HDF5_HAS_Fortran2003-test.log")
  else()
    set(USE_HDF5 ON)
  endif()
endif()

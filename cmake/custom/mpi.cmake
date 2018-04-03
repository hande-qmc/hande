#.rst:
#
# Enables MPI support and usage of ScaLAPACK.
# ScaLAPACK will be searched for in system paths if the Intel MKL BLAS/LAPACK
# implementation is found. The user will have to specify the link line for
# ScaLAPACK if any other implementation of the libraries should be used.
#
# Variables used::
#
#   ENABLE_MPI
#   ENABLE_SCALAPACK
#
# Variables defined::
#
#   USE_MPI
#   USE_ScaLAPACK
#
# Variables modified::
#
#   CMAKE_Fortran_COMPILER
#   CMAKE_C_COMPILER
#
# autocmake.yml configuration::
#
#   docopt:
#     - "--mpi Enable MPI parallelization [default: False]."
#     - "--mpi-with-scalapack Enable ScaLAPACK usage with MPI [default: False]."
#   define:
#     - "'-DENABLE_MPI=\"{0}\"'.format(arguments['--mpi'])"
#     - "'-DENABLE_SCALAPACK=\"{0}\"'.format(arguments['--mpi-with-scalapack'])"

option_with_print(ENABLE_MPI "Enable MPI parallelization" OFF)
set(USE_MPI OFF)
if(ENABLE_MPI)
  find_package(MPI REQUIRED)
  if(MPI_Fortran_FOUND AND MPI_C_FOUND)
    set(USE_MPI ON)
    set(scratch_directory ${CMAKE_CURRENT_BINARY_DIR}/feature_tests/HAS_MPI3_features-test)
    try_compile(HAS_MPI3_features
      ${scratch_directory}
      SOURCES
        ${CMAKE_CURRENT_LIST_DIR}/test_MPI-3_features_simple.f90
      CMAKE_FLAGS
        -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
        -DINCLUDE_DIRECTORIES=${MPI_Fortran_INCLUDE_PATH}
        -DCMAKE_Fortran_FLAGS=${MPI_Fortran_COMPILE_FLAGS}
      LINK_LIBRARIES
        ${MPI_Fortran_LINK_FLAGS}
        ${MPI_Fortran_LIBRARIES}
      )
    if(HAS_MPI3_features)
      message(STATUS "MPI-3 available!")
    endif()
  else()
    message(FATAL_ERROR "You asked for MPI, but CMake could not find any MPI installation, check $PATH")
  endif()
endif()

include(CMakeDependentOption)
cmake_dependent_option(
  ENABLE_SCALAPACK "Enable usage of ScaLAPACK" OFF
  "NOT ENABLE_MPI" ON
  )
if(ENABLE_SCALAPACK)
  set(USE_ScaLAPACK ON)
endif()

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
    # Check for MPI-3
    set(HAS_MPI3 OFF)
    # Lazily enough, the check for MPI-3 is done only on Fortran
    # Version of CMake before 3.10 were not able to find the MPI version implemented,
    # so we do this by ourselves adapting from:
    # https://github.com/Kitware/CMake/blob/v3.11.1/Modules/FindMPI.cmake
    if(NOT DEFINED MPI_Fortran_VERSION)
      set(WORK_DIR "${CMAKE_CURRENT_BINARY_DIR}/feature_tests/MPI_version-test")
      set(SRC_DIR  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/custom")
      set(BIN_FILE "${WORK_DIR}/mpiver_Fortran.bin")
      set(MPI_Fortran_INCLUDE_LINE "use mpi_f08\n      implicit none")
      configure_file("${SRC_DIR}/mpiver.f90.in" "${WORK_DIR}/mpiver.f90" @ONLY)
      set(MPI_TEST_SOURCE_FILE "${WORK_DIR}/mpiver.f90")
      try_compile(MPI_RESULT_Fortran_mpiver
        "${CMAKE_CURRENT_BINARY_DIR}" SOURCES "${MPI_TEST_SOURCE_FILE}"
        CMAKE_FLAGS
          -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
          -DINCLUDE_DIRECTORIES=${MPI_Fortran_INCLUDE_PATH}
          -DCMAKE_Fortran_FLAGS=${MPI_Fortran_COMPILE_FLAGS}
        LINK_LIBRARIES
          ${MPI_Fortran_LINK_FLAGS}
          ${MPI_Fortran_LIBRARIES}
        COPY_FILE "${BIN_FILE}")
      if(MPI_RESULT_Fortran_mpiver)
        file(STRINGS ${WORK_DIR}/mpiver_Fortran.bin _MPI_VERSION_STRING LIMIT_COUNT 1 REGEX "INFO:MPI-VER")
        if("${_MPI_VERSION_STRING}" MATCHES ".*INFO:MPI-VER\\[([0-9]+)\\.([0-9]+)\\].*")
          set(MPI_Fortran_VERSION_MAJOR "${CMAKE_MATCH_1}")
          set(MPI_Fortran_VERSION_MINOR "${CMAKE_MATCH_2}")
          set(MPI_Fortran_VERSION "${MPI_Fortran_VERSION_MAJOR}.${MPI_Fortran_VERSION_MINOR}")
        endif()
      endif()
      unset(WORK_DIR)
      unset(SRC_DIR)
      unset(BIN_FILE)
      unset(MPI_Fortran_INCLUDE_LINE)
    endif()
    if(MPI_Fortran_VERSION VERSION_GREATER 3.0)
      set(HAS_MPI3 ON)
    endif()
    if(HAS_MPI3)
      message(STATUS "MPI-3 available!")
    endif()
  else()
    message(FATAL_ERROR "You asked for MPI, but CMake could not find any MPI installation, check $PATH")
  endif()
endif()

option_with_print(ENABLE_SCALAPACK "Enable usage of ScaLAPACK" OFF)
# Usage of ScaLAPACK is conditional to:
#   1. Having enabled and detected a working MPI implementation
#   2. Having explicitly requested to link against ScaLAPACK
# Thus, we introduced an internal _enable_scalapack dependent option
# as an helper here.
include(CMakeDependentOption)
cmake_dependent_option(
  _enable_scalapack "Enable usage of ScaLAPACK" OFF
  "ENABLE_MPI;USE_MPI;NOT ENABLE_SCALAPACK" ON
  )
set(USE_ScaLAPACK OFF)
if(_enable_scalapack)
  set(USE_ScaLAPACK ON)
endif()

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
  else()
    message(FATAL_ERROR "You asked for MPI, but CMake could not find any MPI installation, check $PATH")
  endif()
endif()

include(CMakeDependentOption)
cmake_dependent_option(
  ENABLE_SCALAPACK "Enable usage of ScaLAPACK" OFF
  "ENABLE_MPI" ON
  )
if(ENABLE_SCALAPACK)
  set(USE_ScaLAPACK ON)
endif()

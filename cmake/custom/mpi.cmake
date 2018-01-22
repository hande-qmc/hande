#.rst:
#
# Enables MPI support and usage of ScaLAPACK.
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
#   define:
#     - "'-DENABLE_MPI=\"{0}\"'.format(arguments['--mpi'])"
#     - "'-DENABLE_SCALAPACK=\"{0}\"'.format(arguments['--mpi'])"

option_with_print(ENABLE_MPI "Enable MPI parallelization" OFF)
set(USE_MPI OFF)
if(ENABLE_MPI)
  find_package(MPI REQUIRED)
  if(MPI_Fortran_FOUND AND MPI_C_FOUND)
    set(USE_MPI ON)
    set(USE_ScaLAPACK ON)
    set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
    set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
  else()
    message(FATAL_ERROR "You asked for MPI, but CMake could not find any MPI installation, check $PATH")
  endif()
endif()

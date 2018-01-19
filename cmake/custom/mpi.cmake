#.rst:
#
# Enables MPI support and usage of ScaLAPACK.
#
# Variables used::
#
#   ENABLE_MPI
#   ENABLE_ScaLAPACK
#
# Variables defined::
#
#   USE_MPI
#   USE_ScaLAPACK
#   ScaLAPACK_LIBRARIES
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
#     - "--scalapack=<ScaLAPACK_LIBRARIES> Set ScaLAPACK libraries to be linked in [default: ]."
#   define:
#     - "'-DENABLE_MPI=\"{0}\"'.format(arguments['--mpi'])"
#     - "'-DENABLE_ScaLAPACK=\"{0}\"'.format(arguments['--mpi'])"
#     - "'-DScaLAPACK_LIBRARIES=\"{0}\"'.format(arguments['--scalapack'])"

option_with_print(ENABLE_MPI "Enable MPI parallelization" OFF)
set(USE_MPI OFF)
if(ENABLE_MPI)
  find_package(MPI REQUIRED)
  if(MPI_Fortran_FOUND AND MPI_C_FOUND)
    set(USE_MPI ON)
    set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
    set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
  else()
    message(FATAL_ERROR "You asked for MPI, but CMake could not find any MPI installation, check $PATH")
  endif()
endif()

option_with_print(ENABLE_ScaLAPACK "Use ScaLAPACK (parallel compilation only)" OFF)
set(USE_ScaLAPACK OFF)
if(ENABLE_ScaLAPACK)
  if(NOT ScaLAPACK_LIBRARIES)
    # List of ScaLAPACK_LIBRARIES not set, abort
    message(FATAL_ERROR "You have to explicitly provide the list of ScaLAPACK libraries to link in!")
  elseif(NOT USE_MPI)
    # Valid MPI implementation not found
    message(FATAL_ERROR "ScaLAPACK needed for MPI build, but MPI was not enabled")
  else()
    set(USE_ScaLAPACK ON)
  endif()
endif()

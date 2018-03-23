#.rst:
#
# Enables OpenMP support.
#
# Variables used::
#
#   ENABLE_OPENMP
#
# Variables defined::
#
#   USE_OPENMP
#
# autocmake.yml configuration::
#
#   docopt: "--omp Enable OpenMP parallelization [default: False]."
#   define: "'-DENABLE_OPENMP=\"{0}\"'.format(arguments['--omp'])"

option_with_print(ENABLE_OPENMP "Enable OpenMP parallelization" OFF)
set(USE_OPENMP OFF)
if(ENABLE_OPENMP)
  find_package(OpenMP REQUIRED)
  set(USE_OPENMP ON)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
endif()

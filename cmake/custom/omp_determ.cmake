#.rst:
#
# Enables OpenMP support.
#
# Variables used::
#
#   OPENMP_DETERM
#
# Variables defined::
#
# autocmake.yml configuration::
#
#   docopt: "--omp-determ Force OpenMP parallelization to be deterministic [default: False]."
#   define: "'-DOPENMP_DETERM=\"{0}\"'.format(arguments['--omp-determ'])"

option_with_print(OPENMP_DETERM "Deterministic OpenMP parallelization" OFF)
if(OPENMP_DETERM)
  find_package(OpenMP REQUIRED)
endif()

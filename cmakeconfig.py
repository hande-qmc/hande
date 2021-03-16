#!/usr/bin/env python

# This file is autogenerated by Autocmake v1.0.0 http://autocmake.org
# Copyright (c) 2015-2021 by Radovan Bast, Roberto Di Remigio, Jonas Juselius, and contributors.

import os
import sys
assert sys.version_info >= (2, 6), 'Python >= 2.6 is required'

sys.path.insert(0, 'cmake')
from autocmake import configure
from autocmake.external import docopt


options = """
Usage:
  ./cmakeconfig.py [options] [<builddir>]
  ./cmakeconfig.py (-h | --help)

Options:
  --fc=<FC>                              Fortran compiler [default: gfortran].
  --extra-fc-flags=<EXTRA_FCFLAGS>       Extra Fortran compiler flags [default: ''].
  --cc=<CC>                              C compiler [default: gcc].
  --extra-cc-flags=<EXTRA_CFLAGS>        Extra C compiler flags [default: ''].
  --python=<PYTHON_INTERPRETER>          The Python interpreter (development version) to use. [default: ''].
  --add-definitions=<STRING>             Add preprocesor definitions [default: ''].
  --lua=<LUA_ROOT>                       Specify the path to the Lua installation to use [default: ''].
  --mpi                                  Enable MPI parallelization [default: False].
  --mpi-with-scalapack                   Enable ScaLAPACK usage with MPI [default: False].
  --omp                                  Enable OpenMP parallelization [default: False].
  --blas=<BLAS>                          Detect and link BLAS library (auto or off) [default: auto].
  --lapack=<LAPACK>                      Detect and link LAPACK library (auto or off) [default: auto].
  --mkl=<MKL>                            Pass MKL flag to the Intel compiler and linker and skip BLAS/LAPACK detection (sequential, parallel, cluster, or off) [default: off].
  --scalapack=<SCALAPACK_LIBRARIES>      Link line for ScaLAPACK libraries [default: '']
  --blacs=<BLACS_IMPLEMENTATION>         Implementation of BLACS for MKL ScaLAPACK (openmpi, intelmpi, sgimpt) [default: openmpi]
  --explicit-libs=<EXPLICIT_LIBS>        Explicit linker options for extra libraries to be linked in [default: ''].
  --dsfmt-mexp=<HANDE_DSFMT_MEXP>        An integer among 521, 1279, 2203, 4253, 11213, 19937, 44497, 86243, 1322049, 216091 [default: 19937].
  --det-size=<HANDE_DET_SIZE>            An integer among 32 or 64 [default: 32].
  --pop-size=<HANDE_POP_SIZE>            An integer among 32 or 64 [default: 32].
  --exe-name=<HANDE_EXE_NAME>            [default: "hande.cmake.x"].
  --hdf5=<HDF5_ROOT>                     Specify the path to the HDF5 installation to use [default: ''].
  --uuid=<UUID>                          Whether to activate UUID generation [default: True].
  --single                               Enable usage of single precision, where appropriate [default: False].
  --backtrace                            Enable backtrace functionality [default: False].
  --popcnt                               Enable use of intrinsic popcnt [default: False].
  --type=<TYPE>                          Set the CMake build type (debug, release, relwithdebinfo, minsizerel) [default: release].
  --generator=<STRING>                   Set the CMake build system generator [default: Unix Makefiles].
  --show                                 Show CMake command and exit.
  --cmake-executable=<CMAKE_EXECUTABLE>  Set the CMake executable [default: cmake].
  --cmake-options=<STRING>               Define options to CMake [default: ''].
  --prefix=<PATH>                        Set the install path for make install.
  <builddir>                             Build directory.
  -h --help                              Show this screen.
"""


def gen_cmake_command(options, arguments):
    """
    Generate CMake command based on options and arguments.
    """
    command = []
    command.append(arguments['--cmake-executable'])
    command.append('-DCMAKE_Fortran_COMPILER={0} -DEXTRA_FCFLAGS="{1}"'.format(arguments['--fc'], arguments['--extra-fc-flags']))
    command.append('-DCMAKE_C_COMPILER={0} -DEXTRA_CFLAGS="{1}"'.format(arguments['--cc'], arguments['--extra-cc-flags']))
    command.append('-DPYTHON_INTERPRETER="{0}"'.format(arguments['--python']))
    command.append('-DPREPROCESSOR_DEFINITIONS="{0}"'.format(arguments['--add-definitions']))
    command.append('-DLUA_ROOT="{0}"'.format(arguments['--lua']))
    command.append('-DENABLE_MPI="{0}"'.format(arguments['--mpi']))
    command.append('-DENABLE_SCALAPACK="{0}"'.format(arguments['--mpi-with-scalapack']))
    command.append('-DENABLE_OPENMP="{0}"'.format(arguments['--omp']))
    command.append('-DENABLE_BLAS={0}'.format(arguments['--blas']))
    command.append('-DENABLE_LAPACK={0}'.format(arguments['--lapack']))
    command.append('-DMKL_FLAG={0}'.format(arguments['--mkl']))
    command.append('-DMATH_LIB_SEARCH_ORDER="MKL;ESSL;OPENBLAS;ATLAS;ACML;SYSTEM_NATIVE"')
    command.append('-DBLAS_LANG=Fortran')
    command.append('-DLAPACK_LANG=Fortran')
    command.append('-DSCALAPACK_LIBRARIES="{0}"'.format(arguments['--scalapack']))
    command.append('-DBLACS_IMPLEMENTATION="{0}"'.format(arguments['--blacs']))
    command.append('-DEXPLICIT_LIBS="{0}"'.format(arguments['--explicit-libs']))
    command.append('-DHANDE_DSFMT_MEXP="{0}"'.format(arguments['--dsfmt-mexp']))
    command.append('-DHANDE_DET_SIZE="{0}"'.format(arguments['--det-size']))
    command.append('-DHANDE_POP_SIZE="{0}"'.format(arguments['--pop-size']))
    command.append('-DHANDE_EXE_NAME="{0}"'.format(arguments['--exe-name']))
    command.append('-DENABLE_HDF5="{0}"'.format('FALSE' if arguments['--hdf5'] in ['False', 'false', 'OFF', 'off'] else 'TRUE'))
    command.append('-DHDF5_ROOT="{0}"'.format(arguments['--hdf5']))
    command.append('-DENABLE_UUID="{0}"'.format(arguments['--uuid']))
    command.append('-DENABLE_SINGLE_PRECISION="{0}"'.format(arguments['--single']))
    command.append('-DENABLE_BACKTRACE="{0}"'.format(arguments['--backtrace']))
    command.append('-DENABLE_INTRINSIC_POPCNT="{0}"'.format(arguments['--popcnt']))
    command.append('-DCMAKE_BUILD_TYPE={0}'.format(arguments['--type']))
    command.append('-G"{0}"'.format(arguments['--generator']))
    if arguments['--cmake-options'] != "''":
        command.append(arguments['--cmake-options'])
    if arguments['--prefix']:
        command.append('-DCMAKE_INSTALL_PREFIX="{0}"'.format(arguments['--prefix']))

    return ' '.join(command)


# parse command line args
try:
    arguments = docopt.docopt(options, argv=None)
except docopt.DocoptExit:
    sys.stderr.write('ERROR: bad input to {0}\n'.format(sys.argv[0]))
    sys.stderr.write(options)
    sys.exit(-1)


# use extensions to validate/post-process args
if configure.module_exists('extensions'):
    import extensions
    arguments = extensions.postprocess_args(sys.argv, arguments)


root_directory = os.path.dirname(os.path.realpath(__file__))


build_path = arguments['<builddir>']


# create cmake command
cmake_command = '{0} -H{1}'.format(gen_cmake_command(options, arguments), root_directory)


# run cmake
configure.configure(root_directory, build_path, cmake_command, arguments)

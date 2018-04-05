.. _compilation-with-cmake:

Compilation using CMake
=======================

It is possible to configure and build HANDE using CMake. At least version 3.6
of CMake is required. You can get CMake either *via* your package manager or by
downloading an executable tarball from `here <https://cmake.org/download/>`_
Unpacking and adding to your ``PATH`` will do the trick:

.. code-block:: bash

   $ curl -L https://cmake.org/files/v3.10/cmake-3.10.2-Linux-x86_64.tar.gz | tar -xz
   $ export PATH=$HOME/Software/cmake-3.10.2-Linux-x86_64/bin${PATH:+:$PATH}

where we have assumed the tarball was downloaded in the ``$HOME/Software``
directory.

Once dependencies are installed, you can configure HANDE either by running the
``cmake`` executable directly:

.. code-block:: bash

   $ cmake -H. -Bbuild

or by using the frontend Python script ``cmakeconfig.py``:

.. code-block:: bash

   $ ./cmakeconfig.py build

The result of using the two methods is exactly the same: a subdirectory
``build`` will be created containing the build system.
Using the frontend script however results in more compact configure lines.

After ensuring HANDE's dependencies are installed, produce a makefile by running the
``cmakeconfig.py`` (residing in the tools subdirectory) script in the root directory:

Configuration options
---------------------

Building of HANDE can be tweaked in various ways passing options to the
frontend script (or CMake directly).
The help menu for the ``cmakeconfig.py`` script shows the available options:

.. code-block:: bash

   Usage:
     ./cmakeconfig.py [options] [<builddir>]
     ./cmakeconfig.py (-h | --help)

   Options:
     --fc=<FC>                              Fortran compiler [default: gfortran].
     --extra-fc-flags=<EXTRA_FCFLAGS>       Extra Fortran compiler flags [default: ''].
     --cc=<CC>                              C compiler [default: gcc].
     --extra-cc-flags=<EXTRA_CFLAGS>        Extra C compiler flags [default: ''].
     --python=<PYTHON_INTERPRETER>          The Python interpreter (development version) to use. [default: ''].
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
     --hdf5=<HDF5>                          Enable HDF5 [default: True].
     --uuid=<UUID>                          Whether to activate UUID generation [default: True].
     --lanczos=<TRLan_LIBRARIES>            Set TRLan libraries to be linked in [default: ''].
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

These options are translated to CMake native options. For more detailed information on
HANDE-specific compile-time settings, see :ref:`compilation-settings`. The following list
is a translation guide between the frontend script and "bare" CMake:

- ``--fc=FC``/``-DCMAKE_Fortran_COMPILER=FC``. To set the Fortran compiler. Default
  is ``gfortran``.
- ``--extra-fc-flags="list-of-flags"``/``-DEXTRA_FCFLAGS="list-of-flags"``. To set additional flags
  for the Fortran compiler.
- ``--cc=CC``/``-DCMAKE_C_COMPILER=CC``. To set the C compiler. Default is ``gcc``.
- ``--extra-cc-flags="list-of-flags"``/``-DEXTRA_CFLAGS="list-of-flags"``. To set additional flags
  for the C compiler.
- ``--python=INTERP``/``-DPYTHON_INTERPRETER=INTERP``. To set the Python interpreter. The
  default is empty, so that CMake will attempt to find a suitable version.
- ``--lua=LUA``/``-DLUA_ROOT=LUA``. To set the Lua installation to use. Minimum
  required version of Lua is 5.3. The default is empty, so that CMake will attempt to
  find a suitable version.
  See below for Lua detection issues.

  .. warning::

     CMake will not pick up Lua from a nonstandard location, even though it is on
     path (any or all of ``CPATH``, ``LIBRARY_PATH``, ``LD_LIBRARY_PATH``,
     ``PATH``)

- ``--mpi``/``-DENABLE_MPI=ON``. Enables MPI parallelization. CMake will
  attempt to find a suitable implementation of MPI and set the compilers
  accordingly.

  .. warning::

     To use a specific MPI implementation, pass the appropriate MPI compiler
     wrappers as arguments to ``--fc`` (``-DCMAKE_Fortran_COMPILER``) and
     ``--cc`` (``-DCMAKE_C_COMPILER``)

- ``--mpi-with-scalapack``/``-DENABLE_SCALAPACK=OFF``. Enables linking to
  ScaLAPACK. This requires that MPI is enabled and that a ScaLAPACK
  implementation is available.
- ``--omp``/``-DENABLE_OPENMP=ON``. Enables OpenMP parallelization. CMake will
  check which flags are supported by your choice of compilers and add them to
  the compiler flags.
- ``--blas=auto``/``-DENABLE_BLAS=auto``. Triggers autodetection of BLAS libraries.
  See below for math libraries detection issues.
- ``--lapack=auto``/``-DENABLE_LAPACK=auto``. Triggers autodetection of BLAS libraries.
  See below for math libraries detection issues.
- ``--mkl=VALUE``/``-DMKL_FLAG=VALUE``. Sets the ``-mkl=VALUE`` flag for the Intel
  compiler and linker. Valid values are ``sequential``, ``parallel``, ``cluster``, or
  ``off``, with ``off`` being the default.

  .. warning::

     Passing this option overrides automatic math detection

- ``--scalapack="link-line"``/``-DSCALAPACK_LIBRARIES="link-line"``. Link line for ScaLAPACK libraries.
  If using Intel MKL, CMake will be able to correctly locate and set these for
  you. Use this option in case you run into trouble with detecting ScaLAPACK
  and prefer setting the link line explictly.
- ``--blacs=openmpi``/``-DBLACS_IMPLEMENTATION=openmpi``. Sets the implementation of
  BLACS for the Intel MKL ScaLAPACK libraries. Valid values are ``openmpi``,
  ``intelmpi`` and ``sgimpt``, with ``openmpi`` being the default.
- ``--explicit-libs="link-line"``/``-DEXPLICIT_LIBS="link-line"``. Sets explicit linker options for
  extra libraries to be linked in.
  See below for math libraries detection issues.
- ``--dsfmt-mexp=VALUE``/``-DHANDE_DSFMT_MEXP=VALUE``. Set exponent for the period of the
  Mersenne Twister (MT) random number generator (RNG). Valid values are 521,
  1279, 2203, 4253, 11213, 19937, 44497, 86243, 1322049, and 216091. with 19937
  being the default.
- ``--det-size=VALUE``/``-DHANDE_DET_SIZE=VALUE``. Set the integer length for representing
  Slater determinants as bit strings. Valid values are 32 and 64, with 32
  being the default.
- ``--pop-size=VALUE``/``-DHANDE_POP_SIZE=VALUE``. Set the integer length for storing
  walker populations. Valid values are 32 and 64, with 32
  being the default.
- ``--exe-name=NAME``/``-DHANDE_EXE_NAME=NAME``. Set the name for the generated HANDE executable.
  The default is ``hande.cmake.x``. The executable is copied to the ``bin``
  directory in the root of the project and symlinked to ``hande.x``. Passing
  the executable name will let you preserve executables generated with
  different configuration settings.
- ``--hdf5=<ON/OFF>``/``-DENABLE_HDF5=<ON/OFF>``. Enables use of HDF5. By
  default, this is turned on. At least HDF5 1.8.15 is required and with Fortran
  2003 bindings enabled. CMake will search for a suitable version of HDF5 and
  check that all necessary components are available.
  See below for HDF5 detection issues.
- ``--uuid=<ON/OFF>``/``-DENABLE_UUID=<ON/OFF>``. Enables use of the UUID library.
  By default, this is turned on.
- ``--lanczos="link-line"``/``-DTRLan_LIBRARIES="link-line"``. Set the TRLan
  libraries to be linked in. By default empty, thus disabling use of TRLan.
- ``--single``/``-DENABLE_SINGLE_PRECISION=ON``. Enables use of single
  precision, where appropriate.
- ``--backtrace``/``-DENABLE_BACKTRACE=ON``. Enables backtrace functionality.
- ``--popcnt``/``-DENABLE_INTRINSIC_POPCNT=ON``. Enables usage of popcnt
  intrinsic (requires hardware support)
- ``--type=debug``/``-DCMAKE_BUILD_TYPE=Debug``. Switches build type. Valid
  values are ``debug``, ``release``, ``releasewithdebinfo`` and ``minsizerel``.
  The default is a debug build.
- ``--cmake-options="-DTHIS -DTHAT"``. Sets options to be forwarded as-is to
  CMake.

CMake compilation issues
------------------------

When dependencies are not in standard search paths, CMake needs to be nudged
and pointed in the right direction. This can be done directly using either ``cmake`` or
``cmakeconfig``; the equivalent commands for both are given below but only one should be
used.

- Detection of math libraries is usually the trickiest part. The CMake math
  detection scripts shipped with HANDE rely on the ``MATH_ROOT`` environment
  variable being set to point to the root of the math libraries installation
  you want to use.
  The detection scripts will attempt to provide a link line for math libraries
  based on the search order in the CMake variable ``MATH_LIB_SEARCH_ORDER``.
  By default, Intel MKL is searched for first, using the ``MKLROOT``
  environment variable.
  If math detection fails, libraries can be set manually:

  .. code-block:: bash

     $ ./cmakeconfig.py --blas=off --lapack=off --explicit-libs="-L/usr/lib -lblas -llapack"
     $ cmake -H. -DENABLE_BLAS=OFF -DENABLE_LAPACK=OFF -DEXPLICIT_LIBS="-L/usr/lib -lblas -llapack"

- Lua in a non-standard directory. Exporting the root directory of the Lua
  installation as ``LUA_ROOT`` (or ``LUA_DIR``) or directly passing it as an option:

  .. code-block:: bash

     $ ./cmakeconfig.py --lua=/install/dir/for/Lua build
     $ cmake -H. -Bbuild -DLUA_ROOT=/install/dir/for/Lua

- HDF5 in a non-standard directory. Exporting the root directory of the HDF5
  installation as ``HDF5_ROOT`` os directly passing it as an option:

  .. code-block:: bash

     $ ./cmakeconfig.py --hdf5 --cmake-options="-DHDF5_ROOT=/install/dir/for/HDF5" build
     $ cmake -H. -Bbuild -DENABLE_HDF5=ON -DHDF5_ROOT=/install/dir/for/HDF5

For compiler- and library-specific issues, see :ref:`compiler-issues`.

Compiling with MPI
------------------

To compile with MPI it is necessary to pass **both** the ``--mpi`` option
**and** the correct compiler wrappers with the ``--cc`` and ``--fc``:

.. code-block:: bash

   $ ./cmakeconfig.py --mpi --fc=mpif90 --cc=mpicc
   $ cmake -H. -Bbuild -DCMAKE_Fortran_COMPILER=mpif90 -DCMAKE_C_COMPILER=mpicc -DENABLE_MPI=ON

CMake can in fact botch the identification of the compiler wrappers and MPI
libraries, a mismatch that will result in linker errors.
Here are some examples of configuration lines. In all cases, remember to set
the ``MATH_ROOT`` variable to point to the location of the math libraries:

- OpenMPI with GNU compilers.

  .. code-block:: bash

     $ ./cmakeconfig.py --mpi --fc=mpif90 --cc=mpicc
     $ cmake -H. -Bbuild -DCMAKE_Fortran_COMPILER=mpif90 -DCMAKE_C_COMPILER=mpicc -DENABLE_MPI=ON

- OpenMPI with Intel compilers.

  .. code-block:: bash

     $ ./cmakeconfig.py --mpi --fc=mpif90 --cc=mpicc
     $ cmake -H. -Bbuild -DCMAKE_Fortran_COMPILER=mpif90 -DCMAKE_C_COMPILER=mpicc -DENABLE_MPI=ON

- IntelMPI with Intel compiler.

  .. code-block:: bash

     $ ./cmakeconfig.py --mpi --fc=mpiifort --cc=mpiicc
     $ cmake -H. -Bbuild -DCMAKE_Fortran_COMPILER=mpiifort -DCMAKE_C_COMPILER=mpiicc -DENABLE_MPI=ON

- OpenMPI with GNU compilers and OpenBLAS ScaLAPACK.

  .. code-block:: bash

     $ ./cmakeconfig.py --mpi --fc=mpif90 --cc=mpicc --mpi-with-scalapack --scalapack="-L/location/of/scalapack -lscalapack"
     $ cmake -H. -Bbuild -DCMAKE_Fortran_COMPILER=mpif90 -DCMAKE_C_COMPILER=mpicc -DENABLE_MPI=ON -DENABLE_SCALAPACK=ON -DSCALAPACK_LIBRARIES="-L/location/of/scalapack -lscalapack"

- OpenMPI with Intel compilers and MKL ScaLAPACK. The math detection script
  will use the OpenMPI implementation of BLACS by default.

  .. code-block:: bash

     $ ./cmakeconfig.py --mpi --fc=mpif90 --cc=mpicc --mpi-with-scalapack
     $ cmake -H. -Bbuild -DCMAKE_Fortran_COMPILER=mpif90 -DCMAKE_C_COMPILER=mpicc -DENABLE_MPI=ON -DENABLE_SCALAPACK=ON

- IntelMPI with Intel compiler and MKL ScaLAPACK. In this case we need to tell
  CMake what BLACS implementation to use with ScaLAPACK.

  .. code-block:: bash

     $ ./cmakeconfig.py --mpi --fc=mpiifort --cc=mpiicc --mpi-with-scalapack --blacs=intelmpi
     $ cmake -H. -Bbuild -DCMAKE_Fortran_COMPILER=mpiifort -DCMAKE_C_COMPILER=mpiicc -DENABLE_MPI=ON -DENABLE_SCALAPACK=ON -DBLACS_IMPLEMENTATION=intelmpi

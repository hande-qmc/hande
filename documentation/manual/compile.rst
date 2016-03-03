Compilation
===========

After ensuring HANDE's dependencies are installed, produce a makefile by running the
``mkconfig.py`` (residing in the tools subdirectory) script in the root directory:

.. code-block:: bash

    tools/mkconfig.py config/conf

where ``conf`` is one of the platforms available and is simply the name of the relevant
file residing in the ``config/`` directory.  Various configurations are provided and it is
simple to adapt one to the local environment (e.g. changing compiler or library paths).

Run

.. code-block:: bash

    tools/mkconfig.py --help

to see the options available, including inspecting available configurations.

A configuration is defined using a simple ini file, consisting of three sections:
main, opt and dbg.  For instance::

    [main]
    fc = gfortran
    ld = gfortran
    libs = -llapack -lblas

    [opt]
    fflags = -O3

    [dbg]
    fflags = -g

Any option not specified in the 'opt' and 'dbg' sections is inherited from the
'main' section.  The settings in 'opt' are used by default; the debug options
can be selected by passing the -g option to mkconfig.

All options are strings unless otherwise specified.  Available options are:

fc
    Set the fortran compiler.
fflags
    Set flags to be passed to the fortran compiler during compilation.
f90_module_flag
    Set the flag used by the compiler which is used to specify the directory
    where module (.mod) files are placed when created and where they should be
    searched for.
f90_module_flag_pad [boolean]
    True if a space needs to be inserted between the defined f90_module_flag
    and the corresponding directory argument.  Default: true.
cc
    Set the C compiler.
cflags
    Set flags to be passed to the C compiler during compilation.
ccd
    Set the C compiler used to generate the C dependency files.  Only required
    if cc doesn't support -MM and -MT flags.  Default: use cc.
cdflags
    Set the flags for the c++ compiler used to generate the C++ dependency files.
    Default: $CFLAGS -MM -MT
cxx
    Set the C++ compiler.
cxxflags
    Set flags to be passed to the C++ compiler during compilation.
cxxd
    Set the C compiler used to generate the C++ dependency files.  Only required
    if cc doesn't support -MM and -MT flags.  Default: use cxx.
cxxdflags
    Set the flags for the c++ compiler used to generate the C++ dependency files.
    Default: $CXXFLAGS -MM -MT
cpp
    Set the C preprocessor to be used on Fortran source files.  If not defined
    then the Fortran compiler is used to do the preprocessing.
cppflags
    Set flags to be used in the C preprocessing step.
    C preprocessing is applied to .F90, .F, .c and .cpp files (and not .f90
    files).
ld
    Set the linker program.
ldflags
    Set flags to be passed to the linker during linking of the compiled objects.
libs
    Set libraries to be used during the linking step.
ar
    Set the archive program.  Default: ar.
arflags
    Set the flags to be passed to the archive program.  Default: -rcs.

To compile the code run 

.. code-block:: bash

    make
    
HANDE's build system uses the ``sfmakedepend`` script (http://people.arsc.edu/~kate/Perl/,
supplied in ``tools/``) by Kate Hedstrom to generate the list of dependencies for each
Fortran source file.  These are generated automatically when make is run if the dependency
files do not exist.

The executable, ``hande.x``, is placed in the ``bin`` subdirectory.  Note that this is
actually a symbolic link: a unique executable is produced for each platform and
optimisation level and ``hande.x`` merely points to the most recently compiled executable
for convenience.  This makes testing against multiple platforms particularly easy.

There are various goals in the makefile.  Run

.. code-block:: bash

    make help

to see the available goals.

Compile-time settings
---------------------

The behaviour of the program can be changed in various ways by some choices at
compile-time by using C pre-processing.  These choices largely influence the
speed, memory usage, inclusion of parallel code and workarounds for certain
compilers.

The pre-processing options which accept a value are set by::

    -DOPTION=VAL

which defines the pre-processing definition OPTION to have value VAL.
Similarly, the options which just need to be defined to be used are set by::

    -DOPTION

These should be added to the cppflags or cppdefs lines in the configuration
files or in the Makefile, as desired.


.. warning::

    Certain options, for technical reasons, change the Markov chain of QMC calculations.
    Results should be in statistical agreement but the precise data produced by the
    calculation (even using the same random number seed) may well be changed.

    This currently applies to the following options: POP_SIZE and
    SINGLE_PRECISION.

DET_SIZE
    Default: 32.

    HANDE uses bit strings to store Slater determinants, where each bit
    corresponds to an occupied spin-orbital if the bit is set and an unoccupied
    spin orbital otherwise.  As Fortran does not include a type for a single
    bit, integers are used.  Note that this does lead to some wasted memory when
    the number of spin-orbitals is not a multiple of the size of the integer used.
    An array of integers is used to store the determinant bit string if
    a single integer is not sufficient.

    This option sets the integer length to be used.  Allowed values are 32 and
    64, corresponding to using 32-bit and 64-bit integers respectively.  As bit
    operations on a 64-bit integer are faster than those on two 32-bit integers,
    using DET_SIZE=64 is recommended for production calculations.  (Note,
    however, that this will use more memory than DET_SIZE=32 if the number of
    basis functions is closer to a multiple of 32 rather than 64.  This is
    rarely a concern in practice.)
POP_SIZE
    Default: 32

    This option is used to specify whether 32 or 64-bit integers are used to
    store walker populations in HANDE. It is unlikely that 64-bit integers will
    be needed when using the integer code but this option is more critical
    when the **real_amplitudes** option is being used. When using the
    **real_amplitudes** option with POP_SIZE=32, the largest walker amplitude
    that can be stored is :math:`2^{20}=1048576`, while the smallest fractional part that
    can be represented is :math:`2^{-11}=0.00049`. When using this option and POP_SIZE=64
    the largest amplitude is :math:`2^{32}=4.3\times10^9` and the smallest fractional part
    is :math:`2^{-31}=4.66\times10^{-10}`.
DEBUG
    Default: not defined.

    If defined then add additional information in output (e.g. stack traces) that might be
    useful for debugging.  Recommended for developers only.  The format and content of the
    additional debug output should not be relied upon.
DISABLE_LANCZOS
    Default: not defined.

    If defined then Lanczos diagonalisation is disabled.  This removes the dependency on the TRLan library.
DISABLE_HDF5
    Default: not defined.

    If defined then the QMC restart functionality is disabled and the dependency on HDF5
    (which can be tricky to compile on some machines) is removed.  Note that restart
    functionality is extremely useful in production simulations so this option should
    only be used during initial porting efforts.
DISABLE_UUID
    Default: not defined.

    If defined then each calculation will not print universally unique identifier. This removes the
    dependency on libuuid.
DSFMT_MEXP 
    Default: 19937.

    HANDE uses the dSFMT random number generator (RNG).  It is based on
    a Mersenne Twister algorithm, is extremely fast and produces high quality
    random numbers.  See http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/SFMT/index.html 
    for more details. 

    DSFMT_EXP sets the exponent of the period of the RNG.  Allowed values are
    521, 1279, 2203, 4253, 11213, 19937, 44497, 86243,
    132049 and 216091 and lead to, for example, random numbers with a period of
    a Mersenne Prime such as :math:`2^{512}-1`.
NAGF95  
    Default: not defined.

    If defined then code specific to, and necessary for compilation using, the
    NAG Fortran compiler is included.
PARALLEL  
    Default: not defined.

    Include source code required for running in parallel.
SINGLE_PRECISION  
    Default: not defined.

    Set the precision (where possible) to be single precision.  The default is
    double precision.  This is faster, but (of course) can change results
    significantly.  Use with care.
USE_POPCNT
    Default: not defined.

    Use the intrinsic popcnt function instead of the version implemented in HANDE.

    An important procedure involves counting the number of set bits in an integer.  HANDE
    includes a very efficient, branchless procedure to do this.  However, the Fortran
    2008 standard includes an intrinsic function, popcnt, for this exact operation.
    The performance of this intrinsic will be implementation-dependent and, with
    standard compilation flags, we expect the HANDE version to be competitive or more
    performant, based upon some simple tests.  The key difference is on modern
    processors containing the popcnt instruction: the popcnt intrinsic can then
    make use of this instruction and will be much faster than the implementation
    in HANDE.  The existence of the popcnt instruction can be found, on Unix
    and Linux platforms, by inspecting the flags field in ``/proc/cpuinfo``: if
    it contains ``popcnt``, then the processor contains the popcnt instruction.

    Using the popcnt instruction often involves a compiler-specific flag to
    tell the compiler to use that instruction set; often compilers include the
    popcnt instruction with the flag that specifies the use of the SSE4.2
    instruction set.  The use of the popcnt instruction can be tested using
    objdump.  For example:

    .. code-block:: bash

        $ objdump -d bin/hande.x | grep popc
        0000000000400790 <__popcountdi2@plt>:
          400931:e8 5a fe ff ff         callq  400790 <__popcountdi2@plt>

    indicates that HANDE is using a compiler-supplied function for popcnt.  Exact output
    (especially the function name) is compiler dependent.  In contrast:

    .. code-block:: bash

        $ objdump -d bin/hande.x | grep popc
          4008ac:f3 0f b8 c0            popcnt %eax,%eax

    indicates HANDE is using the popcnt instruction.  If the above command does not give
    any output, then USE_POPCNT has most likely not been defined.

Compilation issues
------------------

We attempt to work round any compiler and library issues we encounter but sometimes this
is not possible.  Issues and, where known, workarounds we have found are:

* HDF5 1.8.14 (and possibly 1.8.13) has a bug revealed by Intel compilers v15 onwards.
  This results in unusual error messages and/or segmentation faults when writing out
  restart files.  Possibly workarounds:

  * use HDF5 1.8.15 (best).
  * recompile HDF5 with ``-assume nostd_value``.
  * recompile HDF5 with an earlier version of the Intel compilers.
  * recompile HANDE with HDF5 support disabled.

* Compiling with GCC and linking the Intel MKL library leads to segmentation faults or
  incorrect answers for FCI calculation on systems with complex-valued integrals when
  run in parallel.  Either use a different ScaLAPACK library, or use the Intel compilers.

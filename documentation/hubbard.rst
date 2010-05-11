hubbard
=======

Introduction
------------

hubbard can currently perform Full Configuration Interaction (FCI) calculations
of the Hubbard model using either the real space or momentum space formulation.
Full and Lanczos diagonalisation methods are implemented using external
libraries (lapack/scalapack and TRLan respectively) and can be performed in
both serial and parallel.  Lanczos diagonalisation can be performed with or
without precomputing the Hamiltonian matrix.

Recently the Full Configuration Interaction Quantum Monte Carlo method has been
implemented.

Directory structure
--------------------

./
    Root directory of the program.
bin/
  Directory containing the compiled program, hubbard.x.  Created during
  compilation.
config/
  Directory containing the configuration input files used to generate makefiles.
dest/
  Directory containing the compiled object files and dependency files.  Created
  during compilation.
documentation/
   Directory containing documentation on the hubbard program.  The
   documentation is written in reStructured Text and can be converted
   into a wide range of output formats.
src/
    Directory containing the main source files.
lib/
   Directory containing "library" source files.  These a procedures which are
   not specific to the hubbard code but are generally useful.  Some are written
   by the authors, some are freely available (as noted in the source files).
tools/
    Directory containing scripts and tools for compiling, running and analysing
    output from hubbard.

Compilation
-----------

hubbard requires the lapack (http://www.netlib.org/lapack/), blas
(http://www.netlib.org/blas) and TRLan
(http://crd.lbl.gov/~kewu/ps/trlan_.html) libaries.

After meeting these requirements, produce a makefile by running the mkconfig.py
(residing in the tools subdirectory) script in the root directory:

.. code-block:: bash

    tools/mkconfig.py config

where config is one of the platforms available.  The config name is simply the
name of the relevant file residing in the config/ directory.  Various configurations
are provided and it is simple to adapt one to the local environment (e.g. changing
compiler or library paths).

Run

.. code-block:: bash

    tools/mkconfig.py --help

to see the options available, including inspecting available platforms.
A platform is defined using a simple ini file, consisting of three sections:
main, opt and dbg.  For instance::

    [main]
    cc = gfortran
    ld = gfortran
    libs = -llapack -lblas

    [opt]
    cflags = -O3

    [dbg]
    cflags = -g

Any option not specified in the 'opt' and 'dbg' sections is inherited from the
'main' section.  The optimised settings in 'opt' are used by default; the debug
options can be selected by passing the -g option to mkconfig.py.

Available options are:

fc
    Set the fortran compiler.
fflags
    Set flags to be passed to the fortran compiler during compilation.
cppdefs
    Set definitions to be used in the C pre-processing step.
cppflags
    Set flags to be used in the C pre-processing step.
ld
    Set the linker program.
ldflags
    Set flags to be passed to the linker during linking of the compiled objects.
libs
    Set libraries to be used during the linking step.
module_flag
    Set the compiler-specific flag which specifies the directory where module
    (.mod) files are placed when created and where they should be searched for.

To compile the code run 

.. code-block:: bash

    make
    
hubbard.x uses the sfmakedepend script (http://www.arsc.edu/~kate/Perl/,
supplied in tools/) by Kate Hedstrom to generate the dependencies.  These are
generated automatically when make is run if the dependency files don't exist.

The executable, hubbard.x, is placed in the bin subdirectory.  Note that this is
actually a symbolic link: a unique executable is produced for each platform and
optimisation level and hubbard.x merely points to the most recently compiled executable
for convenience.  This makes testing against multiple platforms particularly easy.

There are various goals in the makefile.  Run

.. code-block:: bash

    make help

to see the available goals.

Compile-time settings
^^^^^^^^^^^^^^^^^^^^^

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

DET_SIZE
    Default: 32.

    hubbard uses bit strings to store Slater determinants, where each bit
    corresponds to an occupied spin-orbital if the bit is set and an unoccupied
    spin orbital otherwise.  As fortran does not include a type for a single
    bit, integers are used.  Note that this does lead to some wasted memory when
    the number of spin-orbitals is not a multiple of the size of the integer used.
    An array of integers is used to store the determinant bit string if
    a single integer is not sufficient.

    This option sets the integer length to be used.  Allowed values are 8, 16,
    32 and 64, corresponding to using 8-bit, 16-bit, 32-bit and 64-bit integers
    respectively.  Note that using 8-bit or 16-bit integers is much slower on
    modern platforms.  The recommended value is 32 unless more than 32 basis
    functions are used, in which case 64 is also a good choice.  The parallel FCIQMC
    algorithm requires the determinant bit-strings to be made up of either 32- or 64-bit
    integers.
32BIT
    Default: not defined.

    Must be defined if using 64-bit integers as the determinant bit-strings
    with a 32-bit compiler for performing parallel FCIQMC calculations.
DSFMT_MEXP 
    Default: 19937.

    hubbard uses the dSFMT random number generator (RNG).  It is based on
    a Mersenne Twister algorithm, is extremely fast and produces high quality
    random numbers.  See http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/SFMT/index.html 
    for more details. 

    DSFMT_EXP sets the exponent of the period of the RNG.  Allowed values are
    521, 1279, 2203, 4253, 11213, 19937, 44497, 86243,
    132049 and 216091 and lead to, for example, random numbers with a period of
    a Mersenne Prime such as 2^512-1.
NAGF95  
    Default: not defined.

    If defined then code specific to, and necessary for compilation using, the
    NAG Fortran compiler is included.
PGI  
    Default: not defined.

    If defined then code required to work around a bug in the PGI compiler (only 
    version 10.1 was tested) is included.  This is required for successful
    compilation if DET_SIZE is set to be 8 or 16.
PARALLEL  
    Default: not defined.

    Include source code required for running in parallel.
SINGLE_PRECISION  
    Default: not defined.

    Set the precision (where possible) to be single precision.  The default is
    double precision.  This is faster, but (of course) can change results
    significantly.  Use with care.

Usage
-----

.. code-block:: bash

    hubbard.x [input_filename]

If no input filename is provided then the input options are read from STDIN.
Note that this feature is not guaranteed to work when run in parallel!

Output is sent to STDOUT and can be redirected as desired.

hubbard.x only performs i/o operations on the root processor when run on
multiple processors.

Input options
-------------

Input options are case insensitive and can be given in any order.  A new line
is required for each keyword.  Keywords are given in **bold** text.  Items
following a keyword that are in *italics* are given as input values to that
keyword.  Optional arguments are enclosed in square brackets.

With the exception of the **lattice** keyword, all values associated with
a specific keyword should appear on the same line as that keyword.

Items enclosed in parentheses are treated as comments.  All input options are
echoed in the output and so comments allow for notes on the calculation to be
made in the input which are then automatically included in the output.

The current input options allow the formulation of the Hubbard model, the
system parameters and the nature of the calculation to be given.

System type
^^^^^^^^^^^

These options select the type of system to use.

**k_space**
    Default system type.

    Use the momentum space formulation of the Hubbard model.  Slater
    determinants are formed in the basis of Bloch functions :math:`\psi_k`:

    .. math::

        \psi_k(r) = e^{ik.r} \sum_i \phi_i(r)

    where :math:`\phi_i(r)` is the basis function centred on site :math:`i`.
**momentum_space**
    Synonym for **k_space**.
**real_space**
    Use the real space formulation of the Hubbard model.  Slater determinants
    are formed from the basis functions, :math:`\phi_i`, which are each centred
    on a lattice site.  Periodic boundary conditions are imposed through the
    kinetic 'hopping' term in the Hamiltonian.

System
^^^^^^

These options describe the system which is to be investigated.

**electrons** *nel*
    Integer.

    Required.

    Set the number of electrons in the system to be *nel*.
**lattice** *lattice vectors*
    Integer matrix.

    Required.

    Set the lattice vectors (and as a result the dimensionality) of the system.
    The lines immediately after **lattice** are assumed to be the :math:`n
    \times n` matrix containing the lattice vectors of the crystal cell (i.e.
    one lattice vector per line).  1D, 2D and 3D systems can be specified using
    vectors of the appropriate dimensionality.
**nel** *nel*
    Synonym for **electrons**.
**T** *t*
    Real.

    Default: 1.

    Set the kinetic term in the Hamiltonian to be *t*, i.e. the kinetic operator is:

    .. math::

        T = -t \sum_{i,j,\sigma} a_{i\sigma}^{\dagger} a_{j\sigma}

**U** *U*
    Real.

    Default: 1.

    Set the Coulomb term in the Hamiltonian to be *U*.
**twist** *t1 [t2 [t3]]*
    Real.

    Default: 0.0.

    Apply a twist to the wavevector grid.  The twist is an *ndim*-dimensional
    vector in units of :math:`2\pi`.  The twist angle should be within the
    first Brillouin zone, and hence the components should be between -0.5 and
    +0.5.

    Applicable only in the momentum space formulation of the Hubbard model.

Calculation type
^^^^^^^^^^^^^^^^

The following options select which kind of calculation(s) are performed on the
chosen system.  If no calculation type is given, then only the calculation
initialisation (mainly the enumeration of the basis) is performed.

**exact**
    Perform an full diagonalisation of the Hamiltonian matrix.
**fci**
    Synonym for **exact**.
**simple_fciqmc**
    Perform an FCIQMC calculation using an extremely simple (but wasteful, in
    terms of CPU and memory resources) algorithm.  This should be used for testing only.
**fciqmc**
    Perform an FCIQMC calculation.
**ifciqmc**
    Perform an initiator-FCIQMC calculation.
**lanczos**
    Perform a Lanczos diagonalisation of the Hamiltonian matrix.
**lanczos_direct**
    Perform a Lanczos diagonalisation of the Hamiltonian matrix but calculate
    the required Hamiltonian matrix elements on the fly rather than
    pre-computing the entire Hamiltonian matrix (as is done with **lanczos**).
    This is slower but requires much less memory.  This is currently only
    implemented in serial.

Calculation options: symmetry options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

FCI calculations consider the full Hamiltonian matrix.  This is automatically
contructed in a block diagonal form via the use of symmetry, allowing for the
Hamiltonian matrix to be considered a block at a time.  This results in
a substantial reduction in CPU and memory demands.  The default behaviour is to
diagonalise all blocks of the Hamiltonian matrix but this can be controlled by
the following options.

In contrast, an FCIQMC calculation can only consider a single block of the
Hamiltonian matrix.  The spin polarisation must be specified and the symmetry
of the determinant is currently hard-coded.

**ms** *ms*
    Integer.

    Diagonalise only blocks containing determinants with the specified value of Ms,
    in units of electron spin (i.e. 1/2).
**symmetry** *isym*
    Integer.

    Only relevant for the momentum space formulation.  Diagonalise only blocks containing
    determinants of the same symmetry as the specified symmetry block *isym*.  *isym* refers
    to a wavevector label (as given in the output).  To see the symmetry labels for a specific
    crystal cell, run the calculation without any calculation type specified.  The :math:`\Gamma`
    wavevector is always given by *isym*:math:`=1`.
**sym** *isym*
    Synonmym for **symmetry**.

Calculation options: diagonalisation options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These options are only valid when a diagalisation (either full or Lanczos)
calculation is performed.

**eigenvalues**
    Default behaviour.

    Find only the eigenvalues of the Hamiltonian matrix.
**eigenvectors**
    Find the eigenvectors and eigenvalues of the Hamiltonian matrix.  This is
    much slower.  Currently the eigenvectors are not used or even outputted.

Calculation options: Lanczos options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These options are only valid when a Lanczos diagonalisation calculation is
performed.

**lanczos_basis** *nbasis*
    Integer.

    Default: 40.

    Set the number of Lanczos vectors to be used.  This determines the main
    memory requirements of the Lanczos routine.  The size of the basis can have
    an impact on the performance of the Lanczos diagonalisation and which
    excited eigensolutions are found.  See the TRLan documentation,
    http://crd.lbl.gov/~kewu/ps/trlan_.html, for more details.
**lanczos_solutions** *nsolns*
    Integer.

    Default: 5.  

    Set the number of eigenvalues (and eigenvectors, if required) to be found
    via Lanczos diagonlisation.  The Hamiltonian matrix is constructed in block
    diagonal form using spin and crystal momentum conservation rules.  nsolns
    is the number of solutions found per block.
**lanczos_solns** *nsolns*
    Synonym for **lanczos_solutions**.

Calculation options: FCIQMC options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following options are valid for FCIQMC calculations.

**mc_cycles** *mc_cycles*
    Integer.

    Number of Monte Carlo cycles to perform per "report loop".
**nreports** *nreports*
    Integer.

    Number of "report loops" to perform.  Each report loop consists of 
    *mc_cycles* cycles of the FCIQMC algorithm followed by updating the shift
    and output of information on the current state of the walker populations, in
    particular the instantaneous energy estimators.

    The total number of Monte Carlo cycles performed in an FCIQMC calculation
    is *nreports* x *mc_cycles*.
**seed** *seed*
    Integer.

    Default: 7.

    Set the seed used to initialise the dSFMT random number generator.
    In parallel the seed on each processor is *seed* + iproc, where iproc is
    the processor index (as supplied by MPI) and ranges from 0 to nprocs-1.
**tau** *tau*
    Real.

    Set the timestep to be used.  Each Monte Carlo cycle amounts to propogating
    the walker population by the *tau* in units of imaginary time.

    A small timestep causes the walker population to evolve very slowly.  Too
    large a timestep, on the other hand, leads to a rapid particle growth which
    takes a long time to stabilise, even once the shift begins to vary, and
    coarse population dynamics.
**initial_shift** *initial_shift*
    Real.

    Default: 0.

    Set the value of the shift to use during the period before the shift is
    allowed to vary.  Positive values lead to faster growth in the number of
    walkers due to cloning.  Using too large a value can lead to poor sampling
    as large numbers of walkers reside on the same small number of determinants
    rather than diffusing appropriately through the determinant space.
**varyshift_target** *varyshift_target*
    Integer.

    Default: 10000.

    Set the target number of particles to be reached before the shift is
    allowed to vary.  This is only checked at the end of each report loop.
**shift_damping** *xi*
    Real.

    Default: 0.05.

    Once the *varyshift_target* has been reached, the shift is updated according to:

    .. math::

        S(\beta) = S(\beta-A*\tau) - \xi*log(N_w(\tau)/N_w(\beta-A*\tau))/(A*\tau)

    where :math:`\beta` is the current imaginary time, :math:`A\tau` is the
    amount of imaginary time between shift updates, :math:`N_w` is the number of
    walkers at the given time and :math:`xi` is a damping factor to prevent
    wild fluctations in the population dynamics and can be set using the
    **shift_damping** keyword.
**reference_det** *electron_1 electron_2 ... electron_nel*
    Integer list.

    Default: use the first nalpha alpha spin-orbitals and first nbeta beta
    spin-orbitals, where nalpha and nbeta are the number of alpha and beta
    electrons respectively, as defined by the **ms** input option.  Note that
    this can lead to using a 'bad' reference determinant which is a long way
    from the ground state energy. This is particularly true when using the real
    space formulation of the Hubbard model, as it causes as many sites as
    possible to be doubly occupied.

    Set the reference determinant to occupy the specified spin-orbitals.
    The index of each spin-orbital is printed out in the basis functions
    section of the output.  This will be overridden by a restart file.
**reference_det_population** *pop*
    Integer.

    Default: 10.

    Set the initial walker population on the reference determinant.  This will
    be overridden by a restart file.
**walker_length** *walker_length*
    Integer.

    Size of walker array.  This is allocated at the start of the calculation
    and is used to store the population of walkers on determinants with
    a non-zero population and the associated energy of the determinant.

    Care: this needs to be large enough to hold the number of unique
    determinants with a non-zero population of walkers in the simulation.  The
    code does not currently check whether this size is exceeded and so setting
    **walker_length** to be too small can lead to memory problems and
    segmentation faults.  For large calculations this should be substantial
    smaller than the full size of determinant space.

    Not valid for simple_fciqmc calculations, where the population of walkers
    on each determinant is stored.
**spawned_walker_length** *spawned_walker_length*
    Integer.

    Size of the spawned walker array.  This is allocated at the start of the
    calculation and is used to store the population of spawned walkers on child
    determinants.

    Care: this needs to be large enough to store all the particles which are spawned
    during a Monte Carlo cycle and so needs to be a reasonable fraction of the 
    targetted number of total number of walkers.  The code does not currently
    check whether this size is exceeded and so setting
    **spawned_walker_length** to be too small can lead to memory problems and
    segmentation faults.

    Not valid for simple_fciqmc calculations, where the population of spawned
    walkers on each determinant is stored.
**dump_restart** [*id*]
    Optional integer.

    Write out information required for restarting an FCIQMC calculation to
    a file called restart.x, where x is *id* if *id* is given.  Otherwise 
    x is chosen to be the smallest integer possible such that restart.x does
    not exist in the calculation directory.

    Restart is currently only implemented in serial.

    Warning: these files can become very large, so care should be taken when
    not re-using the same filenames.
**restart** [*id*]
    Optional integer.

    Restart an FCIQMC calculation using a previous restart file, restart.x,
    where x is a non-negative integer.  If *id* is given, then the file
    restart.id is used, otherwise x is chosen to be the largest integer such
    that restart.x exists and restart.x+1 does not.

    The restart file does not contain system information such as the U and
    T parameter, lattice vectors, number of electrons or if the walker
    population were evolved using standard FCIQMC or initiator-FCIQMC. Thus it
    is important use the same system parameters when restarting a calculation.
    The consistency of the restart file with the input options supplied is not
    checked.
    
    Please note that the RNG is not stored in the restart file, so running two
    shorter calculations via the restart facility is not completely identical
    to running a single calculation for the same number of Monte Carlo cycles.

    Furthermore, the current implementation does not allow restart files
    produced with one value of DET_SIZE to be used with binaries produced with
    a different value of DET_SIZE.  However, this is not checked!

Calculation options: initiator-FCIQMC options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In addition to the options for general FCIQMC calculations, the following
options are also valid in initiator-FCIQMC calculations:

**initiator_population** *population*
    Integer.

    Default: 3.

    Set the population at which a determinant is considered to be an initiator
    determinant.  Setting this value to 0 retrieves the FCIQMC result.
**cas** *N* *M*
    Integers.

    Default: 0 0.

    Set the complete active space (CAS) to be (*N*, *M*), which defines the CAS
    such that the lowest *nel* - *N* spin-orbitals are core (occupied)
    spin-orbitals; precisely *N* electrons occupy the next 2 *M* "active"
    spin-orbitals and the remaining spin-orbitals form the "external" space and
    are unoccupied.  Any determinant within the CAS is considered to be an
    initiator determinant, no matter what the population of walkers on that
    determinant.

    A CAS of (0,0) contains only the determinant with the *nel* lowest energy
    spin-orbitals occupied and a CAS of (*nel*, *norbs*) contains the full
    space of determinants, where *norbs* is the number of spin-orbitals used in
    the simulation (i.e. twice the number of sites in the crystal cell in the
    case of the Hubbard model).

    Note that the CAS is somewhat meaningless when using the real space
    formulation of the Hubbard model (as the spin-orbitals used as the basis do
    not have an associated energy) and so great care should be used.

Calculation options: parallel options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These options control the behaviour when run in parallel.  They do not affect
the result but can have a significant impact on performance.

**block_size** *block_size*
    Integer.

    Default: 64.

    Set the block size used to distribute the Hamiltonian matrix across the
    processors.  The Hamiltonian matrix is divided into :math:`n \times n`
    sub-matrices, where :math:`n` is the block size, which are the distributed
    over the processors in a cyclic fashion.  Applicable only to FCI
    calculations.

output options
^^^^^^^^^^^^^^

These options increase the verbosity but can be useful for debugging.  Note that
the filesizes scale factorially with system size.  These should not currently
be used in parallel.

**determinants** [*filename*]
    Optional character string.

    Default: off.  Default filename: DETS.

    Write out the enumerated list of determinants to the given *filename* or
    to the default filename if no filename is give.
**det** [*filename*]
    Synonym for **determinants**.
**hamiltonian** [*filename*]
    Optional character string.

    Default: off.  Default filename: HAMIL.

    Write out the diagonal and the non-zero off-diagonal elements of the
    Hamiltonian matrix to the given *filename*, or to the default filename if
    not filename is given.
**hamil** [*filename*]
    Synonym for **hamiltonian**.

other options
^^^^^^^^^^^^^

**end**
    End of input.  Any subsequent lines in an input file are ignored.  It is
    only strictly required if the input is given via STDIN.

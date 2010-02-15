hubbard
=======

Introduction
------------

hubbard can currently perform FCI calculations of the Hubbard model using
either the real space or momentum space formulation.  Full and Lanczos
diagonalisation methods are implemented using external libraries
(lapack/scalapack and TRLan respectively) and can be performed in both serial
and parallel.  Lanczos diagonalisation can be performed with or without
precomputing the Hamiltonian matrix.

The target is to implement an efficient FCIQMC algorithm for the Hubbard model.

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
'main' section.  The settings in 'opt' are used by default; the debug options
can be selected by passing the -g option to mkconfig.py.

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
    Synonym for **electrons**
**T** *t*
    Integer.

    Default: 1.

    Set the kinetic term in the Hamiltonian to be *t*.
**U** *U*
    Integer.

    Default: 1.

    Set the Coulomb term in the Hamiltonian to be *U*.

Calculation type
^^^^^^^^^^^^^^^^

The following options select which kind of calculation(s) are performed on the
chosen system.  If no calculation type is given, then only the calculation
initialisation (mainly the enumeration of the basis) is performed.

**exact**
    Perform an full diagonalisation of the Hamiltonian matrix.
**fci**
    Synonym for **exact**.
**lanczos**
    Perform a Lanczos diagonalisation of the Hamiltonian matrix.
**lanczos_direct**
    Perform a Lanczos diagonalisation of the Hamiltonian matrix but calculate
    the required Hamiltonian matrix elements on the fly rather than
    pre-computing the entire Hamiltonian matrix (as is done with **lanczos**).
    This is slower but requires much less memory.  Currently only implemented
    in serial.

Calculation options: symmetry options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

FCI calculations consider the full Hamiltonian matrix.  This is automatically
contructed in a block diagonal form via the use of symmetry, allowing for the
Hamiltonian matrix to be considered a block at a time.  This results in
a substantial reduction in CPU and memory demands.  The default behaviour is to
diagonalise all blocks of the Hamiltonian matrix but this can be controlled by
the following options.

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
    over the processors in a cyclic fashion.

Output options
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

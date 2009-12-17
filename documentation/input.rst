Hubbard
=======

Hubbard can currently perform FCI calculations of the Hubbard model using either
the real space or momentum space formulation of the Hubbard Hamiltonian.  Full and
Lanczos diagonalisation methods are implemented using external libraries (lapack/scalapack
and trlan respectively) and can be performed in both serial and parallel.  Lanczos diagonalisation can be performed with or without precomputing the Hamiltonian matrix.

Input options
-------------

Input options are case insensitive.

Items enclosed in parentheses are treated as comments.  All input options are
echoed in the output and so comments allow for notes on the calculation to be
made in the input which are then automatically included in the output.

The current input options allow the formulation of the Hubbard model, the system
parameters and the nature of the calculation to be given.

System type
^^^^^^^^^^^

REAL_SPACE
K_SPACE, MOMENTUM_SPACE

System
^^^^^^

LATTICE
NEL
ELECTRONS
T
U

Calculation type
^^^^^^^^^^^^^^^^

EXACT, FCI
LANCZOS
LANCZOS_DIRECT

Calculation options
^^^^^^^^^^^^^^^^^^^

LANCZOS_BASIS
LANCZOS_SOLUTIONS
LANCZOS_SOLNS
BLOCK_SIZE
EIGENVALUES
EIGENVECTORS

Output options
^^^^^^^^^^^^^^

HAMIL, HAMILTONIAN
DET, DETERMINANTS

Other options
^^^^^^^^^^^^^

END

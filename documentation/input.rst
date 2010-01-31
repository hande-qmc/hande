Hubbard
=======

Hubbard can currently perform FCI calculations of the Hubbard model using
either the real space or momentum space formulation.  Full and Lanczos
diagonalisation methods are implemented using external libraries
(lapack/scalapack and trlan respectively) and can be performed in both serial
and parallel.  Lanczos diagonalisation can be performed with or without
precomputing the Hamiltonian matrix.

Input options
-------------

Input options are case insensitive and can be given in any order.  A new line
is required for each keyword.  Keywords are given in **bold** text.  Items
following a keyword that are in *italics* are given as input values to that
keyword.  Optional items are enclosed in square brackets.

With the exception of the **lattice** keyword, all values associated with
a specific keyword should appear on the same line as that keyword.

Items enclosed in parentheses are treated as comments.  All input options are
echoed in the output and so comments allow for notes on the calculation to be
made in the input which are then automatically included in the output.

The current input options allow the formulation of the Hubbard model, the
system parameters and the nature of the calculation to be given.

System type
^^^^^^^^^^^

**k_space**
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

**electrons** *nel*
    Set the number of electrons in the system to be *nel*.
**lattice** *lattice vectors*
    Set the lattice vectors (and as a result the dimensionality) of the system.
    The lines immediately after **lattice** are assumed to be the :math:`n
    \times n` matrix containing the lattice vectors of the crystal cell (i.e.
    one lattice vector per line).  1D, 2D and 3D systems can be specified using
    vectors of the appropriate dimensionality.
**nel** *nel*
    Synonym for **electrons**
**T** *t*
    Set the kinetic term in the Hamiltonian to be *t*.
**U** *U*
    Set the Coulomb term in the Hamiltonian to be *U*.

Calculation type
^^^^^^^^^^^^^^^^

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

Calculation options
^^^^^^^^^^^^^^^^^^^

**block_size** *block_size*
    Set the block size used to distribute the Hamiltonian matrix across the
    processors.  The Hamiltonian matrix is divided into :math:`n \times n`
    sub-matrices, where :math:`n` is the block size, which are the distributed
    over the processors in a cyclic fashion.  Only relevant in parallel
    execution, where it can have a significant impact on performance.  The
    default value is 64.
**eigenvalues**
    Find only the eigenvalues of the Hamiltonian matrix from either diagonalisation
    method.  This is the default behaviour.
**eigenvectors**
    Find the eigenvectors and eigenvalues of the Hamiltonian matrix from either
    diagonalisation method.  This is much slower.  Currently the eigenvectors
    are not used.
**lanczos_basis** *nbasis*
    Set the number of Lanczos vectors to be used.  This determines the main
    memory requirements of the Lanczos routine.  The size of the basis can have
    an impact on the performance of the Lanczos diagonalisation and which
    excited eigensolutions are found.  See the TRLan documentation,
    http://crd.lbl.gov/~kewu/ps/trlan_.html, for more details.  The default is 40.
**lanczos_solutions** *nsolns*
    Set the number of eigenvalues (and eigenvectors, if required) to be found
    via Lanczos diagonlisation.  The default is 5.  The Hamiltonian matrix
    is constructed in block diagonal form using spin and crystal momentum
    conservation rules.  nsolns is the number of solutions found per block.
**lanczos_solns** *nsolns*
    Synonym for **lanczos_solutions**.

Output options
^^^^^^^^^^^^^^

**determinants** [*filename*]
    Write out the enumerated list of determinants to the given *filename*.  The
    default filename used is DETS.
**det** [*filename*]
    Synonym for **determinants**.
**hamiltonian** [*filename*]
    Write out the diagonal and the non-zero off-diagonal elements of the
    Hamiltonian matrix to the given *filename*.  The default filename used is
    HAMIL.
**hamil** [*filename*]
    Synonym for **hamiltonian**.

other options
^^^^^^^^^^^^^

**end**
    End of input.  Any subsequent lines in an input file are ignored.  It is
    only strictly required if the input is given via STDIN.

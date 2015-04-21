Full Configuration Interaction
==============================

.. code-block:: lua

    fci {
        sys = system,
        fci = { ... },
        lanczos = { ... },
        reference = { ... },
    }

.. note::

    The FCI engine in HANDE is particularly simple (i.e. slow, dumb, memory hungry) and is
    designed mainly for testing.  A conventional quantum chemistry package, such as
    MOLPRO, Q-Chem or PSI4, is highly recomended for production FCI calculations as these
    implement substantially more efficient algorithms.

Options
-------

sys
    The system on which to perform the calculation.  Must be created via a system
    function.
fci
    type: lua table.

    Optional.  No default.

    Further FCI options.  See below.
lanczos
    type: lua table.

    Optional.  No default.

    Table containing Lanczos-specific options; see below.  If present, the diagonalisation
    is performed via an iterative Lanczos algorithm.  Otherwise diagonalisation is
    performed using LAPACK or ScaLAPACK.
reference
    type: lua table.

    Optional.  No default.

    If not specified, the entire Hilbert space is used.

fci options
-----------

The ``fci`` table can take the following options:

write_hamiltonian
hamiltonian_file
write_determinants
determinant_file
write_nwfns
wfn_file
nanalyse
blacs_block_size
    type: integer.

    Optional.  Default: 64.

    Set the block size used by BLACS to distribute the Hamiltonian matrix across the
    processors with MPI parallelism.  The Hamiltonian matrix is divided into :math:`n
    \times n` sub-matrices, where :math:`n` is the block size, which are the distributed
    over the processors in a cyclic fashion.
rdm 
    type: table of integers.

    Optional.  No default.

    If present, calculate the eigenvalues for the reduced density matrix consisting of the
    specified list of sites, with a trace performed over all other sites.

    .. note::

        The ``rdm`` option is only currently available for Heisenberg systems and cannot
        be used with the Lanczos algorithm.

lanczos options
---------------

The ``lanczos`` table can take the following options:

neigv
    type: integer.

    Optional.  Default: 5.

    Number of lowest eigenstates to be found.
nbasis
    type: integer.

    Optional.  Default: 40.

    Number of Lanczos vectors used.   The size of the basis can have an impact on the
    performance of the Lanczos diagonalisation and which excited eigensolutions are found.
    See the `TRLan documentation, <http://crd.lbl.gov/~kewu/ps/trlan_.html>`_, for more
    details.
direct
    type: boolean.

    Optional.  Default: false.

    If true, generate the Hamiltonian matrix on the fly (very slow).  Otherwise generate
    the Hamiltonian once and store it for use at each Lanczos iteration.  Not implemented
    with MPI parallelism.
sparse
    type: boolean.

    Optional.  Default: true.

    If true store the Hamiltonian in a sparse matrix format.  The generation of the
    Hamiltonian matrix takes longer but requires consequently *much* less memory.  Not
    implemented with MPI parallelism.

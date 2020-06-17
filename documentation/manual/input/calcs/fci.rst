Full Configuration Interaction
==============================

Calculate the ground state of a system via a full diagonalisation of the Hamiltonian matrix [Knowles89]_.

.. code-block:: lua

    fci {
        sys = system,
        fci = { ... },
        reference = { ... },
    }

.. note::

    The FCI engine in HANDE is particularly simple (i.e. slow, dumb, memory hungry) and is
    designed mainly for testing.  A conventional quantum chemistry package, such as
    MOLPRO, or PSI4, is highly recommended for production FCI calculations as these
    implement substantially more efficient algorithms.

Options
-------

``sys``
    type: system object.

    Required.

    The system on which to perform the calculation.  Must be created via a system
    function.
``fci``
    type: lua table.

    Optional.  No default.

    Further FCI options.  See below.
``reference``
    type: lua table.

    Optional.  No default.

    If not specified, the entire Hilbert space is used.  See :ref:`reference_table`.

fci options
-----------

The ``fci`` table can take the following options:

``write_hamiltonian``
    type: boolean.

    Optional.  Default: false.

    Write out the diagonal and the non-zero off-diagonal elements of the Hamiltonian
    matrix.
``hamiltonian_file``
    type: string.

    Optional. Default: 'HAMIL'.

    Filename to which the Hamiltonian matrix is written.
``write_determinants``
    type: boolean.

    Optional.  Default: false.

    Write out the enumerated list of determinants in the FCI Hilbert space.
``determinant_file``
    type: string.

    Optional. Default: 'DETS'.

    Filename to which the list of determinants (or, more generally, many-body
    basis functions) is written.
``write_nwfns``
    type: integer.

    Optional.  Default: 0.

    Number of wavefunctions to write out (in the basis of Slater determinants).
    A negative value indicates all wavefunctions are to be written out.
``wfn_file``
    type: string.

    Optional. Default: 'FCI_WFN'.

    Filename to which the wavefunctions are written.
``nanalyse``
    type: integer.

    Optional.  Default: 0.

    Calculate properties of the first *nwfn* FCI wavefunctions from each spin and
    symmetry block.  If nwfn is negative (default) then all wavefunctions are
    analysed.  This is slow, and uses a very simple algorithm.  It is only
    designed for debugging purposes.  The properties evaluated depend upon the system
    and are liable to change without warning.
``blacs_block_size``
    type: integer.

    Optional.  Default: 64.

    The block size used by BLACS to distribute the Hamiltonian matrix across the
    processors with MPI parallelism.  The Hamiltonian matrix is divided into :math:`n
    \times n` sub-matrices, where :math:`n` is the block size, which are the distributed
    over the processors in a cyclic fashion.
``rdm``
    type: table of integers.

    Optional.  No default.

    If present, calculate the eigenvalues for the reduced density matrix consisting of the
    specified list of sites, with a trace performed over all other sites.

    .. note::

        The ``rdm`` option is only currently available for Heisenberg systems.

.. note::

    The ``write_wfn``, ``nanalyse`` and ``rdm`` options require the eigenvectors to be
    calculated in addition to the eigenvalues, which requires additional computational
    time.


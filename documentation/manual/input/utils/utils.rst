Utilities
=========

Redistribution of restart files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: lua

    redistribute {
        -- options
    }

Returns:
    nil

For speed in reading in restart files and for simplicity, HANDE produces restart files
specific to the number of MPI ranks used in the calculation and hence by default
calculations can only be restarted on the same number of MPI ranks the original
calculation ran on.  The ``redistribute`` function reads in a set of restart files and
produces a new set to be used on a different number of processors.

.. note::

   * It is convenient to place this before the QMC calculation call in the input file.
     However, the process of redistributing particles is a somewhat serial task and hence
     ``redistribute`` may not scale well to large numbers of processors.  Hence it may be
     more computationally efficient to do the redistribution targetting a large (ie 100s
     or 1000s) of processors using a much smaller number of processors in a separate run
     of HANDE.
   * Load balancing settings are reset to their default values.

HANDE uses one restart file per MPI rank with a filename of the form ``HANDE.RS.X.pY.H5``,
where ``X`` is the restart index and ``Y`` is the MPI rank.

Options:

nprocs
    type: integer.

    Optional.  Default: number of processors the calculation is running on.

    Set the number of processors that the new set of restart files are to be used on.
read
    type: integer.

    Optional.  Default: highest non-negative integer for which a set of restart files
    exists.

    Set the index, ``X`` of the set of restart files to be read in.
write
    type: integer.

    Optional.  Default: highest non-negative integer for which a set of restart files does
    not yet exist.

    Set the index, ``X`` of the set of restart files to be written out.

.. warning::

   Each processor must be able to access the entire set of existing restart files, which
   are assumed to be in the working directory.

MPI information
^^^^^^^^^^^^^^^

.. code-block:: lua

    mpi_root()

Returns:
    true if the processor is the MPI root processor and false otherwise.

The input file is processed and run by each processor.  It is occasionally useful to
perform (for example) additional I/O from lua but only on one processor.  Testing if
the procesor is the MPI root processor is a safe way to do this, e.g.

.. code-block:: lua

    if mpi_root() then
        print('root says hello from lua!')
    end

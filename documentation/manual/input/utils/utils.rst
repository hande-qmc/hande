.. _utils:


Utilities
=========

Redistribution of restart files
-------------------------------

.. code-block:: lua

    redistribute {
        -- options
    }

For speed in reading in restart files and for simplicity, HANDE produces restart files
specific to the number of MPI ranks used in the calculation and hence by default
calculations can only be restarted on the same number of MPI ranks the original
calculation ran on.  The ``redistribute`` function reads in a set of restart files and
produces a new set to be used on a different number of processors.

.. note::

   * It is convenient to place this before the QMC calculation call in the input file.
     However, the process of redistributing particles is a somewhat serial task and hence
     ``redistribute`` may not scale well to large numbers of processors.  Hence it may be
     more computationally efficient to do the redistribution targeting a large (ie 100s
     or 1000s) of processors using a much smaller number of processors in a separate run
     of HANDE.
   * Load balancing settings are reset to their default values.

HANDE uses one restart file per MPI rank with a filename of the form ``HANDE.RS.X.pY.H5``,
where ``X`` is the restart index and ``Y`` is the MPI rank.

Options:

``nprocs``
    type: integer.

    Optional.  Default: number of processors the calculation is running on.

    Set the number of processors that the new set of restart files are to be used on.
``read``
    type: integer.

    Optional.  Default: highest non-negative integer for which a set of restart files
    exists.

    Set the index, ``X`` of the set of restart files to be read in.
``write``
    type: integer.

    Optional.  Default: highest non-negative integer for which a set of restart files does
    not yet exist.

    Set the index, ``X`` of the set of restart files to be written out.
``sys``
    type: system object.

    Optional.

    Only used to determine the number of basis functions, if changing the value of DET_SIZE
    for the restart files.

.. warning::

   Each processor must be able to access the entire set of existing restart files, which
   are assumed to be in the working directory.

MPI information
---------------

.. code-block:: lua

    mpi_root()

Returns:
    true if the processor is the MPI root processor and false otherwise.

The input file is processed and run by each processor.  It is occasionally useful to
perform (for example) additional I/O from lua but only on one processor.  Testing if
the processor is the MPI root processor is a safe way to do this, e.g.

.. code-block:: lua

    if mpi_root() then
        print('root says hello from lua!')
    end

Memory management
-----------------

Objects returned from functions (e.g. system and qmc_state objects) are deallocated by
Lua's garbage collector when they are no longer required.  This can either be because the
variable goes out of scope or is set to :code:`nil`.  This level of memory management is
sufficient in most calculations.  However, there may be a substantial memory overhead when
running multiple separate calculations in the same input file as the garbage collection
need not take place immediately.  As such, objects which are no longer required can be
explicitly freed using :code:`free` methods on all objects returned by HANDE's functions.
For example, for :code:`qmc_state` objects:

.. code-block:: lua

    system = hubbard_k {
        lattice = { { 10 } },
        electrons = 6,
        ms = 0,
        sym = 1,
        U = 1,
    }

    qs1 = fciqmc {
        sys = system,
        qmc = {
            tau = 0.01,
            init_pop = 10,
            mc_cycles = 20,
            nreports = 100,
            target_population = 50000,
            state_size = 5000,
            spawned_state_size = 500,
        },
    }

    -- Deallocate all memory associated with qs1 produced by the first FCIQMC calculation.
    qs1:free()

    qs2 = fciqmc {
        sys = system,
        qmc = {
            tau = 0.02,
            init_pop = 10,
            mc_cycles = 10,
            nreports = 100,
            target_population = 50000,
            state_size = 5000,
            spawned_state_size = 500,
        },
    }

and similarly for :code:`system` objects.

.. _utils_hdf5_system_dump:

Write HDF5 system file
---------------------

.. code-block:: lua

    write_read_in_system {
        sys = system,
        filename = filename,
    }

Options:

``sys``
    type: system object.

    Required.

    The system on which to perform the calculation.  Must be created via the read_in
    function.

``filename``
    type: string. Optional.

    Filename to dump system hdf5 file to. If unset will generate a filename to dump to
    based on the template: int_file + CAS_information + .H5, where ``int_file`` and the
    CAS information are set in the call to ``read_in`` which create the ``system`` object.

Returns:
    name of HDF5 file created.  This is currently only available on the root processor and
    can be passed into subsequent calls to ``read_in`` safely as only the root processor
    reads from integral and system files.

When running a calculation using a system generated from a FCIDUMP, the `:`system`` object
created by lua_read_in can be dumped in HDF5 format for reuse in subsequent calculations;
this speeds initialisation by a factor of ~100x and reduces the required file size by ~16x
for large FCIDUMPs.  When running in parallel on a large number of cores this is
particularly important to utilise as it overcomes an inherent serialisation point in the
calculation initialisation.

For example:

.. code-block:: lua

     sys = read_in {
         int_file = "FCIDUMP",
         nel = 24,
         ms = 0,
         sym = 0,
     }

     hdf5_name = write_read_in_system {
         sys = sys,
     }

produces an HDF5 file entitled "FCIDUMP.H5" and return this value to the variable
``hdf5_name``.  Passing this as the argument to ``int_file`` within ``read_in`` will use
it in future calculations -- the HDF5 format of the file is automatically detected.

If a CAS is used to produce the system object used to produce such a file it will be
labelled as such and only information for basis functions within the CAS will be stored;
conversion between different CAS within this functionality is not currently supported.

.. important::

    When using a HDF5 file to initialise a system either both of nel and ms must be
    specified or neither; if neither are specified the values stored within the system
    HDF5 file will be used and otherwise the given values override those stored.

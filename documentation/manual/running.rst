Usage
=====

.. code-block:: bash

    $ hande.x [input_filename]

Output is sent to STDOUT and can be redirected as desired.

Parallel Usage
--------------

Using MPI only:

.. code-block:: bash

    $ mpirun -np n hande.x [input_filename]

where ``n`` is the number of processors to run on in parallel. On an HPC system this may
differ (for example mpirun -np n may be replaced with mpiexec), depending on how the
environment has been set up.

Using OpenMP parallelism:

.. code-block:: bash

    $ export OMP_NUM_THREADS=n
    $ hande.x [input_filename]

OpenMP parallelism is currently only implemented for CCMC.

Using OpenMP and MPI parallelism:

.. code-block:: bash

    $ export OMP_NUM_THREADS=n
    $ mpirun -np m hande.x [input_filename]

where ``m`` is the number of MPI processes and ``n`` is the number of OpenMP threads per
MPI process.  HANDE prints this information at the top of the output, so one can easily
check there environment is set up correctly.

HANDE only performs I/O operations on the root processor when run on multiple processors.

.. note::

    Due to the implementation of efficient Monte Carlo algorithms, running the Monte Carlo
    algorithms in HANDE on different numbers of processors (or using OpenMP) results in
    different Markov chains and hence such calculations will not agree exactly but instead
    statistically.

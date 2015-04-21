Usage
=====

.. code-block:: bash

    hande.x [input_filename]

If no input filename is provided then the input options are read from STDIN.
Note that this feature is not guaranteed to work when run in parallel!

Output is sent to STDOUT and can be redirected as desired.

Parallel Usage
--------------

Using MPI only:

.. code-block:: bash

    mpirun -np n hande.x [input_filename]

Where n is the number of processors to run on in parallel. On a HPC this may
differ (for example mpirun -np n may be replaced with mpiexec), depending on 
how the environment has been set up.

Using OpenMP parallelism: 

.. code-block:: bash

    export OMP_NUM_THREADS=n
    hande.x [input_filename]

OpenMP parallelism is currently only implemented for CCMC.

Using OpenMP parallelism: 

.. code-block:: bash

    export OMP_NUM_THREADS=n
    mpirun -np m hande.x [input_filename]

Where m is the number of mpi processes and n is the number of OpenMP threads
per mpi process.
Hande prints this information at the top of the output file, so one can easily
check there environment is set up correctly:

| Number of MPI processes running on: m
| Running with n threads per MPI process.

hande.x only performs i/o operations on the root processor when run on
multiple processors.

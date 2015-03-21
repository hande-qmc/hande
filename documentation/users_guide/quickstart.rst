Quick Start Guide
=================

Documentation
-------------
You've found the quick start guide it seems.  Other documentation is found in 
documentation/users_guide/manual.rst and more detailed developer documentation at
documentation/developers_guide/developers_guide.rst .

Basic compilation
-----------------
For a machine which has gfortran, but no special extra libraries, there is a very simple
configuration to get you started:

.. code-block:: bash

    u@h:~/code$ git clone hande@tycpc15.cmth.ph.ic.ac.uk:hande.git
    u@h:~/code$ cd hande 
    u@h:~/code/hande$ tools/mkconfig.py gfortran.basic
    u@h:~/code/hande$ make -j2

The -j2 runs the compilation on two threads, but you can use as many as are available.
This simple configuration links to the following: libstdc++ libz libuuid liblapack libblas
which should all be available on systems configured for scientific computing, but others
are available which enable OpenMP, MPI, HDF5 and other features.

You should checkout and run testcode to make sure that nothing untoward has happened
in the compilation:
.. code-block:: bash

    u@h:~/code/hande$ cd ..
    u@h:~/code$ git clone https://github.com/jsspencer/testcode.git 
    u@h:~/code$ cd hande/test_suite
    u@h:~/code/hande/test_suite$ ../../testcode/bin/testcode.py --total-processors=2

You can modify the --total-processors argument to the number of available cores on which
to run concurrent tests for faster runs.  For reasonably fast processors this will take
about 10 core-minutes, and will result in a printout of something like
.. code-block:: bash

67 out of 67 tests passed (60 not checked).

If any tests have failed, consult the user guide for some hints as to what you might need
to do do change configuration.

A simple run
------------
We'll do a simple FCIQMC calculation to show how to run and analyse data.
We'll base this on one of the test suite examples.

Pick a temporary directory to do some runs in.  I'll call this wdir, and the directory HANDE
is in as $hdir.  We'll copy the input file first.

.. code-block:: bash

    u@h:wdir$ cp $hdir/test_suite/fciqmc/np1/H2O-RHF-cc-pVTZ/h2o.in .

Looking at h2o.in, there a few types of line 
.. code-block:: bash

    varyshift_target 10000000

    (    SCF calculation produced by Q-Chem using:  )

    #    $molecule

The first is a standard line which contains a keyword followed by a value.
Keywords can take numbers, text, and combinations of both.  They are all documented in the
user guide.
The final two lines are two styles of comment.

More on this when the lua input file is done.

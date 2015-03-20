Quick Start Guide
=================

Basic compilation
-----------------
For a machine which has gfortran, but no extra libraries, there is a very simple
configuration to get you started:

.. code-block:: bash

    u@h:~/code$ git clone hande@tycpc15.cmth.ph.ic.ac.uk:hande.git
    u@h:~/code$ cd hande 
    u@h:~/code/hande$ tools/mkconfig.py gfortran.basic

This links to the following: libstdc++ libz libuuid liblapack libblas
which should all be available on systems configured for scientific computing.

You should checkout and run testcode to make sure that nothing untoward has happened:
.. code-block:: bash

    u@h:~/code/hande$ cd ..
    u@h:~/code$ git clone https://github.com/jsspencer/testcode.git 
    u@h:~/code$ cd hande/test_suite
    u@h:~/code/hande/test_suite$ ../../testcode/bin/testcode.py --total-processors=2

You can modify the --total-processors argument to the number of available cores on which
to run concurrent tests for faster runs.  For reasonably fast processors this will take about 10 core-minutes, and will result in a printout of something like
.. code-block:: bash

67 out of 67 tests passed (60 not checked).

If any tests have failed, consult the user guide for some hints as to what you might need
to do do change configuration.
 

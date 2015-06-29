Test suite
==========

HANDE has an extensive test suite covering all core functionality.
The tests are run using the ``testcode`` package (https://github.com/jsspencer/testcode).
Note that the data extraction scripts for HANDE require the ``pandas`` python library.

testcode can be run from the test_suite subdirectory:

.. code-block:: bash

    testcode.py

Note that the default set of tests are serial only.  The entire test suite is
run every night using buildbot (http://www.cmth.ph.ic.ac.uk/buildbot/hande/).

Selected data from the HANDE output is compared to known 'good' results
('benchmarks').

testcode is quite flexible and it's easy to run subsets of tests, check against
different benchmarks, compare previously run tests, run tests concurrently for
speed, etc.  Please see the testcode documentation for more details.

.. note::

    For algorithmic reasons, certain compilation and runtime options (principally
    POP_SIZE and processor/thread count) result in different Markov chains
    and hence different exact results (but same results on average).  The tests
    should therefore be run using the same compilation options and the same
    parallel distribution as was used for the benchmarks.  The latter for MPI
    parallelisation is done automatically by testcode.  Separate tests exist
    for both POP_SIZE=32 and POP_SIZE=64.

    Currently there are no QMC tests suitable for OpenMP parallelisation due to
    difficulties with making the scheduler behave deterministically without
    affecting performance of production simulations.
    It is advised that you make sure to set the shell varialble OMP_NUM_THREADS
    to 1 when running the test suite - otherwise these will all be marked SKIPPED.

What if the tests fail?
-----------------------

A common cause for tests failing is that the configuration causes a different Markov
Chain to be run, or part of the code has been disabled in your build.
testcode should determine that some tests are inappropriate and skip them.
To force testcode to skip some tests, see below.

A second cause of failure is that some floating point values have rounded differently on
different architectures.
The tolerances used for the tests can also be adjusted as specified below:

Skipping Tests
--------------

If there is a unique line printed out in the output for jobs which are to be skipped, 
this can be used to tell testcode this, by modifying the skip_args line in the 
test_suite/userconfig file.  See the testcode documentation for more details

Adjusting Test Tolerances
-------------------------

The tolerance for an individual job can be modified as specified in the testcode documentation.
As an example, to modify the tolerance because of the following failure:

::

    dmqmc/np1/heisenberg_1d - replica.in: **FAILED**.
    \sum\rho_{ij}M2{ji}
        ERROR: absolute error 1.00e-06 greater than 1.00e-10. (Test: 17.378583.  Benchmark: 17.378584.)

The follow section can be inserted into test_suite/jobconfig.  Note the backslash-quoting of the 
backslashes, as the tolerance value is interpreted as a python tuple containing a python string.

::

    #Job specific tolerances:                                                                 
    [dmqmc/np1/heisenberg_1d/]                                                                
    tolerance = (1e-5,1e-5,'\\sum\\rho_{ij}M2{ji}')          

Test suite
==========

The tests are run using the ``testcode`` package (https://github.com/jsspencer/testcode).
The data extraction scripts for HANDE require the ``pandas`` python library (http://pandas.pydata.org/),
which could for example be installed by

.. code-block:: bash
    
    pip install pandas

or alternatively where pip is not available, one can install it locally:

.. code-block:: bash

    wget https://github.com/pydata/pandas/archive/v0.15.2.tar.gz
    tar -xzvf v0.15.2.tar.gz
    cd pandas-0.15.2
    python setup.py build
    python setup.py install
    
If you do not have root access, you can install the library locally with:

.. code-block:: bash

    python setup.py install --user
    echo 'export PYTHONPATH=$HOME/.local/lib/python2.7/site-packages:$PYTHONPATH' >>.bashrc



testcode can be run from the test_suite subdirectory:

.. code-block:: bash

    testcode.py

Note that the default set of tests are serial only.  The entire test suite is
run every night using buildbot (http://www.cmth.ph.ic.ac.uk/buildbot/hande/).

Selected data from the HANDE output is compared to known 'good' results
('benchmarks').  The python script which extracts this data uses the pandas
module and, unfortunately, importing pandas is actually the time-consuming step
in the data analysis.  To help alleviate this, the data extraction script, can
be run in a server-client mode.  The server can be launched using:

.. code-block:: bash

    tools/tests/extract_test_data.py --socket &

If a server (on the default port) is running, the data extraction script used
by testcode will automatically use it, greatly speeding up the data analysis
step.

testcode is quite flexible and it's easy to run subsets of tests, check against
different benchmarks, compare previously run tests, run tests concurrently for
speed, etc.  Please see the testcode documentation for more details.

.. note::

    For algorithmic reasons, certain compilation options (principally POP_SIZE
    and DET_SIZE and processor/thread count) result in different Markov chains
    and hence different exact results (but same results on average).  The tests
    should therefore be run using the same compilatition options and the same
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
.. code-block::
    dmqmc/np1/heisenberg_1d - replica.in: **FAILED**.
    \sum\rho_{ij}M2{ji}
        ERROR: absolute error 1.00e-06 greater than 1.00e-10. (Test: 17.378583.  Benchmark: 17.378584.)

The follow section can be inserted into test_suite/jobconfig.  Note the backslash-quoting of the 
backslashes, as the tolerance value is interpreted as a python tuple containing a python string.
.. code-block::
    #Job specific tolerances:                                                                 
    [dmqmc/np1/heisenberg_1d/]                                                                
    tolerance = (1e-5,1e-5,'\\sum\\rho_{ij}M2{ji}')          



.. _analysis:

Analysis
=====================

The following provides a brief overview for the most common analysis required for each
type of Monte Carlo calculation.  The guides in :ref:`tutorials` provide a step-by-step
guide to analysing HANDE calculations and explain the reasoning behind the required
analysis parameters.

HANDE includes a variety of scripts and utilities in the ``tools`` subdirectory.  However,
these only provide a simple, command-line interface.  A comprehensive python module,
:ref:`pyhande`, drives all the analysis.  :ref:`pyhande` is extremely powerful for dealing
with complex analysis, data-driven investigation or bulk data analysis.

FCIQMC and CCMC
---------------------
Fundamental Usage
^^^^^^^^^^^^^^^^^^^^^

QMC calculations print out data from a block of iterations (a 'report loop'), the length
of which is controlled by the **mc_cycles** input option.  Care should be taken analysing
this data and, in particular, producing accurate estimates of the errors in the means of
the energy estimators.  Almost all data is averaged over the report loop (see output for
further details).

Note that no data is lost when quantities are summed over report loops, as the
correlation length in the data is substantially longer than the length of the
report loop (typically 10-20 iterations).

As the particle distribution at one iteration is not independent from the distribution at
the previous iteration, estimators at each iteration are not independent.  This
correlation in the data needs to be taken into account when estimating standard errors.
A simple and effective way of doing this is to use a blocking analysis
[Flyvbjerg89]_.

The ``reblock_hande.py`` script (in the ``tools`` subdirectory) does this.  Run

.. code-block:: bash

    $ reblock_hande.py --help

to see the available options.  Estimates for the shift and projected energy are
typically obtained using

.. code-block:: bash

    $ reblock_hande.py --start N out

respectively, where ``N`` is the iteration from which data should be blocked (i.e.
after the calculation has equilibrated) and ``out`` is the file to which the
calculation output was saved. Without --start option, this script automatically 
estimates the appropriate ``N``, so you don't have to its value basically.

Note that reblock_hande.py can accept multiple output files for the case when
a calculation is restarted as follows:

.. code-block:: bash

    $ reblock_hande.py -m out1 out2 out3

More complicated analysis can be performed in python by
using the ``pyhande`` library --- ``reblock_hande.py`` simply provides a convenient
interface for the most common analysis tasks.

Hybrid method
^^^^^^^^^^^^^^^^^^^^^

Hybrid method is another choice to estimate 
errors, which is available as 

.. code-block:: bash

    $ reblock_hande.py -a hybrid out

It has been shown by our experiment that
hybrid method gives more reliable estimation
of errors than blocking method: In the experiment,
1000 different CCMC-SD energy time-series were 
obtained for Nitrogen atom, using the same calculation 
setting but with different random seeds.
The mean and the errors were obtained by hybrid method
and blocking method, and it is examined how many means
coincide with the CCSD energy within the calculated
standard errors with respect to both methods.
The expected coincidence rate for standard errors 
is 68.27%. Thus, when the actual coincidence rate 
is closer to this value, the using post-analysis 
method is more reliable.

Two types of coincidence rate were used to evaluate 
reliabilities, conditional coincidence rate (CCR) 
and unconditional coincidence rate (UCR).
They are given as 

    CCR = Hit / ( Total - Failed ) * 100 改行 
    UCR = Hit /       Total        * 100

Here, 'Total' is the total number of post-analyses (=1000),
'Failed' is the number of post-analyses which fails to
make an estimation of the error(*), and 
'Hit' is the number of post-analyses which make an estimation
of the error and the estimated error coincided with 
the CCMC-SD energy within the standard error.
(*: e.g. 'Shift is not started yet' in the case of blocking method)

The below figures compare the CCRs and UCRs obtained for 
different lengths of time-series using hybrid method and
blocking method, respectively. Both figures shows that
hybrid method is more reliable than blocking method
especially for long lengths of time-series. 

.. plot::

    import pandas
    import matplotlib.pyplot as plt
    
    filename='calcs/hybrid.result'
    res_hyb = pandas.read_csv(filename, delim_whitespace=True)
    plt.plot(res_hyb['data_size'], res_hyb['ccr'])
    
    filename='calcs/blocking.result'
    res_blk = pandas.read_csv(filename, delim_whitespace=True)
    plt.plot(res_blk['data_size'], res_blk['ccr'])

    plt.axhline(y=68.27, linestyle='-')


.. plot::

    import pandas
    import matplotlib.pyplot as plt
    
    filename='calcs/hybrid.result'
    res_hyb = pandas.read_csv(filename, delim_whitespace=True)
    plt.plot(res_hyb['data_size'], res_hyb['ucr'])
    
    filename='calcs/blocking.result'
    res_blk = pandas.read_csv(filename, delim_whitespace=True)
    plt.plot(res_blk['data_size'], res_blk['ucr'])

    plt.axhline(y=68.27, linestyle='-')



MSER minimization
^^^^^^^^^^^^^^^^^^^^^

MSER minimization method is another choice 
to estimate starting iterations. 
One can use the method as

.. code-block:: bash

    $ reblock_hande.py -b mser_min out

Applying MSER min. and WREE min. to 1000
different CCMC-SD energy time-series,
MSER min. gave almost the same starting
iteration for different length of time-series 
and, meanwhile, the starting iteration predicted 
by WREE min. increased according to the length of time-series.


Canonical Total Energy MC
---------------------------

The configurations and resulting estimates in a canonical total energy
calculation are statistically independent and therefore no blocking analysis is
required. The ``analyse_canonical.py`` script is available in ``tools/canonical_energy/`` which
performs the appropriate averaging and standard error analysis on the output file
using the pyhande suite.

DMQMC
-----

No blocking analysis is required for the error analysis of DMQMC calculations
as estimates are averaged over statistically independent runs. The
``finite_temp_analysis.py`` script in ``tools/dmqmc`` can be used to perform a
standard error analysis of the Monte Carlo data for a number of different observables.

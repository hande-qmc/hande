.. _analysis:

Analysis
========

The following provides a brief overview for the most common analysis required for each
type of Monte Carlo calculation.  The guides in :ref:`tutorials` provide a step-by-step
guide to analysing HANDE calculations and explain the reasoning behind the required
analysis parameters.

HANDE includes a variety of scripts and utilities in the ``tools`` subdirectory.  However,
these only provide a simple, command-line interface.  A comprehensive python module,
:ref:`pyhande`, drives all the analysis.  :ref:`pyhande` is extremely powerful for dealing
with complex analysis, data-driven investigation or bulk data analysis.

FCIQMC and CCMC
---------------

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
calculation output was saved.

Note that reblock_hande.py can accept multiple output files for the case when
a calculation is restarted.  More complicated analysis can be performed in python by
using the ``pyhande`` library --- ``reblock_hande.py`` simply provides a convenient
interface for the most common analysis tasks.

Canonical Total Energy MC
---------------------------

The configurations and resulting estimates in a canonical total energy
calculation are statistically independent and therefore no blocking analysis is
required. The ``analyse_canonical.py`` script is available in ``tools/dmqmc/`` which
performs the appropriate averaging and standard error analysis on the output file
using the pyhande suite.

DMQMC
-----

No blocking analysis is required for the error analysis of DMQMC calculations
as estimates are averaged over statistically independent runs. The
``finite_temp_analysis.py`` script in ``tools/dmqmc`` can be used to perform a
standard error analysis of the Monte Carlo data for a number of different observables.

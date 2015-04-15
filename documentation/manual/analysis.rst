Analysis
========

FCIQMC and iFCIQMC calculations print out data from a block of iterations (a
'report loop'), the length of which is controlled by the **mc_cycles** input
option.  Care should be taken analysing this data and, in particular, producing
accurate estimates of the errors in the means of the energy estimators.

Users are encouraged to read the notes in
documentation/theory/projected_energy/proje.tex.  As the psip distribution at
one iteration is not independent from the distribution at the previous
iteration, the energy at each iteration is not independent.  This correlation in
the data needs to be taken into account when estimating standard errors.
A simple and effective way of doing this is to use a blocking analysis
[FlyvbjergPetersen89]_.

Each report loop prints out the following data:

iterations
    The number of completed iterations.
Instant shift
    The value of the shift (growth estimator, in DMC language) based upon the
    current psip distribution.
Av. shift
    The running average of the shift.  This is accumulated from the first
    iteration that the shift is allowed to vary within the current calculation
    (i.e. it is not preserved when a calculation is restarted).  As such, it
    does not exclude an equilibration period and is not always a good estimate
    of the true mean as a result.  
\sum H_0j Nj
    The numerator of the projected energy summed over the iterations in the
    report loop; the sum over the determinants connected to the reference
    determinant multiplied by the psip population on the determinant, in term
    summed over iterations in the report loop.
Av. Proj. E
    The running average of the projected energy.  This is accumulated from the 
    start of the current calculation (i.e. it is not preserved when
    a calculation is restarted).  As such, it does not exclude an
    equilibration period and is not always a good estimate of the true mean as
    a result.

    Note that the numerator and denominator are accumulated separately and the
    ratio printed out to avoid a bias caused by the ratio of means being
    different from the mean of a ratio.
# D0
    The denominator of the projected energy summed over the iterations in the
    report loop; the psip population on the reference determinant summed over
    the iterations in the report loop.
# particles
    The total psip population at the end of the report loop.
R_spawn
    The average rate of spawning for each iteration in the report loop; the
    fraction of spawning attempts which were successful.
time
    The average time each iteration took between report loops.

Note that no data is lost when quantities are summed over report loops, as the
correlation length in the data is substantially longer than the length of the
report loop (typically 20 iterations).

The running averages of the shift and projected energy can be reset using the
**zero_means** option with :ref:`HANDE.COMM`.  However, it should be
emphasised that the best estimates of the energy and associated standard error
are obtained via re-blocking the data as a post-processing step.  Often the
averaged values printed out are only adequate for (at best) monitoring
convergence and stability.  The reblock_hande.py script (in the tools subdirectory)
does this.  Run

.. code-block:: bash

    reblock_hande.py --help

to see the available options.  Estimates for the shift and projected energy are
typically obtained using

.. code-block:: bash

    reblock_hande.py --start N out

respectively, where N is the iteration from which data should be blocked (i.e.
after the calculation has equilibrated) and out is the file to which the
calculation output was saved.

Note that reblock_hande.py can accept multiple output files for the case when
a calculation is restarted.  More complicated analysis can be performed in python by
using the pyhande library (which reblock_hande is a simple wrapper around).

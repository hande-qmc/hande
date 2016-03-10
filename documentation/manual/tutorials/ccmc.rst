.. _ccmc_tutorial:

Coupled Cluster Monte Carlo
===========================

In this tutorial we will run CCMC on the carbon monoxide molecule in a cc-pVDZ
basis.  For details of the theory see [Thom10]_ and [Spencer15]_.

This tutorial only presents the basic options available in a CCMC calculation;
for the full range of options see the main :ref:`documentation <ccmc>`.

To perform calculations on a molecular system in HANDE, we need the one- and
two- electron :ref:`integrals <generating_integrals>` in some appropriate basis
from an external source.  For the calculations in this tutorial, the integrals
were calculated using Psi4; input and output files can be found with the files
from the calculations herein in the ``documentation/manual/tutorials/calcs/ccmc``
subdirectory.

The system definition is exactly the same as for FCIQMC:

.. code-block:: lua

    sys = read_in {
        int_file = "CO.CCPVDZ.FCIDUMP",
        nel = 14,
        ms = 0,
    }

Note that we have not specified an overall symmetry.  In this case HANDE uses
the Aufbau principle to select a reference determinant.

A CCMC calculation can be run in a very similar way to FCIQMC.  As for FCIQMC we
can substantially reduce stochastic error by using real amplitudes, which we do
for all calculations presented here.  The most significant difference from an
FCIQMC input is that it is standard to use truncation with CCMC, specified by the
``ex_level`` option, (i.e. 2 for CCSD, 3 for CCSDT, etc.).  The determination of
a plateau and hence a suitable value for ``target_population`` is
exactly analogous to :ref:`FCIQMC <fciqmc_tutorial>`, as the sign problem is
similar between the two methods; we will not discuss it
further here.  The CCSDTMC calculation can be run using an input file such as:

.. literalinclude:: calcs/ccmc/co_ccmc.lua
    :language: lua

Note the much larger initial population compared to an FCIQMC calculation; if
this is too low the correct wavefunction will not be obtained.

Looking at the :download:`output <calcs/ccmc/co_ccmc.out>`, we see the evolution of
the population has a similar form to FCIQMC:

.. plot::

    import pyhande
    import matplotlib.pyplot as plt
    (metadata, qmc_data) = pyhande.extract.extract_data('calcs/ccmc/co_ccmc.out')[0]
    plt.plot(qmc_data['iterations'], qmc_data['# H psips'])
    plt.xlabel('iteration')
    plt.ylabel('# particles')

and the shift and projected energy vary about the correlation energy:

.. plot::

    import pyhande
    import matplotlib.pyplot as plt
    (metadata, qmc_data) = pyhande.extract.extract_data('calcs/ccmc/co_ccmc.out')[0]
    plt.plot(qmc_data['iterations'], qmc_data['Shift'], label=r'$S(\tau)$')
    plt.plot(qmc_data['iterations'], qmc_data['\sum H_0j N_j']/qmc_data['N_0'], label=r'$E(\tau) = \sum_j H_{0j} N_j(\tau)/N_0(\tau)$')
    plt.legend()
    plt.xlabel('iteration')
    plt.ylabel('Correlation energy / $t$')

The output of the calculation can be analysed in exactly the same way as for
FCIQMC:

.. code-block:: bash

    $ reblock_hande.py --quiet --start 100000 co_ccmc.out

giving

.. literalinclude:: calcs/ccmc/co_ccmc.block

Due to the sampling of the wavefunction in CCMC, it is more prone to "blooming"
events where many particles are created in a single spawning event than is
FCIQMC.  Details of blooming during a calculation are reported at the end of the
output.  It can be seen that significant blooming occurred.  This substantially
increases the stochastic error, and in particularly severe cases can cause the
calculation to not give a correct result due to the instability.  These
events can be avoided by reducing the timestep, but the timestep
required to eliminate them entirely is often prohibitively small.  Another way
of reducing them is the use of the ``cluster_multispawn_threshold`` keyword,
whereby large spawning attempts are divided into a number of smaller spawns:

.. literalinclude:: calcs/ccmc/co_ccmc_multispawn.lua
    :language: lua

Running as before, and inspecting the :download:`output <calcs/ccmc/co_ccmc_multispawn.out>`,
it can be seen that, despite using a much larger timestep, there are now no blooms.  This
substantially reduces the stochastic error:

.. literalinclude:: calcs/ccmc/co_ccmc_multispawn.block

The extra spawning causes the calculation to run more slowly, but the reduction in
error bars can often more than make up for this.

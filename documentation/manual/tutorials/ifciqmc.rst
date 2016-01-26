.. _ifciqmc_tutorial:

Initiator Approximation to FCIQMC
=================================

We shall again calculate the ground state energy of the 18-site 2D Hubbard model at
half-filling and with :math:`U/t=1.3`, as in the :ref:`FCIQMC tutorial <fciqmc_tutorial>`.
The initiator approximation [Cleland10]_ greatly reduces the number of particles required
to sample the wavefunction.  The drawback, however, is that the approximation must be
carefully controlled to obtain an accurate estimate of the FCI energy by running multiple
calculations with increasing populations.

It is efficient (both computationally and in terms of elapsed time) to treat each
calculation separately.  For compactness, we shall simply run multiple calculations with
different ``target_population`` values one after the other in the same HANDE calculation.
This is trivial to do by using a lua loop as ``fciqmc`` is simply a function call:

.. literalinclude:: calcs/ifciqmc/hubbard_ifciqmc.lua

The only difference between the above input and an FCIQMC calculation is the setting
``initiator = true``.  As before, this can be run using:

.. code-block:: bash

    $ mpiexec hande.x hubbard_ifciqmc.lua >  hubbard_ifciqmc.out

Again, the exact command to launch MPI will vary with MPI implementation and local
configurations.

Inspecting the :download:`output <calcs/ifciqmc/hubbard_ifciqmc.out>`, we see one iFCIQMC
calculation was run for each call to the :ref:`fciqmc <fciqmc>` function.  :ref:`pyhande`
(and, by extension, ``reblock_hande.py``) can handle such cases, so we easily extract and
inspect the data for each calculation.

Let's start by inspecting instantaneous projected energy estimator for the three smallest
populations:

.. plot::

    import pyhande
    import matplotlib.pyplot as plt
    ifciqmc_data = pyhande.extract.extract_data('calcs/ifciqmc/hubbard_ifciqmc.out')
    for (md, data) in ifciqmc_data[:3]:
        plt.plot(data['iterations'], data['\sum H_0j N_j']/data['N_0'], label=r'$N_{\mathrm{target}} = %s$' % int(md['qmc']['target_particles']))
    plt.ylim(-0.40, -0.25)
    plt.xlabel('iteration')
    plt.ylabel(r'$E(\tau) = \sum_j H_{0j} N_j(\tau)/N_0(\tau)$ / $t$')
    plt.legend()
    plt.tight_layout()

Whilst the difference is small on this scale, it is evident that the calculation with the
smallest population has a slightly higher mean than calculations with larger populations.
To confirm this, we will plot the energy as a function of population.  As
``target_population`` is the population at which the population **starts** to be
controlled, we should consider the average population (which is somewhat higher).
We can also compare directly to the FCIQMC energy in this case, as the population required
for the FCIQMC calculation is sufficiently small:

.. plot::

    import pyhande
    import matplotlib.pyplot as plt
    import pandas as pd
    # Start averaging after the population has stabilised in all calculations.
    calcs = pyhande.lazy.std_analysis(['calcs/ifciqmc/hubbard_ifciqmc.out'], 30000, extract_psips=True)
    opt_blocks = [calc.opt_block for calc in calcs]
    estimates = pd.concat([opt.stack() for opt in opt_blocks], keys=range(len(opt_blocks)), axis=1).T
    plt.errorbar(estimates[('# H psips', 'mean')],estimates[('Proj. Energy', 'mean')], 
                 xerr=estimates[('# H psips', 'standard error')],
                 yerr=estimates[('Proj. Energy', 'standard error')], label='iFCIQMC')

    fciqmc_calc = pyhande.lazy.std_analysis(['calcs/fciqmc/hubbard_fciqmc.out'], 30000)[0]
    fciqmc_energy = fciqmc_calc.opt_block['mean']['Proj. Energy']
    fciqmc_err = fciqmc_calc.opt_block['standard error']['Proj. Energy']
    plt.axhline(fciqmc_energy, label='FCIQMC', color='b')
    plt.axhspan(fciqmc_energy+fciqmc_err, fciqmc_energy-fciqmc_err, color='b', alpha=0.25)

    plt.xscale('log')
    plt.xlabel('# particles')
    plt.ylabel('Projected energy / $t$')
    plt.legend()


The light blue region indicates the extent of the FCIQMC stochastic error, as
calculated in the :ref:`FCIQMC tutorial <fciqmc_tutorial>`.  In this case, the initiator
approximation reduces the population required by a factor of :math:`\sim 2`.  However,
many studies (including on the electron gas and molecular systems) have demonstrated the
initiator approximation can reduce the population required by many orders of magnitude.

The estimates for each calculation can be found directly by using ``reblock_hande.py``:

.. code-block:: bash

    $ reblock_hande.py --quiet --start 30000 hubbard_ifciqmc.out

where again we chose the start point from inspecting the population growth.  This gives:

.. literalinclude:: calcs/ifciqmc/hubbard_ifciqmc.block

``reblock_hande.py`` can also handle the case where each calculation is run separately and
each separate file is passed in as a separate argument on the command line.

.. note::

    We highly recommend a visual inspection of the convergence plot of the initiator error
    as the convergence can be non-monotonic and, as a result, requiring more than two
    populations give statistically equivalent results.

Finally, using real populations can, as with :ref:`fciqmc <fciqmc_tutorial>`, have
a significant impact on the stochastic error.  Again, this is done by setting
``real_amplitudes = true`` in the input file (see
:download:`hubbard_ifciqmc_real.lua <calcs/ifciqmc/hubbard_ifciqmc_real.lua>`) and running:

.. code-block:: bash

    $ mpiexec hande.x hubbard_ifciqmc_real.lua >  hubbard_ifciqmc_real.out

followed by the blocking analysis on the :download:`output <calcs/ifciqmc/hubbard_ifciqmc_real.out>`:

.. code-block:: bash

    $ reblock_hande.py --quiet --start 30000 hubbard_ifciqmc_real.out

results in

.. literalinclude:: calcs/ifciqmc/hubbard_ifciqmc_real.block

.. [review] - AJWT: Plotting both the reals and non-reals on the same graph is quite instructive.
.. [review] - AJWT: As such you need some more data up to 10^6 to get to a convergerd number by the look of it.

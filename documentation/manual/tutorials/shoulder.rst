.. _shoulder_tutorial:

Shoulder Plots 
==============

This tutorial looks further into finding the optimal target particle population
in more detail. It is advisable to have read the :ref:`FCIQMC <fciqmc_tutorial>` and :ref:`CCMC <ccmc_tutorial>` tutorials
before this one. More information and details on shoulder plots can be found in [Spencer15]_.

The example used here is a CCSDT calculation on water in a cc-pVDZ basis [Dunning89]_.
As for the :ref:`CCMC tutorial <ccmc_tutorial>`, the
integrals were calculated with PSI4 (see :ref:`generating_integrals` for details).
Input and output files are in ``documentation/manual/tutorials/calcs/shoulder/``.

The first calculation was run using

.. literalinclude:: calcs/shoulder/h2o_plat.lua
	:language: lua

As in :ref:`FCIQMC <fciqmc_tutorial>`, a plateau can be seen in the total population vs iteration
plot, which indicates roughly the minimum particle number to make the calculation
stable.

.. plot::

    import pyhande
    import matplotlib.pyplot as plt
    (metadata, qmc_data) = pyhande.extract.extract_data('calcs/shoulder/h2o_plat.out')[0]
    fig, axisY1 = plt.subplots()
    axisY1.set_yscale('log')
    axisY1.plot(qmc_data['iterations'], qmc_data['# H psips'], 'r-')
    axisY1.set_xlabel('iteration')
    axisY1.set_ylabel('# particles', color='r')
    axisY2 = axisY1.twinx()
    axisY2.set_yscale('log')
    axisY2.plot(qmc_data['iterations'], qmc_data['N_0'], 'b-')
    axisY2.set_ylabel('# particles on reference determinant', color='b')
    for tickY1 in axisY1.get_yticklabels() :
	tickY1.set_color('r')
    for tickY2 in axisY2.get_yticklabels() :
	tickY2.set_color('b')
    plt.show()

The plateau is clearly visible at around 20000 particles.  This is one technique but the
plateau is frequently not so easy to observe by visual inspection, especially for CCMC.

.. [review] - JSS:

    Alex has some thoughts about terminology (determinant vs excitor, reference
    determinant/excitor/reference) etc.  It would be good to be consistent.

    Isn't the key observation that prior to the plateau the total population grows faster
    than the reference, after the plateau the populations grow at the same rate (because
    we have a converged wavefunction) and during the population on the reference increases
    whilst the total remains approximately constant?  Hence the plateau corresponds to
    a peak in the N_tot/N_0 vs N_tot.  This isn't clear from your text, which implies that
    the total population grows faster than the reference after the plateau.

In the beginning, the only occupied determinant is the reference determinant. These
particles then spawn onto other available determinants, making the total population grow
at a greater pace than the population on the reference determinant does. At the plateau
point, death and annihilation cancel out the spawning until the total population grows
again. However, the new rate of growing is less than the growth rate of particles on the
reference determinant. See [Spencer15]_ for details. This implies that it is informative
to consider the ratio of total population to population on the reference determinant which
we do in the form of "shoulder" plots as shown below:

.. plot::

    import pyhande
    import matplotlib.pyplot as plt
    (metadata, qmc_data) = pyhande.extract.extract_data('calcs/shoulder/h2o_plat.out')[0]
    plt.loglog(qmc_data['# H psips'], (qmc_data['# H psips'] / qmc_data['N_0']))
    plt.xlabel('# particles')
    plt.ylabel('# particles / # particles on reference determinant')

The position of the shoulder is at about 20000 which corresponds to the position
of the plateau. 

.. note::
   
    :ref:`pyhande` contains two functions to estimate the position of the
    plateau/shoulder: :func:`pyhande.analysis.plateau_estimator`, which looks
    for the peak in the shoulder plot [Spencer15]_, and :func:`pyhande.analysis.plateau_estimator_hist`,
    which uses a histogram approach to identifying the plateau [Shepherd14]_.  As a result
    of the difference in approaches, the former tends to pick up the population at the
    start of the plateau whilst the latter favours the end of the plateau and is less well
    suited to cases without a clear plateau.

    In this case, ``plateau_estimator`` gave 18480 with an estimated standard error of 40
    for the shoulder height and ``plateau_estimator_hist`` gave 20155.  The difference is
    not important as the plateau is not exactly constant; its value to a few significant
    values is the important quantity.

The position of the plateau/shoulder is somewhat sensitive to input parameters and can be
varied with changing the time step ``tau`` or the ``cluster_multispawn_threshold`` (if
applicable) for example, more details below. A large initial population ``init_pop`` can
also lead to overshooting of the shoulder.

Effects of the Time Step
------------------------

Now we will run another calculation with a higher time step, see input file below:

.. literalinclude:: calcs/shoulder/h2o_plat_bigtau.lua
        :language: lua

The two resulting shoulders are shown in the following graph:

.. plot::

    import pyhande
    import matplotlib.pyplot as plt
    (metadata, qmc_data) = pyhande.extract.extract_data('calcs/shoulder/h2o_plat.out')[0]
    (metadata, qmc_data2) = pyhande.extract.extract_data('calcs/shoulder/h2o_plat_bigtau.out')[0]
    plt.loglog(qmc_data['# H psips'], (qmc_data['# H psips'] / qmc_data['N_0']) , label=r'$\tau = 10^{-4}$')
    plt.loglog(qmc_data2['# H psips'], (qmc_data2['# H psips'] / qmc_data2['N_0']), label=r'$\tau = 10^{-3}$')
    plt.legend()
    plt.xlabel('# particles')
    plt.ylabel('# particles / # particles on reference determinant')

A smaller time step leads to fewer particles at the shoulder position, as described in [Booth09]_, [Vigor16]_.

Effects of Cluster Multispawn Threshold
---------------------------------------

.. review - JSS: some explanation of why this helps would provide welcome insight into the method.

This part looks at changing the multispawn threshold. This is another feature which can change the number of particles at the shoulder. 
Positive effects of that have already been shown in :ref:`CCMC <ccmc_tutorial>`.
Note that while changing the time step changes the position of the plateau for
FCIQMC for example as well, cluster multispawn threshold is specific to CCMC. 

To demonstrate the effects of decreasing the multispawn threshold, we will run
the following calculation with a low multispawn threshold:

.. literalinclude:: calcs/shoulder/h2o_plat_smallmsc.lua
        :language: lua

The plot below compares the shoulder plot of this and the first calculation on top of this tutorial:

.. plot::

    import pyhande
    import matplotlib.pyplot as plt
    (metadata, qmc_data) = pyhande.extract.extract_data('calcs/shoulder/h2o_plat.out')[0]
    (metadata, qmc_data2) = pyhande.extract.extract_data('calcs/shoulder/h2o_plat_smallmsc.out')[0]
    plt.loglog(qmc_data['# H psips'], (qmc_data['# H psips'] / qmc_data['N_0']) , label='multispawn threshold = none')
    plt.loglog(qmc_data2['# H psips'], (qmc_data2['# H psips'] / qmc_data2['N_0']), label='multispawn threshold = 0.1')
    plt.legend(loc="best")
    plt.xlabel('# particles')
    plt.ylabel('# particles / # particles on reference determinant')

Note that "multispawn threshold = none" means that there is no threshold within
computer number representation limits.

Clearly, setting a low multispawn threshold lowers the total number of particles at
the shoulder.

Effects of Initial Population
-----------------------------

.. review - JSS: again, some explanation of why this overshoot happens would provide welcome insight into the method.

In this part of the tutorial we will see that a large initial population can
lead to overshooting the shoulder. 

As a demonstration, we look at almost the same calculation as the first one but
with a larger initial population. 

.. literalinclude:: calcs/shoulder/h2o_plat_bigstartpop.lua
	:language: lua

The following plot compares the original with the calculation starting with a
large initial population:

.. plot::

    import pyhande
    import matplotlib.pyplot as plt
    (metadata, qmc_data) = pyhande.extract.extract_data('calcs/shoulder/h2o_plat.out')[0]
    (metadata, qmc_data2) = pyhande.extract.extract_data('calcs/shoulder/h2o_plat_bigstartpop.out')[0]
    plt.loglog(qmc_data['# H psips'], (qmc_data['# H psips'] / qmc_data['N_0']), label='initial population = 200')
    plt.loglog(qmc_data2['# H psips'], (qmc_data2['# H psips'] / qmc_data2['N_0']), label='initial population = 800')
    plt.legend(loc="best")
    plt.xlabel('# particles')
    plt.ylabel('# particles / # particles on reference determinant')

.. review - JSS: what do we mean by stable?

We see that the calculation with a larger initial population overshoots the
shoulder. However, it still forms a shoulder and is stable.

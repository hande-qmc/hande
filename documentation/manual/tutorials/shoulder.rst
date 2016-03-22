.. _shoulder_tutorial:

Shoulder Plots 
==============

This tutorial looks further into finding the optimal target particle population
in more detail. It is advisable to have read the :ref:`FCIQMC <fciqmc_tutorial>` and :ref:`CCMC <ccmc_tutorial>` tutorials
before this one. More information and details on shoulder plots can be found in [Spencer15]_.
One technique for finding the minimum number of required particles is looking
for a plateau in the total population vs iteration plot, as described in the
:ref:`FCIQMC <fciqmc_tutorial>` tutorial.
However, this plateau might not be easy to find by visual inspection and in fact might not even be there. Looking
at the ratio of total population to population on reference determinant vs total
population plot might be easier as there often is a "shoulder" at the plateau point.

.. review:: RSTF - What about the plateau detection functions in pyhande (pyhande.analysis.plateau_estimator{,_hist})?  Surely they should be at the very least mentioned.
.. reviewanswer:: VAN - added here as a note. Could add python script if I get pyhande path to work.

.. review:: RSTF - you keep referring to the CCMC tutorial, but there is minimal discussion of plateaus/shoulders there (and certainly no relevant graph) - better to refer to the FCIQMC one.
.. reviewanswer:: VAN - it seems you have fixed that already. Oh and thanks for grammar fixes

The example here is water in cc-pVDZ basis [Dunning89]_ [correctcite???]. Similarly to :ref:`CCMC <ccmc_tutorial>`, the
integrals were calculated with PSI4 (see :ref:`<generating_integrals>` for details).
Input and output files are in ``documentation/manual/tutorials/calcs/shoulder/``.
The first calculation was started as

.. literalinclude:: calcs/shoulder/h2o_plat.lua
	:language: lua

As in :ref:`CCMC <ccmc_tutorial>`, a plateau can be seen in the total population vs iteration
plot which indicates roughly the minimum particle number to make the calculation
stable.

.. review:: RSTF - maybe useful to also plot N_0 on here?
.. reviewanswer:: VAN - done.

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

The plateau is found at around 20000 particles. In the beginning, the only
occupied determinant is the reference determinant. These particles then spawn
onto other available determinants, making the total population grow at a greater
pace than the population on the reference determinant does. At the plateau
point, death and annihilation cancel out the spawning until the total population
grows again. However, the new rate of growing is less than the growth rate of
particles on the reference determinant. See [Spencer15]_ for details. This
implies that it is informative to consider the ratio of total population to
population on the reference determinant which we do in the form of "shoulder"
plots as shown below:

.. review:: RSTF - Explain why it is a relevant/meaningful thing to plot and the relationship with the plateau
.. reviewanswer:: VAN - done.

.. plot::

    import pyhande
    import matplotlib.pyplot as plt
    (metadata, qmc_data) = pyhande.extract.extract_data('calcs/shoulder/h2o_plat.out')[0]
    plt.loglog(qmc_data['# H psips'], (qmc_data['# H psips'] / qmc_data['N_0']))
    plt.xlabel('# particles')
    plt.ylabel('# particles / # particles on reference determinant')

The position of the shoulder is at about 20000 which corresponds to the position
of the plateau. 
.. review:: RSTF - I'm not sure what you mean by that last sentence - you seem to be making a different distinction between shoulder and plateau than usual.
.. reviewanswer:: VAN - I have deleted that last sentence.

.. note:: :ref:`pyhande` contains two functions to estimate the position of the
	plateau/shoulder. These are ``plateau_estimator`` and ``plateau_estimator_hist`` and they are described in
	detail in ``pyhande.analysis``, see :ref:`pyhande`. 
	``plateau_estimator`` gave 18480 with an estimated standard error of 40 for the shoulder height and ``plateau_estimator_hist`` gave 20155. 

The position of the shoulder can be varied with changing the time step ``tau`` or the
``cluster_multispawn_threshold`` (if applicable) for example, more details below. A large
initial population ``init_pop`` can lead to overshooting of the shoulder. 

Effects of the Time Step
------------------------

.. review:: RSTF - The other tutorials don't use the passive voice, it's a bit jarring.
.. reviewanswer:: VAN - changed to active future/present.

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

A smaller time step leads to less particles at the shoulder position.


Effects of Cluster Multispawn Threshold
---------------------------------------

.. review:: RSTF - mention this is specific to CCMC (whereas everything previously applies equally much to FCIQMC (and DMQMC?))
.. reviewanswer:: VAN - modified. What about DMQMC?

This part looks at changing the multispawn threshold. This is another feature which can change the number of particles at the shoulder. 
Positive effects of that have already been shown in :ref:`CCMC <ccmc_tutorial>`.
Note that while changing the time step changes the position of the plateau for
FCIQMC for example as well, cluster multispawn threshold is specific to CCMC. 

To demonstrate the effects of decreasing the multispawn threshold, we will run
the following calculation with a low multispawn threshold:

.. literalinclude:: calcs/shoulder/h2o_plat_smallmsc.lua
        :language: lua

The plot below compares the shoulder plot of this and the first calculation on top of this tutorial:

.. review:: RSTF - Probably better not to refer to the multispawn threshold as a 'cutoff' to avoid confusion with the spawn_cutoff option for real amplitudes.
.. reviewanswer:: VAN - Changed.


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

We see that the calculation with a larger initial population overshoots the
shoulder. However, it still forms a shoulder and is probably stable.  

.. _shoulder_tutorial:

Shoulder Plots 
==============

This tutorial looks further into finding the optimal target particle population
in more detail. It is advisable to have read the ccmc and/or fciqmc tutorials
before this one. More information and details on shoulder plots can be found in [Spencer15]_. 
One technique for finding the minimum number of required particles is looking
for a plateau in the total population vs iteration plot, as described in the :ref:`CCMC <ccmc_tutorial>`. 
However, this plateau might not be easy to find and in fact might not even be there. Looking 
at the ratio of total population to population on reference determinant vs total
population plot might be easier as there often is a "shoulder" at the plateau point. 

The example here is water in cc-pVDZ basis [Dunning89]_ [correctcite???]. Similarly to :ref:`CCMC <ccmc_tutorial>`, the 
integrals were calculated with Psi4 [cite???] (see :ref:`CCMC <ccmc_tutorial>` for details). 
Input and output files are in ``documentation/manual/tutorials/calcs/shoulder/``. 
The first calculation was started as 

.. literalinclude:: calcs/shoulder/h2o_plat.lua
	:language: lua

As in :ref:`CCMC <ccmc_tutorial>`, a plateau can be seen in the total population vs iteration 
plot which indicates roughly the minimum particle number to make the calculation 
stable. 

.. plot::

    import pyhande
    import matplotlib.pyplot as plt
    (metadata, qmc_data) = pyhande.extract.extract_data('calcs/shoulder/h2o_plat.out')[0]
    plt.semilogy(qmc_data['iterations'], qmc_data['# H psips'])
    plt.xlabel('iteration')
    plt.ylabel('# particles')

The plateau is found just below 20000 particles. Another possible plot to look at is a
"shoulder plot" as shown below:

.. plot::

    import pyhande
    import matplotlib.pyplot as plt
    (metadata, qmc_data) = pyhande.extract.extract_data('calcs/shoulder/h2o_plat.out')[0]
    plt.loglog(qmc_data['# H psips'], (qmc_data['# H psips'] / qmc_data['N_0']))
    plt.xlabel('# particles')
    plt.ylabel('# particles / # particles on reference determinant')

The position of the shoulder is just under 20000 which corresponds to the position
of the plateau. Shoulders are often easier to find than plateaus. 

The position of the shoulder can be varied with changing the time step ``tau`` or the
``cluster_multispawn_threshold`` for example, more details below. 

Effects of the Time Step
------------------------

Another run was done with a higher time step, see input file below:

.. literalinclude:: calcs/shoulder/h2o_plat_bigtau.lua
        :language: lua

The two resulting shoulders are shown in the following graph:

.. plot::

    import pyhande
    import matplotlib.pyplot as plt
    (metadata, qmc_data) = pyhande.extract.extract_data('calcs/shoulder/h2o_plat.out')[0]
    (metadata, qmc_data2) = pyhande.extract.extract_data('calcs/shoulder/h2o_plat_bigtau.out')[0]
    plt.loglog(qmc_data['# H psips'], (qmc_data['# H psips'] / qmc_data['N_0']) , label=r'$(\tau) = 0.0001$')
    plt.loglog(qmc_data2['# H psips'], (qmc_data2['# H psips'] / qmc_data2['N_0']), label=r'$(\tau) = 0.001$')
    plt.legend()
    plt.xlabel('# particles')
    plt.ylabel('# particles / # particles on reference determinant')

A smaller time step leads to less particles at the shoulder position.


Effects of Cluster Multispawn Threshold
---------------------------------------

Another feature that can decrease the number of particles at the shoulder is
setting a low multispawn threshold. Positive effects of that have already been
shown in :ref:`CCMC <ccmc_tutorial>`.   

The following calculation with a low multispawn threshold was run:

.. literalinclude:: calcs/shoulder/h2o_plat_smallmsc.lua
        :language: lua

The shoulder plot of this and the first calculation on top of this tutorial are
compared in the plot below:

.. plot::

    import pyhande
    import matplotlib.pyplot as plt
    (metadata, qmc_data) = pyhande.extract.extract_data('calcs/shoulder/h2o_plat.out')[0]
    (metadata, qmc_data2) = pyhande.extract.extract_data('calcs/shoulder/h2o_plat_smallmsc.out')[0]
    plt.loglog(qmc_data['# H psips'], (qmc_data['# H psips'] / qmc_data['N_0']) , label='multispawn cutoff = none')
    plt.loglog(qmc_data2['# H psips'], (qmc_data2['# H psips'] / qmc_data2['N_0']), label=r'multispawn cutoff = $0.1$')
    plt.legend(loc="best")
    plt.xlabel('# particles')
    plt.ylabel('# particles / # particles on reference determinant')

Note that "multispawn cutoff = none" means that there is no cutoff within
computer number representation limits. 

Clearly, setting a low multispawn cutoff lowers the total number of particles at
the shoulder.  

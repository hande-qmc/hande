.. _dmqmc_tutorial:

Density Matrix Quantum Monte Carlo
==================================

In this tutorial we will run DMQMC on the 2D Heisenberg model and the uniform electron gas.
The input and output files can be found under the ``documentation/manual/tutorials/calcs/dmqmc``
subdirectory of the source distribution.  Knowledge of the terminology and theory given in
[Booth09]_, [Blunt14]_ and [Malone15]_ is assumed.

To begin we will focus on the 6x6 antiferromagnetic Heisenberg model on a square lattice with periodic
boundary conditions. The input file for this system is given as

.. literalinclude:: calcs/dmqmc/heisenberg_dmqmc.lua

and is largely analagous to that found in the :ref:`FCIQMC tutorial <fciqmc_tutorial>`. We
refer the reader to the discussion there and the manual for system specific input options.
Note that ``init_pop`` here controls the population with which the density matrix at
:math:`\tau=0` is sampled. Typically the shift is allowed to vary from the beginning of
a simulation by setting ``target_pop`` equal to ``init_pop``. Here we will attempt to run to
a final temperature of :math:`\tau=5\beta J`.
The ``beta_loops`` option determines the number of independent simulations over which
observables are averaged, see :ref:`dmqmc_table` for more options. The operators table
specifies which observables are to be evaluated in a given simulation. Here only the total
energy is considered, a full list is available in :ref:`operators_table`.

An issue encountered when applying DMQMC to larger systems is that the population on the diagonal
(denoted Trace in the output file) decays with increasing :math:`\beta` which results in
poor estimates for observables. The seriousness of this problem needs to be assessed on
a system by system basis and should be tested for as a first step, which we'll do now.

To do this we set ``beta_loops`` to 1 in the input file and run the code as:

.. code-block:: bash

    $ aprun -N 1 hande.x heisenberg_dmqmc_excit_dist.lua > heisenberg_dmqmc.out

We find that for this system the population on the diagonal does indeed decay to quite
a low value:

.. plot::

    import pyhande
    import matplotlib.pyplot as plt
    dmqmc_data = pyhande.extract.extract_data('calcs/dmqmc/heisenberg_dmqmc.out
    for (md, data) in 


The source of and solution to this problem is outlined in [Blunt14]_. The source of this
problem can be investigated by analysing the distribution of psips on different excitation
levels of the density matrix where we see the weight redistributes from the diagonal to
highly excited determinants. This was calculated in aticipation of this result using the
excit_dist option in the operators table.

.. plot::
    excitation level

To overcome this [Blunt]_ invented an importance sampling scheme to encourage psips to
stay on or near the diagonal by penalising spawning moves away from excitation levels.
This is justified as typically the majority of the weight contributing to most physically
significant observables originates from the determinants at lower excitation levels which
we wish to sample more regularly.

Practically this amounts to first running a calculation with the ref:`find_weights`
option. This will output the importance sampling weights necessary as input for the
production calculation. It is worthwhile to run the calculation for a few ref:`beta_loops`
to ensure the weights are fluctuating too much, and also check they don't fluctuate too
much with the ref:`target_population`. The algorithm currently tries to ensure that the
number of walkers on each excitation level is roughly constant once the ground state is
thought to have  been to be reached. The iteration number where this is deemed to have
been reached is controlled by the ref:`find_weights_start` option.

Running the calculation with the chosen weights results in a much improved distribution of
psips among excitation levels:

.. plot::
    importance sampling

and we are now in a position to run the production calculation by setting the number of
beta loops to 100.
.. code-block:: bash
    aprun -N 1 heisenberg_dmqmc_production.lua > heisenberg_dmqmc_production.out 

The results can be analysed using the finite_temperature_analysis.py script provided in
tool/dmqmc:

.. code-block:: bash
    finite_temp_analysis.py heisenberg_dmqmc_production.out > heisenberg_dmqmc_production_analysis.out

Finally we can plot the results of the total energy as a function of temperature:

Interaction Picture Density Matrix Quantum Monte Carlo
======================================================

It turns out that the original formulation of DMQMC can run into problems for moderately
weakly interacting systems which are relatively well described by Hartree--Fock theory. An
extreme example of this is the uniform electron gas (UEG) especially at higher densities
(low :math:`r_s`). This issue is largely overcome by switching to the interaction picture
which enables us to start from a (temperature dependent) mean-field distribution at
:math:`tau=0` ensuring low energy determinants are initially sampled. See [Malone15]_ for
details. 

Most of the running details for IP-DMQMC are the same as for DMQMC, however there are some
additional considerations. This is best demonstrated by running a simulation. We will
focus on a 7-electron, spin polarised system in 319 plane waves at :math:`r_s=1`.

Looking at the input file 

.. literalinclude:: calcs/dmqmc/ipdmqmc_ueg.lua

we see most of the same options are present as for dmqmc. Note that unlike DMQMC where
estimates for the whole temperature range are gathered in a single simulation, in IP-DMQMC
only one temperature value is (directly) accessible, specified by the ref:`init_beta`
option. We've also set the energy scale to be determined by the Fermi energy of the
corresponding (thermodynamic limit) free electron gas so that the temperatures are
interpreted as fractions of the Fermi temperature (here :math:`\Theta = 0.5`.
ref:`all_sym_sectors` ensures all momentum symmetry sectors are averaged over. To average
over spin polarisation the ref:`all_spin_sectors` option must be specified.

Moving on through the ipdmqmc table we've set the the ref:`initial_matrix` to be the free
electron density matrix, i.e., Fermi-Dirac like. Additionally we're using the
ref:`grand_canonical_initialisation` option to initialise this density matrix (see
[Malone15]_). This is the recommended method to initialise the density matrix, the
Metropolis algorithm should only be used for testing. This options requires an additional
level of input in the form of the ref:`chem_pot` option in the system table. The chemical
potential for this system can be determined using /tools/dmqmc/chem_pot.py script as
follows:

.. code-block:: bash

    ./chem_pot.py ueg -p 7 2 1 10 

Finally we will use the asymmetric form of the original ip-dmqmc algorithm by specifying
symmetric=false. The symmetric algorithm is somewhat experimental but can lead to better
estimates for quantities other that the internal energy especially at lower temperatures. 

Running the code

.. code-block:: bash

    aprun -N 1 hande.x ipdmqmc_ueg.lua > ueg.out

and analysing the output we find:

.. code-block:: bash

    ./finite_temp_analysis.py ueg.out | tail -n 1 

where again only estimates at the final iteration are physical, i.e., when
:math:`tau=\beta`.

The initiator approximation can significantly extend the range of applicability of dmqmc
but is somewhat experimental. See the options, inparticular :ref:`initiator_level` in the
manual for more discussion. The user should ensure results are meaningful by comparing
answers at various walker populations.

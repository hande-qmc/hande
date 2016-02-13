.. _canonical_energy_tutorial:

.. [review] - JSS: appropriate name?
.. [reply]  - FDM: ?

Canonical Energies
==================

In this tutorial we will discuss how estimates for various expectation values of
partitions of the Hamiltonian can be evaluated in the canonical ensemble at finite
temperatures. These estimates are useful for basis set extrapolation as well as
comparison to the fully interacting results and are non-trivial to evaluate analytically.
See [Malone15]_ for details.

The input file is fairly simple:

.. literalinclude:: calcs/canonical_energy/canonical_energy.lua
    :language: lua

Much like in FCIQMC we bin data into blocks of ``nattempts`` and then run the simulation for
``ncycles*nattempts`` iterations in total. The only other options available are the inverse
temperature desired, which can be scaled by the Fermi temperature (where appropriate).
Here we restrict ourself to the fully spin polarised UEG in M=389 plane waves, which can be
compared to the IP-DMQMC simulation in the :ref:`DMQMC tutorial <dmqmc_tutorial>`.

Running the input file we find

.. [review] - JSS: is this ever expensive enough to require running in parallel?
.. [reply]  - FDM: mmm not necessarily, it is essentially embarrassingly parallel
.. [reply]  - (although not how I've implemented it) so you can get answers way quicker, but theres's no requirement.

.. code-block:: bash

    $ hande.x canonical_energy.lua >  canonical_energy.out

Inspecting the :download:`output <calcs/canonical_energy/canonical_energy.out>`, we see
a number of columns for various estimates including the kinetic, potential and internal
energy - precise definitions of everything can be found in the output file. Analysing the
data

.. [review] - JSS: any issue with correlated data?  Each iteration is independent right?  (So multiple calculations can be concatenated?)
.. [reply] - FDM: Completely independent you can combine multiple calculations.

.. code-block:: bash

    $ ./tools/dmqmc/analyse_canonical.py canonical_energy.out

Inspecting the output we see

.. [review] - JSS: not really for the tutorial, but is it worth adding options to analyse_canonical.py (and other DMQMC tools) for pretty-printing the statistics, or CSV output, as in reblock_hande.py?
.. [reply] - FDM: Probably, I don't like reading csv files though.

.. literalinclude:: calcs/canonical_energy/canonical_energy_res.out

and in particular, we can compare the values of :math:`\langle H \rangle_0` and
:math:`\langle H\rangle_{HF}` to the value of 32.92(2) Ha from the IP-DMQMC tutorial.

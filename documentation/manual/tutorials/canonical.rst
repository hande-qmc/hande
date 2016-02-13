.. _canonical_energy_tutorial:

Canonical Energies
==================

In this tutorial we will discuss how estimates for various expectation values of
partitions of the Hamiltonian can be evaluated in the canonical ensemble at finite
temperatures. These estimates are useful for basis set extrapolation as well as
comparison to the fully interacting results and are non-trivial to evaluate analytically.

The input file is fairly simple:

.. literalinclude:: calcs/canonical_energy/canonical_energy.lua
    :language: lua

Much like in FCIQMC we bin data into blocks of nattempts and then run the simulation for
ncycles*nattempts iterations in total. The only other options available are the inverse
temperature desired, which can be scaled by the fermi_temperature (where appropriate).
Here we restrict ourself to the fully spin polarised UEG in M=389 plane waves which can be
compared to the ipdmqmc simulation in the :ref:`DMQMC tutorial <dmqmc_tutorial>`.

Running the input file we find

.. code-block:: bash

    $ ./hande.x canonical_energy.lua >  canonical_energy.out

Inspecting the :download:`output <calcs/canonical_energy/canonical_energy.out>`, we see
a number of columns for various estimates including the kinetic, potential and internal
energy - precise definitions of everything can be found in the output file. Analysing the
data

.. code-block:: bash

    $ ./tools/dmqmc/analyse_canonical.py canonical_energy.out > canonical_energy_res.out

Inspecting the output we see

.. literalinclude:: calcs/canonical_energy/canonical_energy_res.out

and in particular, we can compare the values of :math:`\langle H \rangle_0` and
:math:`\langle H\rangle_{HF}` to the value of 32.92(2) Ha from the ip-dmqmc tutorial.

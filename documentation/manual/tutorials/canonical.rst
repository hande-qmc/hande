.. _canonical_estimates_tutorial:

Canonical Estimates
===================

In this tutorial we will discuss how estimates for various mean-field properties of
a system can be evaluated in the canonical ensemble at finite temperatures. These
estimates are useful for basis set extrapolation as well as comparison to the fully
interacting results and are non-trivial to evaluate analytically.  See [Malone15]_ for
details.

The input file is fairly simple:

.. literalinclude:: calcs/canonical_estimates/canonical_estimates.lua
    :language: lua

Here we attempt to generate N particle states making ``nattempts`` attempts and then run
the simulation for ``ncycles*nattempts`` iterations in total. The only other options
available are the inverse temperature desired, which can be scaled by the Fermi
temperature (where appropriate).  Here we restrict ourself to the fully spin polarised UEG
in M=389 plane waves, which can be compared to the IP-DMQMC simulation in the :ref:`DMQMC
tutorial <dmqmc_tutorial>`.

Running the input file we find

.. code-block:: bash

    $ hande.x canonical_estimates.lua >  canonical_estimates.out

Inspecting the :download:`output <calcs/canonical_estimates/canonical_estimates.out>`, we
see a number of columns for various estimates including the kinetic, potential, internal,
free energy and entropy - precise definitions of everything can be found in the output
file.  The data can be analysed to find the mean and standard error using the
``analyse_canonical.py`` script in the ``tools/dmqmc`` subdirectory:

.. code-block:: bash

    $ analyse_canonical.py canonical_estimates.out

which gives

.. literalinclude:: calcs/canonical_estimates/canonical_estimates_res.out

In particular, we can compare the values of :math:`U_0` and
:math:`U_{\mathrm{HF}}` to the value of 32.92(2) Ha from the IP-DMQMC tutorial.

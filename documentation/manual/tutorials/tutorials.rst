.. _tutorials:

Tutorials
=========

The tutorials below demonstrate how to set up and run Monte Carlo calculations in HANDE.
The input files in the test suite also demonstrate how calculations can be performed.

The input and output files from the calculations performed in the tutorials can be
found under the ``documentation/manual/tutorials/calcs/`` directory and can be downloaded
from [tutorial_dataset]_.

.. todo - create dataset on Zenodo and provide DOI link.

Note that none of the tutorials fix a random number seed (as this is the best approach for
running multiple production calculations on the same system) so results will not be
exactly identical (but should agree statistically) from those in the above dataset unless
the same seeds (which can be found in the output files) are used.

The test suite contains example input files for a wide variety of scenarios.

We recommend working through the FCIQMC tutorial before the iFCIQMC, CCMC or DMQMC tutorials.

.. toctree::
    :maxdepth: 1

    fciqmc
    ifciqmc
    ccmc
    dmqmc
    kinetic

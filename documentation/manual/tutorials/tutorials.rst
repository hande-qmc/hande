.. _tutorials:

Tutorials
=========

The tutorials below demonstrate how to set up and run Monte Carlo calculations in HANDE.
The input files in the test suite also demonstrate how calculations can be performed.  The
aim here is to provide an introduction to setting up, running and analysing calculations
and only basic input options are considered; for advanced options please consult the
appropriate section of the manual.

The tutorials assume that HANDE has been successfully compiled and the test suite has been
sucessfully run.  Any reference to ``hande.x`` should be replaced with the full path to
the HANDE executable and similarly for the ``reblock_hande.py`` script.

The input and output files from the calculations performed in the tutorials can be
found under the ``documentation/manual/tutorials/calcs/`` directory.  The example
calculations are deliberately not trivial and may require up to a few hundred core hours
to run as shown.  Smaller calculations can be performed by reducing the system size (e.g.
using fewer electrons or orbitals) or running for fewer iterations.

.. note::

    None of the tutorials fix a random number seed (as this is the best approach for
    running multiple production calculations on the same system) so results will not be
    exactly identical (but should agree statistically) from those in the above dataset
    unless the same seeds (which can be found in the output files) are used.

We recommend working through the FCIQMC tutorial before the iFCIQMC, CCMC or DMQMC tutorials.

.. toctree::
    :maxdepth: 1

    fciqmc
    ifciqmc
    semi_stoch
    ccmc
    dmqmc
    canonical
    shoulder

All calculations were analysed using :ref:`pyhande` and all graphs were plotted using
`matplotlib <http://matplotlib.org/>`_.  Parts of the plot generation code were
adapted from the matplotlib tutorials.

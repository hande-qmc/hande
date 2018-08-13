.. _solids_tutorial:

Solid state calculations
==================================================

In this tutorial we will run periodic boundary conditions CCMC calculation on a diamond
crystal in STO-3G basis with 2x1x1 sampling of the Brillouin zone. Familiarity with :ref:`CCMC <ccmc_tutorial>` and :ref:`FCIQMC <fciqmc_tutorial>` tutorials is assumed. The input and output files can be found under the ``documentation/manual/tutorials/calcs/ccmc_solids`` subdirectory of the source distribution.

First of all, we need the one- and two- electron integrals from an external source. We will
use PySCF_ software package to perform preliminary Hartree-Fock calculation
and generate the integrals. PySCFDump script can be used to save the integrals in the FCIDUMP
format readable by HANDE.

.. _PySCF: http://www.pyscf.org/

To correctly address exchange divergence, additional exchange integrals are needed.
Those are written in a FCIDUMP_X file by PySCFDump. Theoretical details of this procedure
will be elaborated on as part of an upcoming paper.

Cell object can be conveniently prepared using ``build_cell`` function of the pyscfdump.helpers
module and ASE_ library. Set of basis functions, kinetic energy cutoff and pseudopotential
used are specified.

.. _ASE: http://wiki.fysik.dtu.dk/ase/

The ``run_khf`` function of the pyscfdump.scf module is used to run the HF calculation.
Number of k-points in each dimension is specified. By default, the Monkhorst-Pack grid is used. If
``gamma`` is set to true, the grid will be shifted to include the Gamma point. Exchange divergence
treatment scheme should also be chosen. The function returns the converged HF calculation object
and a list of scaled k-points (i.e. (0.5,0.5,0.5) is at the very corner of the Brillouin zone).

Finally, ``fcidump`` function of the ``pyscfdump.pbcfcidump`` module is used to dump the integrals
to a file. Name of the resulting file, SCF calculation object, number of k points, list of
scaled k points and boolean variable ``MP`` (true if the grid is Monkhorst-Pack - i.e. not shifted to
include the gamma point for even grids) must be provided.

.. literalinclude:: calcs/ccmc_solids/diamond_HF.py
    :language: python

Now we are ready to run the CCMC calculation. System definition is read in from the FCIDUMP files.
There are two points to notice

    - path to files with additional exchange integrals is specified as ``ex_int_file``

    - to properly exploit translational symmetry in the crystal lattice, the orbitals and hence the integrals must be complex. The complex mode of HANDE is enabled by setting ``complex = true``. Consequently, numbers of particles on excitors also have both real and imaginary part, as well as projected energy.

Choice of parameters for the calculation is beyond the scope of this tutorial and generally
requires some trial and error. A shoulder plot should be used to determine target population.
The input script used in this tutorial is:

.. literalinclude:: calcs/ccmc_solids/ccmc.lua
    :language: lua

The Hilbert space of the system we are dealing with is quite small and so is the target population.
The general rule is that in order to use MPI parallelism, each process should contain at least :math:`1 \times 10^{5}`
excitors. Having less excitors on each process is both inefficient and in extreme cases can lead
to biased results. This is why we will use only one process here. However, use of openMP shared 
memory threads is recommended in order to make full use of the available resources.

.. code-block:: bash

    $ hande.x ccmc.lua > diamond_ccmc.out

We can now plot the population

.. plot::

    import pyhande
    import matplotlib.pyplot as plt
    (metadata, qmc_data) = pyhande.extract.extract_data('calcs/ccmc_solids/diamond_ccmc.out')[0]
    plt.plot(qmc_data['iterations'], qmc_data['# H psips'])
    plt.xlabel('iteration')
    plt.ylabel('# particles')

and correlation energy:

.. plot::

    import pyhande
    import matplotlib.pyplot as plt
    (metadata, qmc_data) = pyhande.extract.extract_data('calcs/ccmc_solids/diamond_ccmc.out')[0]
    plt.plot(qmc_data['iterations'], qmc_data['Shift'], label=r'$S(\tau)$')
    plt.plot(qmc_data['iterations'], qmc_data['Re{\sum H_0j N_j}']/qmc_data['Re{N_0}'], label=r'$E(\tau) = \Re(\sum_j H_{0j} N_j(\tau))/\Re(N_0(\tau))$')
    plt.legend()
    plt.ylim(-0.4,0)
    plt.xlabel('iteration')
    plt.ylabel('Correlation energy/' + r'$E_h$')

It is worth noting that the projected energy is in fact a complex quantity, whose imaginary part evaluates to zero
in a non-trivial way. For the plot (as well as the following analysis) we calculate the result by only
using real parts of both :math:`\sum_j H_{0j} N_j(\tau)` and :math:`N_0(\tau)`, using the fact that

:math:`E(\tau) = \Re(\sum_j H_{0j} N_j(\tau)/(N_0(\tau)) = \Re(\sum_j H_{0j} N_j(\tau))/\Re(N_0(\tau))`

where the second equality holds provided that imaginary part of :math:`E(\tau)` is zero.
Anyway, the shift remains a strictly real measure of the correlation energy [#]_.

To analyse the calculation we can use reblock_hande.py script:

.. code-block:: bash

   $ reblock_hande.py --quiet diamond_ccmc.out

which results in:

.. literalinclude:: calcs/ccmc_solids/diamond.block


.. rubric:: Footnotes

.. [#] As discussed  in [Booth13]_ for FCIQMC - the CCMC case is exactly analogous.

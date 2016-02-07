.. _semi_stoch_tutorial:

Semi-Stochastic FCIQMC
======================

In this tutorial we will explain how to run FCIQMC calculations using the
semi-stochastic adaptation to reduce stochastic errors [Petruzielo12]_. We
will consider the half-filled 18-site 2D Hubbard model at :math:`U/t=1.3`,
as previously considered in the basic :ref:`fciqmc_tutorial` tutorial. In
particular, we shall begin from the input file presented at the end of the
FCIQMC tutorial, which introduces the use of non-integer psip amplitudes
through the ``real_amplitudes`` keyword:

.. literalinclude:: calcs/fciqmc/hubbard_fciqmc_real.lua
    :language: lua

which results in the following simulation:

.. plot::

    import pyhande
    import matplotlib.pyplot as plt
    (metadata, qmc_data) = pyhande.extract.extract_data('calcs/fciqmc/hubbard_fciqmc_real.out')[0]
    plt.ylim(-0.3315,-0.328)
    plt.plot(qmc_data['iterations'], qmc_data['\sum H_0j N_j']/qmc_data['N_0'], label='$E(\tau) = \sum_j H_{0j} N_j(\tau)/N_0(\tau)$')
    plt.legend()
    plt.xlabel('iteration')
    plt.ylabel('Correlation energy / $t$')

The semi-stochastic adaptation provides a way to reduce the stochastic noise
in such simulations. It does so by choosing a certain subspace, which is
deemed to be most important (in that most of the wave function amplitude
resdies in this subspace), and performing projection exactly within it.
Projection outside the subspace is performed by usual FCIQMC spawning rules.
Thus, we simply need to specify what subspace to use for the exact projection.
This is done in HANDE using the scheme of [Blunt15]_, where the subspace is
formed from the determinants on which the largest number of psips reside. We
therefore simply need to tell HANDE what iteration to start using the
semi-stochastic adaptation, and how many determinants to form the deterministic
subspace from.

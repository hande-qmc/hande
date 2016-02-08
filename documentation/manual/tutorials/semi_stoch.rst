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
resides in this subspace), and performing projection exactly within it.
Projection outside the subspace is performed by usual FCIQMC spawning rules.
Thus, we simply need to specify what subspace to use for the exact projection.
This is done in HANDE using the scheme of [Blunt15]_, where the subspace is
formed from the determinants on which the largest number of psips reside. We
therefore simply need to tell HANDE what iteration to start using the
semi-stochastic adaptation, and how many determinants to form the deterministic
subspace from.

Looking at the above simulation, it appears that the energy has converged by
iteration 20000. This is not a guarantee that the wave function is also
fully converged, but full convergence is not critical. A reasonable deterministic
space size is 10000. So, to start using a deterministic space of size 20000 at
iteration 10000, we modify the above input to the following:

.. literalinclude:: calcs/semi_stoch/hubbard_semi_stoch_high.lua
    :language: lua

Here, the ``semi-stoch`` table contains four keywords. The use of ``size`` and
``start_iteration`` keywords is hopefully clear. The ``space`` keyword determines
which method is used to generate the deterministic space - in this case it done
by choosing the determiniants with the highest weights. The ``separate_annihilation``
keyword is somewhat technical - setting it to false avoids the use of an extra MPI
communication per iteration when using semi-stochastic. While this can be a benefit,
it can slow other areas of annihilation down substantially in certain cases, and so
we advise the inexperienced user against this.

This results in the following simulation:

.. plot::

    import pyhande
    import matplotlib.pyplot as plt
    (metadata, qmc_data) = pyhande.extract.extract_data('calcs/semi_stoch/hubbard_semi_stoch_high.out')[0]
    plt.ylim(-0.3315,-0.328)
    plt.plot(qmc_data['iterations'], qmc_data['\sum H_0j N_j']/qmc_data['N_0'], label='$E(\tau) = \sum_j H_{0j} N_j(\tau)/N_0(\tau)$')
    plt.legend()
    plt.xlabel('iteration')
    plt.ylabel('Correlation energy / $t$')

As can be seen, at iteration 20000 there is a large reduction in stochastic error.

When performing a blocking analysis, the user should not begin averaging data until
after the semi-stochastic adaptation has been turned on, since there is a
significant change in the probability distributions of data beyond this point. This
is particularly true in initiator FCIQMC simulations, where the use of semi-stochastic
can alter the initiator error. We can therefore analyse the above simulation using

.. code-block:: bash

    $ reblock_hande.py --quiet --start 30000 hubbard_semi_stoch_high.out

which results in:

.. literalinclude:: calcs/semi_stoch/hubbard_semi_stoch_high.block

Compared to the equivalent non-semi-stochastic simulation performed in the
:ref:`fciqmc_tutorial` tutorial, the error bars on the shift and projected
energy estimators have reduced from :math:`4 \times 10^{-5}` and
:math:`3 \times 10^{-6}` to :math:`2 \times 10^{-5}` and :math:`8 \times 10^{-7}`,
respectively.

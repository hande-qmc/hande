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

and which results in the following simulation:

.. plot::

    import pyhande
    import matplotlib.pyplot as plt
    (metadata, qmc_data) = pyhande.extract.extract_data('calcs/fciqmc/hubbard_fciqmc_real.out')[0]
    plt.ylim(-0.3315,-0.328)
    plt.plot(qmc_data['iterations'], qmc_data['\sum H_0j N_j']/qmc_data['N_0'], label='$E(\\tau) = \sum_j H_{0j} N_j(\\tau)/N_0(\\tau)$')
    plt.legend()
    plt.xlabel('iteration')
    plt.ylabel('Correlation energy / $t$')

The semi-stochastic adaptation provides a way to reduce the stochastic noise
in such simulations. It does so by choosing a certain subspace, which is
deemed to be most important (in that most of the wave function amplitude
resides in this subspace), and performing projection exactly within it.
Projection outside the subspace is performed by the usual FCIQMC stochastic spawning.

.. [review] - JSS: insert some discussion about the memory implications and, as a result, the size of the deterministic space that is realistic?
.. [review] - JSS: is it always best to use as large a subspace as possible or is there a speed/memory/etc tradeoff?
.. [review] - JSS: This is done -> One way of doing this.  Hopefully CAS will be available soon!

Thus, we simply need to specify what subspace to use for the exact projection.
This is done in HANDE using the scheme of [Blunt15]_, where the subspace is
formed from the determinants on which the largest number of psips reside. We
therefore simply need to tell HANDE what iteration to start using the
semi-stochastic adaptation, and how many determinants to form the deterministic
subspace from.

..

    [review] -  JSS: I think the above could do with some expansion:

       + semi-stochastic does not (seem to?) improve initiator convergence in many cases
       + semi-stochastic is slower, so may as well run without it whilst converging.

   Is this a fair summmary?

.. [review] - JSS: state why not critical and why 10000 is appropriate size of subspace.
.. [review] - JSS: would an option to turn the semi-stochastic on after the shift is turned on be useful?

Looking at the above simulation, it appears that the energy has converged by
iteration 20000. This is not a guarantee that the wave function is also
fully converged, but full convergence is not critical. A reasonable deterministic
space size is 10000. So, to start using a deterministic space of size 20000 at
iteration 10000, we modify the above input to the following:

.. literalinclude:: calcs/semi_stoch/hubbard_semi_stoch_high.lua
    :language: lua

Here, the ``semi-stoch`` table contains three keywords. The use of ``size`` and
``start_iteration`` keywords is hopefully clear. The ``space`` keyword determines
which method is used to generate the deterministic space - in this case by choosing the
determiniants with the highest weights.

This results in the following simulation:

.. plot::

    import pyhande
    import matplotlib.pyplot as plt
    (metadata, qmc_data) = pyhande.extract.extract_data('calcs/semi_stoch/hubbard_semi_stoch_high.out')[0]
    plt.ylim(-0.3315,-0.328)
    plt.plot(qmc_data['iterations'], qmc_data['\sum H_0j N_j']/qmc_data['N_0'], label='$E(\\tau) = \sum_j H_{0j} N_j(\\tau)/N_0(\\tau)$')
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
:ref:`FCIQMC tutorial <fciqmc_tutorial>` tutorial, the error bars on the shift and projected
energy estimators have reduced from :math:`4 \times 10^{-5}` and
:math:`3 \times 10^{-6}` to :math:`2 \times 10^{-5}` and :math:`8 \times 10^{-7}`,
respectively.

.. [review] - JSS: But simulation took much longer.  The reduction means it is still worth it though compared to using stochastic propagation, right?

Note that if you do not specify a ``start_iteration`` value in the ``semi_stoch``
table of the input file, then the semi-stochastic adaptation will be turned
on from the first iteration. This should not be done when starting a new
simulation, because wave functions in HANDE are initialised as single determinants.
However, if restarting a simulation from an HDF5 file then this is a sensible approach -
the simulation will begin from the wave function stored in the HDF5 file, and the
deterministic space will be chosen from the most populated determinants in this
wave function. An input file for such a restarted simulation would contain the
following ``semi_stoch`` and ``restart`` tables within the ``fciqmc`` table:


.. [review] - JSS: just place partial input files directly in the source and use literalinclude just for files which are to be run or are outputs.

.. literalinclude:: calcs/semi_stoch/hubbard_semi_stoch_restart.lua
    :language: lua

(see the :ref:`restart_table` entry in the documentation for more options relating to
restarting simulations).

Finally, when restarting simulations which were already using the semi-stochastic
adaptation, it is important to use exactly the same deterministic space to ensure
that estimators are statistically consistent before and after restarting. However,
the approach in HANDE uses the instantaneous FCIQMC wave function to generate the
deterministic space, which changes during the simulation. Using the above approach
would therefore lead to a slightly different space being generated after restarting.
One can get around this by outputting the deterministic space in use to a file, and
reading it back in for the restarted calculation. For example, to generate a
deterministic space from the 10000 most populated determinants at iteration 20000,
and to then print this space to a file, one should use the ``write`` keyword in the
``semi-stoch`` table:

.. [review] - JSS: just place partial input files directly in the source and use literalinclude just for files which are to be run or are outputs.

.. literalinclude:: calcs/semi_stoch/hubbard_semi_stoch_write.lua
    :language: lua

Here, the value of the ``write`` keyword, :math:`0`, is an index used in the
name of the resulting file.

.. [review] - JSS: note write can also be a boolean.

When restarting the simulation, one can then specify the ``space`` option to
read a semi-stochastic HDF5 file, using:

.. [review] - JSS: just place partial input files directly in the source and use literalinclude just for files which are to be run or are outputs.

.. literalinclude:: calcs/semi_stoch/hubbard_semi_stoch_read.lua
    :language: lua

The deterministic space file is an HDF5 file. As such, both writing and reading
of such files requires :ref:`compilation <compilation>` of HANDE with HDF5 enabled, which is the
default compilation behaviour.

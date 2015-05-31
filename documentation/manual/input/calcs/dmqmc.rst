.. _dmqmc:

Density Matrix Quantum Monte Carlo
==================================

.. code-block:: lua

    dmqmc {
        sys = system,
        qmc = { ... },
        dmqmc = { ... },
        ipdmqmc = { ... },
        operators = { ... },
        rdm = { ... },
        restart = { ... },
        reference = { ... },
    }

``dmqmc`` performs a density matrix quantum Monte Carlo (DMQMC) calculation on a system.

Unlike :ref:`ccmc` and :ref:`fciqmc`, where quantites are averaged inside each report
loop, any quantities in DMQMC are evaluated at the **first** iteration of the report loop
only due to the explcit temperature dependence in the algorithm.

.. note::

    Density Matrix Quantum Monte Carlo is currently rather experimental.  In particular,
    it is not implemented for all systems yet and some options are only implemented for
    specific systems.

.. todo - note any options which are only implemented for certain systems.

Options
-------

``sys``
    type: system object.

    Required.

    The system on which to perform the calculation.  Must be created via a system
    function.
``qmc``
    type: lua table.

    Required.

    Further options that are common to all implemented QMC algorithms.  See
    :ref:`qmc_table`.
``dmqmc``
    type: lua table.

    Optional.

    Further options to control the DMQMC algorithm.  See :ref:`dmqmc_table`.
``ipdmqmc``
    type: lua table.

    Optional.

    If set, even to an empty table, then interaction picture DMQMC is performed.  The
    table can contain further options to control the IP-DMQMC algorithm.  See
    :ref:`ipdmqmc_table`.
``operators``
    type: lua table.

    Optional.

    Further options to select the operators for which expectation values are evaluated.
    See :ref:`operators_table`.
``rdm``
    type: lua table.

    Optional.

    Further options to select which (if any) reduced density matrices and corresponding
    operators are to be evaluated.  See :ref:`rdm_table`.
``restart``
    type: lua table.

    Optional.

    Further options to control restarting the calculation from a previous calculation.
    See :ref:`restart_table`.
``reference``
    type: lua table.

    Optional.

    Further options to select the reference state used.  See :ref:`reference_table`.

.. _dmqmc_table:

dmqmc options
-------------

``replica_tricks``
    type: boolean.

    Optional.  Default: false.

    Perform replica simulations (i.e. evolve two density matrices concurrently) if true.
    This allows calculation of unbiased estimators that are quadratic in the density
    matrix.
``fermi_temperature``
    type: boolean.

    Optional.  Default: false.

    .. todo - why is this in dmqmc rather than ipdmqmc?

     Interpret ``initial_beta`` as the inverse reduced temperature, :math:`\beta
     = 1/\Theta = T_F/T`, where :math:`T_F` is the Fermi temperature, instead of the
     inverse temperature in atomic units, :math:`\beta = 1/T`.
``all_sym_sectors``
    type: boolean.

    Optional.  Default: false.

    Sample states with all symmetris of the system instead of just those which conserve
    the symmetry of the reference state.
``all_spin_sectors``
    type: boolean.

    Optional.  Default: false.

    Sample states with all spin polarisations of the system instead of just those which
    conserve the spin polarisation of the reference state.
``beta_loops``
    type: integer.

    Optional.  Default: 100.

    The number of loops over the desired temperature range (each starting from
    :math:`T=\infty` and performing the desired number of iterations) to perform.  Each
    beta loop samples the initial conditions independently.

    .. note::

        Estimators must be averaged at each temperature from different beta loops.  As
        each beta loop is independent, this can be done in separate calculations in an
        embararassingly parallel fashion.

``sampling_weights``
    type: vector of floats.
    .. todo - How many floats need to be specified?

    Optional.  Default: none.

    .. todo - how to generate these?  (Describe in DMQMC tutorial?)

    .. todo - description
``vary_weights``
    type: integer.

    Optional.  Default: 0

    .. todo - description
``find_weights``
    type: boolean.

    Optional.  Default: false.

    .. todo - description
``symmetrize``
    type: boolean.

    Optional.  Default: false.

    Explicitly symmetrize the density matrix, thus only sampling one triangle of the
    matrix.  This can yield significant improvements in stochastic error in some cases.

.. _ipdmqmc_table:

ipdmqmc options
---------------

``initial_beta``
    type: float.

    Optional.  Default: 1.0.

    The inverse temperature to propogate the density matrix to.
``initial_matrix``
    type: string.

    Optional.  Default: 'hartree_fock'.

    Possible values: 'free_electron', 'hartree_fock'.

    Initialisation of the density matrix at :math:`\beta=0`.  'free_electron' samples the
    free electron density matrix, i.e. :math:`\rho = \sum_i e^{-\beta \sum_j \varepsilon_j
    \hat{n}_j} |D_i><D_i|`, where :math:`\varepsilon_j` is the single-particle eigenvalue
    and :math:`hat{n}_j` the corresponding number operator.  'hartree_fock' samples
    a 'Hartree--Fock' density matrix defined by :math:`\rho = \sum e^{-\beta H_ii} |D_i><D_i|`,
    where `:math:`H_ii = <D_i|H|D_i>` and is more efficient than 'free_electron'.

    .. todo - recommendations?  What is 'efficient' measured by?
``grand_canonical_initialsation``
    type: boolean.

    Optional.  Default: false.

    Use the grand canonical partition function to initialise the psip distribution.
    .. todo - how does this differ from the default behaviour?
``metropolis_attempts``
    type: integer.

    Optional.  Default: 0.

    .. todo - Ok, it seems I don't understand how to use these options in combination (or not).  Clearly a tutorial is required!

.. _operators_table:

operators options
-----------------

``renyi2``
    type: boolean.

    Optional.  Default: false.

    Calculate the Renyi-2 entropy of the entire system.  Requires ``replica_tricks`` to be
    enabled.
``energy``
    type: boolean.

    Optional.  Default: false.

    Calculate the thermal expectation value of the Hamiltonian operator.
``energy2``
    type: boolean.

    Optional.  Default: false.

    Calculate the thermal expectation value of the Hamiltonian operator squared.
``staggered_magnetisation``
    type: boolean.

    Optional.  Default: false.

    Calculate the thermal expectation value of the staggered magnetisation operator.
    Only available for bipartite Heisenberg lattices.
``excit_dist``
    type: boolean.

    Optional.  Default: false.

    Calculate the fraction of psips at each excitation level, where the excitation level
    is the number of excitation separating the two states labelling a given density matrix
    element.  Accumulated from ``excit_dist_start`` iterations onwards.
``excit_dist_start``
    type: integer.

    Optional.  Default: 0.

    The iteration number from which ``excit_dist`` is accumulated.
``correlation``
    type: 2D vector of integers.

    Optional.  Default: false.

    Calculate the spin-spin correlation function between the two specified lattice sites,
    :math:`i` and :math:`j`, which is defined as the thermal expectation value of:

    .. math::

    	\hat{C}_{ij} = S_{xi}S_{xj} + S_{yi}S_{yj} + S_{zi}S_{zj}.

.. _rdm_table:

rdm options
-----------

``rdms``
    type: table of 1D vectors.

    Required.

    Each vector corresponds to the subsystem of a reduced density matrix as a list of the
    basis function indices in the subsystem.  For example:

    .. code-block: lua

        rdms = { { 1, 2 } }

    specifies one RDM containing basis functions with indices 1 and 2, and

    .. code-block: lua

        rdms = { { 1, 2 }, { 3, 4} }

    specifies two RDMs, with the first containing basis functions with indices 1 and 2,
    and the second basis functions 3 and 4.

    Either ``instantaneous`` or ``ground_state`` must be enable to set the desired mode of
    evaluating the RDM.

    .. todo - Why do we have ``instantaneous`` and ``ground_state`` when they appear to be
              mutually incompatible?  [Nick?]

``instantaneous``
    type: boolean.

    Optional.  Default: false.

    Calculate the RDMs at each temperature based upon the instantaneous psip distribution.
``ground_state``
    type: boolean.

    Optional.  Default: false.

    Accumulate the RDM once the ground state (as specified by ``ground_state_start``)
    is reached.  This has two limitations: only one RDM can be accumulated in
    a calculation and the subsystem should be at most half the size of the system (which
    is always sufficient for ground-state calculations).
``spawned_state_size``
    type: integer.

    Required if ``instantaneous`` is true.  Ignored otherwise.

    Maximum number of states (ie reduced density matrix elements) to store in the
    "spawned" list, which limits the number of unique RDM elements that each processor can
    set.  Should be a sizable fraction of ``state_size`` (see :ref:`qmc_table`) and
    depends on the size of the subsystem compared to the full space.

    .. todo - should allow -ve numbers to specify the MB usage instead (as in the main
              spawned_state_size?)

    .. note::

        This is a **per processor** quantity.  It is usually safe to assume that each
        processor has approximately the same number of states.

``ground_state_start``
    type: integer.

    Optional.  Default: 0.

    Monte Carlo cycle from which the RDM is to be accumulated in each beta loop.  Relevant
    only if ``ground_state`` is set to true and, as such, should be set to a iteration
    (which is a measure of temperature) such that the system has reached the ground state.
``concurrence``
    type: boolean.

    Optional.  Default: false.

    Calculate the unnormalised concurrence and the trace of the reduced density matrix at
    the end of each beta loop.  The concurrence can be calculated from this using the
    ``average_entropy.py`` script.

    Valid for ``ground_state`` only; temperature-dependent concurrence is not currently
    implemented.
``renyi2``
    type: boolean.

    Optional.  Default: false.

    Calculate the Renyi-2 entropy of each subsystem, :math:`S^n_2 = \sum_{ij} (\rho^n_{ij})^2`,
    where :math:`\rho^n` is the reduced density matrix of the :math:`n`-th subsystem.  The
    temperature-dependent estimate of the Renyi-2 entropy can then be obtained using the 
    ``finite_temp_analysis.py`` script.

    Valid for ``instantaneous`` only; ground-state Renyi-2 averaged over a single beta
    loop is not currently implemented.  Requires ``replica_tricks`` to be enabled in order
    to obtained unbiased estimates.
``von_neumann``
    type: boolean.

    Optional.  Default: false.

    Calculate the unnormalised von Neumann entropy and the trace of the reduced density
    matrix at the end of each beta loop.  The von Neumann entropy can be calculated from
    this using the ``average_entropy.py`` script.

    Valid for ``ground_state`` only; temperature-dependent von Neumann entropy is not
    currently implemented.
``write``
    type: boolean.

    Optional.  Default: false.

    Print out the ground-state RDM to file at the end of each beta loop.  The file
    contains the trace of the RDM in the first line followed by elements of the upper
    triangle of the RDM labelled by their index.

    Valid for ``ground_state`` only.

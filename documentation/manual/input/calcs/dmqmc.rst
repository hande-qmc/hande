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
        qmc_state = qmc_state,
    }

Returns:
    a qmc_state object.

``dmqmc`` performs a density matrix quantum Monte Carlo (DMQMC) calculation on a system.

Unlike :ref:`ccmc` and :ref:`fciqmc`, where quantities are averaged inside each report
loop, any quantities in DMQMC are evaluated at the **first** iteration of the report loop
only. This is because different iterations represent different temperatures in DMQMC,
and so averaging over a report loop would average over different temperatures, which is
not the desired behaviour.

.. note::

    Density Matrix Quantum Monte Carlo is currently rather experimental.  In particular,
    it is not implemented for all systems yet and some options are only implemented for
    specific systems. In particular, DMQMC is only implemented for the Heisenberg model, the UEG,
    the real and momentum-space Hubbard model, and for molecular systems. The evaluation of operators
    other than the total energy, such as correlation functions and entanglement measures,
    is currently only possible for the Heisenberg model. The calculation of the reduced
    density matrices from DMQMC is also only supported for the Heisenberg model (for both
    temperature-dependent and ground state RDM calculations).

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

    If set, even to an empty table, then interaction picture DMQMC [Malone15]_ is
    performed.  The table can contain further options to control the IP-DMQMC algorithm.
    See :ref:`ipdmqmc_table`.
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
``qmc_state``
    type: qmc_state object.

    Optional.

    Output of a previous calculation to resume.

    .. warning::

        The qmc_state object must have been returned by a previous DMQMC calculation.
        The validity of this is not checked.  The system must also be unchanged.

.. _dmqmc_table:

dmqmc options
-------------

``replica_tricks``
    type: boolean.

    Optional.  Default: false.

    Perform replica simulations (i.e. evolve two independent DMQMC simulations
    concurrently) if true. This allows calculation of unbiased estimators that are
    quadratic in the density matrix.
``fermi_temperature``
    type: boolean.

    Optional.  Default: false.

    Rescale tau so that the simulation runs in timesteps of :math:`\Delta\tau / T_F` where :math:`T_F`
    is the Fermi temperature. This is so results are at dimensionless inverse temperatures of :math:`\Theta^{-1}
    =T_F/T`. This option is only valid for systems with a well defined Fermi energy.
``all_sym_sectors``
    type: boolean.

    Optional.  Default: false.

    Sample states with all symmetries of the system instead of just those which conserve
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

    Optional.  Default: none.

    Specify factors used to alter the spawning probabilities in the DMQMC importance
    sampling procedure. See PRB, 89, 245124 (2014) for an explanation, in particular
    section IV and appendix B.

    The length of the vector of floats should be equal to the maximum number of
    excitations from any determinant in the space. For a chemical system with :math:`N`
    electrons and more than :math:`2N` spin orbitals, this would be equal to
    :math:`N`. For a Heisenberg model with :math:`N` spins in the :math:`M_s=0` spin
    sector, this should be equal to :math:`N/2` (each pair of opposite spins flipped is
    one excitation).
``vary_weights``
    type: integer.

    Optional.  Default: 0

    The number of iterations over which to introduce the weights in the importance
    sampling scheme (see PRB, 89, 245124 (2014)). If not set then the full weights
    will be used from the first iteration. Otherwise, the weights will be increased
    by a factor of :math:`(W_{\gamma})^{\beta/\beta_{target}}` each iteration, where
    :math:`W_{\gamma}` is the final weight of excitation level :math:`\gamma` and
    :math:`\beta_{target}` is the beta value to vary the weights until (equal to
    the value specified by this option, multiplied by the time step size).
``find_weights``
    type: boolean.

    Optional.  Default: false.

    Run a simulation to attempt to find appropriate weights for use in the DMQMC
    importance sampling procedure. This algorithm will attempt to find weights such
    that the population of psips is evenly distributed among the various excitation
    levels when the ground state is reached (at large beta values). The algorithm
    should be run for several beta loops until the weights settle down to a roughly
    constant value.

    The weights are output at the end of each beta loop.

    This option should be used together with the **find_weights_start** option,
    which is used to specify at which iteration the ground state is reached
    and therefore when averaging of the excitation distribution begins.

    This option cannot be used together with the **excit_dist** option. The
    **find_weights** option averages the excitation distribution in the ground
    state, whereas the **excit_dist** option accumulates and prints out the
    excitation distribution at every report loop.

    .. warning::
    
        This feature is found to be unsuccessful for some larger lattices (for example,
        6x6x6, for the Heisenberg model). The weights output should be checked. Increasing
        the number of psips used may improve the weights calculated.
``find_weights_start``
    type: integer.

    Optional.  Default: 0.

    The iteration number at which averaging of the excitation distribution begins,
    when using the **find_weights** option.
``symmetrize``
    type: boolean.

    Optional.  Default: false.

    Explicitly symmetrize the density matrix, thus only sampling one triangle of the
    matrix.  This can yield significant improvements in stochastic error in some cases.
``initiator_level``
    type: integer.

    Optional.  Default: -1.

    Set all density matrix elements at excitation level **initiator_level** and
    below to be initiator determinants. An **initiator_level** of -1 indicates
    that no preferential treatment is given to density matrix elements and the
    usual initiator approximation is imposed, 0 indicates that the diagonal
    elements are initiators, etc.

    This is experimental and the user should identity when convergence has been
    reached.

.. _ipdmqmc_table:

ipdmqmc options
---------------

``target_beta``
    type: float.

    Optional.  Default: 1.0.

    The inverse temperature to propagate the density matrix to.
    If fermi_temperature is set to True then target_beta is interpreted as the inverse reduced temperature
    :math:`\tilde{\beta} = 1/\Theta = T_F/T`, where :math:`T_F` is the Fermi temperature. Otherwise target_beta is taken
    to be in atomic units.
``initial_matrix``
    type: string.

    Optional.  Default: 'hartree_fock'.

    Possible values: 'free_electron', 'hartree_fock'.

    Initialisation of the density matrix at :math:`\tau=0`.  'free_electron' samples the
    free electron density matrix, i.e. :math:`\hat{\rho} = \sum_i e^{-\beta \sum_j \varepsilon_j
    \hat{n}_j} |D_i\rangle\langle D_i|`, where :math:`\varepsilon_j` is the single-particle eigenvalue
    and :math:`\hat{n}_j` the corresponding number operator.  'hartree_fock' samples
    a 'Hartree--Fock' density matrix defined by :math:`\hat{\rho} = \sum e^{-\beta H_{ii}} |D_i\rangle\langle D_i|`,
    where :math:`H_{ii} = \langle D_i|\hat{H}|D_i\rangle`.

    It is normally best to use the hartree-fock option as this removes cloning/death on the diagonal if the shift
    is fixed at zero. This requires slightly more work when also using the grand_canonical_initialisation, but this
    is negligeable.

``grand_canonical_initialisation``
    type: boolean.

    Optional.  Default: false.

    Use the grand canonical partition function to initialise the psip distribution.
    The default behaviour will randomly distribute particles among the determinants
    requiring a non-zero value of metropolis_attempts to be set for the correct
    distribution to be reached.

``metropolis_attempts``
    type: integer.

    Optional.  Default: 0.

    Number of Metropolis moves to perform (per particle) on the initial distribution.
    It is up to the user to determine if the desired distribution has been reached,
    i.e. by checking if results are independent of metropolis_attempts.

``symmetric``
    type: boolean.

    Optional. Default: false.

    Use symmetric version of ip-dmqmc where now :math:`\hat{f}(\tau) =
    e^{-\frac{1}{2}(\beta-\tau)\hat{H}^0}e^{-\tau\hat{H}}e^{-\frac{1}{2}(\beta-\tau)\hat{H}^0}`.

    .. warning::

        This feature is experimental and only tested for the 3D uniform electron
        gas.

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
    Only available for the Heisenberg model.
``staggered_magnetisation``
    type: boolean.

    Optional.  Default: false.

    Calculate the thermal expectation value of the staggered magnetisation operator.
    Only available for the Heisenberg model and with bipartite lattices.
``excit_dist``
    type: boolean.

    Optional.  Default: false.

    Calculate the fraction of psips at each excitation level, where the excitation level
    is the number of excitations separating the two states labelling a given density matrix
    element. This fraction is then output to the data table at each report loop, and so the
    temperature-dependent excitation distribution is printed out.

    This option should not be used with the **find_weights** option, which averages the
    excitation distribution within the ground state.
``correlation``
    type: 2D vector of integers.

    Optional.  Default: false.

    Calculate the spin-spin correlation function between the two specified lattice sites,
    :math:`i` and :math:`j`, which is defined as the thermal expectation value of:

    .. math::

    	\hat{C}_{ij} = \hat{S}_{xi}\hat{S}_{xj} + \hat{S}_{yi}\hat{S}_{yj} + \hat{S}_{zi}\hat{S}_{zj}.

    Only available for the Heisenberg model.
``potential_energy``
    type: boolean

    Optional. Default: false

    Evaluate the bare Coulomb energy. Only available for the UEG.
``kinetic_energy``
    type: boolean

    Optional. Default: false

    Evaluate the kinetic energy. Only available for the UEG.
``H0_energy``
    type: boolean

    Optional. Default: false

    Evaluate the thermal expectation value of the zeroth order Hamiltonian
    where :math:`\hat{H} = \hat{H}^0 + \hat{V}`. See **initial_matrix**
    option. Only available when using the ip-dmqmc algorithm.
``HI_energy``
    Evaluate the expectation value of the interaction picture Hamiltonian where

    .. math::

        \hat{H}_I(\frac{1}{2}(\beta-\tau)) =
            e^{\frac{1}{2}(\beta-\tau)\hat{H}^0}\hat{H}e^{-\frac{1}{2}(\beta-\tau)\hat{H}^0}.



.. _rdm_table:

rdm options
-----------

Note that the use of RDMs is currently only available with the Heisenberg model.

``rdms``
    type: table of 1D vectors.

    Required.

    Each vector corresponds to the subsystem of a reduced density matrix as a list of the
    basis function indices in the subsystem.  For example:

    .. code-block:: lua

        rdms = { { 1, 2 } }

    specifies one RDM containing basis functions with indices 1 and 2, and

    .. code-block:: lua

        rdms = { { 1, 2 }, { 3, 4} }

    specifies two RDMs, with the first containing basis functions with indices 1 and 2,
    and the second basis functions 3 and 4.

    Either ``instantaneous`` or ``ground_state`` must be enabled to set the desired mode of
    evaluating the RDM (but both options cannot be used together).
``instantaneous``
    type: boolean.

    Optional.  Default: false.

    Calculate the RDMs at each temperature based upon the instantaneous psip distribution.

    Cannot be used with the ground_state option (either ground_state or instantaneous RDMs
    can be calculated, but not both concurrently).
``ground_state``
    type: boolean.

    Optional.  Default: false.

    Accumulate the RDM once the ground state (as specified by ``ground_state_start``)
    is reached.  This has two limitations: only one RDM can be accumulated in
    a calculation and the subsystem should be at most half the size of the system (which
    is always sufficient for ground-state calculations).

    Cannot be used with the instantaneous option (either ground_state or instantaneous RDMs
    can be calculated, but not both concurrently).
``spawned_state_size``
    type: integer.

    Required if ``instantaneous`` is true.  Ignored otherwise.

    Maximum number of states (i.e. reduced density matrix elements) to store in the
    "spawned" list, which limits the number of unique RDM elements that each processor can
    set.  Should be a sizeable fraction of ``state_size`` (see :ref:`qmc_table`) and
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
    only if ``ground_state`` is set to true and, as such, should be set to an iteration
    (which is a measure of temperature) such that the system has reached the ground state.
``concurrence``
    type: boolean.

    Optional.  Default: false.

    Calculate the unnormalised concurrence and the trace of the reduced density matrix at
    the end of each beta loop.  The normalised concurrence can be calculated from this using
    the ``average_entropy.py`` script.

    Valid for ``ground_state`` only; temperature-dependent concurrence is not currently
    implemented.
``renyi2``
    type: boolean.

    Optional.  Default: false.

    Calculate the Renyi-2 entropy of each subsystem. More accurately, the quantity output
    to the data table is :math:`S^n_2 = \sum_{ij} (\rho^n_{ij})^2`, (which differs from the
    Renyi-2 entropy by a minus sign and a logarithm) where :math:`\rho^n` is the reduced
    density matrix of the :math:`n`-th subsystem. The temperature-dependent estimate of
    the Renyi-2 entropy can then be obtained using the ``finite_temp_analysis.py`` script.

    Valid for ``instantaneous`` only; ground-state Renyi-2 averaged over a single beta
    loop is not currently implemented.  Requires ``replica_tricks`` to be enabled in order
    to obtained unbiased estimates.
``von_neumann``
    type: boolean.

    Optional.  Default: false.

    Calculate the unnormalised von Neumann entropy and the trace of the reduced density
    matrix at the end of each beta loop.  The normalised von Neumann entropy can be
    calculated from this using the ``average_entropy.py`` script.

    Valid for ``ground_state`` only; temperature-dependent von Neumann entropy is not
    currently implemented.
``write``
    type: boolean.

    Optional.  Default: false.

    Print out the ground-state RDM to a file at the end of each beta loop.  The file
    contains the trace of the RDM in the first line followed by elements of the upper
    triangle of the RDM labelled by their index.

    Valid for ``ground_state`` only.

.. _fciqmc:

Full Configuration Interaction Quantum Monte Carlo
==================================================

.. code-block:: lua

    fciqmc {
        sys = system,
        qmc = { ... },
        fciqmc = { ... },
        semi_stoch = { ... },
        restart = { ... },
        reference = { ... },
        load_bal = { ... },
        logging = { ... },
        output = { ... },
        blocking = { ... },
        qmc_state = qmc_state,
    }

Returns:
    a qmc_state object.

``fciqmc`` performs a full configuration interaction quantum Monte Carlo (FCIQMC)
calculation [Booth09]_ on a system.

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
``fciqmc``
    type: lua table.

    Optional.

    Further options to control the FCIQMC algorithm.  See :ref:`fciqmc_table`.
``semi_stoch``
    type: lua table.

    Optional.

    Further options to control using a semi-stochastic projection of the Hamiltonian
    operator instead of a purely stochastic projection.  Note that some options in the
    ``semi_stoch`` table are required to be set if the table is given.  See
    :ref:`semi_stoch_table`.
``restart``
    type: lua table.

    Optional.

    Further options to control restarting the calculation from a previous calculation.
    See :ref:`restart_table`.
``reference``
    type: lua table.

    Optional.

    Further options to select the reference state used.  See :ref:`reference_table`.
``load_bal``
    type: lua table.

    Optional.

    Further options to improve the parallel load balancing of an FCIQMC simulation.  If
    present (even if empty) an advanced load-balancing algorithm is used
    [Malone16a]_.  See :ref:`load_bal_table` for more details.
``logging``
    type: lua table.

    Optional.

    Further options to enable various logging outputs from a FCIQMC simulation. Only
    available when compiled in debug mode. See :ref:`logging_table` for information
    on current options.
``output``
    type: lua_table.

    Optional.

    Further options to enable direction of calculation output to a different file.
    See :ref:`output_table` for more information.
``blocking``
    type: lua table.

    Optional.

    Further options to switch on and control blocking on the fly. See :ref:`blocking_table`.
``qmc_state``
    type: qmc_state object.

    Optional.

    Output of a previous calculation to resume.

    .. warning::

        The qmc_state object must have been returned by a previous FCIQMC calculation.
        The validity of this is not checked.  The system must also be unchanged.

    .. warning::

        This destroys the qmc_state object and so it cannot be re-used in subsequent
        QMC calculations.

.. _fciqmc_table:

fciqmc options
--------------

``select_reference_det``
    type: boolean or Lua table.

    Optional.  Default: false.

    If true or if a lua table is provided, attempt to automatically set the reference
    state to be the state with the greatest population.  A lua table can contain the
    following options and need only be provided in order to modify the defaults.

    .. note::

        Care should be take when analysing the projected estimator to ensure that
        all quantities averaged have the same reference state.

    .. warning::

        Excitation levels are relative to the reference state and hence this should
        **not** be used with a truncated CI calculation.

    ``update_every``
        type: integer

        Optional.  Default: 20.

        The number of report loops between attempts to update the reference state.
    ``pop_factor``
        type: float.
        
        Optional.  Default: 1.5.

        The factor of the reference population another state must have in order for the
        reference to be changed.  This helps prevent continually switching between states
        with similar or degenerate populations.

``non_blocking_comm``
    type: boolean.

    Optional.  Default: false.

    Use non-blocking MPI communications instead of blocking MPI communications.

    .. note::

        This is an experimental option and may or may not improve performance.  In
        particular, its efficiency is highly dependent upon architecture and MPI
        implementation.  For expert use only!

``load_balancing``
    type: boolean.

    Optional.  Default: false.

    Enable dynamic load balancing of determinants among processors. This will move
    determinants to try and keep the number of walkers on each processor roughly
    constant. See :ref:`load_bal_table` for more details.

``init_spin_inverse_reference_det``
    type: boolean.

    Optional.  Default: false.

    In addition to initialising the reference determinant with an initial
    population, initialise the spin-inversed determinant (if different) with
    the same population.  Overridden by a restart file.
``trial_function``
    type: string.

    Optional.  Default: 'single_basis'.

    Possible values: 'single_basis', 'neel_singlet' (Heisenberg model only).

    The trial function to use in the projected energy estimator.  'single_basis'
    uses the single reference state as the trial function.  'neel_singlet' uses the Neel
    singlet state, :math:`|NS \rangle = \sum_{i} a_i |D_i \rangle`, where the amplitudes
    :math:`a_i` are defined in K. Runge, Phys. Rev. B 45, 7229 (1992).

    Using a multi-reference trial function can substantially reduce stochastic noise.
 
``guiding_function``
    type: string.

    Optional.  Default: 'none'.

    Possible values: 'none', 'neel_singlet' (Heisenberg model only).

    The importance sampling transformation to apply to the Hamiltonian.

    'neel_singlet' uses the Neel singlet state (K. Runge, Phys. Rev. B 45, 7229 (1992))
    to transform the Hamiltonian such that the matrix elements, :math:`H_{ij}`, are
    replaced with :math:`a_i H_{ij} / a_j`. Using 'neel_singlet' automatically sets
    ``trial_function`` to 'neel_singlet'.
``replica_tricks``
    type: boolean.

    Optional.  Default: false.

    Perform replica simulations (i.e. evolve two independent FCIQMC simulations
    concurrently).

``density_matrices``
    type: boolean.

    Optional.  Default: false.

    Samples the 2-RDM and reports the trace and numerator of the energy estimate.
    Prints the Hermitian part of the 2-RDM once the calculation has completed.

    Requires ``replica_tricks`` to be enabled.  
    
    .. note::

        Both replicas must be in variable shift mode in order for the 2-RDM 
        to be calculated. 

``density_matrix_file``
    type: string.

    Optional. Default: RDM

    The name of the output file that the final 2-RDM is printed to.
    The (normalised, Hermitian part of) the 2-RDM is given in physical notation.

``density_matrix_report``
    type: int.

    Optional. Default: 3000

    The minimum report cycle the 2-RDM statistics are collected after.

.. _load_bal_table:

load_bal options
----------------

The default values are usually sufficient if load balancing is enabled.  It is highly
recommended to only attempt to improve load balancing for large calculations and once the
population has been stabilised by the shift.  It may be easiest to do this by monitoring
a calculation carefully until this condition is reached, producing a restart file and then
running a production calculation with load balancing enabled.

``nslots``
    type: integer.

    Optional.  Default: 20.

    The average number of slots per processor used to distribute the list of occupied
    states via a hashing of the states.  A large value will affect performance but could
    potentially result in a better distribution of walkers.
``min_pop``
    type: integer.

    Optional.  Default: 1000.

    The minimum total population required before load balancing is attempted.  This is
    a system dependent value and, in order to maximise performance improvements, should be
    set such that the population is roughly stable.
``target``
    type: float.

    Optional.  Default: 0.05.

    Desired imbalance (as a percentage of the average population per processor) between
    the most and least populated processors.  Note that the workload on a processor is not
    entirely determined by its population and that, due to the algorithms used, an
    arbitrary small population imbalance is not usually possible.
``max_attempts``
    type: integer.

    Optional.  Default: 2.

    The number of attempts to make to improve load balancing.  Often multiple attempts can
    improve the balancing but each attempt may be non-negligible and there are usually
    diminishing returns.
``write``
    type: boolean.

    Optional.  Default: false.

    Write out the population of the most and least heavily populated processor
    before and after load balancing is carried out. Also print out the
    minimum slot population on the most populated processor which will
    indicate if load balancing is possible.

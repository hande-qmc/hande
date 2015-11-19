.. _qmc_table:

qmc options
===========

The following options in the ``qmc`` table are common to the FCIQMC, CCMC and DMQMC
algorithms and control the core settings in the algorithms.

``tau``
    type: float.

    Required.

    The timestep to use.

    A small timestep causes the particles sampling the wavefunction/matrix to evolve very
    slowly.  Too large a timestep, on the other hand, leads to a rapid particle growth
    which takes a long time to stabilise, even once the shift begins to vary, and coarse
    population dynamics.
``init_pop``
    type: float.

    Required.  Overridden if the calculations is initialised from a restart file.

    Set the initial population on the reference determinant.  For DMQMC calculations this
    option sets the number of psips which will be randomly distributed along the diagonal
    at the start of each beta loop.
``mc_cycles``
    type: integer.

    Required.

    Number of Monte Carlo cycles to perform per "report loop".
``nreports``
    type: integer.

    Required.

    Number of "report loops" to perform.  Each report loop consists of ``mc_cycles``
    cycles of the QMC algorithm followed by updating the shift (if appropriate) 
    and output of information on the current state of the particle populations, including
    terms in the energy estimators.
``state_size``
    type: integer.

    Required.

    Maximum number of states (i.e. determinants, excitors or density matrix elements) to
    store in the "main" list, which holds the number of particles on the state and related
    information such as the diagonal Hamiltonian matrix element.  The number of elements
    that can be stored usually should be of the same order as the target population.

    If negative, then the absolute value is used as the maximum amount of memory in MB to
    use for this information.

    .. note::

        This is a **per processor** quantity.  It is usually safe to assume that each
        processor has approximately the same number of states.

``spawned_state_size``
    type: integer.

    Required.

    Maximum number of states (i.e. determinants, excitors or density matrix elements) to
    store in the "spawned" list, i.e. the maximum number of states which can be spawned onto
    at a given timestep.  The amount of memory required for this is usually a small
    fraction of that required for ``state_size``, unless ``real_amplitudes`` is in use,
    in which case this should be a sizeable fraction (or potentially even greater than the
    memory for ``state_size``, if load balancing of states across processors is poor).
    The amount of memory required is also dependent on the value of ``tau``.

    If negative, then the absolute value is used as the maximum amount of memory in MB to
    use for this information.

    .. note::

        This is a **per processor** quantity.  It is recommended that a short trial
        calculation is run and the spawning rate for the desired timestep examined in
        order to estimate a reasonable value for ``spawned_state_size``.

``rng_seed``
    type: integer.

    Optional.  Default: generate a seed from a hash of the time and calculation UUID.

    The seed used to initialise the random number generator.
``target_population``
    type: float.

    Optional.  Default: none.

    Set the target number of particles to be reached before the shift is allowed to vary.
    This is only checked at the end of each report loop.  Once the ``target_population`` is reached, the shift is varied according to 

    .. math::

        S(t) = S(t-A\tau) - \frac{\xi}{A\tau} log\left( \frac{N_p(t)} {N_p(t-A\tau)} \right)

    where :math:`S` is the shift, :math:`t` the current imaginary time, :math:`\tau` the
    timestep, :math:`A` ``mc_cycles``, :math:`\xi` ``shift_damping``, and :math:`N_p` the
    number of particles.
``real_amplitudes``
    type: boolean.

    Optional.  Default: false.

    Allow amplitudes to take non-integer weights.  This will often significantly reduce
    the stochastic noise in the Monte Carlo estimates.

    Automatically enabled if semi-stochastic is used.

    .. note::

        Real amplitudes are handled using fixed precision and so numbers which can not be
        exactly represented are stochastically rounded to values that can be stored.

        The preprocessor option POP_SIZE=32 (default) uses 32-bit integers to store the
        amplitudes and stores amplitudes to within a precision/resolution of
        :math:`2^{-11}` and to a maximum absolute population of :math:`2^{20}`.

        Consider using the preprocessor option POP_SIZE=64 to allow a greater range of
        amplitudes to be encoded (precision of :math:`2^{-31}` and maximum absolute
        population of :math:`2^{32}` at the cost of doubling the memory required to store
        the amplitudes.

        By default uses integer weights, i.e. with the minimum resolution of 1.

``real_amplitude_force_32``
    type: boolean.

    Optional.  Default: false.

    Force the precision of the real amplitudes to that used for POP_SIZE=32 irrespective
    of the actual POP_SIZE compile-time parameter.

    .. note::

        The main use-case for this is reproducing results produced by binaries compiled
        using POP_SIZE=32 with binaries compiled using POP_SIZE=64; it is not intended for
        use in production calculations.

``spawn_cutoff``
    type: float.

    Optional.  Default: 0.01 if ``real_amplitudes`` is used, 0 otherwise.

    The minimum absolute value for the amplitude of a spawning event. If a spawning event
    with a smaller amplitude occurs then its amplitude will probabilistically be rounded
    up to the cutoff or down to zero in an unbiased manner.  A spawning event with an
    amplitude above the cutoff is stochastically rounded such that it can be stored in a
    fixed precision value.  If ``real_amplitudes`` is not in use, the fixed precision
    corresponds to unit values.

    Only relevant when using ``real_amplitudes``.
``excit_gen``
    type: string

    Optional.

    Possible values: 'renorm', 'no_renorm'.

    ============  =================     =========
    System        Implemented           Default
    ============  =================     =========
    chung_landau  renorm, no_renorm     renorm
    heisenberg    renorm, no_renorm     renorm
    hubbard_k     renorm, no_renorm     renorm
    hubbard_real  renorm, no_renorm     renorm
    ueg           no_renorm             no_renorm
    ringium       no_renorm             no_renorm
    read_in       renorm, no_renorm     renorm
    ============  =================     =========

    The type of excitation generator to use.  Note that not all types are implemented for
    all systems, usually because a specific type is not suitable for (large) production
    calculations or not feasible or useful.

    The 'renorm' generator requires an orbitals to be selected such that a valid
    excitation is possible, e.g. for a double excitation :math:`(i,j)\rightarrow(a,b)`,
    the combination :math:`i,j,a` is only selected if there exists at least one unoccupied
    orbital for :math:`b` which conserves any symmetry and spin quantum numbers.  This is
    efficient in terms of generating allowed excitations but involves an expensive
    renormalisation step.  The 'no_renorm' generator lifts this restriction at the cost of
    generating (and subsequently rejecting) such excitations; the excitation generation is
    consequently much faster.  In general, 'renorm' is a good choice for small basis sets
    and 'no_renorm' is a good choice for large basis sets, especially with a small number
    of electrons (such that forbidden excitations are rarely generated).

``pattempt_single``
    type: float.

    Optional.  Default: use the fraction of symmetry-allowed excitations from the
    reference determinant that correspond to single excitations.

    The probability of generating a single excitation.
``pattempt_double``
    type: float.

    Optional.  Default: use the fraction of symmetry-allowed excitations from the
    reference determinant that correspond to double excitations.

    The probability of generating a double excitation.
``initial_shift``
    type: float.

    Optional.  Default: 0.0.

    The initial value of the shift.
``shift_damping``
    type: float.

    Optional.  Default: 0.05.

    The shift damping factor, :math:`\xi`.
``vary_shift_from``
    type: float or string.

    Optional.  Default: ``initial_shift``.

    Specify a value to set the shift to when ``target_population`` is reached.  If the
    string 'proje' is specified then the instantaneous projected energy is used.  By
    instantly setting the shift to a value closer to the correlation energy, the total
    population can be stabilised substantially faster.

   There is no guarantee that the instantaneous projected energy is a good
   estimate of the ground state (particularly in the real-space formulation of
   the Hubbard model), but it is likely to be closer to it than the default
   shift value of 0.

``initiator``
    type: boolean.

    Optional.  Default: false.

    Enable the initiator approximation [Cleland10]_, in which spawned particles are only kept if they
    are created onto states which already have a non-zero population, or were produced by
    states which are already highly occupied (see ``initiator_threshold``), or multiple
    spawning events onto a previously unoccupied state occurred in the same timestep.

    .. note::

        The initiator approximation should be considered experimental for CCMC and DMQMC (see
        ``initiator_level`` option for DMQMC).

    .. warning::

        The initiator approximation is non-variational (due to the non-variational
        energy estimator used) and the error should be carefully converged by
        performing repeated calculations with increasing ``target_population`` values.

``initiator_threshold``
    type: float.

    Optional.  Default: 3.0.

    Set the (absolute) population above which a state is considered to be an initiator
    state.  A value of 0 is equivalent to disabling the initiator approximation.
``tau_search``
    type: boolean.

    Optional.  Default: false.  Not currently implemented in DMQMC.

    Update the timestep, ``tau``, automatically if by scaling it by 0.95 if a bloom event
    is detected.  A bloom event is defined as one which spawns more than three particles
    in a single spawning event in FCIQMC and one which spawns more than 5% of the total
    current population in a single spawning event in CCMC.

    .. note::

        Experimental option.  Feedback on required flexibility or alternative approaches
        is most welcome.

``use_mpi_barriers``
    type: boolean.

    Optional.  Default: false.

    Perform MPI_Barrier calls before the main MPI communication calls (both
    for communication of the spawned list, and any semi-stochastic
    communication). These are timed, and the total time spent in these calls
    is reported at the end of a simulation.  This is useful for assessing
    issues in load balancing, as it will allow you to see when certain
    processors take longer to perform their work than others. This is turned
    off by default because such calls may have an initialisation time which
    scales badly to many processors.

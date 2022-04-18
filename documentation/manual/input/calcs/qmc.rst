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

    Required unless the calculations is initialised from a restart file or qmc_state.

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

    Required unless qmc_state is given.

    Maximum number of states (i.e. determinants, excitors or density matrix elements) to
    store in the "main" list, which holds the number of particles on the state and related
    information such as the diagonal Hamiltonian matrix element.  The number of elements
    that can be stored usually should be of the same order as the target population.

    If negative, then the absolute value is used as the maximum amount of memory in MB to
    use for this information.

    Ignored if qmc_state is given.

    .. note::

        This is a **per processor** quantity.  It is usually safe to assume that each
        processor has approximately the same number of states.

``spawned_state_size``
    type: integer.

    Required unless qmc_state is given.

    Maximum number of states (i.e. determinants, excitors or density matrix elements) to
    store in the "spawned" list, i.e. the maximum number of states which can be spawned onto
    at a given timestep.  The amount of memory required for this is usually a small
    fraction of that required for ``state_size``, unless ``real_amplitudes`` is in use,
    in which case this should be a sizeable fraction (or potentially even greater than the
    memory for ``state_size``, if load balancing of states across processors is poor).
    The amount of memory required is also dependent on the value of ``tau``.

    If negative, then the absolute value is used as the maximum amount of memory in MB to
    use for this information.

    Ignored if qmc_state is given.

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
``reference_target``
    type: float.

    Optional.  Default: none.

    Set a target reference population to be reached before the shift is allowed to vary.
    Cannot be used in conjunction with ``target_population``.
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

    Optional. Default: system dependent.

    Possible values are system dependent (alternative, deprecated names in bracket):

    ============  ==========================  =========
    System        Implemented                 Default
    ============  ==========================  =========
    chung_landau  renorm, no_renorm           renorm
    heisenberg    renorm, no_renorm           renorm
    hubbard_k     renorm, no_renorm           renorm
    hubbard_real  renorm, no_renorm           renorm
    read_in       renorm, no_renorm,          renorm
                  renorm_spin,
                  no_renorm_spin,
                  heat_bath,
                  heat_bath_uniform_singles
                  (heat_bath_uniform),
                  heat_bath_exact_singles
                  (heat_bath_single),
                  uniform_power_pitzer
                  (power_pitzer_orderM),
                  heat_bath_power_pitzer
                  (power_pitzer_orderM_ij),
                  heat_bath_power_pitzer_ref
                  (power_pitzer_orderN),
                  uniform_cauchy_schwarz
                  (cauchy_schwarz_orderM),
                  heat_bath_cauchy_schwarz
                  (cauchy_schwarz_orderM_ij)
    ringium       no_renorm                   no_renorm
    ueg           no_renorm,                  no_renorm
                  power_pitzer
    ============  ==========================  =========

    The type of excitation generator to use.  Note that not all types are implemented for
    all systems, usually because a specific type is not suitable for (large) production
    calculations or not feasible or useful.

    The 'renorm' generator requires orbitals to be selected such that a valid
    excitation is possible, e.g. for a double excitation :math:`(i,j)\rightarrow(a,b)`,
    the combination :math:`i,j,a` is only selected if there exists at least one unoccupied
    orbital for :math:`b` which conserves any symmetry and spin quantum numbers.  This is
    efficient in terms of generating allowed excitations but involves an expensive
    renormalisation step.  The 'no_renorm' generator lifts this restriction at the cost of
    generating (and subsequently rejecting) such excitations; the excitation generation is
    consequently much faster.  In general, 'renorm' is a good choice for small basis sets
    and 'no_renorm' is a good choice for large basis sets, especially with a small number
    of electrons (such that forbidden excitations are rarely generated).
    'renorm_spin' and 'no_renorm_spin' are very similar to 'renorm' and 'no_renorm'
    respectively but when selecting :math:`i` and :math:`j`, they first decide with
    probability ``pattempt_parallel`` whether :math:`i` and :math:`j` should have
    parallel spins or not. The idea is by Alavi and co-workers, see [Booth09]_ and [Booth14]_
    for example for more details on these excitation generators.

    Note that the implementations of the weighted excitation generators here are all
    described in [Neufeld19]_.

    The 'heat_bath' excitation generator is very similar to the "original" heat bath
    excitation generator described by Holmes et al. [Holmes16]_. :math:`i,j,a,b` are chosen
    with weighted, precalculated probabilities that aim to make :math:`|H_{ij}|/p_\mathrm{gen}` as constant
    as possible. The difference to Holmes et al. is that we never do a single and a double
    excitation at the same time. When Holmes et al. decide to do both, we do a single
    excitation with probability of 0.5 and a double with 0.5. The 'heat_bath' excitation
    generator can have a bias if for a valid excitation :math:`i` going to :math:`a`,
    there might be no occupied :math:`j` that lets us select :math:`ija`. See Holmes et al.
    for details. We check for the bias in the beginning of a calculation and stop it if
    necessary.
    The Cauchy-Schwarz ([Smartunpub]_, described in [Blunt17]_)
    and Power-Pitzer excitation generators use approximate upper bounds
    for these weights. A version of Cauchy-Schwarz excitation generators is described in [Schwarz]_
    but the weights used here and the implementation differ.
    Here, Cauchy-Scharz uses Coulomb integrals and Power-Pitzer uses
    exchange integrals to approximate weights.
    'heat_bath_uniform_singles' is very similar to 'heat_bath' but samples single excitations
    uniformly (mentioned by Holmes et al.) and 'heat_bath_exact_singles' is also very similar
    but samples single excitations with the correct weighting (following a
    recommendation by Pablo Lopez Rios). 'heat_bath_uniform_singles' and 'heat_bath_exact_singles' do
    not have this potential bias that 'heat_bath' can have.

    Some of the Power-Pitzer excitation generators use elements of the heat-bath excitation
    generators ([Holmes16]_) and their approximations for selecting :math:`a` and :math:`b`
    are inspired by the Cauchy-Schwarz excitation generators by Alavi and co-workers
    [Smartunpub]_. See more details on all these weighted excitations generator in Ref. [Neufeld19]_.

    The 'power_pitzer' excitation generator generates double excitations using a Power-Pitzer
    [Power74]_ upper bound for the value of the Hamiltonian matrix element, 
    :math:`|\langle ij|ab\rangle|^2 => \langle ia|ai\rangle\langle jb|bj\rangle`
    (:math:`|\langle ij|ab\rangle|^2 => \langle ia|ia\rangle\langle jb|jb\rangle` for
    Cauchy-Schwarz excitation generators).
    This involves some precalcalated weights and alias tables, but should reduce both noise
    and shoulder heights. The weights to select a certain excitation are calculated for
    the reference in the beginning of the QMC calculation. Each time the excitation
    generator is called, the weights are mapped from the reference to the actual 
    determinant we attempt a spawn from. Only available for the UEG and read_in systems.
    The time spent in this excitation generator scales as :math:`\mathcal{O}(N)`, where
    :math:`N` is the number of electrons and the memory requirements are :math:`\mathcal{O}(N M)`,
    where :math:`M` is the number of basis functions.  Single excitations are done uniformly.

    'uniform_power_pitzer' uses a more refined upper bound for the Hamiltonian matrix
    elements, where the weights for selecting an excitation are calculated each time the
    excitation is called for the actual determinant we are spawning from. This requires
    :math:`\mathcal{O}(M)` time cost for each particle being spawned from. The 
    memory requirements are of :math:`\mathcal{O}(M)`. 'heat_bath_power_pitzer'
    is similar to 'uniform_power_pitzer' but samples selects :math:`i` and :math:`j`
    similarly to the heat bath excitation generators. The memory cost is
    :math:`\mathcal{O}(M^2)`.
    'uniform_cauchy_schwarz' is similar to 'uniform_power_pitzer' and 'heat_bath_cauchy_schwarz'
    is similar to 'heat_bath_power_pitzer', the distinction being the types of weights used
    to select :math:`ab`.

    The 'heat_bath_power_pitzer_ref' excitation generator [Neufeld19]_ uses precalculated weights and unlike
    'uniform_power_pitzer', it also samples :math:`i` and :math:`j` with weighted probabilities.
    It also samples single excitations in a weighted manner. Its memory cost is
    :math:`\mathcal{O}(M^2)`.
    This excitation generator can be useful in single-referenced systems when doing
    CCMC especially where the basis set size gets too big for 'heat_bath_power_pitzer' and
    'heat_bath_uniform_singles'. The computational scaling is also more favourable than
    with 'heat_bath_power_pitzer'.

    In the case of the UEG, the 'power_pitzer' excitation generator pre-calculates
    Power-Pitzer like weights for the selecting of orbital :math:`a`. :math:`i` and
    :math:`j` are selected like the 'no_renorm' UEG excitation generator.  If :math:`a` is
    occupied, the excitation is forbidden.

    .. note::
        Our current advice for selecting an excitation generator to use with read_in systems [Neufeld19]_:
        First consider the 'heat_bath' excitation generator. A bias test will be run at the beginning of
        the calculation then. If the bias test fails, try 'heat_bath_uniform_singles'.
        If 'heat_bath' and/or 'heat_bath_uniform_singles' fail due to memory constraints,
        try 'heat_bath_power_pitzer_ref'. Note that only 'heat_bath' requires a bias test.
    
    .. note::
        The Cauchy-Schwarz excitation generators are not implemented for complex read_in systems.

    .. note::
        Currently only the no_renorm and renorm excitation generators are available in
        DMQMC.

``power_pitzer_min_weight``
    type: float.

    Optional. Default: 0.01.

    Only used in 'power_pitzer_orderN' excitation generator or in 'read_in' systems if
    the 'power_pitzer' excitation generator is used.
    This number (approximately) sets the minimum value of
    weight(orbital to excite to)/(total weights times number of orbitals to excite to).
    The aim of this is to reduce the number of spawns with larger :math:`|H_{ij}|/p_\mathrm{gen}`
    which can happen if orbital connections with small values of :math:`p_\mathrm{gen}` are mapped to
    orbital connections with large values of :math:`|H_{ij}|`.

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

    .. note::
        If ``pattempt_single`` and ``pattempt_double`` do not sum to 1, we renormalize them.
``pattempt_update``
    type: boolean.

    Optional. Default: False.
    
    If true, then ``pattempt_single`` is varied during the run
    to attempt to align the means of :math:`|H_{ij}|/p_\mathrm{gen}` for single and double excitations.
    Mentioned in [Holmes16]_.
    Update of pattempt_single only happens if shift has not started varying yet.
    Not applicable to "original" heat bath algorithm excitation generator (excit_gen="heat_bath").
    When restarting a calculation, if ``pattempt_update`` is set to true and both ``pattempt_single``
    and ``pattempt_double``
    are specified by the user, previous update information is lost and the update (provided
    shift has not started varying yet) starts from scratch (the information to update
    ``pattempt_single`` from previous runs is lost).
    If ``pattempt_single`` or ``pattempt_double`` are in danger of getting too small, they will
    be set to 1/the number of allowed spawn attempts needed before they are updated again
    which is 10000 currently. A warning will be printed "WARNING: min. pattempt_single/double!" if
    that is the case. Do make sure that before accepting a final ``pattempt_single`` or
    ``pattempt_double``, this warning will have not been printed for a while.

    .. note::
        Currently not available in DMQMC.

    .. note::
        By the way we set the minimum values for ``pattempt_single`` and ``pattempt_double``, the
        minimum value for these is 0.0001. If that is too high, consider setting them manually by
        specifying both (only one is not sufficient) in the input file.
``pattempt_zero_accum_data``
    type: boolean

    Optional. Default: False.

    If true and restarting a calculation, accumulated data needed to update ``pattempt_single``
    and ``pattempt_double`` is reset (set to zero, overflow boolean is set to false).
    Only to be used together with ``pattempt_update``. Only to be used by experienced users.
``pattempt_parallel``
    type: float.

    Optional. Default: Estimate it using :math:`\frac{ \sum_{i_{\Vert}j_{\Vert}ab} |H_{ijab}| }{ \sum_{ijab} |H_{ijab}| }`, where :math:`i_{\Vert} j_{\Vert}` indicates :math:`i, j` are restricted to having parallel spins. 

    Probability that :math:`i, j` have parallel spins.
    Only to be used with ``excit_gen`` == 'no_renorm_spin' and 'renorm_spin'.

    Cannot be bigger than 1 and if negative, the default estimate is applied.
    It is recalculated in the beginning of each (restarted) calculation.
    
``initial_shift``
    type: float.

    Optional.  Default: 0.0.

    The initial value of the shift.
``shift_damping``
    type: float.

    Optional.  Default: 0.05.

    The shift damping factor, :math:`\xi`. This can be optimised using the
    ``auto_shift_damping`` keyword (see :ref:`blocking_table`).
    On restarting the final value in the previous calculation will replace
    the usual default value if ``shift_damping`` is not specified.

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

    Enable the initiator approximation (FCIQMC: [Cleland10]_; CCMC: [Spencer15]_; DMQMC:
    [Malone16]_) in which spawned particles are only kept if they are created onto states
    which already have a non-zero population, or were produced by states which are already
    highly occupied (see ``initiator_threshold``), or multiple spawning events onto
    a previously unoccupied state occurred in the same timestep.

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
``quadrature_initiator``
    type: logical.

    Optional. Default: true.

    The initiator approximation in a complex spaces could be applied in (at least) two different
    ways.
    If this parameter is true, the magnitude of the instantaneous complex coefficient at each site
    is used to determine initiator properties for both real and imaginary parents.

    If this parameter is false, the magnitude of the real and imaginary populations are compared
    separately and initiator flags for real and imaginary set individually.

    .. note::

        The comparative efficacy of these two approaches is currently under investigation.

``quasi_newton``
    type: boolean.

    Optional. Default: False.

    Turn on quasi-Newton steps.  Conventional FCIQMC and related methods take steps which are
    the equivalent of a scaled steepest-descent approach, which results in very long equilibration
    times, and requires smaller values of tau for stability.
    The quasi-Newton approach (partially) scales the steps according to the inverse difference in Fock energy to
    the reference determinant, reducing the contributions from very high-energy determinants.

    For more details see V. A. Neufeld, A. J. W. Thom, JCTC (2020), 16, 3, 1503-1510.

    .. note::

        Not currently available for DMQMC.
        Due to Fock value calculations, only supported for read_in systems and the 3D uniform electron gas.
        For semistochastic FCIQMC, determinants in the deterministic space are given weighting 1.

``quasi_newton_threshold``
    type: float.
    
    Optional. Default: Energy difference between LUMO and HOMO.

    Used when ``quasi_newton`` is true.
    The quasi-Newton approach (partially) scales the steps according to the inverse difference in Fock energy to
    the reference determinant (with Fock energy :math:`F_0`) for each determinant.  Any determinant with energy
    less than :math:`F_0 + \Delta_{\mathrm{QN}}`, where :math:`\Delta_{\mathrm{QN}}` is the value
    given to ``quasi_newton_threshold``, will have weighting :math:`v_{\mathrm{QN}}^{-1}`,
    where :math:`v_{\mathrm{QN}}` is the value given by ``quasi_newton_value``.
    The shift containing term in the death step are scaled by ``quasi_newton_pop_control`` instead. This
    makes sure that that term is scaled by a constant, independent of the determinant/excitor involved,
    so that the energy does not diverge with fluctuations around the true energy.

    For more details see V. A. Neufeld, A. J. W. Thom, JCTC (2020), 16, 3, 1503-1510.

``quasi_newton_value``
    type: float.

    Optional. Default: ``quasi_newton_threshold``.

    See ``quasi_newton_threshold``.

``quasi_newton_pop_control``
    type: float

    Set to 1 for original/non quasi-Newton propagation and otherwise for quasi-Newton,
    the default is 1/``quasi_newton_threshold``.

    See ``quasi_newton_threshold``.

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
``vary_shift``
    type: boolean.

    Optional.

    If present, overrides any value of ``vary_shift`` set by a previous calculation
    contained either in a restart file or a qmc_state object.  If set to true, the shift
    is set to ``initial_shift``.

    .. note::

        The shift will still be varied when ``target_population``, if set, is reached.

``shift_harmonic_forcing``
    type: float.

    Optional.  Default: 0.00

    If present, this sets the restoring force factor value in the harmonic population 
    control algorithm. This differs from the canonical two-step population control by
    an additional term based on the target population, as follows

    .. math::

        S(t) = S(t-A\tau) - \frac{\xi}{A\tau} log\left( \frac{N_p(t)} {N_p(t-A\tau)} \right)
            - \frac{\zeta}{A\tau} log\left( \frac{N_p(t)} {N_t} \right)

    where where :math:`S` is the shift, :math:`t` the current imaginary time, :math:`\tau` the
    timestep, :math:`A` ``mc_cycles``, :math:`\xi` ``shift_damping``, :math:`\zeta` 
    is the restoring force factor described here, :math:`N_p` the number of particles and
    :math:`N_t` is the target population. 

    For more details see M. Yang, E. Pahl and J. Brand, J. Chem. Phys. 153, 174103 (2020) 
    (DOI:10.1063/5.0023088). 

    .. note::
  
        The original population control algorithm is obtained if set equal to zero.
    
    .. note::

        When used, the shift will vary throughout the entire simulation, even if the 
        target population has not been reached. 
        
    .. note::

        The harmonic population control algorithm will not work with target populations
        less than or equal to zero. 

``shift_harmonic_crit_damp``
    type: boolean.

    Optional.  Default: false.

    If set to true, the value of ``shift_harmonic_forcing`` will be set to the square
    of ``shift_damping`` divided by four to obtain critical damping.  

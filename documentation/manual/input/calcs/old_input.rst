.. _old_input:

Old input options
=================

.. warning::

    The following options were used in an old (statement-based) input parser and are
    included here only whilst new documentation is being written.  Some names have changed
    (see comments in src/lua_hande_calc.f90 and examples in test_suite).  Behaviour is,
    however, identical.

Algorithm options
^^^^^^^^^^^^^^^^^

The following are modes which can be used on top of some of the calculation
types below. They are turned off by default.

**real_amplitudes**
    Allow walker amplitudes to have a non-zero fractional part.

    This will often significantly reduce the stochastic noise in the various
    Monte Carlo estimates. One should consider setting the pre-processor option
    POP_SIZE=64 when using this option as this allows a greater range of
    amplitudes to be encoded.

    This option is only implemented with the **fciqmc**, **ccmc** and **dmqmc**
    options currently.
**spawn_cutoff** *cutoff*
    Real.

    Default when using **real_amplitudes**: 0.01.
    Default otherwise: 0.0.

    Set the minimum absolute value for the amplitude of a spawning event. If a
    smaller spawn occurs then its amplitude will probabilistically be rounded up
    to *cutoff* or down to zero in an unbiased manner.

    This parameter is relevant when using the **real_amplitudes** option. When
    not using the **real_amplitudes** option, all spawning occurs in multiples
    of 1.
**semi_stoch_high_pop** *space_size*
    Perform a semi-stochastic calculation. The deterministic space is created
    by choosing the *space_size* most populated determinants in the simulation.
    If there are less than *space_size* determinants in the simulation then all
    determinants will be used in the deterministic space.

    If the **semi_stoch_iteration** option is used then this option will use
    the walker configuration at the specified iteration, else the deterministic
    space will be created using the determinants present before the start of
    the first iteration. Therefore, one should only use this option in
    conjuction with the **restart** option or with the **semi_stoch_iteration**
    option.

    This option is only implemented with the **fciqmc** method.

    If this option is used then the **real_amplitudes** option will be turned on
    automatically.
**semi_stoch_read**
    Perform a semi-stochastic calculation. The deterministic space is created
    by reading in determinants from an HDF5 file produced using the
    **write** option.
    
    The filename will be of the form SEMI.STOCH.x.H5, where x is the file id.
    By default, x will be the smallest available id of all existing files.
    However, a particular id can be chosen using the **read** option.
**semi_stoch_iteration** *iter*
    Turn the semi-stochastic algorithm on at iteration number *iter*.
**semi_stoch_shift_start** *iter*
    Turn the semi-stochastic algorithm on *iter* iterations after the shift
    starts to vary.
**semi_stoch_combine_annihil**
    Default: false.

    This option will allow the semi-stochastic method to be used without an
    extra MPI call per iteration. Instead, deterministic spawnings are added to
    the spawned list and communicated with all other spawnings. One may find
    large variations in the time to perform each iteration, depending on
    whether this option is used or not.

**estimate_hilbert_space** *ncycles*
    Integer.

    Estimate the size of the Hilbert space within the desired symmetry block of
    the Hamiltonian by performing *ncycles* cycles of a Monte Carlo algorithm.
    The overall spin must be set using **ms**.
    Appropriate spatial and momentum symmetries are taken into account.
    The symmetry block can be selected by specifying a reference determinant.
    When run on multiple processors, an estimate of the error in the size is produced.
    This is not available on a single processor, and the user is warned to test the
    value by changing seeds or number of cycles, as not all printed figures may be significant.

    For the real space formulation of the Hubbard model and the Heisenberg
    model, the exact size of the space (at least to the first 8 significant
    figures) is found by simple combinatorics.
**redistribute_restart** *nprocs*
    Redistribute a set of restart files to be used in a calculation parallelised over
    *nprocs* MPI ranks.

    .. warning::

        Currently this is not parallelised and must be done separately from any QMC
        calculations.

**estimate_canonical_kinetic_energy**

    Estimate the free-electron thermal kinetic energy in the canonical ensemble
    in all momentum symmetry sectors by performing nkinetic_cycles*init_pop cycles of a
    Monte Carlo algorithm. Also estimate <H>_0, which is a form of Hartree-Fock energy.
    Estimates for mean and variance are printed out every init_pop cycles.

Calculation options: FCIQMC options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following options are valid for FCIQMC calculations.

**mc_cycles** *mc_cycles*
    Integer.

    Number of Monte Carlo cycles to perform per "report loop".

    Note that *mc_cycles* is set to be 1 for the continuous time algorithm.
**nreports** *nreports*
    Integer.

    Number of "report loops" to perform.  Each report loop consists of 
    *mc_cycles* cycles of the FCIQMC algorithm followed by updating the shift
    and output of information on the current state of the walker populations, in
    particular the instantaneous energy estimators.

    If *nreports* is set to be a negative number, then the FCIQMC algorithm
    will effectively loop indefinitely (strictly speaking: *nreports* is set to
    the largest possible number that can be held in the standard integer type).
    In such cases calculations should be cleanly exited using the
    HANDE.COMM functionality.

    The total number of Monte Carlo cycles performed in an FCIQMC calculation
    is *nreports* x *mc_cycles*.
**seed** *seed*
    Integer.

    Default: random value based upon a hash of the time and (if available) the
    universally unique identifier (UUID) of the calculation.

    Set the seed used to initialise the dSFMT random number generator.
    In parallel the seed on each processor is *seed* + iproc, where iproc is
    the processor index (as supplied by MPI) and ranges from 0 to nprocs-1.
**tau** *tau*
    Real.

    Set the timestep to be used.  Each Monte Carlo cycle amounts to propagating
    the walker population by the *tau* in units of imaginary time.

    A small timestep causes the walker population to evolve very slowly.  Too
    large a timestep, on the other hand, leads to a rapid particle growth which
    takes a long time to stabilise, even once the shift begins to vary, and
    coarse population dynamics.
**tau_search**
    Update the **tau** automatically by scaling it by 0.95 if a bloom event is
    detected.  A bloom event is defined as one which spawns more than three
    particles in a single spawning event in FCIQMC and one which spawns more than 5% of
    the total current population in a single spawning event in CCMC.

    .. note::

        This is an experimental option and feedback on required flexibility or
        alternative approaches is most welcome.

        **tau_search** is currently ignored in DMQMC calculations.

**initial_shift** *initial_shift*
    Real.

    Default: 0.

    Set the value of the shift to use during the period before the shift is
    allowed to vary.  Positive values lead to faster growth in the number of
    walkers due to cloning.  Using too large a value can lead to poor sampling
    as large numbers of walkers reside on the same small number of determinants
    rather than diffusing appropriately through the determinant space.
**vary_shift_from** **proje** | *shift* 
    String or real.

    Default: off.

    Set the shift to be either the instantaneous projected energy or the value
    specified by *shift* when *varyshift_target* is reached.  Most calculations
    start with setting the shift to be 0; by instantly setting the shift to
    a value closer to the true ground state, the simulation can stabilise the
    total walker population substantially faster.

    Note that the last option out of **initial_shift** or **vary_shift_from**
    *shift* is used.  Only use both options if you know what you're doing.

    There is no guarantee that the instantaneous projected energy is a good
    estimate of the ground state (particularly in the real-space formulation of
    the Hubbard model), but it is likely to be closer to it than the default
    shift value of 0.
**varyshift_target** *varyshift_target*
    Long integer.

    Default: 10000.

    Set the target number of particles to be reached before the shift is
    allowed to vary.  This is only checked at the end of each report loop.
**shift_damping** *xi*
    Real.

    Default: 0.05.

    Once the *varyshift_target* has been reached, the shift is updated according to:

    .. math::

        S(\beta) = S(\beta-A\tau) - \frac{\xi}{A\tau} log\left( \frac{N_w(\tau)} {N_w(\beta-A\tau)} \right)

    where :math:`\beta` is the current imaginary time, :math:`A\tau` is the
    amount of imaginary time between shift updates, :math:`N_w` is the number of
    walkers at the given time and :math:`\xi` is a damping factor to prevent
    wild fluctuations in the population dynamics and can be set using the
    **shift_damping** keyword.
**reference_det** *electron_1 electron_2 ... electron_nel*
    Integer list.

    Default: Momentum-space formulation of the Hubbard model
    Uses the Hartree--Fock determinant (ie that formed from occupying the
    nalpha and nbeta spin-orbitals with the lowest kinetic energy); 

    Default: Real-space formulation of the Hubbard model
    Attempt to minimise the number of doubly-occupied sites.  
    Note that this is not guaranteed (especially in the
    real-space formulation) to give a reference determinant which is close to
    the ground state.  Further, the default ignores any value of
    the symmetry as defined by the **sym** input option.
    
    Default: Heisenberg model
    For ferromagnetic cases (J>0) the default will attempt to group the up
    spins together, which often will result in the best reference determinant.
    For antiferromagnetic cases, first it will attempt to choose sites
    which do not neighbour each other. Then, if more spins are required
    it will choose the remaining spins in order of site label.
    This will usually give a good reference determinant, but it is not guaranteed
    always. For bipartite lattices however, the antiferromagnetic determinant 
    chosen should be the best one possible.
    
    Set the reference determinant to occupy the specified spin-orbitals.
    The index of each spin-orbital is printed out in the basis functions
    section of the output.  This will be overridden by a restart file and
    in a simple_fciqmc calculation, where the determinant with the lowest
    energy is set to the reference determinant.
    
    For the Heisenberg model, the electron positions will actually represent the
    positions on the lattice of the up spins in the reference basis vector.
    (Note that the number of up spins is deduced from the ms value specified and the
    total number of sites).
**init_pop** *pop*
    Integer.

    Default: 10.

    Set the initial walker population on the reference determinant.  This will
    be overridden by a restart file.

    For DMQMC calculations this option sets the number of psips which will
    be randomly distributed along the diagonal at the start of each beta loop.
**cluster_multispawn_threshold** *thresh*
    real.
    
    Default: huge  (i.e. off).

    When selecting clusters the generations probabilities can vary over orders of
    magnitude.  If after having selected the cluster, the value of
    cluster%amplitude/cluster%pselect
    is greater than *thresh*, then the number of spawning attempts from that cluster,
    nspawn_attempts, will be set to the smallest number such that
    cluster%amplitude/(cluster%pselect*nspawn_attempts) is less than *thresh*.
    The overall effect will be to reduce population blooms which raise plateau heights.
    The lower this number is the slower a calculation will be, though a larger tau might
    be able to be used.
    To enable, set to a number such as 0.1.
    NB the probability that the spawning is successful is still also dependent on 
    tau*(the spawning matrix element)/(the probability of generating the spawning excitation),
    and so estimates of these might be able to be used to set sensible values of *thresh*.

**init_spin_inverse_reference_det**
    Default: false.

    In addition to initialsing the reference determinant with an initial
    population, initialise the spin-inversed determinant (if different) with
    the same population.  This will be overridden by a restart file.
**select_reference_det** [*N* [*pop_fac*]]
    Default: off, 20 and 1.5.

    This option is only available when using the *fciqmc* method.

    Set the reference determinant to be the determinant with the largest
    population every *N* cycles if that population is greater than the
    population on the current reference determinant by a factor larger than
    *pop_fac*.  *pop_fac* should be greater than 1 to avoid repeated switching
    between degenerate determinants.

    .. warning::

        Care must be taken with averaging quantities when using this option.
        In particular, one should only average the projected estimator over
        imaginary time during which the reference determinant is constant.

**walker_length** *walker_length* [**MB**]
    Integer.

    Size of walker array.  This is allocated at the start of the calculation
    and is used to store the population of walkers on determinants with
    a non-zero population and the associated energy of the determinant.

    If **MB** is specified, then the walker_length is given in terms of MB per
    core rather than number of elements per core in each array
    associated with the parent walkers.

    Care: this needs to be large enough to hold the number of unique
    determinants with a non-zero population of walkers in the simulation.  The
    code does not currently check whether this size is exceeded and so setting
    **walker_length** to be too small can lead to memory problems and
    segmentation faults.  For large calculations this should be substantial
    smaller than the full size of determinant space.

    Not valid for simple_fciqmc calculations, where the population of walkers
    on each determinant is stored.
**spawned_walker_length** *spawned_walker_length* [**MB**]
    Integer.

    Size of the spawned walker array.  This is allocated at the start of the
    calculation and is used to store the population of spawned walkers on child
    determinants.

    If **MB** is specified, then the spawned_walker_length is given in terms of
    MB per core rather than number of elements per core in each array
    associated with the spawned walkers.

    Care: this needs to be large enough to store all the particles which are spawned
    during a Monte Carlo cycle and so needs to be a reasonable fraction of the 
    targeted number of total number of walkers.  The code does not currently
    check whether this size is exceeded and so setting
    **spawned_walker_length** to be too small can lead to memory problems and
    segmentation faults.

    Not valid for simple_fciqmc calculations, where the population of spawned
    walkers on each determinant is stored.
**no_renorm**
    Default (uniform electron gas): On.

    Default (all other systems): Off.

    Generate (and then reject) excitations which involve exciting an electron
    into a spin-orbital which is already occupied.  Whilst this is wasteful, it
    avoids having to renormalise the excitation generation probabilities, which
    can be expensive for large systems.
**dump_restart** [**shift**]  [*id*]
    Optional integer.

    Write out information required for restarting an FCIQMC calculation to
    a file called HANDE.RS.x.py.H5, where x is *id* if *id* is given and y is 
    the processor rank. If x is not given, it is chosen to be the smallest 
    integer possible such that HANDE.RS.x.py.H5 does not exist in the
    calculation directory.

    If **shift** is specified, then the restart information is dumped out before
    the shift turns on. Both dump_restart and dump_restart shift may be specified
    in the input file but the optional *id* (if specified) for both must be different.

    Restarting a parallel run with a different number of processors is not 
    currently supported.

    Warning: these files can become very large, so care should be taken when
    not re-using the same filenames.
**dump_restart_every** *nreport*
    Integer.  Default: off.

    Write out a restart file every *nreport* report cycles.

    .. warning::

         Unless **dump_restart** is specified with a file id, this will create
         a new restart file every *nreport* report cycles.  The disk space used
         with this option can therefore be very large.  Small values of
         *nreport* should only be used for diagnostic purposes and not in
         production calculations on large systems.

         Furthermore, writing to (for instance) a network disk will degrade performance
         substantially.

**write** *id*
    Default: off.

    Write the determinants in any used semi-stochastic deterministic space to a
    file.

    The filename will have the form SEMI.STOCH.x.H5, where x is the file *id*.
**ascii_format_out**
    The default format for restart files is binary, as for the most part the files
    are meant purely for reading by Hubbard, and having the file in human-readable
    ASCII format is both wasteful of space and unnecessary. 

    If the **ascii_format_out** keyword is specified, however, this overrides the default
    and the restart file is written out in ASCII. Beware; these files can become
    very large.
**ascii_format_in**
    Similar behaviour to **ascii_format_out** except that this one specifies that the restart
    file to be read (specified with the **restart** keyword) is in non-standard ASCII format
    as opposed to binary format.
**ascii_format**
    An Alias for both **ascii_format_in** and **ascii_format_out**
**restart** [*id*]
    Optional integer.

    Restart an FCIQMC calculation using a previous restart file,
    HANDE.RS.x.py.H5, where x is a non-negative integer and y is the processor
    rank. If *id* is given, x is set to *id*; otherwise x is chosen to be the
    largest integer such that HANDE.RS.x.py.H5 exists and HANDE.RS.x+1.py.H5
    does not.

    The restart file does not contain system information such as the U and
    T parameter, lattice vectors, number of electrons or if the walker
    population were evolved using standard FCIQMC or initiator-FCIQMC. Thus it
    is important use the same system parameters when restarting a calculation.
    The consistency of the restart file with the input options supplied is not
    checked.
    
    Please note that the RNG is not stored in the restart file, so running two
    shorter calculations via the restart facility is not completely identical
    to running a single calculation for the same number of Monte Carlo cycles.

    Furthermore, the current implementation does not allow restart files
    produced with one value of DET_SIZE to be used with binaries produced with
    a different value of DET_SIZE.  However, this is not checked!
**uniform_combination**
    For the Heisenberg model only. If this keyword is specified then instead of using a
    single reference detereminant to calculate the projected energy, a linear combination
    of all basis functions with amplitudes 1 is used:

    .. math::

    	|\psi \rangle = \sum_{i} |D_i \rangle

    hence the estimator used is


    .. math::

        E_0 = \frac{ \langle \psi|H|\psi_0 \rangle }{ \langle \psi|\psi_0 \rangle }
            = \frac{ \sum_{i,j} \langle D_i|H|D_j \rangle c_j } { \sum_{i} c_i }
                  
    A unitary transformation will be applied to the Hamiltonian so that all the
    off-diagonal elements are multiplied by -1. This has the effect of making
    the transformed ground state have all positive components, and hence the above
    trial function has a large overlap with this transformed ground state.
    
    This can only be used for bipartite lattices.
**neel_singlet_estimator**
    For the Heisenberg model only. If this keyword is specified then instead of
    using a single reference detereminant to calculate the projected energy,
    the Neel singlet state is used. This is a state,
    :math:`|NS \rangle = \sum_{i} a_i |D_i \rangle`, where the amplitudes
    :math:`a_i` are defined in K. Runge, Phys. Rev. B 45, 7229 (1992). For
    further details, see the comments in the subroutine
    update_proj_energy_heisenberg_neel_singlet in heisenberg_estimator.F90.
    
    This can only be used for bipartite lattices.
**neel_singlet_guiding**
    For the Heisenberg model only. If this keyword is specified then the Neel
    singlet state is used as a guiding state for importance sampling. This
    means that the the matrix elements of the Hamiltonian, :math:`H_{ij}`, are
    replaced by new components

    .. math::
    
        H_{ij} \leftarrow (a_i H_{ij})/a_j
    
    where :math:`a_i` is a component of the Neel state, as specified above.
    
    When this guiding function is used, the Neel singlet must be used in the
    projected energy, so the neel_singlet_estimator option is automatically
    applied.

Calculation options: CCMC options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**move_freq** [x]
    Optional integer.  Default: 5.

    Excitors are allowed to move processors every 2^x iterations in order to
    allow all composite excitors to be correctly sampled.  Relevant only when
    performing CCMC calculations with multiple MPI processes.

**ccmc_full_nc**
    Default: off.

    The original CCMC algorithm involves randomly selected a cluster of arbitrary size
    consisting of any set of excitors and then making spawning attempts from it.
    The full non-composite algorithm is a simple modification in which all occupied
    non-composite clusters (i.e. those consisting of the reference or just a single
    excitor) are (deterministically) selected and composite clusters (involving two or
    more excitors) are randomly selected to make spawning attempts.  This has been shown
    to give substantially more stable dynamics and reduce the plateau height in
    several systems.

**ccmc_linked**
    Default: off

    The original CCMC algorithm solves the equations

    .. math::

        \langle D_m | \hat{H} - E | \psi_{CC} \rangle = 0.

    It is possible to instead sample the equivalent equations

    .. math::

        \langle D_m | e^{-\hat{T}} (\hat{H} - E) | \psi_{CC} \rangle = 0.

    Using the Hausdorff expansion of the Hamiltonian and the linked cluster theorem means 
    that the only clusters which contribute are those with at most four excitors and where 
    the exitation sampled from the Hamiltonian has an orbital in common with each excitor 
    in the cluster operator. Using this option can give substantial reductions in the 
    plateau height.

Calculation options: DMQMC options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In addition to the options for FCIQMC calculations, the following options are additional to the 
configuration of a Density Matrix Quantum Monte Carlo (DMQMC) calculation

Note: The DMQMC features have only been coded and tested for the Heisenberg model.

**beta_loops**
    Integer.

    Default: 100.

    Set the number of beta loops. This is the number of times that the complete range of beta values
    will be looped over before the simulation finishes.
**dmqmc_energy**
    Calculate the thermal expectation value of the Hamiltonian operator.

    This value will be calculated from the first iteration of each report loop.
**dmqmc_energy_squared**
    Calculate the thermal expectation value of the Hamiltonian squared operator.

    This value will be calculated from the first iteration of each report loop.
**dmqmc_staggered_magnetisation**
    Calculate the thermal expectation value of the staggered magnetisation operator.

    This value will be calculated from the first iteration of each report loop.

    This option is only available for bipartite lattices.
**dmqmc_correlation_function** *site_1* *site_2*
    Integers.

    Calculate the spin-spin correlation function between the two lattice sites *site_1* and
    *site_2*, defined as the thermal expectation value of the following operator:

    .. math::

    	\hat{C}_{ij} = S_{xi}S_{xj} + S_{yi}S_{yj} + S_{zi}S_{zj}.

    This value will be calculated from the first iteration of each report loop.

    Note: the correlation function can only be calculated for one pair of spins in a single simulation.
**dmqmc_full_renyi_2**
    Calculate the Renyi-2 entropy of the entire system.

    This option must only be used when the **replica_tricks** option is also used.

    The quantity output in the column 'Full S2' is the instantaneous estimate of
    :math:`\sum_{ij}\rho_{ij}^2`. The traces of the two replicas are in the columns
    named 'Trace' and 'Trace 2'. The finite_temp_analysis.py script in the tools
    directory can then be used to obtain a final temperature-dependent estimate of
    the Renyi-2 entropy from these quantities.
**truncation_level** *truncation_level*
    Integer.

    Consider only elements of the density matrix where the determinants differ
    by at most *truncation_level* excitations.
**half_density_matrix**
    Symmetrise the density matrix explicitly. This may slightly improve the efficiency
    of the algorithm in some situations.
**output_excitation_distribution**
    Output the fraction of psips on each excitation level.
**use_all_sym_sectors**
    Run a DMQMC calculation in all symmetry sectors simultaneously. Psips will be
    distributed across all symmetry sectors for the initial density matrix.
**use_all_spin_sectors**
    Run a DMQMC calculation in all spin symmetry sectors simultaneously. Psips will be
    distributed across all spin symmetry sectors for the initial density matrix.
**dmqmc_weighted_sampling** *number_weights* Integer.
                            *w_{01} w_{12} ... w_{n-1,n}* Real list.

    This option will allow a form of importance sampling to be applied to the DMQMC calculation.

    The values of :math:`w_{01}, \ldots, w_{n-1,n}` will define weights which alter the spawning probabilities
    between the various excitation levels. When attempting to spawn from an excitation level
    i to a different excitation level j, the spawning probability will be altered by a factor
    :math:`1/w_{ij}`. Also, :math:`w_{ji} = 1/w_{ij}`. This can be used to help keep psips near the diagonal elements
    and hence improve the quality of sampling when calculating estimators, which typically depend upon
    psips on the diagonal and first one or two excitation levels. This is particularly useful for larger
    lattices where typically no psips will reside on the diagonal elements when the ground state is
    reached.

    To account for the altered spawning probabilities, different weights are given to different
    psips when calculating estimators, such that the same mean values are estimated, but with an
    improved quality of sampling.

    The value *number_weights* must equal the number of weights which have been specified.
    The weights :math:`w_{01}, \ldots, w_{n-1,n}` should be input on the lines directly after
    **dmqmc_weighted_sampling**, and can be input over as many lines as required.
**dmqmc_vary_weights** *N*
    Integer.

    If this option is specified then the importance sampling procedure used with the
    dmqmc_weighted_sampling is applied with weights which are introduced gradually. The weights
    :math:`w_{01}, \ldots, w_{n-1,n}` are altered, from 1 initially, by a factor of :math:`w^{1/N}` at
    the end of each Monte Carlo cycle, so that after N cycles the weights will have reached the values
    specified. They are then held constant until the end of the beta loop, at which point they are
    reset to 1.

    This helps psips to diffuse more appropriately initially.
**dmqmc_find_weights**
    Run a simulation to attempt to find appropriate weights for use in the DMQMC importance sampling
    procedure. This algorithm will attempt to find weights such that the population of psips is
    evenly distributed among the various excitation levels when the ground state is reached (at large
    beta values). The algorithm should be run for several beta loops until the weights settle down to a
    roughly constant value.

    This option should be used with **start_averaging**, to specify when the ground state
    has been reached.

    Warning: This feature is found to be unsuccessful for some larger lattices (for example, 6x6x6).
    The weights output should be checked. Increasing the number of psips used may improve the weights
    calculated.

    The weights are output at the end of each beta loop, in a form which can be copied directly into
    the input file.
**reduced_density_matrix** *nrdm* Integer.
                           *site_1 site_2 ... site_n* Integer list.

    Option to specify which reduced density matrices (RDMs) to obtain results for.
    
    *nrdm* specifies the number of RDMs which will be calculated. Then, on the next *nrdm* lines,
    a list of the sites making up the subsystem(s) to study should be given.

    With this option, one of the two options **ground_state_rdm** or **instantaneous_rdm** should
    also be used. Both options cannot be used together. Only one RDM may be considered (*nrdm*
    must be equal to 1) when using the **ground_state_rdm** option. Moreover, when using the
    **ground_state_rdm** option, the subsystem specified should be at most half the size of the
    system (which will always be sufficient for ground-state calculations).
**ground_state_rdm**
    For the subsystem specified with the **reduced_density_matrix** option, only accumulate the
    RDM when the ground state is reached. This is specified by the user using the
    **start_averaging** option. For each beta loop, the RDM will be averaged from this first
    iterations until the end of the beta loop. Results will then be output before the next loop
    is started.
**instantaneous_rdm**
    For the subsystem(s) specified with the **reduced_density_matrix** option, calculate the RDM(s)
    from the instantaneous psip distribution. This is done on the first iteration of every
    report loop.

    Results will only be output if using an option which makes use of these instantaneous RDM
    estimates, for example, **renyi_entropy_2**.
**output_rdm**
    Only available with the **ground_state_rdm** option.

    At the end of each beta loop, output the ground-state RDM accumulated to a file. This
    file will contain the RDM trace on the first line, followed by all RDM elements above and
    including the diagonal (labelled by their index).
**start_averaging** *N*
    Integer.

    If this option is specified then averaging of the ground-state reduced density matrix only begins at Monte
    Carlo cycle *N*. Hence, when only ground state properties are desired, the cycle at which the ground
    state is deemed to have been reached should be decided, and averaging should be started from this point.
    Thus, this feature should be used when calculating values which depend on the ground-state reduced
    density matrix (using **ground_state_rdm**).

    Futhermore, this option should also used when using **dmqmc_find_weights**, again, to specify
    when the ground state is reached.
**renyi_entropy_2**
    For all the subsystems specified with the **reduced_density_matrix** option, calculate the
    Renyi-2 entropy.

    The quantity output in the 'RDM(n) S2' columns is the instantaneous estimate of
    :math:`\sum_{ij}(\rho^n_{ij})^2`, where :math:`\rho^n` is the reduced density
    matrix for the nth subsystem specified by the user. The traces of the two replicas
    are in the columns named 'RDM(n) Trace 1' and 'RDM(n) Trace 2'. The finite_temp_analysis.py
    script in the tools directory can then be used to obtain a final temperature-dependent
    estimate of the Renyi-2 entropy from these quantities.

    This option cannot be used with **ground_state_rdm**.
**concurrence**
    At the end of each beta loop, the unnormalised concurrence and the trace of the reduced density matrix
    are output. The concurrence can then be calculated by running the average_entropy.py script in the tools
    subdirectory.

    This option should be used with the **ground_state_rdm** option. Temperature-dependent concurrence is
    not implemented in HANDE.
**von_neumann_entropy**
    At the end of each beta loop, the unnormalised von Neumann entropy and the trace of the reduced density matrix
    are output. The von Neumann entropy can then be calculated by running the average_entropy.py script in the tools
    subdirectory.

    This option should be used with the **ground_state_rdm** option. Temperature-dependent von Neumann entropy
    is not implemented in HANDE.
**exact_rdm_eigenvalues**
    When performing an **exact** calculaton, using this option will cause the eigenvalues of the RDM specified
    with the **reduced_density_matrix** option to be calculated and output.

    Note that the **ground_state_rdm** option must also be used. RDM eigenvalues can only be calculated for
    one subsystem in one simulation.

    The **use_all_sym_sectors** option is not implemented with **exact** calculations, and so cannot be used
    here.
**propagate_to_beta**
    Default False.

    Propagate a particular trial density matrix to a specific value of :math:`\beta` so that in the last step we are sampling the
    actual density matrix at this :math:`\beta`.
    To see this consider the function

    .. math::

        f(\tau) = \rho^{T}(\beta-\tau)\rho(\tau),

    where :math:`\rho^{T}` is a "trial" density matrix and :math:`\rho` is our usual density matrix.
    Note that

    .. math::

        f(0) = \rho^{T}(\beta) \rho(0) = \rho^{T}(\beta)

    and

    .. math::

        f(\tau=\beta) = \rho(\beta).

    Thus by propagating :math:`f` using the (appropriately modified) DMQMC algorithm we can sample the density matrix at a particular beta.
    This removes the difficulty of sampling the infinite temperature density matrix for systems with strong reference components, as typically
    the reference will be highly populated in the trial density matrix at any :math:`\beta > 0`.
    Currently only implemented for the UEG and k-space Hubbard model.
**init_beta** *beta*
    Real.

    Beta value the (trial) density matrix will initially be sampled at when using propagate_to_beta option.
    When analysing observables using the finite_temp_analysis.py script it is only this temperature value
    which has any meaning (in terms of averages with respect to the thermal density matrix) and is the last
    iteration in the simulation.
**metropolis_attempts** *nattempts*
    Integer.

    Default 0.

    Number of metropolis iterations per psip to be carried out when attempting to sample a trial density matrix.
**max_metropolis_moves** *max_move*
    Integer.

    Default: 2.

    A metropolis move is defined as a nfold excitation of the determiant under consideration.
    max_metropolis_move gives the maximum n considered in that nfold excitation.
**free_electron_trial**
    Default use "Hartree-Fock" trial density matrix.

    Use the non-interacting Hamiltonian in our trial density matrix. This is not as efficient as the default "Hartree-Fock" density matrix.
    If using the grand_canonical_intialisation option then metropolis_attempts can be set to zero as the canonical free-electron trial density
    matrix is already being sampled.
**grand_canonical_initialisation**
    Default False.

    Use the grand canonical partition function to guide the initialisation of the trial density matrix.
    This is usually a good starting point for the Metropolis algorithm and *is* also the starting point when using the free-electron trial
    density matrix.
**chem_pot** *chemical potential*
    Real.

    Chemical potential to be used to initialise the density matrix in the grand canonical ensemble. This can be calculated using chem_pot.py in tools/dmqmc/.
    If using the free-electron trial density matrix this chemical potential will produce the correct single particle occupancies, :math:`p_i`, so that the probability of occupying
    a given determinant is given by

    .. math::

        p(i_1, i_2, \dots, i_N) = \prod_i^N p_i,

    where,

    .. math::
        p_i = \frac{1}{e^{\beta (\varepsilon_i - \mu)} + 1}

    is the usual Fermi factor.
**fermi_temperature**
    Default: False.

    Rescale time step to be a multiple of :math:`1/T_F`, where :math:`T_F = E_F/k_B` is the Fermi Temperature, :math:`E_F` is the Fermi energy and :math:`k_B` is the Boltzman constant.
    This allows results to be output in terms of :math:`\Theta=T/T_F` which a useful quantity when comparing energy scales.

Calculation options: initiator-FCIQMC options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In addition to the options for general FCIQMC calculations, the following
options are also valid in initiator-FCIQMC calculations:

**initiator_population** *population*
    Integer.

    Default: 3.

    Set the (unsigned) population at which a determinant is considered to be an
    initiator determinant.  Setting this value to 0 retrieves the FCIQMC
    result.

Calculation options: parallel options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These options control the behaviour when run in parallel.  They do not affect
the result but can have a significant impact on performance.

**doing_load_balancing**
    Attempt to dynamically modify the hashing of determinants to processors
    so as to get a more even distribution of walkers across processors.
    See top-level comments in load_balancing.F90 for details.
**load_balancing_slots**
    Integer.

    Default: 20.

    Set the number of slots the walker list hash range is divided into.
    proc_map then contains N_p*load_balancing_slots number of slots.
    Setting this to too large a value will affect performance but could
    potentially result in a better distribution of walkers.
**load_balancing_pop**
    Long integer.

    Default 1000.

    Attempt to perform load balancing after the total number of walkers
    across processors is greater than load_balancing_pop. This is a
    system dependent variable and should be set so that the population
    is roughly stable at this value.
**percent_imbal**
    Real.

    Default 0.05.

    Desired percentage imbalance between the most/least populated processor
    and the average population. So, min_pop ~ (1-percent_imbal)*av_pop and
    max_pop ~ (1+percent_imbal)*av_pop.
**max_load_attempts**
    Integer.

    Default 2.

    Load balancing will be attempted once per report loop until max_load_attempts
    is reached.
**write_load_info**
    Default: false.

    Write out the population of the most and least heavily populated processor
    before and after load balancing is carried out. Also print out the
    minimum slot population on the most populated processor which will
    indicate if load balancing is possible.

**use_mpi_barriers**
    Default: false.

    Perform MPI_Barrier calls before the main MPI communication calls (both
    for communication of the spawned list, and any semi-stochastic
    communication). These are timed, and the total time spent in these calls
    is output at the end of a simulation. This option is useful for assessing
    issues in load balancing, as it will allow you to see when certain
    processors take longer to perform their work than others. This is turned
    off by default because such calls may have an initialisation time which
    scales badly to many processors.


Calculation options: estimate canonical kinetic energy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    **nkinetic_cycles**
    Integer.

    Default 1.

    Perform nkinetic_cycles * init_pop Mote Carlo iterations for estimating the
    canonical kinetic energy.

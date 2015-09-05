.. _semi_stoch_table:

semi_stoch options
==================

The semi-stochastic approach divides the Hilbert space into two regions: a small region in
which the action of the Hamiltonian is applied exactly, and the remainder of the Hilbert
space, in which the action is applied stochastically.  This can substantially reduce the
stochastic error in many cases.

``space``
    type: string.

    Required.

    Possible values: 'read', 'high'.

    The type of deterministic space to use.  Using 'read' uses a deterministic space
    produced from a previous calculation and saved to file using the semi_stoch 
    ``write`` option (the ``write_determ_space`` can be used but is now deprecated).
    Using 'high' sets the deterministic space to consist of the states with
    the highest population when the semi-stochastic projection is enabled.
``size``
    type: integer.

    Required if ``space`` is 'high', otherwise ignored.

    The number of states to include in the semi-stochastic space.
``start_iteration``
    type: integer.

    Optional.  Default: 1.

    The number of iterations to perform, during which the action of the Hamiltonian is
    applied entirely stochastically, before semi-stochastic projection is enabled.  This
    allows for a period for the population to grow and the ground-state wavefunction to
    emerge before the deterministic space is selected.
``shift_start_iteration``
    type: integer.

    Optional.  Default: None.  Overrides ``start_iteration``.

    The number of iterations to perform after the shift is varied (i.e. after the
    ``target_population`` is reached) before the semi-stochastic projection is enabled.
``separate_annihilation``
    type: boolean.

    Optional.  Default: true.

    If true, the deterministic amplitudes are communicated separately at the cost of an
    additional MPI call.  If false, the annihilation of particles created from
    deterministic and stochastic projections are performed separately, which removes the
    need for an additional MPI call at the cost of communicating an additional
    :math:`\mathcal{O}(N_p N_D)` more amplitudes, where :math:`N_p` is the number of
    processors and :math:`N_D` the size of the deterministic space.  If the deterministic
    space is small and communication latency high, setting ``separate_annihilation`` to
    false might improve performance.  For most systems and computer architectures, the
    default value is faster.
``write``
    type: boolean or integer.

    Optional.  Default: false.

    Write out the deterministic space to file of form ``SEMI.STOCH.X.H5``, where ``X`` is
    the file id.  If set to true, ``X`` will be the smallest non-negative id such that
    ``SEMI.STOCH.X.H5`` does not already exist, otherwise the value provided is used as
    the file id.
``read``
    type: integer.

    Optional.  Default: largest value of ``X`` such that the file ``SEMI.STOCH.X.H5`` exists.

    Index of the file containing the deterministic space produced from a previous
    calculation.

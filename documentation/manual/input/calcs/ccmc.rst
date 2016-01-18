.. _ccmc:

Coupled Cluster Monte Carlo
===========================

.. code-block:: lua

    ccmc {
        sys = system,
        qmc = { ... },
        ccmc = { ... },
        restart = { ... },
        reference = { ... },
        qmc_state = qmc_state,
    }

Returns:
    a qmc_state object.

``ccmc`` performs a coupled cluster Monte Carlo (CCMC) calculation [Thom10]_ on a system.

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
``ccmc``
    type: lua table.

    Required.

    Further options to control the CCMC algorithm.  See :ref:`ccmc_table`.
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

        The calculation must be of the same type to succesfully resume, but this is not checked.

.. _ccmc_table:

ccmc options
------------

``move_frequency``
    type: integer

    Optional.  Default: 5.

    Allow excitors to move processors every :math:`2^x` iterations, where :math:`x` is the
    value of ``move_frequency``, in order to allow all composite excitors to be correctly
    sampled.  Relevant only when performing CCMC with MPI parallelisation.  A large value
    may introduce a bias.  Modify with caution.
``cluster_multispawn_threshold``
    type: float.

    Optional.  Default: :math:`2^{31}-1`.

    Set the maximum value of :math:`A_C/p_C`, where :math:`A_C` is the cluster amplitude
    and :math:`p_C` is the probability of selecting the cluster.  A cluster with a value
    above this is split into multiple spawning attempts.  The default value essentially
    disables this but a smaller option can substantially reduce population blooms, albeit
    potentially at a significant computational cost.

    .. note::

        This is an experimental option and feedback is most welcome.  The current
        recommendation is to use the smallest setting such that large blooms do not occur.

``full_non_composite``
    type: boolean.

    Optional.  Default: false.

    If true, allow all non-composite clusters to attempt to spawn each iteration.  The
    original CCMC algorithm involves randomly selecting a cluster of arbitrary size
    consisting of any set of excitors and then making spawning attempts from it.  The full
    non-composite algorithm is a simple modification in which all occupied non-composite
    clusters (i.e. those consisting of the reference or just a single excitor) are
    (deterministically) selected and composite clusters (involving two or more excitors)
    are randomly selected to make spawning attempts.  This has been shown to give
    substantially more stable dynamics and reduce the plateau height in several systems.
``linked``
    type: boolean.

    Optional.  Default: false.

    If true, sample the linked coupled cluster equations instead of the unlinked coupled
    cluster equations.  The original CCMC algorithm solves the equations

    .. math::

        \langle D_m | \hat{H} - E | \psi_{CC} \rangle = 0.

    It is possible to instead sample the equivalent equations

    .. math::

        \langle D_m | e^{-\hat{T}} (\hat{H} - E) | \psi_{CC} \rangle = 0.

    Using the Hausdorff expansion of the Hamiltonian and the linked cluster theorem means 
    that the only clusters which contribute are those with at most four excitors and where 
    the excitation sampled from the Hamiltonian has an orbital in common with each excitor 
    in the cluster operator. Using this option can give substantial reductions in the 
    plateau height.

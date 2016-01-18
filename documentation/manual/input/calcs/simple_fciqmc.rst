.. _simple_fciqmc:

Full Configuration Interaction Quantum Monte Carlo (simple)
===========================================================

Find the ground state of a system via FCIQMC [Booth09]_.

.. code-block:: lua

    simple_fciqmc {
        sys = system,
        sparse = true/false,
        qmc = { ... },
        restart = { ... },
        reference = { ... },
        qmc_state = qmc_state,
    }

Returns:
    a qmc_state object.

``simple_fciqmc`` performs a full configuration interaction quantum Monte Carlo (FCIQMC)
calculation on a system using an explicitly calculated and stored Hamiltonian matrix.

.. warning::

    This is an **extremely** simple implementation of FCIQMC.  In particular it makes no
    effort to be efficient (in time or memory), is not parallelised, and does not include
    any advanced features.  It is, however, useful for educational purposes and
    (occasionally) hacking experimental ideas quickly.  Do **not** use for production
    calculations.

Options
-------

``sys``
    type: system object.

    Required.

    The system on which to perform the calculation.  Must be created via a system
    function.
``sparse``
    type: boolean.

    Optional.  Default: true.

    Store the Hamiltonian matrix in a sparse matrix format.
``qmc``
    type: lua table.

    Required.

    Further options that are common to all implemented QMC algorithms.  Note that 
    options relating to memory usage, excitation generation and real amplitudes are not
    implemented for ``simple_fciqmc``.  See :ref:`qmc_table`.
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

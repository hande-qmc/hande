.. _uccmc:

Unitary Coupled Cluster Monte Carlo
===================================

.. code-block:: lua

    uccmc {
        sys = system,
        qmc = { ... },
        ccmc = { ... },
        uccmc = { ... },
        restart = { ... },
        reference = { ... },
        logging = { ... },
        output = { ... },
        blocking = { ... },
        qmc_state = qmc_state,
    }

Returns:
    a qmc_state object.

``uccmc`` performs a coupled cluster Monte Carlo (CCMC) calculation [Filip20]_ on a system.

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
``uccmc``
    type: lua table.

    Required.

    Further options to control the UCCMC algorithm.  See :ref:`uccmc_table`.
``restart``
    type: lua table.

    Optional.

    Further options to control restarting the calculation from a previous calculation.
    See :ref:`restart_table`.
``reference``
    type: lua table.

    Optional.

    Further options to select the reference state used.  See :ref:`reference_table`.
``logging``
    type: lua table.

    Optional.

    Further options to enable various logging outputs from a UCCMC simulation. Only
    available when compiled in debug mode. See :ref:`logging_table` for information
    on current options.
``output``
    type: lua_table.

    Optional.

    Further options to enable direction of calculation output to a different file.
    See :ref:`output_table` for more information.

.. _uccmc_table:

uccmc options
-------------
``pow_trunc``
    type: integer.

    Optional.  Default: 12.

    Polynomial order at which to truncate UCCMC Taylor expansion.
``average_wfn``
    type: logical.

    Optional. Default: false.

    Whether to print out the average wavefunction sampled over the course of the UCCMC calculation.
``trot``
    type: logical.

    Optional. Default: false.

    Whether to run a trotterized UCCMC calculation, rather than using the full exponential ansatz.
``threshold``
    type: float.

    Optional. Default: -1.0.

    Threshold in the ratio of probability selection and cluster amplitude below which selections are discarded.
    Can be used to stabilise UCCMC calculations.

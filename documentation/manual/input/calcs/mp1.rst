MP1 wavefunction
================

.. code-block:: lua

    psip_list = mp1_mc {
        sys = system,
        mp1 = {...},
        qmc = {...},
        ccmc = {...},
    }

Returns:
    a ``particle_t`` object.


``mp1_mc`` creates a deterministic MP1 wavefunction and stochastically coarse-grains (rounds down small amplitudes) 
it into a ``particle_t`` object that can be used to initialise a subsequent CCMC calculation. The MP1 wavefunction will be ignored 
if restarting from a restart file.

Options
-------

``sys``
    type: system object.

    Required.

    The system on which to perform the calculation.  Must be created via a system
    function.
``mp1``
    type: lua table.

    Required.

    Further options controlling the calculation.

``qmc``
    type: lua table.

    Optional.

    If provided, will override the ``mp1`` input options where relevant. 
``ccmc``
    type: lua table.

    Optional.

    If provided, will override the ``mp1`` input options where relevant. 

.. note::

    If you wish to chain together a MP1 calculation and a CCMC calculation, one way to do it would be:

    .. code-block:: lua

        sys = {...}
        qmc_opt = {...}
        ccmc_opt = {...}

        psip_list = mp1_mc{
            sys = sys,
            qmc = qmc_opts,
            ccmc = ccmc_opts,
        }

        ccmc {
            sys = sys,
            qmc = qmc_opts,
            ccmc = ccmc_opts,
            psip_list = psip_list,
            reference = {...},
        }

    Note the lack of commas after the main tables. This makes sure the MP1 wavefunction (the ``psip_list`` object) is compatible 
    with the subsequent CCMC calculation.


MP1 options
-----------
``D0_population``
    type: integer.

    Required.

    Set the initial population on the reference determinant. 

    .. note::

        This will overwrite the value given in the ``qmc`` table.

``state_size``
    type: integer.

    Maximum number of excitors to store in the “main” list, 
    which holds the number of particles on the state and related information such as the diagonal Hamiltonian matrix element. 
    The number of elements that can be stored usually should be of the same order as the target population. 

    If negative, then the absolute value is used as the maximum amount of memory in MB to use for this information.

    .. note::

        This will overwrite the value given in the ``qmc`` table.

``real_amplitudes``
    type: boolean.

    Optional. Default: false.

    Allow amplitudes to take non-integer weights. This will often significantly reduce the stochastic noise 
    in the Monte Carlo estimates.

    .. note::

        This should be the same as the subsequent calculation.

``spawn_cutoff``
    type: float.

    Optional. Default: 0.01.

    The threshold for stochastic rounding.

``rng_seed``
    type: integer.

    Optional. Default: generate a seed from a hash of the time and calculation UUID.

    The seed used to initialise the random number generator.

``even_selection``
    type: boolean.

    Optional. Default: false.

    .. note::

        Must be true if true in the subsequent CCMC calculation.

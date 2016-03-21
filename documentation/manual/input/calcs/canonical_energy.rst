Canonical total energy
======================

.. code-block:: lua

    canonical_estimates {
        sys = system,
        canonical_estimates = { ... },
    }


``canonical_estimates`` calculates various estimates for properties of a system in the
canonical ensemble at a given temperature, using knowledge of the grand canonical ensemble
and the single-particle eigenvalues of the underlying non-interacting system.  See
[Malone15]_ for details. Currently only implemented for the UEG and read_in.

Options
-------

``sys``
    type: system object.

    Required.

    The system on which to perform the calculation.  Must be created via a system
    function.
``canonical_estimates``
    type: lua table.

    Required.

    Further options controlling the calculation.

kinetic options
---------------

``ncycles``
    type: integer.

    Required.

    The number of Monte Carlo iterations to perform.  Each iteration produces
    independent estimates based upon the ``nattempts`` made.
``nattempts``
    type: integer.

    Required.

    Number of determinants within the canonical ensemble we attempt to generate each Monte
    Carlo cycle.
``beta``
    type:  float.

    Required.

    The temperature of the system.
``fermi_temperature``
    type: boolean.

    Optional.  Default: false.

    If true, rescale ``beta`` as the inverse reduced temperature: :math:`\tilde{\beta} = 1/\Theta = T_F/T`,
    where :math:`T_F` is the Fermi temperature.  If false, ``beta`` is taken to be in
    atomic units.
``rng_seed``
    type: integer.

    Optional.  Default: generate a seed from a hash of the time and calculation UUID.

    The seed used to initialise the random number generator.

Electron gases
==============

An electron gas contains interacting electrons in some geometry with a constant
compensating positive charge.

Uniform electron gas
--------------------

.. code-block:: lua

    ueg {
        -- options,
    }

Returns:
    a system object.

``ueg`` creates a system object for the (conventional) electron gas:

.. math::

    H = -\frac{1}{2} \sum_i \nabla_i^2 + \sum_{i<j} \frac{1}{r_{ij}}

(including an appropriate uniform background potential to counteract the charge),
using a single-particle basis of plane waves, :math:`\psi_{\bf k} = e^{i {\bf k}.{\bf r}}`.

Options
^^^^^^^

``sys``
    type: system object produced by a previous call.

    Optional.

    If provided, a previously created system object is updated with the new settings
    supplied, otherwise a new system object is created.
``electrons``
    type: integer.

    Required.

    Number of electrons in the unit cell.
``ms``
    type: integer.

    Required.

    Set the spin polarisation of the system in units of electron spin (i.e. a single      
    electron can take values 1 or -1).
``sym``
    type: integer.

    Required.

    Set the symmetry (i.e. crystal momentum) of the system.  This is the index of
    a specific wavevector; see the output produced by creating a system for possible
    values and their corresponding wavevectors.
``rs``
    type: float.

    Optional.  Default: 1.

    Set the density, :math:`r_s`, of the UEG.
``cutoff``
    type: float.

    Optional.  Default: 3.

    Set the maximum kinetic energy of the orbitals included in the basis set.

    Note that this is in scaled units of :math:`(2\pi/L)^2`, where :math:`L` is the
    dimension of simulation cell defined by *electrons* and *rs* and is compared to
    the kinetic energy of each plane-wave without the twist angle included.  In
    this way the cutoff can be kept constant whilst the twist is varied and the
    basis set used will remain consistent.
``dim``
    type: integer.

    Optional.  Default: 3.

    Set the dimension of the electron gas.  2- and 3-dimensional gases are implemented.
``twist``
    type: :math:`N`-dimensional vector

    Optional.  Default: 0 in each dimension.

    Apply a twist to the wavevector grid.  The twist is an :math:`N`-dimensional vector in
    units of :math:`2\pi`.  The twist angle should be within the first Brillouin zone, and
    hence the components should be between -0.5 and +0.5.
``chem_pot``
    type: float.

    Optional.  Used only in interaction-picture DMQMC calculations.

    Set the chemical potential of the system.

Ringium
-------

.. code-block:: lua

    ringium {
        -- options,
    }

Returns:
    a system object.

Ringium [Loos13]_, is a 1D system of electrons confined to a ring of radius :math:`R`:

.. math::

    H = -\frac{1}{2R^2} \sum_i \frac{\partial^2}{\partial\theta_i^2} + \sum_{i<j} \frac{1}{r_{ij}}

where :math:`r_{ij} = R\sqrt{2-2\cos(\theta_i-\theta_j)}`, using a single-particle
basis of functions :math:`\psi_n = e^{i n \theta}`.  As it is 1D, the different 
spin polarisations are degenerate, so without loss of generality all electrons
are forced to be spin up.

Options
^^^^^^^

``sys``
    type: system object produced by a previous call.

    Optional.

    If provided, a previously created system object is updated with the new settings
    supplied, otherwise a new system object is created.
``electrons``
    type: integer

    Required.

    Number of electrons in the system.
``radius``
    type: float

    Required.

    The radius of the ring.
``maxlz``
    type: integer

    Required.

    The maximum angular momentum of the orbitals used in the basis set.

    Note that this is in units of :math:`\frac{\hbar}{2}` and must have opposite
    parity to the number of electrons.

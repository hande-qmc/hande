Model lattice systems
=====================

Hubbard model (momentum-space)
------------------------------

.. code-block:: lua

    hubbard_k {
        -- options,
    }

Returns:
    a system object.

``hubbard_k`` creates a system object for the Hubbard model:

.. math::

    H = -t \sum_{<r,r'>,\sigma} c^\dagger_{r,\sigma} c_{r',\sigma} + U \sum_r n_{r,\uparrow} n_{r,\downarrow}

using a single-particle basis of Bloch functions, :math:`\psi_k`:

.. math::

    \psi_k(r) = e^{ik.r} \sum_i \phi_i(r)

where :math:`\phi_i(r)` is a single-particle basis function centred on site :math:`i`
in real space.

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
``lattice``
    type: :math:`N\ N`-dimensional vectors of floats.

    Required.

    Unit cell on which periodic boundary conditions are placed.  See below.
``ms``
    type: integer.

    Required.

    Set the spin polarisation of the system in units of electron spin (i.e. a single
    electron can take values 1 or -1).
``sym``
    type: integer.

    Required for deterministic calculations and highly recommended for Monte Carlo calculations.

    Set the symmetry (i.e. crystal momentum) of the system.  This is the index of
    a specific wavevector; see the output produced by creating a system for possible
    values and their corresponding wavevectors.  If not specified (and no reference
    determinant supplied for a calculation) the symmetry used is that of a determinant
    selected using the Aufbau principle.
``U``
    type: float.

    Optional.  Default: 1.

    Specifies the :math:`U` parameter in the Hamiltonian.
``t``
    type: float.

    Optional.  Default: 1.

    Specifies the :math:`t` parameter in the Hamiltonian.
``twist``
    type: :math:`N`-dimensional vector

    Optional.  Default: 0 in each dimension.

    Apply a twist to the wavevector grid.  The twist is an *ndim*-dimensional vector in
    units of :math:`2\pi`.  The twist angle should be within the first Brillouin zone, and
    hence the components should be between -0.5 and +0.5.
``verbose``
    type: boolean.

    Optional.  Default: true.

    Print out the single-particle basis set.

Hubbard model (real-space)
--------------------------

.. code-block:: lua

    hubbard_real {
        -- options,
    }

Returns:
    a system object.

``hubbard_real`` creates a system object for the Hubbard model:

.. math::

    H = -t \sum_{<r,r'>,\sigma} c^\dagger_{r,\sigma} c_{r',\sigma} + U \sum_r n_{r,\uparrow} n_{r,\downarrow}

using a single-particle basis of functions in real-space.

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
``lattice``
    type: :math:`N\ N`-dimensional vectors of floats.

    Required.

    Unit cell on which periodic boundary conditions are placed.  See below.
``ms``
    type: integer.

    Required.

    Set the spin polarisation of the system in units of electron spin.
``U``
    type: float.

    Optional.  Default: 1.

    Specifies the :math:`U` parameter in the Hamiltonian.
``t``
    type: float.

    Optional.  Default: 1.

    Specifies the :math:`t` parameter in the Hamiltonian.
``finite``
    type: boolean.

    Optional.  Default: false.

    If false then periodic boundary conditions are applied to the unit cell, otherwise the
    system specified by the lattice is treated as an isolated set of sites.
``verbose``
    type: boolean.

    Optional.  Default: true.

    Print out the single-particle basis set.

Heisenberg model
----------------

.. code-block:: lua

    heisenberg {
        -- options,
    }

Returns:
    a system object.

``heisenberg`` creates a system object for the Heisenberg model, which models a set of
spin 1/2 particles on a lattice:

.. math::

    \hat{H} = -J \sum_{\langle i,j \rangle} \hat{\boldsymbol{S}}_i \cdot \hat{\boldsymbol{S}}_j  - h_z \sum_i \hat{S}_{iz} - h_z' \sum_i \hat{S}_{iz}^{\xi},

where :math:`h_z` and :math:`h_z'` denote the magnetic field strength and
staggered magnetic field strength, respectively, and :math:`\xi`
is equal to +1 for sites on sublattice 1 and is equal to -1 for sites on
sublattice 2.

Options
^^^^^^^

``sys``
    type: system object produced by a previous call.

    Optional.

    If provided, a previously created system object is updated with the new settings
    supplied, otherwise a new system object is created.
``lattice``
    type: :math:`N\ N`-dimensional vectors of floats.

    Required.

    Unit cell on which periodic boundary conditions are placed.  See below.

    .. warning::

        For efficiency reasons it is assumed that the smallest dimension lattice vector is
        greater than 2 if periodic boundary conditions are used.

``ms``
    type: integer.

    Required.

    Set the spin polarisation of the system in units of 1/2.
``J``
    type: float.

    Optional.  Default: 1.

    Set the coupling constant for the Heisenberg model.
``magnetic_field``
    type: float.

    Optional.  Default: 0.
``staggered_magnetic_field``
    type: float.

    Optional.  Default: 0.

    .. note:: 

        Specifying non-zero values for both ``magnetic_field`` and ``staggered_magnetic_field``
        is not currently possible.

``finite``
    type: boolean.

    Optional.  Default: false.

    If false then periodic boundary conditions are applied to the unit cell, otherwise the
    system specified by the lattice is treated as an isolated set of sites.
``triangular``
    type: boolean.

    Optional.  Default: false.

    If true, then a triangular lattice of sites on which the spins reside is used,
    requiring a 2D lattice.  The default is to use a :math:`N`-dimensional cubic
    arrangement of sites.
``verbose``
    type: boolean.

    Optional.  Default: true.

    Print out the single-particle basis set.

Chung-Landau model
------------------

.. code-block:: lua

    chung_landau {
        -- options,
    }

Returns:
    a system object.

``chung_landau`` creates a system object for the system of spinless fermions proposed by
Chung and Landau:

.. math::

    H = -t \sum_{\langle r,r' \rangle} c^\dagger_{r} c_{r'} + U \sum_{\langle r,r' \rangle} n_{r} n_{r'}

using a single-particle basis of functions in real-space.

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

    Number of fermions in the unit cell.
``lattice``
    type: :math:`N\ N`-dimensional vectors of floats.

    Required.

    Unit cell on which periodic boundary conditions are placed.  See below.
``U``
    type: float.

    Optional.  Default: 1.

    Specifies the :math:`U` parameter in the Hamiltonian.
``t``
    type: float.

    Optional.  Default: 1.

    Specifies the :math:`t` parameter in the Hamiltonian.
``finite``
    type: boolean.

    Optional.  Default: false.

    If false then periodic boundary conditions are applied to the unit cell, otherwise the
    system specified by the lattice is treated as an isolated set of sites.
``verbose``
    type: boolean.

    Optional.  Default: true.

    Print out the single-particle basis set.

Specifying the lattice
----------------------

The lattice is specified as a table of vectors.  Sites (on which a spin or electron
resides) are at unit locations on the grid.  The unit cell (or, if periodic boundary
conditions are not used, the geometry of the 'flake' essentially cut out of the infinite
lattice) are given in this basis.  The lattice variable hence requires :math:`N` vectors,
each of dimension :math:`N`.  This is specified in lua by a nested table.  For example:

.. code-block:: lua

    lattice = { { 10 } }

sets a 1D system, with the unit cell containing 10 sites;

.. code-block:: lua

    lattice = { { 2, 0 }, { 0, 2 } }

sets a 2D system, with the unit cell containing 4 sites; and

.. code-block:: lua

    lattice = { { 3, 3 }, { 3, -3 } }

sets a 2D system, with the (square) unit cell containing 18 sites and rotated by
:math:`45^\circ` relative to the primitive lattice.

HANDE supports 1-, 2- and 3-dimensional lattices.  Lattice vectors must be orthogonal.

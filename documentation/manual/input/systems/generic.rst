Generic systems
===============

.. code-block:: lua

    read_in {
        -- options
    }

Returns:
    a system object.

A generic system, including atoms and molecules, can be specified by providing a file
containing information about the single-particle basis set and the one- and two-body
integrals between these basis functions.  This file is in FCIDUMP format
[Knowles89]_, which can be produced by several quantum chemistry packages including
MOLPRO, Q-Chem (via additions from Alex Thom) and PSI4 (via a plugin from James Spencer).

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

    Optional.

    Set the symmetry of the system.  This is the index of a specific irreducible
    representation from the FCIDUMP file; see the output produced by creating a system for
    possible values.  If not specified (and no reference determinant supplied for a calculation)
    the symmetry used is that of a determinant selected using the Aufbau principle.
``Lz``
    type: boolean.

    Optional.  Default: false.

    If true, enable :math:`L_z` symmetry.  See below for details.
``int_file``
    type: string.

    Optional.  Default: 'FCIDUMP'.

    Specify the FCIDUMP file containing the integrals and information relating to the
    single-particle basis.
    Can also be an HDF5 file previously produced by HANDE from a FCIDUMP for considerable gains in
    initialisation speed and read_in file size. To produce such a file you can use the 
    dump_hdf5_generic_system lua function. For further information see :ref:`utils_hdf5_system_dump`.
``dipole_int_file``
    type: string.

    Optional.  No default.

    Specify a FCIDUMP-like file containing the dipole integrals, i.e. :math:`\langle i | x | i \rangle`, in a given direction.
    
    Not currently used. 
``CAS``
    type: 2D-vector of integers.

    Optional.  No default.

    If specified, then the basis set is restricted to a given complete active space,
    whereby ``CAS = {N,M}`` corresponds to allowing only :math:`N` electrons to be distributed
    among :math:`M` spin orbitals.  Any additional electrons are 'frozen' (i.e. forced to
    be in the lowest spin orbitals) and any additional high-energy spin orbitals are
    removed from the basis set.
    Functionality not compatible with reading from an hdf5 file; to use this in combination with
    hdf5 initialisation use setting when generating hdf5 file.
``verbose``
    type: boolean.

    Optional.  Default: true.

    Print out the single-particle basis set.

``complex``
    type: boolean.

    Optional. Default: false.

    Specify if the calculation should use complex dynamics in any calculation performed, 
    and if the FCIDUMP supplied is complex-formatted. Currently only compatible with
    fci calculations.


:math:`L_z` symmetry
--------------------

For cylindrically symmetrical systems, the :math:`L_z` (z-component of orbital angular momentum)
operator commutes with the Hamiltonian, and this can be a convenient symmetry to conserve.
:math:`L_z` is measured in units of :math:`\hbar`.  Normal FCIDUMP files do not contain orbitals which are
eigenfunctions of the :math:`L_z` operator, so they must be transformed using post-processing.  The
TransLz  script from the `NECI <https://github.com/ghb24/NECI_STABLE>`_ project can be
used for this purpose. The FCIDUMP file header format has been modified to include
additional parameters: SYML, and SYMLZ which have a list of values, one for each orbital.
SYML gives the magnitude of L for the orbital if known (or -20 if not) but is not used.
SYMLZ give the eigenvalue of :math:`L_z` (the :math:`m_l` value).  Orbitals with defined values of :math:`L_z` are
likely to be complex-valued, but luckily the integrals involving them are not, so although
the FCIDUMP file must be translated, it still retains the same format (see comments in
``src/read_in.F90`` and ``src/molecular_integrals.F90`` for details if you wish to create
FCIDUMP files by other means).  

.. warning::

    These transformed integral files require you to enforce :math:`L_z` symmetry and will produce
    incorrect results if you do not.

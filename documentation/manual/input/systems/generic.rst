.. _generic_systems:

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
See :ref:`generating_integrals` for more details.

Options
^^^^^^^

``sys``
    type: system object produced by a previous call.

    Optional.

    If provided, a previously created system object is updated with the new settings
    supplied, otherwise a new system object is created.
``electrons``
    type: integer.

    Optional.  If specified, then ``ms`` **must** be specified.

    Number of electrons in the unit cell.  If not provided, the value in the FCIDUMP file is used.
``ms``
    type: integer.

    Optional.  If specified, then ``electrons`` **must** be specified.

    Set the spin polarisation of the system in units of electron spin (i.e. a single
    electron can take values 1 or -1).  If not provided, the value in the FCIDUMP file is used.
``sym``
    type: integer or string.

    Optional. Default: ``aufbau``.

    Set the symmetry of the system if a reference determinant is not provided. This can be
    set to:

    - An integer specifying the index of a specific irreducible representation from the FCIDUMP
      file; see the output produced by creating a system for possible values.
    - ``aufbau``. Uses the symmetry of a determinant selected using the Aufbau principle.
    - ``tot_sym``. Uses the totally symmetric representation, whatever its index may be.;
``Lz``
    type: boolean.

    Optional.  Default: false.

    If true, enable :math:`L_z` symmetry.  See below for details.
``int_file``
    type: string.

    Optional.  Default: 'FCIDUMP'.

    Specify the FCIDUMP file containing the integrals and information relating to the
    single-particle basis. For details of the format see :ref:`fcidump_format`.
    This can also be an HDF5 file previously produced by HANDE from a FCIDUMP via the
    ``write_read_in_system`` function (see :ref:`utils_hdf5_system_dump`), which is
    both more compact in size and considerably faster to process.
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
    among :math:`2M` spin orbitals.  Any additional electrons are 'frozen' (i.e. forced to
    be in the lowest spin orbitals) and any additional high-energy spin orbitals are
    removed from the basis set.

    .. warning::

        This functionality is not compatible with reading from an HDF5 file; to use a CAS
        in combination with HDF5 initialisation, create the HDF5 file using a system with
        the desired CAS.

``verbose``
    type: boolean.

    Optional.  Default: true.

    Print out the single-particle basis set.

``complex``
    type: boolean.

    Optional. Default: false.

    Specify if the calculation should use complex dynamics in any calculation performed, 
    and if the FCIDUMP supplied is complex-formatted. Currently compatible with
    fci, fciqmc, ccmc and dmqmc (including ip-dmqmc) calculations.

``max_integral_chunk``
    type: integer

    Optional. Default :math:`2^{31} - 1`.

    Maximum number of MPI objects to broadcast in a single call for two body integrals.
    Above this value a contiguous MPI type is used instead.

    .. warning::

        This functionality is included only for ease of testing. It should not be used
        for production calculations.

:math:`L_z` symmetry
--------------------

For cylindrically symmetrical systems, the :math:`L_z` (z-component of orbital angular momentum)
operator commutes with the Hamiltonian, and this can be a convenient symmetry to conserve.
:math:`L_z` is measured in units of :math:`\hbar`.  Normal FCIDUMP files do not contain orbitals
which are eigenfunctions of the :math:`L_z` operator, so they must be transformed using
post-processing. 

SYMLZ give the eigenvalue of :math:`L_z` (the :math:`m_l` value).  Orbitals with defined values of 
:math:`L_z` are likely to be complex-valued, but luckily the integrals involving them are not, so
althoughthe FCIDUMP file must be translated, it still retains the same format (see comments in
``src/read_in.F90``, ``src/molecular_integrals.F90`` and :ref:`fcidump_format` for details if you wish to create
FCIDUMP files by other means).  

The FCIDUMP file header format has been modified to include
additional parameters: SYML, and SYMLZ which have a list of values, one for each orbital.

SYML gives the magnitude of L for the orbital if known (or -20 if not) but is not used.

.. note::

    There is a tool provided in ``tools/fcidump/lz_fcidump.py`` that can generate :math:`L_z`-transformed FCIDUMPs 
    from PySCF calculations. To run this script, you need to:

    - compile the ``lz_transform.f90`` Fortran script,
    - have PySCF, 
    - install a Python package called ``f90nml``.

    See the comments within the Python script for further help. 

.. warning::

    These transformed integral files require you to enforce :math:`L_z` symmetry and will produce
    incorrect results if you do not.

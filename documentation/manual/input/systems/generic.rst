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
[KnowlesHandy1989]_, which can be produced by several quantum chemistry packages including
MOLPRO, Q-Chem (via additions from Alex Thom) and PSI4 (via a plugin from James Spencer).

Options:

sys
    type: system object produced by a previous call.

    Optional.

    If provided, a previously created system object is updated with the new settings
    supplied, otherwise a new system object is created.
electrons
    type: integer.

    Required.

    Number of electrons in the unit cell.
ms
    type: integer.

    Required.

    Set the spin polarisation of the system in units of electron spin.
sym
    type: integer.

    Required.

    Set the symmetry of the system.  This is the index of a specific irreducible
    representation from the FCIDUMP file; see the output produced by creating a system for
    possible values.
Lz
    type: boolean.

    Optional.  Default: false.

    If true, enable Lz symmetry.  See below for details.
int_file
    type: string.

    Optional.  Default: 'FCIDUMP'.

    Specify the FCIDUMP file containing the integrals and information relating to the
    single-particle basis.
dipole_int_file
    type: string.

    Optional.  No default.

    Specify a FCIDUMP-like file containing the dipole integrals, i.e. `:math:\langle i | x | \rangle`, in a given direction.
    
    Not currently used. 
CAS
    type: 2D-vector of integers.

    Optional.  No default.

    If specified, then the basis set is restricted to a given complete active space,
    whereby ``CAS = {M,N}`` corresponds to allowing only N electrons to be distributed
    among :math:`2M` spin orbitals.  Any additional electrons are 'frozen' (i.e. forced to
    be in the lowest spin orbitals) and any additional high-energy spin orbitals are
    removed from the basis set.

Lz symmetry
-----------

For cylindrically symmetrical systems, the Lz (z-component of orbital angular momentum)
operator commutes with the Hamiltonian, and this can be a convenient symmetry to conserve.
Lz is measured in units of hbar.  Normal FCIDUMP files do not contain orbitals which are
eigenfunctions of the Lz operator, so they must be transformed using post-processing.  The
TransLz  script from the `NECI <https://github.com/ghb24/NECI_STABLE>`_ project can be
used for this purpose. The FCIDUMP file header format has been modified to include
additional parameters: SYML, and SYMLZ which have a list of values , one for each orbital.
SYML gives the magnitude of L for the orbital if known (or -20 if not) but is not used.
SYMLZ give the eigenvalue of Lz (the m_l value).  Orbitals with defined values of Lz are
likely to be complex-valued, but luckily the integrals involving them are not, so although
the FCIDUMP file must be translated, it still retains the same format (see comments in
``src/read_in.F90`` and ``src/molecular_integrals.F90`` for details if you wish to create
FCIDUMP files by other means).  

.. warning::

    These transformed integral files require you to enforce Lz symmetry and will produce
    incorrect results if you do not.

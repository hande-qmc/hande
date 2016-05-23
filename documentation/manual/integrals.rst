.. _generating_integrals:

Generating integrals
====================

HANDE can treat :ref:`systems <generic_systems>` other than model Hamiltonians by reading in the necessary
integrals in the FCIDUMP format [Knowles89]_.  Many quantum chemistry packages can
generate them following Hartree-Fock calculations, including:

HORTON
   https://theochem.github.io/horton/
MOLPRO
   https://www.molpro.net/
PSI4
    http://psicode.org.  Requires the fcidump plugin (https://github.com/hande-qmc/fcidump).
Q-Chem
   http://www.q-chem.com.  FCIDUMP code contributed by Alex Thom.

We most frequently use PSI4 and Q-Chem and so these tend to be better tested.  Note that
the computational cost of the calculations in HANDE vastly outweighs the cost of the
underlying SCF calculations and so the efficiency of the code used to generate the
integrals is usually not a key factor.  Please consult the documentation of the code of
interest regarding how to run SCF calculations and generate the integrals in the FCIDUMP
format.

.. _fcidump_format:

FCIDUMP format
==============

The format of FCIDUMP files used by HANDE is partially defined in [Knowles89]_. It consists
of a namelist header, containing various pieces of information about the system, and a body containing
all integral values.

``&FCI``
    Starts FCI namelist.

``/``
    Terminates a namelist.  Most compilers also
    implement the extension where ``&END`` is used to
    terminate the namelist instead.

``x``  ``i``  ``a``  ``j``  ``b``
    Format for integral values within body of the FCIDUMP. 
    ``x`` is a float or complex value as appropriate for the system.
    ``i``, ``j``, ``a`` and ``b`` are integers.

&FCI namelist
^^^^^^^^^^^^^

``NORB``
    Number of orbitals in the basis.  See note on basis indices below.
    Must be provided in FCIDUMP namelist.
``NELEC``
    Number of electrons in system.
    Must be provided either in FCIDUMP namelist or input file.
``MS2``
    Spin polarisation.
    Must be provided either in FCIDUMP namelist or input file.
``ORBSYM``
    Array containing symmetry label of each orbital.  See
    symmetry notes below.
    If not provided in FCIDUMP namelist we assume the system has no symmetry.
``UHF``
    True if FCIDUMP file was produced from an unrestricted
    Hartree-Fock calculation.  See note on basis indices below.
    If not provided in FCIDUMP namelist RHF calculation is assumed.

    .. note::

         We assume that in UHF calculations the number of spin-up basis
         functions is equal to the number of spin-down basis functions.

``ISYM``
    Currently unused.  Defined solely for compatibility with NECI
    FCIDUMP files.  Gives the symmetry of the wavefunction formed by
    occupied the NELEC lowest energy spin-orbitals.

``SYML``
    Currently unused.  Defined solely for compatibility with NECI
    FCIDUMP files.  Array containing L (angular momentum) for each orbital.
    Set to :math:`-1` if L is not a good quantum number.

``SYMLZ``
    Array containing :math:`L_z` (angular momentum along the z-axis) for each orbital.
    For example :math:`d_xz` would have :math:`L=2` and :math:`L_z=1`, and
    :math:`d_yz L=2`, :math:`L_z=-1`.
    If not provided in FCIDUMP assume no :math:`L_z` symmetry in system.

``NPROP``
    Currently unused.  Defined solely for compatibility with NECI
    FCIDUMP files. Dimensions of the supercell used in translationally
    symmetric systems.

``PROPBITLEN``
    Currently unused.  Defined solely for compatibility with NECI
    FCIDUMP files. Length in bits of each kpoint index dimension in
    translationally symmetric systems.

Integrals
^^^^^^^^^

if :math:`i = j = a = b = 0`, :math:`E_{core} = x` , where :math:`E_{core}` contains the
nuclear-nuclear and other non-electron contributions to the
Hamiltonian.

if :math:`a = j = b = 0`, :math:`\epsilon_i = x`, the single-particle eigenvalue
of the i-th orbital.

if :math:`j = b = 0`, :math:`< i | h | a > = x`, the one-body Hamiltonian matrix element
between the i-th and a-th orbitals, where :math:`h = T+V_{ext}`.

otherwise :math:`< i j | 1/r_{12} | a b > = x`, the Coulomb integral between
the i-a co-density and the j-b codensity.  Note the Coulomb
integrals are given in Chemists' notation.

Basis indices
-------------
``RHF``
    All indices are in terms of spatial orbitals.  NORB is the
    number of spatial orbitals.

``UHF``
    All indices are in terms of spin orbitals.  NORB is the
    number of spin orbitals.

    .. note::

        Basis functions (as stored by basis_fns) are always stored as spin
        orbitals (the memory saving involved in storing only spatial orbitals
        is not worth the additional overhead/headache, as FCIQMC involves
        working in spin orbitals).  Integrals are expensive to store, so we
        store them in as compressed format as possible.

.. warning::

    The single-particle basis is assumed to be orthonormal.

Symmetry
--------

Molecular orbitals are defined by the D2h point group (or a subgroup
thereof)by the quantum chemistry packages (QChem, MOLPRO) used to
produce FCIDUMP files , so we need only concern ourselves with Abelian
symmetries.

ORBSYM(i) = S+1, where S is the symmetry label defining the
irreducible representation spanned by the i-th orbital.
See notes in pg_symmetry about the symmetry label for Abelian point
groups.

If ORBSYM(i) = 0, then the symmetry of the i-th orbital is not
well-defined.  In this case, we can only resort to turning off all
symmetry (i.e. set all orbitals to be totally symmetric).

.. warning::

    Note that this has memory implications for the integral storage.

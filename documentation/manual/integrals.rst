.. _generating_integrals:

Generating integrals
====================

HANDE can treat :ref:`systems <generic_systems>` other than model Hamiltonians by reading in the necessary
integrals in the FCIDUMP format [Knowles89]_.  Many quantum chemistry packages can
generate them following Hartree-Fock calculations, including:

HORTON
    https://theochem.github.io/horton/
MOLPRO
    https://www.molpro.net
PSI4
    http://psicode.org

    A short tutorial can be found :ref:`below <psi4_tutorial>`.
Q-Chem
   http://www.q-chem.com  FCIDUMP code contributed by Alex Thom.

We most frequently use PSI4 and Q-Chem and so these tend to be better tested.  Note that
the computational cost of the calculations in HANDE vastly outweighs the cost of the
underlying SCF calculations and so the efficiency of the code used to generate the
integrals is usually not a key factor.  Please consult the documentation of the code of
interest regarding how to run SCF calculations and generate the integrals in the FCIDUMP
format.

Please note that not all programs use exactly identical FCIDUMP formats and some may not
be compatible with HANDE. The differences are typically in the namelist header. It may be
possible to resolve these differences by hand or a script.

.. _psi4_tutorial:

A short tutorial for generating FCIDUMPs by Psi4
------------------------------------------------

Installation
^^^^^^^^^^^^
The Psi4 package can be simply installed in an Anaconda environment via

.. code-block:: bash
   
   $ conda install psi4 python=3.8 -c psi4

Psi4 provides two APIs to access it, the Psithon interface (:code:`$ psi4 H2.in > H2.out`), or by using it like a Python module. The latter is covered by this tutorial as it is more intuitive to use and renders post-calculation analysis easier. 

Typical usage
^^^^^^^^^^^^^
The following examples runs a RHF calculation on carbon dimer at cc-pVDZ at a clamped orbital occupancy, producing a FCIDUMP file at the end.

.. code-block:: python

    import psi4

    psi4.core.clean() # Clean local scratch files
    psi4.set_memory('2000 MB')

    dump = 'c2_1.200.FCIDUMP'
    outfile = 'c2_1.200.psi4out'
    psi4.core.set_output_file(outfile, False)
    # There are many ways to specify geometry, see documentation
    mol = psi4.geometry("""
    C
    C 1 r
    r = 1.200
    """)
    # Default SCF_TYPE is density fitting, which produces slightly different results than other packages like PySCF.
    # Provide 'SCF_TYPE':'DIRECT' in the option dictionary if it's concerning
    psi4.set_options({'basis':'cc-pvdz','docc':[2,0,0,0,0,2,1,1]})
    E, wfn = psi4.energy('scf', return_wfn=True)
    # oe_ints has to be specified exactly like this
    psi4.fcidump(wfn, fname=dump, oe_ints=['EIGENVALUES'])

Custom basis set
^^^^^^^^^^^^^^^^
Sometimes a custom basis set is needed (for example many chromium dimer benchmarks are done with an 'Ahlrich's SV' basis, which requires deleting a diffuse p basis from the def2-SV(P) basis set for Cr). This can be done either in the input file (https://psicode.org/psi4manual/master/basissets.html#user-defined-basis-sets) or by writing a basis set definition in the python module folder (which should be :code:`/home/{user}/anaconda3/envs/{env-name}/share/psi4/basis/`, see the link above for basis set file naming convention).

Post-calculation analysis
^^^^^^^^^^^^^^^^^^^^^^^^^
The returned wavefunction object :code:`wfn` can be inspected and analysed. For methods summary see `here <https://psicode.org/psi4manual/master/api/psi4.core.Wavefunction.html#psi4.core.Wavefunction>`_. But most basically you can call

.. code-block:: python

    wfn.epsilon_a().nph

which returns all the HF eigenvalues grouped by symmetry.

Symmetry
^^^^^^^^
Psi4 uses 'Cotton ordering' for the irreps of :math:`D_{2h}`, albeit inconsistently (e.g. the :code:`DOCC` option takes in a list of irrep occupation with normal ordering, i.e., :math:`A_{g},\ B_{1g},\ B_{2g},\dots`). But in the &FCI namelist, the symmetry labels are Cotton-ordered, i.e. :math:`[1,2,3,4,5,6,7,8]` means :math:`[A_{1g},B_{3u},B_{2u},B_{1g},B_{1u},B_{2g},B_{3g},A_u]`.

Freezing orbitals
^^^^^^^^^^^^^^^^^
For large systems, if you're already planning on freezing electrons in the HANDE calculation, it might be sensible to freeze them in the FCIDUMP. 
Psi4 can do this for you by just adding :code:`'freeze_core':True` in the options dictionary above (more precise control is available, see `here <https://psicode.org/psi4manual/1.5.0/autodir_options_c/globals__freeze_core.html>`_), and export FCIDUMP in exactly the same way.

.. _fcidump_format:

FCIDUMP format
--------------

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
    For example :math:`d_{xz}` would have :math:`L=2` and :math:`L_z=1`, and
    :math:`d_{yz} L=2`, :math:`L_z=-1`.
    If not provided in FCIDUMP assume no :math:`L_z` symmetry in system. 
    See :ref:`generic_systems` for more details, and also on how to generate :math:`L_z`-transformed FCIDUMPs.

``NPROP``
    Dimensions of the supercell used in translationally symmetric systems.

``PROPBITLEN``
    Length in bits of each kpoint index dimension in translationally symmetric systems.

Integrals
^^^^^^^^^

if :math:`i = j = a = b = 0`, :math:`E_{\text{core}} = x` , where :math:`E_{\text{core}}` contains the
nuclear-nuclear and other non-electron contributions to the
Hamiltonian.

if :math:`a = j = b = 0`, :math:`\epsilon_i = x`, the single-particle eigenvalue
of the i-th orbital.

if :math:`j = b = 0`, :math:`\langle i | h | a \rangle = x`, the one-body Hamiltonian matrix element
between the i-th and a-th orbitals, where :math:`h = T+V_{\text{ext}}`.

otherwise :math:`\langle i j | 1/r_{12} | a b \rangle = x`, the Coulomb integral between
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

For periodic systems symmetries are defined by their kpoint vector.
ORBSYM(i) contains this vector in a format defined by PROPBITLEN,
which is decoded within HANDE.

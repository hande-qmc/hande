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

The format of the FCIDUMP file expected by HANDE is documented in the comments to
``src/read_in.F90``; these might be useful for users wishing to generate integrals from
upon a specific Hamiltonian to feed into HANDE.

.. warning::

    The single-particle basis is assumed to be orthonormal.

.. _reference_table:

reference options
=================

The ``reference`` table contains options used to control the Hilbert space used in the
calculation and trial function for the projected estimator.

``det``
    type: vector of integers.

    Optional.  Default: a simple (but potentially not optimal) guess which satisfies the spin
    and, if provided, symmetry options using the Aufbau principle.  In most cases the
    default (which for molecules typically corresponds to the Hartree--Fock determinant)
    is sufficient.

    Specify the determinant (as a list of indices corresponding to occupied
    single-particle orbitals) to be used as the reference determinant, which is used in
    the trial function for calculating the projected energy estimator.  Typically this
    should be the determinant expected to have the greatest overlap with the
    desired wavefunction.
``hilbert_space_det``
    type: vector of integers.

    Optional.  Default: set to ``det``.

    Specify the determinant (as a list of indices corresponding to occupied single-particle 
    orbitals) used to generate the Hilbert space.  Using different determinants to control
    the Hilbert space and the trial function allows, for example, spin-flip calculations
    to be performed.

    .. note::

        Only relevant if the Hilbert space is not equivalent to the FCI space, i.e.
        ``ex_level`` is smaller than the number of electrons in the system.

``ex_level``
    type: integer.

    Optional.  Default: set to the number of electrons in the system (i.e. consider all
    determinants in the FCI space).

    Maximum excitation level to consider relative to the determinant given by
    ``hilbert_space_det``.

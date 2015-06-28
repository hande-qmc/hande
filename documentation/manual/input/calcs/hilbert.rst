Monte Carlo estimate of size of the Hilbert space
=================================================

Whilst calculating the size of an entire Hilbert space is straightforward via
combinatorics, calculating the size of a specific part of the Hilbert space meeting
a given set of quantum numbers (e.g. spin and symmetry) is more challenging.  Instead,
the size of this subspace can be estimated via a simple Monte Carlo approach[BoothPhD]_.

.. code-block:: lua

    hilbert_space {
        sys = system,
        hilbert = { ... },
    }

Options
-------

All options should be in the hilbert table bar the sys option.

``sys``
    type: system object.

    Required.

    The system on which to perform the calculation.  Must be created via a system
    function.
``hilbert``
    type: lua table.

    Required.

    Further options to control the Monte Carlo estimation of the Hilbert space.  See
    below.

hilbert options
---------------

The ``hilbert`` table can take the following options:

``ncycles``
    type: integer.

    Required.

    Number of cycles  to perform (i.e. number of random determinants to generate).
``rng_seed``
    type: integer.

    Optional.  Default: generate a seed based upon the time and UUID (if available).

    Seed for initialising the random number generator.
``reference``
    type: vector of integers.

    Optional.  Default: attempt to make a good guess based upon the spin and symmetry
    quantum numbers of the system.

    The reference determinant as a list of occupied spin-orbitals.  The reference
    determinant is used in the generation of truncated Hilbert spaces only.
``ex_level``
    type: integer.

    Optional.  Default: set to the number of electrons in the system (i.e. generate the
    FCI space).

    Maximum excitation level to consider relative to the reference determinant.

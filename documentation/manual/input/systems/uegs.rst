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
    type: integer or string.

    Optional. Default: ``aufbau``.

    Set the symmetry (i.e. crystal momentum) of the system if a reference determinant is
    not provided. This can be set to:

    - An integer specifying the index of a specific wavevector; see the output produced by creating a system for possible
      values and their corresponding wavevectors.
    - ``aufbau``. Uses the symmetry of a determinant selected using the Aufbau principle.
    - ``tot_sym``. Uses the totally symmetric representation, whatever its index may be.
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

    Additional note, when ``cubic_lattice`` is set to body-centerd or
    face-centered the scaled units are modified to account for
    1) a normalization for their respective reciprocal lattice vectors and
    2) include an additional factor of 1/2 to make the sequence behavior
    more consistent. For additional details see ``cubic_lattice`` below.
    For determining the actual kinetic energy defining a given shell please
    refer to the basis functions table reported in the output for a calculation.

    For convenience the first 100 basis sets are documented in
    :ref:`cutoff_table` below.
``dim``
    type: integer.

    Optional.  Default: 3.

    Set the dimension of the electron gas.  2- and 3-dimensional gases are implemented.
``twist``
    type: :math:`N`-dimensional vector

    Optional.  Default: 0 in each dimension.

    Apply a twist to the wavevector grid.  The twist is an *ndim*-dimensional vector in
    terms of the reciprocal lattice vectors.  The twist angle is chosen to be within the
    primitive reciprocal cell, and hence the components should be between -0.5 and +0.5.
``verbose``
    type: boolean.

    Optional.  Default: true.

    Print out the single-particle basis set.
``use_mom_vec``
    type: boolean.

    Optional. Default: false.

    For full configuration interaction calculations only, indicates that the
    Hilbert space should be partitioned according to the determinants momentum
    vector rather than mapping the determinants to the orbital basis. This must
    be true when supplying a ``sym`` outside the basis spatial indicies.
    To determine valid values for ``sym``, use the ``determinant_table_only``
    keyword in the :ref:`fci` table. See the folder 'ueg_n6_rs5_fci' folder in
    the fci test_suite for an example.
``trunc_table``
    type: boolean.

    Optional. Default: false.

    Only used when ``use_mom_vec`` is true. This flag is used to truncate
    the determinant table reported during full configuration interaction
    so that only the reference symmetry is reported. Useful when a system
    has a large determinant table which can bloat the file sizes.
``cubic_lattice``
    type: string.

    Optional. Default: 'simple'.

    Possible values: 'simple', 'body_centered', and 'face_centered'.

    Only used for the three-dimensional jellium. Defines the underlying
    cubic lattice type to simulate. This is achieved by changing the
    reciprocal lattice vectors used to generate the wavevectors/sites
    for the orbital basis.

    Note that for the non-simple lattice types the energy cutoff defining
    the basis will progress in slightly different ways. Normally the
    simple basis proceeds in shells with an energy gap between the
    shells typically being 1/2 or 1 (in the scaled units).
    Whereas for the body-centered lattice this remains mostly true but
    the sequence for these gaps will be different.
    For the face-centered lattice the gaps between shells typically take on
    the values 1/4, 3/4, 1, 5/4, and 7/4.
    As a general rule of thumb for including additional shells in the system
    one should increase the ``cutoff`` in steps of 0.5 for the simple and
    body-centered types and steps of 0.25 for the face-centered type.

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
``verbose``
    type: boolean.

    Optional.  Default: true.

    Print out the single-particle basis set.

.. _cutoff_table:
Three-dimensional UEG basis table
---------------------------------

Example code for generating the results found within the table below
for the simple cubic lattice type.

.. code-block:: lua

    for ecut = 0.0, 100.0, 0.5 do
        fci {
            sys = ueg {
                nel = 1,
                ms = 1,
                sym = 1,
                dim = 3,
                cutoff = ecut + 0.001,
                rs = 1.0,
                use_mom_vec = true,
                cubic_lattice = 'simple',
            },
            fci = {
                determinant_table_only = true,
            },
        }
    end

A table for convenient access to the first 100 basis sets that
can be found in the three-dimensional jellium for various lattice types.
For finding additional basis sets we recommend using the example above with
the appropriate modifications to target the relevant outcome.

.. table::
    A table for the the number of spin orbitals and energy cutoff
    defining a basis set containing the first 100 shells for a given
    three-dimensional jellium lattice type. Note that an offset of 0.001
    is included in the cutoff value to account for differences related to
    numerical precision.

    +------+----------------+----------------+----------------+
    |      |     Simple     | Body-centered  | Face-centered  |
    +======+=======+========+=======+========+=======+========+
    |Shells|   M   | Cutoff |   M   | Cutoff |   M   | Cutoff |
    +------+-------+--------+-------+--------+-------+--------+
    |  1   |   2   | 0.001  |   2   | 0.001  |   2   | 0.001  |
    +------+-------+--------+-------+--------+-------+--------+
    |  2   |  14   | 0.501  |  26   | 0.501  |  18   | 0.751  |
    +------+-------+--------+-------+--------+-------+--------+
    |  3   |  38   | 1.001  |  38   | 1.001  |  30   | 1.001  |
    +------+-------+--------+-------+--------+-------+--------+
    |  4   |  54   | 1.501  |  86   | 1.501  |  54   | 2.001  |
    +------+-------+--------+-------+--------+-------+--------+
    |  5   |  66   | 2.001  |  110  | 2.001  |  102  | 2.751  |
    +------+-------+--------+-------+--------+-------+--------+
    |  6   |  114  | 2.501  |  158  | 2.501  |  118  | 3.001  |
    +------+-------+--------+-------+--------+-------+--------+
    |  7   |  162  | 3.001  |  174  | 3.001  |  130  | 4.001  |
    +------+-------+--------+-------+--------+-------+--------+
    |  8   |  186  | 4.001  |  270  | 3.501  |  178  | 4.751  |
    +------+-------+--------+-------+--------+-------+--------+
    |  9   |  246  | 4.501  |  282  | 4.001  |  226  | 5.001  |
    +------+-------+--------+-------+--------+-------+--------+
    |  10  |  294  | 5.001  |  354  | 4.501  |  274  | 6.001  |
    +------+-------+--------+-------+--------+-------+--------+
    |  11  |  342  | 5.501  |  402  | 5.001  |  338  | 6.751  |
    +------+-------+--------+-------+--------+-------+--------+
    |  12  |  358  | 6.001  |  450  | 5.501  |  362  | 8.001  |
    +------+-------+--------+-------+--------+-------+--------+
    |  13  |  406  | 6.501  |  498  | 6.001  |  458  | 8.751  |
    +------+-------+--------+-------+--------+-------+--------+
    |  14  |  502  | 7.001  |  642  | 6.501  |  518  | 9.001  |
    +------+-------+--------+-------+--------+-------+--------+
    |  15  |  514  | 8.001  |  738  | 7.501  |  566  | 10.001 |
    +------+-------+--------+-------+--------+-------+--------+
    |  16  |  610  | 8.501  |  762  | 8.001  |  614  | 10.751 |
    +------+-------+--------+-------+--------+-------+--------+
    |  17  |  682  | 9.001  |  858  | 8.501  |  662  | 11.001 |
    +------+-------+--------+-------+--------+-------+--------+
    |  18  |  730  | 9.501  |  918  | 9.001  |  678  | 12.001 |
    +------+-------+--------+-------+--------+-------+--------+
    |  19  |  778  | 10.001 | 1062  | 9.501  |  774  | 12.751 |
    +------+-------+--------+-------+--------+-------+--------+
    |  20  |  874  | 10.501 | 1110  | 10.001 |  822  | 13.001 |
    +------+-------+--------+-------+--------+-------+--------+
    |  21  |  922  | 11.001 | 1206  | 10.501 |  918  | 14.001 |
    +------+-------+--------+-------+--------+-------+--------+
    |  22  |  970  | 12.001 | 1254  | 11.001 | 1062  | 14.751 |
    +------+-------+--------+-------+--------+-------+--------+
    |  23  | 1030  | 12.501 | 1350  | 11.501 | 1074  | 16.001 |
    +------+-------+--------+-------+--------+-------+--------+
    |  24  | 1174  | 13.001 | 1354  | 12.001 | 1122  | 16.751 |
    +------+-------+--------+-------+--------+-------+--------+
    |  25  | 1238  | 13.501 | 1534  | 12.501 | 1218  | 17.001 |
    +------+-------+--------+-------+--------+-------+--------+
    |  26  | 1382  | 14.501 | 1582  | 13.001 | 1290  | 18.001 |
    +------+-------+--------+-------+--------+-------+--------+
    |  27  | 1478  | 15.001 | 1774  | 13.501 | 1402  | 18.751 |
    +------+-------+--------+-------+--------+-------+--------+
    |  28  | 1502  | 16.001 | 1870  | 14.001 | 1450  | 19.001 |
    +------+-------+--------+-------+--------+-------+--------+
    |  29  | 1598  | 16.501 | 1918  | 14.501 | 1498  | 20.001 |
    +------+-------+--------+-------+--------+-------+--------+
    |  30  | 1694  | 17.001 | 2110  | 15.501 | 1642  | 20.751 |
    +------+-------+--------+-------+--------+-------+--------+
    |  31  | 1790  | 17.501 | 2122  | 16.001 | 1738  | 21.001 |
    +------+-------+--------+-------+--------+-------+--------+
    |  32  | 1850  | 18.001 | 2278  | 16.501 | 1786  | 22.001 |
    +------+-------+--------+-------+--------+-------+--------+
    |  33  | 1898  | 18.501 | 2338  | 17.001 | 1882  | 22.751 |
    +------+-------+--------+-------+--------+-------+--------+
    |  34  | 2042  | 19.001 | 2362  | 17.501 | 1930  | 24.001 |
    +------+-------+--------+-------+--------+-------+--------+
    |  35  | 2090  | 20.001 | 2578  | 18.001 | 2074  | 24.751 |
    +------+-------+--------+-------+--------+-------+--------+
    |  36  | 2282  | 20.501 | 2818  | 18.501 | 2134  | 25.001 |
    +------+-------+--------+-------+--------+-------+--------+
    |  37  | 2378  | 21.001 | 2866  | 19.001 | 2278  | 26.001 |
    +------+-------+--------+-------+--------+-------+--------+
    |  38  | 2426  | 21.501 | 2962  | 19.501 | 2422  | 26.751 |
    +------+-------+--------+-------+--------+-------+--------+
    |  39  | 2474  | 22.001 | 3010  | 20.001 | 2486  | 27.001 |
    +------+-------+--------+-------+--------+-------+--------+
    |  40  | 2618  | 22.501 | 3106  | 20.501 | 2582  | 28.751 |
    +------+-------+--------+-------+--------+-------+--------+
    |  41  | 2714  | 23.001 | 3202  | 21.001 | 2726  | 29.001 |
    +------+-------+--------+-------+--------+-------+--------+
    |  42  | 2730  | 24.001 | 3406  | 21.501 | 2822  | 30.001 |
    +------+-------+--------+-------+--------+-------+--------+
    |  43  | 2838  | 24.501 | 3418  | 22.001 | 2918  | 30.751 |
    +------+-------+--------+-------+--------+-------+--------+
    |  44  | 3006  | 25.001 | 3586  | 22.501 | 2942  | 32.001 |
    +------+-------+--------+-------+--------+-------+--------+
    |  45  | 3102  | 25.501 | 3706  | 23.501 | 3182  | 32.751 |
    +------+-------+--------+-------+--------+-------+--------+
    |  46  | 3150  | 26.001 | 3718  | 24.001 | 3278  | 33.001 |
    +------+-------+--------+-------+--------+-------+--------+
    |  47  | 3294  | 26.501 | 4186  | 24.501 | 3374  | 34.001 |
    +------+-------+--------+-------+--------+-------+--------+
    |  48  | 3486  | 27.001 | 4246  | 25.001 | 3518  | 34.751 |
    +------+-------+--------+-------+--------+-------+--------+
    |  49  | 3582  | 28.001 | 4342  | 25.501 | 3614  | 35.001 |
    +------+-------+--------+-------+--------+-------+--------+
    |  50  | 3678  | 28.501 | 4486  | 26.001 | 3674  | 36.001 |
    +------+-------+--------+-------+--------+-------+--------+
    |  51  | 3726  | 29.001 | 4630  | 26.501 | 3786  | 36.751 |
    +------+-------+--------+-------+--------+-------+--------+
    |  52  | 3870  | 29.501 | 4682  | 27.001 | 3834  | 37.001 |
    +------+-------+--------+-------+--------+-------+--------+
    |  53  | 4014  | 30.501 | 4898  | 27.501 | 3978  | 38.001 |
    +------+-------+--------+-------+--------+-------+--------+
    |  54  | 4206  | 31.001 | 5018  | 28.501 | 4170  | 38.751 |
    +------+-------+--------+-------+--------+-------+--------+
    |  55  | 4218  | 32.001 | 5090  | 29.001 | 4218  | 40.001 |
    +------+-------+--------+-------+--------+-------+--------+
    |  56  | 4410  | 32.501 | 5234  | 29.501 | 4266  | 40.751 |
    +------+-------+--------+-------+--------+-------+--------+
    |  57  | 4602  | 33.001 | 5330  | 30.001 | 4458  | 41.001 |
    +------+-------+--------+-------+--------+-------+--------+
    |  58  | 4650  | 33.501 | 5426  | 30.501 | 4554  | 42.001 |
    +------+-------+--------+-------+--------+-------+--------+
    |  59  | 4746  | 34.001 | 5642  | 31.501 | 4794  | 42.751 |
    +------+-------+--------+-------+--------+-------+--------+
    |  60  | 4938  | 34.501 | 6110  | 32.001 | 4842  | 43.001 |
    +------+-------+--------+-------+--------+-------+--------+
    |  61  | 5034  | 35.001 | 6206  | 32.501 | 4890  | 44.001 |
    +------+-------+--------+-------+--------+-------+--------+
    |  62  | 5106  | 36.001 | 6302  | 33.001 | 5130  | 44.751 |
    +------+-------+--------+-------+--------+-------+--------+
    |  63  | 5202  | 36.501 | 6602  | 33.501 | 5274  | 45.001 |
    +------+-------+--------+-------+--------+-------+--------+
    |  64  | 5442  | 37.001 | 6662  | 34.001 | 5370  | 46.001 |
    +------+-------+--------+-------+--------+-------+--------+
    |  65  | 5554  | 37.501 | 6782  | 34.501 | 5466  | 46.751 |
    +------+-------+--------+-------+--------+-------+--------+
    |  66  | 5602  | 38.001 | 6878  | 35.001 | 5482  | 48.001 |
    +------+-------+--------+-------+--------+-------+--------+
    |  67  | 5794  | 38.501 | 6902  | 35.501 | 5674  | 48.751 |
    +------+-------+--------+-------+--------+-------+--------+
    |  68  | 5890  | 39.001 | 6926  | 36.001 | 5782  | 49.001 |
    +------+-------+--------+-------+--------+-------+--------+
    |  69  | 5938  | 40.001 | 7238  | 36.501 | 5950  | 50.001 |
    +------+-------+--------+-------+--------+-------+--------+
    |  70  | 6142  | 40.501 | 7286  | 37.001 | 6142  | 50.751 |
    +------+-------+--------+-------+--------+-------+--------+
    |  71  | 6238  | 41.001 | 7490  | 37.501 | 6238  | 51.001 |
    +------+-------+--------+-------+--------+-------+--------+
    |  72  | 6382  | 41.501 | 7562  | 38.001 | 6286  | 52.001 |
    +------+-------+--------+-------+--------+-------+--------+
    |  73  | 6478  | 42.001 | 7682  | 38.501 | 6430  | 52.751 |
    +------+-------+--------+-------+--------+-------+--------+
    |  74  | 6574  | 42.501 | 7802  | 39.501 | 6574  | 53.001 |
    +------+-------+--------+-------+--------+-------+--------+
    |  75  | 6814  | 43.001 | 7850  | 40.001 | 6766  | 54.001 |
    +------+-------+--------+-------+--------+-------+--------+
    |  76  | 6862  | 44.001 | 8606  | 40.501 | 6958  | 54.751 |
    +------+-------+--------+-------+--------+-------+--------+
    |  77  | 7150  | 44.501 | 8762  | 41.001 | 7054  | 56.001 |
    +------+-------+--------+-------+--------+-------+--------+
    |  78  | 7390  | 45.001 | 8930  | 41.501 | 7294  | 56.751 |
    +------+-------+--------+-------+--------+-------+--------+
    |  79  | 7486  | 45.501 | 9026  | 42.001 | 7390  | 57.001 |
    +------+-------+--------+-------+--------+-------+--------+
    |  80  | 7582  | 46.501 | 9242  | 42.501 | 7438  | 58.001 |
    +------+-------+--------+-------+--------+-------+--------+
    |  81  | 7774  | 47.001 | 9254  | 43.001 | 7534  | 58.751 |
    +------+-------+--------+-------+--------+-------+--------+
    |  82  | 7822  | 48.001 | 9470  | 43.501 | 7678  | 59.001 |
    +------+-------+--------+-------+--------+-------+--------+
    |  83  | 7918  | 48.501 | 9518  | 44.001 | 7886  | 60.751 |
    +------+-------+--------+-------+--------+-------+--------+
    |  84  | 8134  | 49.001 | 9674  | 44.501 | 8030  | 61.001 |
    +------+-------+--------+-------+--------+-------+--------+
    |  85  | 8278  | 49.501 | 9746  | 45.001 | 8222  | 62.001 |
    +------+-------+--------+-------+--------+-------+--------+
    |  86  | 8338  | 50.001 | 9962  | 45.501 | 8558  | 62.751 |
    +------+-------+--------+-------+--------+-------+--------+
    |  87  | 8674  | 50.501 | 10058 | 46.001 | 8570  | 64.001 |
    +------+-------+--------+-------+--------+-------+--------+
    |  88  | 8770  | 51.001 | 10274 | 46.501 | 8762  | 64.751 |
    +------+-------+--------+-------+--------+-------+--------+
    |  89  | 8914  | 52.001 | 10298 | 47.501 | 8954  | 65.001 |
    +------+-------+--------+-------+--------+-------+--------+
    |  90  | 9106  | 52.501 | 10302 | 48.001 | 9146  | 66.001 |
    +------+-------+--------+-------+--------+-------+--------+
    |  91  | 9250  | 53.001 | 10602 | 48.501 | 9242  | 66.751 |
    +------+-------+--------+-------+--------+-------+--------+
    |  92  | 9394  | 53.501 | 10638 | 49.001 | 9290  | 67.001 |
    +------+-------+--------+-------+--------+-------+--------+
    |  93  | 9458  | 54.001 | 10806 | 49.501 | 9386  | 68.001 |
    +------+-------+--------+-------+--------+-------+--------+
    |  94  | 9602  | 54.501 | 11730 | 50.001 | 9626  | 68.751 |
    +------+-------+--------+-------+--------+-------+--------+
    |  95  | 9890  | 55.001 | 11874 | 50.501 | 9818  | 69.001 |
    +------+-------+--------+-------+--------+-------+--------+
    |  96  | 10082 | 56.501 | 11970 | 51.001 | 9914  | 70.001 |
    +------+-------+--------+-------+--------+-------+--------+
    |  97  | 10274 | 57.001 | 12306 | 51.501 | 10058 | 70.751 |
    +------+-------+--------+-------+--------+-------+--------+
    |  98  | 10370 | 57.501 | 12354 | 52.001 | 10130 | 72.001 |
    +------+-------+--------+-------+--------+-------+--------+
    |  99  | 10514 | 58.001 | 12474 | 52.501 | 10322 | 72.751 |
    +------+-------+--------+-------+--------+-------+--------+
    | 100  | 10754 | 58.501 | 12618 | 53.001 | 10418 | 73.001 |
    +------+-------+--------+-------+--------+-------+--------+

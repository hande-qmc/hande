.. _restart_table:

restart options
===============

The ``restart`` table contains options relating to checkpointing within QMC calculations.

HANDE currently uses one restart file per MPI rank with a filename of the form
``HANDE.RS.X.pY.H5``, where ``X`` is the restart index and ``Y`` is the MPI rank.

``read``
    type: boolean or integer.

    Optional.  Default: false.

    Start a QMC calculation from a previous calculation if ``true`` or an integer.  If
    ``true``, then the highest value of ``X`` is used for which a set of restart files
    exists, otherwise specifies the value of ``X`` to use.
``write``
    type: boolean or integer.

    Optional.  Default: false.

    Write out checkpointing files at the end of the calculation if ``true`` or an
    integer.  If ``true``, then the highest value of ``X`` is used for which a set of
    restart files doesn't exist, otherwise specifies the value of ``X`` to use.
``write_shift``
    type: boolean or integer.

    Optional.  Default: false.

    Write out checkpointing files when the shift is allowed to vary (i.e. once
    ``target_population`` is reached) if ``true`` or an integer.  If ``true``, then the
    highest value of ``X`` is used for which a set of restart files doesn't exist,
    otherwise specifies the value of ``X`` to use.
``write_frequency``
    type: integer.

    Optional: :math:`2^{31}-1`.

    Write out checkpointing files every `N` iterations, where `N` is the
    specified value.

    .. note::

        The index used for the restart files created with this option is the next
        unused index.  Depending upon the frequency used, a large number of restart files
        may be created.  As such, this option is typically only relevant for debugging or
        explicitly examining the evolution of the stochastic representation of the
        wavefunction.

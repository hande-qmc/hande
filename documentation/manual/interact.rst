Interacting with running calculations
=====================================

It is possible to interact with running calculations.

After each report loop, HANDE checks for the existence of the file HANDE.COMM in the
current working directory for all processors. If HANDE.COMM exists, then the file is read
and any modified parameters are then used for the rest of the calculation.  HANDE.COMM is
deleted after it is read in to prevent it from being detected on subsequent report loops
and to enable multiple interactions with a running calculation.

HANDE.COMM is a lua script, in a similar fashion to the input file, but has a much more
restricted range of options.  Options which can be set or modified are:

softexit
    type: boolean.

    End the calculation immediately but still perform any post-processing (e.g. dumping
    out a restart file).  This is useful for cleanly terminating a converged calculation
    or cleanly stopping a calculation before the walltime is reached to allow it to be
    restarted.

    The send_softexit.py script in the tools subdirectory is useful for running
    HANDE on a queueing system as it writes **softexit = true** to HANDE.COMM a certain amount
    of time before the walltime is reached.
tau
    type: float.

    Change the timestep to be used.
varyshift_target
    type: integer.

    Change the number of particles to be reached before the calculation starts varying the
    shift.  Meaningless if the calculation has already started varying the shift.  If
    negative then the shift is immediately allowed to vary.
shift
    type: float or 1D vector of floats.

    Adjust the current value of the shift.  If the calculation has already entered
    variable shift mode then the shift will still be updated every report cycle, otherwise
    this is equivalent to changing the **initial_shift** value.

    Passing a single value such as:

    .. code-block:: lua

        shift = -1

    sets the shift in **all** spaces to the specified value.  Different spaces can be
    modified separately by passing in a vector.  For example:

    .. code-block:: lua

        shift = { -1, -2 }

    sets the shift in the first space to -1, in the second space to -2 and leaves it
    unmodified in all other spaces.

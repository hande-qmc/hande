Interacting with running calculations
=====================================

It is possible to interact with running calculations.

After each report loop, HANDE checks for the existence of the file HANDE.COMM in the
current working directory for all processors. If HANDE.COMM exists, then the file is read
and any modified parameters are then used for the rest of the calculation.  HANDE.COMM is
deleted after it is read in to prevent it from being detected on subsequent report loops
and to enable multiple interactions with a running calculation.

HANDE.COMM is a lua script, in a similar fashion to the input file, but has a much more
restricted range of options.  Settings which can be changed are:

**softexit**
    Boolean.

    End the calculation immediately but still perform any post-processing (e.g. dumping
    out a restart file).  This is useful for cleanly terminating a converged calculation
    or cleanly stopping a calculation before the walltime is reached to allow it to be
    restarted.

    The send_softexit.py script  scripts in the tools subdirectory are useful for running
    HANDE on a queueing system as they write **softexit** to HANDE.COMM a certain amount
    of time before the walltime is reached.
**tau**
    Real.

    Change the timestep to be used.
**varyshift_target**
    Integer.

    Change the number of particles to be reached before the calculation starts varying the
    shift.  Meaningless if the calculation has already started varying the shift.  If
    negative then the shift is immediately allowed to vary.
**shift**
    Array of Reals.

    Adjust the current value of the shift.  If the calculation has already entered
    variable shift mode then the shift will still be updated every report cycle, otherwise
    this is equivalent to changing the **initial_shift** value.
    You should use the syntax
    e.g. shift = { -1 } to set the value of the shift to -1.  If additional spaces are being used
    their shifts may be set by having more values in the array.

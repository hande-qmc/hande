Interacting with running calculations
====================================

It is possible to interact with running calculations.

After each update cycle, HANDE checks for the existence of the file
HANDE.COMM in the current working directory for all processors. If HANDE.COMM
exists, then the file is read and any modified parameters are then used for the
rest of the calculation.  HANDE.COMM is deleted after it is read in to prevent
it from being detected on subsequent update cycles and to enable multiple
interactions with a running calculation.

HANDE.COMM has the same syntax as the input file.  Available options are:

**softexit**
    End the calculation immediately but still perform any
    post-processing (e.g. dumping out a restart file).  This is useful for
    cleanly terminating a converged calculation or cleanly stopping
    a calculation before the walltime is reached to allow it to be restarted.

    The watchdog.py (for PBS queue systems) and send_softexit.py (for other
    queue systems) scripts in the tools subdirectory are useful for running
    HANDE on a queueing system as they write **softexit** to HANDE.COMM a
    certain amount of time before the walltime is reached.
**varyshift_target** *varyshift_target*
    Long integer.

    Change the number of particles to be reached before the calculation starts
    varying the shift.  Meaningless if the calculation has already started
    varying the shift.  If *varyshift_target* is negative then the shift is
    immediately allowed to vary.
**tau** *tau*
    Real.

    Change the timestep to be used.
**zero_means**
    Reset the running averages of the shift and projected energy to 0.
**shift** *shift*
    Real.

    Adjust the current value of the shift.  Please note the impact this has on
    the mean; if used it is not a bad idea to also use **zero_means**.  If the
    calculation has already entered variable shift mode then the shift will
    still be updated every report cycle, otherwise this is equivalent to
    changing the **initial_shift** value.

.. _blocking_table:

blocking options
================

The ``blocking`` table contains options used to control the options for performing
blocking analysis on the fly.

``blocking_on_the_fly``
    type: boolean.

    Optional. Default: false

    If true, the data for blocking analysis is collected every report loop and blocking
    analysis is performed on the fly while the calculation is running. At the end of the calculation,
    Estimated correlation energy and its error together with reference energy is printed to the HANDE
    output file. Zero is printed if insufficient data are collected for the blocking analysis.

``start_save_frequency``
    type: integer.

    Optional. Default: -1

    Log2 of the frequency at which the potential start points of the blocking analysis is
    saved. When negative, the frequency is the nearest integer to the log2(``nreports``) - 8.

``start_point_number``
    type: integer.

    Optional. Default: -1

    Number of potential start points of the blocking analysis that is to be saved. If
    negative, the integer part of ``nreports``/2^(``start_save_frequency``).

``filename``
    type: string.

    Optional. Default: 'BLOCKING'

    Filename to which the blocking analysis report is written.

``start_point``
    type: integer.

    Optional. Default: -1

    The iteration number from which the data for blocking analysis is collected. When
    negative the data is collected when ``target_population`` is reached.

``error_limit``
    type: real.

    Optinal. Default: 0

    One of two conditions for termination of the calculation together with ``blocks_used``.
    This specifies the upper limit of the sum of standard error and the error in error of projected energy.

``min_blocks_used``
    type: real.

    Optional. Default: 10

    The minimum number of optimal reblock lengths required for a calculation to
    terminate. The calculation will not terminate due to the standard error
    falling below ``error_limit`` until at least this number of optimal
    reblock lengths are included within the calculation. This ensures that
    our error estimate is reliable at termination.
    Larger ``min_blocks_used`` ensures a more reliable blocking analysis but
    increases the minimum length of calculations.

``blocks_used``
    type: real.

    Optional. Default: (huge)

    Irrelevant of the error_limit, if the number of blocks used to estimate the standard error of projected energy
    is more than the ``blocks_used``, the calculation is terminated. Larger ``blocks_used`` ensures a more reliable
    blocking analysis.

``auto_shift_damping``
    type: boolean.

    Optional. Default: false

    Whether to automatically optimise the shift damping using information from blocking on the fly. This optimises
    the shift damping to ensure that the standard deviations of the instantaneous projected energy and shift are
    approximately equal. The allowable range of values is currently set to allow the shift standard deviation to
    be between 50% and 200% of that of the instantaneous projected energy, though this could easily be exposed to
    the user if required.

..
    [review] - VAN: What happens if I do not provide a shift damping value in the qmc table but there is also none
    in the restart file due to legacy?
    [reply] - CJCS: The shift damping defaults to the original default, 0.05.

    .. note::
        This approach will modify the shift damping to ensure a reasonable variation in the shift during a calculation.
        Updates to the shift damping will be printed within the output file, and the final shift damping written into
        restart files to be used in any restarted calculations. If no shift_damping is provided to a restarted
        calculation in the qmc table the final value from the restarted calculation will be used.

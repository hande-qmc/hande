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

    One of two conditions for termination of the calculation together with ``inverse_frational_error``.
    This specifies the upper limit of the sum of standard error and the error in error of projected energy.

``min_blocks_used``
    type: real.

    Optional. Default: 10

..
    [review] - AJWT: It isn't clear to me what exactly this does or how it works.

    One of two conditions for termination of the calculation together with ``error_limit``.
    This specifies the lower limit of the number of blocks used for blocking analysis.
    Larger ``min_blocks_used`` ensures a more reliable blocking analysis.

``blocks_used``
    type: real.

    Optional. Default: (huge)

    Irrelevant of the error_limit, if the number of blocks used to estimate the standard error of projected energy
    is less than the ``blocks_used``, the calculation is terminated. Larger ``blocks_used`` ensures a more reliable
    blocking analysis.

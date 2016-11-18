.. _blocking_table:

blocking options
================

The ``blocking`` table contains options used to control the options for performing
blocking analysis on the fly.

``blocking_on_the_fly``
    type: boolean.

    Optional. Default: false

    If true, the data for blocking analysis is collected every report loop and blocking
    analysis is performed on the fly while the calculation is running.

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

    If the sum of error in error and standard deviation of projected energy is smaller than
    this value and the ``min_ratio`` condition is satisfied, ``soft_exit`` = true is returned
    and the calculation is terminated.

``min_ratio``
    type: integer.

    Optional. Default: 3

..
    [review] - AJWT: It isn't clear to me what exactly this does or how it works.


    The ratio between error in error and standard error of projected energy.
    If the ratio is larger, greater number of blocks are used for reblock analysis. If the ``error_limit``
    and ``min_ratio`` condition is satisfied, ``soft_exit`` = true is returned and calculation is terminated

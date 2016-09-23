.. _logging_table:

logging options
===============

The ``logging`` table contains options relating to outputting additional logs from QMC calculations.

Use of this functionality requires compiling HANDE with debug flags (using the -g option with mkconfig).
This enables implementation of logging without having an appreciable impact upon timings of an optimised
build.

This functionality is recommended for developers only. It should allow easy identification of
the causes of any changes in Markov chain between two calculations.

Additional logging functionality can be added upon request. Current coverage is by no means complete.

``calc``
    type: integer

    Optional. Default: ``0``.

    Determines level of logging output related to high-level behaviour within a calculation.
    Currently implemented levels are:

    - ``0`` returns no extra information.
    - ``1`` returns summary of events within a calculation (currently only for FCIQMC and CCMC).

    Any information is produced in the file CALC.log within the working directory.

``spawn``
    type: integer

    Optional. Default: 0.

    Determines level of logging output related to spawning within a calculation. Current levels are:

    - ``0`` returns no extra information.
    - ``1`` returns information on each spawning event creating at least one particle within a
        calculation (currently only for FCIQMC and generic systems).
    - ``2`` returns information on each spawning event within a calculation, regardless of result
        (currently only for FCIQMC and generic systems).

    Any information is produced in the file SPAWN.log within the working directory.


``death``
    type: integer

    Optional. Default: 0.

    Determines level of logging output related to death within a calculation. Current levels are:

    - ``0`` returns no extra information.
    - ``1`` returns information on each death or cloning event resulting in a change in particle number
        within a calculation (currently only for FCIQMC and generic systems).
    - ``2`` returns information on each death or cloning event within a calculation, regardless of result
        (currently only for FCIQMC and generic systems).

    Any information is produced in the file DEATH.log within the working directory.

``start``
    type: integer

    Optional. Default: 0.

    Defines the iteration from which logging information should be produced.

``finish``
    type: integer

    Optional. Default: :math:`2^{31}-1`.

    Defines the iteration after which logging information should cease to be produced.

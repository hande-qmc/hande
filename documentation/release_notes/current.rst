Current development
===================

Added
-----

* Ability to restart the state of the dSFMT RNG stream, allowing for restarted
  calculations to have the same Markov chain as single calculations. Enabled by default.
  Can be disabled in the :ref:`restart_table`. 

Changed
-------

* ``write_frequency`` now is in units of report loops rather than Monte Carlo cycles.

Removed
-------

n/a

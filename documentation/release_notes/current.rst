Current release notes after v1.2

Bug Fixes
----------

* Even selection weighting had an initialization bug affecting any calculation with more than one MPI process where one process has no single excitations.
  The effect is that for those cycles without singles, the even_selection probability for the singles will be uninitialized, and possibly not zero.
  For a sufficiently large calculation there will be a low probability that this is the case, and if there are ever single excitations on the processor,
  that number will be used, and will continue to be used if there are no singles at some later date.  It is expected that there will be no notable effect
  after equlibration.
  Incorrect even selection weightings only affect the efficiency of the selection, and will not in general introduce a bias.
  Effects on systems without single excitations (e.g. UEG and Hubbard models) are undefined. 
  The bug-fix changes Markov chains.

Scriptable Input
================

Having lua control a HANDE simulation allows for some pretty clean ways to run complicated
simulations. Here we will list some examples.

Twist Averaging
^^^^^^^^^^^^^^^

To aid in the removal of single-particle finite size effects it is often helpful to
perform twist averaging. Here we want to average results over multiple twist vectors
:math:`\mathbf{k}_s` where each component of :math:`\mathbf{k}_s` can be chosen to lie
within the simulation cell Brillouin zone. Normally we would need to run multiple
independent simulations yielding many output files, which can be problematic for file
systems. Lua allows us to run the calculations from a single input file.

In the example below we show this for a twist averaged canonical total energy calculation,
which can be useful for correcting incomplete twist averaged QMC calculations which are
typically much more expensive.

.. literalinclude:: examples/canonical_twist.lua
    :language: lua

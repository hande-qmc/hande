Tips
====

Some suggestions from the HANDE developers for using HANDE...heed our words!

Compilation
-----------

For optimised versions of HANDE, explore using:

* compiler-specific optimisation flags

  In general adding 'high-level' optimisation flags (``-O3``, ``-Ofast``, etc.) makes
  a substantial impact on the calculation speed.

* interprocedural optimisation

  Many compilers can perform interprocedural optimisation, whereby optimisations are
  performed at link-time instead of compile-time.  This allows optimisations to be
  performed (including inlining) on procedures specified in different source files.  On
  some compilers (e.g. GCC, Intel) this can have a substantial benefit; on other compilers
  the difference is less marked.

* popcnt instruction

  If the processor being used includes it, uses the popcnt instruction rather than
  a software implementation to count bits set in an integer.  This can have a impact of
  the order of a few percent for the entire calculation.

* DET_SIZE=64

  Use 64-bit integers rather than 32-bit integers to store the representation of the
  determinant/excitor/tensor labels.  This can make certain calculations quicker (i.e.
  those involving more than 32 single-particle basis functions) by reducing the amount of
  bit operations that need to be performed.

Plotting calculation output using gnuplot
-----------------------------------------

The first section of the output file contains information about the basis functions
used in the calculations. This gives spurious data points when the contents of the file
is plotted using gnuplot. They can be removed by creating an executable file gphande
in the path, containing:

.. code-block:: bash

    #!/bin/sed -nf
    1,/iterations/d
    /^ *[0-9]/p

When plotting in gnuplot, using the command

    plot '<gphande file'

instead of

    plot 'file'

will then remove the extra points.

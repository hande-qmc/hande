Tips
====

Compilation
-----------

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

pyblock
=======

`pyblock` is a python module for performing a reblocking analysis on
serially-correlated data.

The algorithms implemented in `pyblock` are not new; please see the documentation for
references.

pyblock is compatible with (and tested on!) python 2.7 and python 3.3.

Documentation
-------------

Documentation and a simple tutorial can be found in the docs subdirectory and on
`readthedocs <http://pyblock.readthedocs.org>`_.

Installation
------------

`pyblock` can be used simply by adding to `$PYTHONPATH`.  Alternatively, it can be
installed using distutils: 

::

    $ pip install pyblock

or from PyPI:

::

    $ pip install pyblock

`pyblock` requires numpy and (optionally) pandas and matplotlib.  Please see the
documentation for more details.

License
-------

Modified BSD license; see LICENSE for more details.

Please cite `pyblock, James Spencer, http://github.com/jsspencer/pyblock` if used to
analyse data for an academic publication.

Author
------

James Spencer, Imperial College London

Acknowledgments
---------------

Will Vigor pointed out and wrote an early implementation of the algorithm to detect the
optimal reblock length.  Comments and suggestions from the HANDE development team.

pyblock
=======

`pyblock` is a python module for performing a reblocking analysis on
serially-correlated data.

The algorithms implemented in `pyblock` are not new; please see the documentation for
references.

pyblock is compatible with (and tested on!) python 2.7 and python 3.3-3.4 and should work
on any other version supported by `pandas`.

.. image:: https://travis-ci.org/jsspencer/pyblock.svg?branch=master
    :target: https://travis-ci.org/jsspencer/pyblock

Documentation
-------------

Documentation and a simple tutorial can be found in the docs subdirectory and on
`readthedocs <http://pyblock.readthedocs.org>`_.

Installation
------------

`pyblock` can be used simply by adding to `$PYTHONPATH`.  Alternatively, it can be
installed using distutils by running:

::

    $ pip install /path/to/pyblock

where `/path/to/` is the (relative or absolute) path to the directory containing
`pyblock`.  To install an editable version (useful for development work) do:

::

    $ pip install -e /path/to/pyblock

`pyblock` can also be installed from PyPI:

::

    $ pip install pyblock

`pyblock` requires numpy and (optionally) pandas and matplotlib.  Please see the
documentation for more details.

License
-------

Modified BSD license; see LICENSE for more details.

Please cite ``pyblock, James Spencer, http://github.com/jsspencer/pyblock`` if used to
analyse data for an academic publication.

Author
------

James Spencer, Imperial College London

Contributing
------------

Contributions are extremely welcome, either by raising an issue or contributing code.
For code contributions, please try to follow the following points:

#. Divide commits into logical units (e.g. don't mix feature development with
   refactoring).
#. Ensure all existing tests pass.
#. Create tests for new functionality.  I aim for complete test coverage.
   (Currently the only function not tested is one that creates plots.)
#. Write nice git commit messages (see `Tim Pope's advice <http://tbaggery.com/2008/04/19/a-note-about-git-commit-messages.html>`_.)
#. Send a pull request!

Acknowledgments
---------------

Will Vigor (Imperial College London) pointed out and wrote an early implementation of
the algorithm to detect the optimal reblock length.

Tom Poole (Imperial College London) contributed code to handle weighted averages.

The HANDE FCIQMC/CCMC development team made several helpful comments and suggestions.

pyhande
=======

`pyhande` is a python module for performing a reblocking analysis on
serially-correlated data.

`pyhande` is compatible with (and tested on!) python 2.7 and python 3.3-3.4 and should
work with any other version supported by `pandas`.

`HANDE` comes with several command-line utilities which provide a wrapper around common
tasks requiring `pyhande`.

Documentation
-------------

API documentation can be found in the main `HANDE` documentation.

Installation
------------

`pyhande` can be used simply by adding to `$PYTHONPATH`.  Alternatively, it can be
installed using distutils:

::

    $ pip install /path/to/pyhande

where `/path/to/pyhande` is the path to the top-level pyhande directory (e.g.
`tools/pyhande`) and can be absolute or relative.

Developers may wish to install an 'editable' version, such that changes to the `pyhande`
source are immediately available in the installed module.  This can be accomplished using

::

    $ pip install -e /path/to/pyhande

As usual with python packages, installing the module is best done inside a
`virtualenv <https://virtualenv.readthedocs.org/en/latest/>`_.

`pyhande` requires numpy, pandas and `pyblock`.  If `pyhande` is used directly from the
HANDE repository, then it will automatically pick up `pyblock`.

License
-------

Modified BSD license; see LICENSE for more details.

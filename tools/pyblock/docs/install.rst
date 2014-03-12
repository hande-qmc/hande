Installation
============

Dependencies
------------

* numpy
* pandas (0.13 and later)
* matplotlib

pandas is only required for :mod:`pyblock.pd_utils` and
:mod:`pyblock.error` and matplotlib for :mod:`pyblock.pd_utils`.  Hence pandas and/or
matplotlib need not be installed if those submodules are not required, in which case
``pyblock/__init__.py`` must be modified to stop the :mod:`pyblock.pd_utils` and
:mod:`pyblock.error` from being automatically imported.

Installation instructions
-------------------------

pyblock can be installed from PyPI:

.. code-block:: bash

    $ pip install pyblock

or from the source package:

.. code-block:: bash

    $ python setup.py install

Both ``pip`` and ``setup.py`` have options for installing in non-default locations, such
as home directories.  Add ``--help`` to the above commands for details.

Alternatively, pyblock can be used directly from source by adding the location of the
pyblock directory to the ``PYTHONPATH`` environment variable.

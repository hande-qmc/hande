.. pyblock documentation master file, created by
   sphinx-quickstart on Mon Feb 24 18:28:48 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

pyblock
=======

pyblock implements the reblocking analysis (see, for example, the
description given by Flyvbjerg and Petersen [1]_), to remove serial correlation
from a data set and hence obtain an improved estimate of the standard error.
Functions for additional analysis, interpretation and manipulation of the
resultant mean and standard error estimates are also provided.

A command-line interface is currently not provided but the :doc:`API <api>` is
simple to use from either a Python/IPython shell or to create an
application-specific script.

.. toctree::
   :maxdepth: 1

   install
   api
   tutorial

References
----------
.. [1]  "Error estimates on averages of correlated data", H. Flyvbjerg and
   H.G. Petersen, J. Chem. Phys. 91, 461 (1989).

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


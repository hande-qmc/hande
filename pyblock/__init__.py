'''pyblock is a python module for analysis of correlated data.

.. toctree::
   :maxdepth: 1
   :glob:

   pyblock.*

:mod:`pyblock.blocking` implements the reblocking algorithm [1]_ and an
algorithm [2]_, [3]_ for suggesting the most appropriate block size (and thus
estimate of the standard error in the data set) for data contained within
:mod:`numpy` arrays.  :mod:`pyblock.pd_utils` provides a nice wrapper around
this using :mod:`pandas`, and it is highly recommended to use this if possible.

:mod:`pyblock.error` contains functions for simple error propagation and
formatting of output of a value and it's associated error.

References
----------
.. [1]  "Error estimates on averages of correlated data", H. Flyvbjerg and
   H.G. Petersen, J. Chem. Phys. 91, 461 (1989).
.. [2] "Monte Carlo errors with less errors", U. Wolff, Comput. Phys. Commun.
       156, 143 (2004) and arXiv:hep-lat/0306017.
.. [3] "Strategies for improving the efficiency of quantum Monte Carlo
       calculations", R. M. Lee, G. J. Conduit, N. Nemec, P. Lopez Rios, and N.
       D.  Drummond, Phys. Rev. E. 83, 066706 (2011).
'''

# copyright: (c) 2014 James Spencer
# license: modified BSD license; see LICENSE for further details.

import warnings
# For convenience, import all submodules so the user need only import pyblock.
import pyblock.error
import pyblock.blocking
import pyblock.pd_utils
try:
    import pyblock.plot
except ImportError:
    warnings.warn('Plotting disabled: matplotlib not available.')

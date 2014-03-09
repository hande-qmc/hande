'''Analysis helper routines for HANDE data.'''

from os import path
import pandas as pd
import sys

sys.path.append(path.join(path.abspath(path.dirname(__file__)), '../pyblock'))
import pyblock.error
import pyblock.pd_utils

def projected_energy(reblock_data, covariance, data_length):
    '''Calculate the projected energy estimator and associated error.

The projected energy estimator is given by

.. math::

    E = \\frac{\\sum_H_0j N_j}{N_0}

The numerator and denominator are correlated and so their covariance must be
taken into account.

Parameters
----------
reblock_data : :class:`pandas.DataFrame`
    reblock data for (at least) the numerator and denominator in the
    projected energy estimator.
covariance: :class:`pandas.DataFrame`
    covariance at each reblock iteration between (at least) the numerator
    and denominator in the projected energy estimator.
data_length: :class:`pandas.DataFrame`
    number of data points in each reblock iteration.

Returns
-------
proje : :class:`pandas.DataFrame`
    The projected energy estimator at each reblock iteration.

See also
--------
:func:`pyblock.pd_utils.reblock` for producing the input parameters.
'''

    proje_sum = reblock_data.ix[:, '\sum H_0j Nj']
    ref_pop = reblock_data.ix[:, '# D0']
    proje_ref_cov = covariance.xs('# D0', level=1)['\sum H_0j Nj']
    proje = pyblock.error.ratio(proje_sum, ref_pop, proje_ref_cov, data_length)
    proje.columns = pd.MultiIndex.from_tuples(
            [('Proj. Energy', col) for col in proje.columns])
    return proje

def qmc_summary(data, keys=('Instant shift', '\sum H_0j Nj', '# D0',
                            'Proj. Energy')):

    '''Summarise a reblocked data set by the optimal block.

Parameters
----------
data : :class:`pandas.DataFrame`
    reblocked data (i.e. data with the reblock iteration as the index).
keys : list of strings
    columns (by top-level index) of the data table to inspect.  Each top-level
    column must contain an optimal block column.

Returns
-------
opt_data : :class:`pandas.DataFrame`
    Data for each column from the optimal block size of that column.
no_opt : list of strings
    list of columns for which no optimal block size was found.
'''

    opt_data = []
    no_opt = []
    for col in keys:
        if col in data:
            summary = pyblock.pd_utils.reblock_summary(data.ix[:, col])
            if summary.empty:
                no_opt.append(col)
            else:
                summary.index = [col]
            opt_data.append(summary)
    opt_data = pd.concat(opt_data)
    return (opt_data, no_opt)

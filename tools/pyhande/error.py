'''Simple error propogation.

.. note::

    We only implement the functions as we need them...

'''

import numpy
import pandas as pd
import pyhande.pd_utils as pd_utils

def ratio(stats_A, stats_B, cov_AB, data_len):
    '''Calculate the mean and standard error of f = A/B.

Parameters
----------
stats_A, stats_B: pandas.Series or pandas.DataFrame
    Statistics (containing at least the 'mean' and 'standard error' fields) for
    variables A and B respectively.  The rows contain different values of these
    statistics (e.g. from a reblocking analysis) if DataFrame are passed.
cov_AB: float or pandas.Series
    Covariance between variables A and B.  If stats_A and stats_B are
    DataFrames, then this must be a Series, with the same index as stats_A and
    stats_B.
data_len: int or pandas.Series
    Number of data points ('observations') used to obtain the statistics given
    in stats_A and stats_B.  If stats_A and stats_B are DataFrames, then this
    must be a Series, with the same index as stats_A and stats_B.

Returns
-------
pandas.Series or pandas.DataFrame
    Mean and standard error (and, if possible/relevant, optimal reblock
    iteration) for f = A/B.  If stats_A, stats_B are DataFrames, this is
    a DataFrame with the same index, otherwise a Series is returned.
'''

    return _quadratic(stats_A, stats_B, cov_AB, data_len, -1)

def product(stats_A, stats_B, cov_AB, data_len):
    '''Calculate the mean and standard error of f = A*B.

Parameters
----------
stats_A, stats_B: pandas.Series or pandas.DataFrame
    Statistics (containing at least the 'mean' and 'standard error' fields) for
    variables A and B respectively.  The rows contain different values of these
    statistics (e.g. from a reblocking analysis) if DataFrame are passed.
cov_AB: float or pandas.Series
    Covariance between variables A and B.  If stats_A and stats_B are
    DataFrames, then this must be a Series, with the same index as stats_A and
    stats_B.
data_len: int or pandas.Series
    Number of data points ('observations') used to obtain the statistics given
    in stats_A and stats_B.  If stats_A and stats_B are DataFrames, then this
    must be a Series, with the same index as stats_A and stats_B.

Returns
-------
pandas.Series or pandas.DataFrame
    Mean and standard error (and, if possible/relevant, optimal reblock
    iteration) for f = A*B.  If stats_A, stats_B are DataFrames, this is
    a DataFrame with the same index, otherwise a Series is returned.
'''

    return _quadratic(stats_A, stats_B, cov_AB, data_len, 1)

def _quadratic(stats_A, stats_B, cov_AB, data_len, sign):
    '''Calculate the mean and standard error of f = g(A,B).

Parameters
----------
stats_A, stats_B: pandas.Series or pandas.DataFrame
    Statistics (containing at least the 'mean' and 'standard error' fields) for
    variables A and B respectively.  The rows contain different values of these
    statistics (e.g. from a reblocking analysis) if DataFrame are passed.
cov_AB: float or pandas.Series
    Covariance between variables A and B.  If stats_A and stats_B are
    DataFrames, then this must be a Series, with the same index as stats_A and
    stats_B.
data_len: int or pandas.Series
    Number of data points ('observations') used to obtain the statistics given
    in stats_A and stats_B.  If stats_A and stats_B are DataFrames, then this
    must be a Series, with the same index as stats_A and stats_B.
sign: int
    f = A*B for sign=1 and f= A/B for sign=-1.

Returns
-------
pandas.Series or pandas.DataFrame
    Mean and standard error (and, if possible/relevant, optimal reblock
    iteration) for f = A*B.  If stats_A, stats_B are DataFrames, this is
    a DataFrame with the same index, otherwise a Series is returned.

Raises
------
ValueError
    sign is not +1 nor -1.
'''

    if sign not in (1,-1):
        raise ValueError('sign must be either +1 or -1')
    
    mean = stats_A['mean'] / stats_B['mean']

    std_err = mean*numpy.sqrt(
                (stats_A['standard error']/stats_A['mean'])**2 +
                (stats_B['standard error']/stats_B['mean'])**2 +
                2*sign*cov_AB/(data_len*stats_A['mean']*stats_B['mean'])
            )
    std_err = abs(std_err)

    stats_dict = dict( (('mean', mean), ('standard error', std_err)) )
    try:
        stats = pd.DataFrame(stats_dict)
        if all('optimal block' in stat.columns for stat in (stats_A, stats_B)):
            opt_A = pd_utils.optimal_block(stats_A)
            opt_B = pd_utils.optimal_block(stats_B)
            if opt_A > opt_B:
                stats['optimal block'] = stats_A['optimal block']
            else:
                stats['optimal block'] = stats_B['optimal block']
    except ValueError:
        # Was given a single data point rather than a set.
        stats = pd.Series(stats_dict)

    return stats

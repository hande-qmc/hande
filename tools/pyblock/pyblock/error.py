'''Simple error propogation.

.. note::

    We only implement the functions as we need them...

'''

# copyright: (c) 2014 James Spencer
# license: modified BSD license; see LICENSE for further details.

import numpy
import pandas as pd
import pyblock.pd_utils as pd_utils

def ratio(stats_A, stats_B, cov_AB, data_len):
    '''Calculate the mean and standard error of :math:`f = A/B`.

Parameters
----------
stats_A : :class:`pandas.Series` or :class:`pandas.DataFrame`
    Statistics (containing at least the 'mean' and 'standard error' fields) for
    variable :math:`A`.  The rows contain different values of these statistics
    (e.g. from a reblocking analysis) if :class:`pandas.DataFrame` are passed.
stats_B : :class:`pandas.Series` or :class:`pandas.DataFrame`
    Similarly for variable :math:`B`.
cov_AB : float or :class:`pandas.Series`
    Covariance between variables :math:`A` and :math:`B`.  If ``stats_A`` and
    ``stats_B`` are :class:`pandas.DataFrame`, then this must be
    a :class:`pandas.Series`, with the same index as ``stats_A`` and
    ``stats_B``.
data_len : int or :class:`pandas.Series`
    Number of data points ('observations') used to obtain the statistics given
    in ``stats_A`` and ``stats_B``.  If ``stats_A`` and ``stats_B`` are
    :class:`pandas.DataFrame`, then this must be a :class:`pandas.Series`, with
    the same index as ``stats_A`` and ``stats_B``.

Returns
-------
ratio_stats : :class:`pandas.Series` or :class:`pandas.DataFrame`
    Mean and standard error (and, if possible/relevant, optimal reblock
    iteration) for :math:`f = A/B`.  If ``stats_A``, ``stats_B`` are
    :class:`pandas.DataFrame`, this is a :class:`pandas.DataFrame` with the
    same index, otherwise a :class:`pandas.Series` is returned.
'''

    return _quadratic(stats_A, stats_B, cov_AB, data_len, -1)

def product(stats_A, stats_B, cov_AB, data_len):
    '''Calculate the mean and standard error of :math:`f = A \\times B`.

Parameters
----------
stats_A : :class:`pandas.Series` or :class:`pandas.DataFrame`
    Statistics (containing at least the 'mean' and 'standard error' fields) for
    variable :math:`A`.  The rows contain different values of these statistics
    (e.g. from a reblocking analysis) if :class:`pandas.DataFrame` are passed.
stats_B : :class:`pandas.Series` or :class:`pandas.DataFrame`
    Similarly for variable :math:`B`.
cov_AB : float or :class:`pandas.Series`
    Covariance between variables `A and B.  If ``stats_A`` and ``stats_B`` are
    :class:`pandas.DataFrame`, then this must be a :class:`pandas.Series`, with
    the same index as ``stats_A`` and ``stats_B``.
data_len : int or :class:`pandas.Series`
    Number of data points ('observations') used to obtain the statistics given
    in ``stats_A`` and ``stats_B``.  If ``stats_A`` and ``stats_B`` are
    :class:`pandas.DataFrame`, then this must be a :class:`pandas.Series`, with
    the same index as ``stats_A`` and ``stats_B``.

Returns
-------
product_stats : :class:`pandas.Series` or :class:`pandas.DataFrame`
    Mean and standard error (and, if possible/relevant, optimal reblock
    iteration) for :math:`f = A \\times B`.  If ``stats_A``, ``stats_B`` are
    :class:`pandas.DataFrame`, this is a :class:`pandas.DataFrame` with the
    same index, otherwise a :class:`pandas.Series` is returned.
'''

    return _quadratic(stats_A, stats_B, cov_AB, data_len, 1)

def _quadratic(stats_A, stats_B, cov_AB, data_len, sign):
    '''Calculate the mean and standard error of :math:`f = g(A,B)`.

Parameters
----------
stats_A : :class:`pandas.Series` or :class:`pandas.DataFrame`
    Statistics (containing at least the 'mean' and 'standard error' fields) for
    variable :math:`A`.  The rows contain different values of these statistics
    (e.g. from a reblocking analysis) if :class:`pandas.DataFrame` are passed.
stats_B : :class:`pandas.Series` or :class:`pandas.DataFrame`
    Similarly for variable :math:`B`.
cov_AB : float or :class:`pandas.Series`
    Covariance between variables `A and B.  If ``stats_A`` and ``stats_B`` are
    :class:`pandas.DataFrame`, then this must be a :class:`pandas.Series`, with
    the same index as ``stats_A`` and ``stats_B``.
data_len : int or :class:`pandas.Series`
    Number of data points ('observations') used to obtain the statistics given
    in ``stats_A`` and ``stats_B``.  If ``stats_A`` and ``stats_B`` are
    :class:`pandas.DataFrame`, then this must be a :class:`pandas.Series`, with
    the same index as ``stats_A`` and ``stats_B``.
sign : int
    :math:`g(A,B) = A*B` for sign=1 and :math:`f = A/B` for sign=-1.

Returns
-------
func_data : :class:`pandas.Series` or :class:`pandas.DataFrame`
    Mean and standard error (and, if possible/relevant, optimal reblock
    iteration) for :math:`f = g(A,B)`.  If ``stats_A``, ``stats_B`` are
    :class:`pandas.DataFrame`, this is a :class:`pandas.DataFrame` with the
    same index, otherwise a :class:`pandas.Series` is returned.

Raises
------
ValueError
    ``sign`` is not +1 nor -1.
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
        # Were given a single data point rather than a set.
        stats = pd.Series(stats_dict)

    return stats

def pretty_fmt_err(val, err):
    '''Pretty formatting of a value and associated error.

Parameters
----------

val : number
    a (noisy) value.
err: number
    error associated with the value.

Returns
-------
val_str : str
    Value to the number of significant digits known, with the error in the
    last digit in brackets.

Examples
--------
>>> pretty_fmt_err(1.2345, 0.01)
'1.23(1)'
>>> pretty_fmt_err(12331, 40)
'12330(40)'

Notes
-----
Rounding is handled with Python's `round` function, which handles rounding numbers at the midpoint in a range (eg 5 if round to the nearest 10) in a slightly odd way.  As we're normally dealing with noisy data and rounding to remove more than just one significant figure, this is unlikely to impact us.
'''
    # Sig figs of accuracy is:
    sig_figs = lambda x: -int(numpy.floor(numpy.log10(x)))

    # Round value and error to accuracy known:
    sig_err = round(err, sig_figs(err))
    # Use the rounded error in case the error was 'rounded up' (e.g. 9.8->10) to
    # get the number of sif figs the value is actually known to.
    sig_err_figs = sig_figs(sig_err)
    sig_val = round(val, sig_err_figs)

    # Format.
    if sig_err > 1:
        sig_val = int(sig_val)
        sig_err = int(round(sig_err))
        fmt = '%i(%i)'
    else:
        sig_err = int(round(sig_err*10**sig_err_figs))
        fmt = '%%.%if(%%i)' % (sig_err_figs)

    return (fmt % (sig_val, sig_err))

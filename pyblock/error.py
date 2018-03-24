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
    '''Calculate the mean and standard error of :math:`f(A,B) = A/B`.

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
stats : :class:`pandas.Series` or :class:`pandas.DataFrame`
    Mean and standard error (and, if possible/relevant, optimal reblock
    iteration) for :math:`f(A,B)`.  If ``stats_A``, ``stats_B`` are
    :class:`pandas.DataFrame`, this is a :class:`pandas.DataFrame` with the
    same index, otherwise a :class:`pandas.Series` is returned.
'''
    (m_A, m_B) = (stats_A['mean'], stats_B['mean'])
    mean = m_A / m_B
    (se_A, se_B) = (stats_A['standard error'], stats_B['standard error'])
    std_err = abs(mean*numpy.sqrt(
                (se_A/m_A)**2 + (se_B/m_B)**2 - 2*cov_AB/(data_len*m_A*m_B)
              ))
    return _stats_summary(mean, std_err, stats_A, stats_B)

def product(stats_A, stats_B, cov_AB, data_len):
    '''Calculate the mean and standard error of :math:`f(A,B) = A \\times B`.

Parameters
----------
See :func:`ratio`.

Returns
-------
See :func:`ratio`.

'''
    (m_A, m_B) = (stats_A['mean'], stats_B['mean'])
    mean = m_A * m_B
    (se_A, se_B) = (stats_A['standard error'], stats_B['standard error'])
    std_err = abs(mean*numpy.sqrt(
                (se_A/m_A)**2 + (se_B/m_B)**2 + 2*cov_AB/(data_len*m_A*m_B)
              ))
    return _stats_summary(mean, std_err, stats_A, stats_B)

def subtraction(stats_A, stats_B, cov_AB, data_len):
    '''Calculate the mean and standard error of :math:`f(A,B) = A - B`.

Parameters
----------
See :func:`ratio`.

Returns
-------
See :func:`ratio`.
'''
    mean = stats_A['mean'] - stats_B['mean']
    se = 'standard error'
    std_err = stats_A[se]**2 + stats_B[se]**2 - 2*cov_AB/data_len
    std_err = abs(numpy.sqrt(std_err))
    return _stats_summary(mean, std_err, stats_A, stats_B)

def addition(stats_A, stats_B, cov_AB, data_len):
    '''Calculate the mean and standard error of :math:`f(A,B) = A \\plus B`.

Parameters
----------
See :func:`ratio`.

Returns
-------
See :func:`ratio`.
'''
    mean = stats_A['mean'] + stats_B['mean']
    se = 'standard error'
    std_err = stats_A[se]**2 + stats_B[se]**2 + 2*cov_AB/data_len
    std_err = abs(numpy.sqrt(std_err))
    return _stats_summary(mean, std_err, stats_A, stats_B)

def _stats_summary(mean, std_err, stats_A, stats_B):
    '''Summarise error propogation.

Parameters
----------
mean : float or :class:`pandas.Series`
    mean of of :math:`f(A,B)`.
std_err: float or :class:`pandas.Series`
    standard error of :math:`f(A,B)`.
stats_A, stats_B:
    see `func:`ratio`.

Returns
-------

See :func:`ratio`.
'''
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
    if numpy.isinf(val):
        return '%s' % (val,)
    if numpy.isinf(err):
        return '%s(%s)' % (val, err)
    if abs(err) < numpy.finfo(float(err)).eps:
        return '%s(0)' % (val,)

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

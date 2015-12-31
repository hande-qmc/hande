'''Pandas-based wrapper around :mod:`pyblock.blocking`.'''

# copyright: (c) 2014 James Spencer
# license: modified BSD license; see LICENSE for further details.

import numpy
import pandas as pd
import pyblock.blocking

def reblock(data, axis=0, weights=None):
    '''Blocking analysis of correlated data.

Parameters
----------
data : :class:`pandas.Series` or :class:`pandas.DataFrame`
    Data to be blocked.  See ``axis`` for order.
axis : int
    If non-zero, variables in data are in rows with the columns
    corresponding to the observation values.  Blocking is then performed along
    the rows.  Otherwise each column is a variable, the observations are in the
    columns and blocking is performed down the columns.  Only used if data is
    a :class:`pandas.DataFrame`.
weights : :class:`pandas.Series` or :class:`pandas.DataFrame`
    A 1D weighting of the data to be reblocked. For multidimensional data an
    identical weighting is applied to the data for each variable.

Returns
-------
data_len : :class:`pandas.Series`
    Number of data points used in each reblocking iteration.  Note some
    reblocking iterations discard a data point if there were an odd number of
    data points in the previous iteration.
block_info : :class:`pandas.DataFrame`
    Mean, standard error and estimated standard error for each variable at each
    reblock step.
covariance : :class:`pandas.DataFrame`
    Covariance matrix at each reblock step.

See also
--------
:func:`pyblock.blocking.reblock`:
    numpy-based implementation; see for documentation and notes on the
    reblocking procedure.  :func:`pyblock.pd_utils.reblock` is a simple wrapper
    around this.
'''

    try:
        columns = [data.name]
        if data.name is None:
            columns = ['data']
        axis = 0
    except AttributeError:
        # Have DataFrame rather than Series.
        if axis:
            columns = data.index.values
        else:
            columns = data.columns.values

    if weights is not None:
        if isinstance(weights, pd.DataFrame):
            if numpy.min(weights.shape) > 1:
                raise RuntimeError("cannot handle multidimensional weights")
            weights = numpy.array(weights.unstack())
        else:
            weights = weights.values

    block_stats = pyblock.blocking.reblock(data.values,
                                           rowvar=axis,
                                           weights=weights)
    data_size = data.shape[axis]
    optimal_blocks = pyblock.blocking.find_optimal_block(data_size, block_stats)

    # Now nicely package it up into a dict of pandas/built-in objects.

    iblock = []
    data_len = []
    block_info = []
    covariance = []
    keys = ['mean', 'standard error', 'standard error error', 'optimal block']
    multi_keys = [(col,k) for col in columns for k in keys]
    multi_keys = pd.MultiIndex.from_tuples(multi_keys)
    null = numpy.zeros(len(columns))
    for stat in block_stats:
        # Contents of stat:
        #     (iblock, data_len, mean, covariance, standard err,
        #      esimate of error in standard error)
        iblock.append(stat.block)
        data_len.append(stat.ndata)

        pd_stat = [stat.mean, stat.std_err, stat.std_err_err, null]
        pd_stat = numpy.array(pd_stat).T.flatten()
        block_info.append(pd.Series(pd_stat, index=multi_keys))

        # Covariance is a 2D matrix (in general) so can't put it into
        # a DataFrame with everything else, so put it in its own.
        cov = numpy.array(stat.cov, ndmin=2)
        covariance.append(pd.DataFrame(cov, index=columns, columns=columns))

    data_len = pd.Series(data_len, index=iblock, name='data length')
    data_len.index.name = 'reblock'

    block_info = pd.concat(block_info, axis=1, keys=iblock).transpose()
    block_info.index.name = 'reblock'
    loc = block_info.columns.get_level_values(1) == 'optimal block'
    block_info.loc[:,loc] = ''

    covariance = pd.concat(covariance, keys=iblock)
    covariance.index.names = ['reblock', '']

    for (ivar, optimal) in enumerate(optimal_blocks):
        if optimal >= 0:
            block_info.loc[optimal,(columns[ivar], 'optimal block')] = '<---    '

    return (data_len, block_info, covariance)

def optimal_block(block_sub_info):
    '''Get the optimal block value from the reblocking data.

Parameters
----------
block_sub_info: :class:`pandas.DataFrame` or :class:`pandas.Series`
    Reblocking data (i.e. the first item of the tuple returned by ``reblock``),
    or a subset thereof containing the statistics columns for one or more data
    items.

Returns
-------
index : int
    Reblocking index corresponding to the reblocking iteration at which serial
    correlation has been removed (as estimated by the procedure in
    ``pyblock.blocking.find_optimal_block``).  If multiple data sets are passed
    in block_sub_info, this is the maximum index out of all data sets.  Set to
    inf if an optimal block is not found for a data set.

Raises
------
ValueError
    block_sub_info contains no Series or column in DataFrame named 'optimal
    block'.
'''

    # Handle the following cases:
    # * Series with optimal block in it.
    # * block_sub_info DataFrame for one variable (no hierarchical column names)
    # * block_sub_info DataFrame for multiple variables (hierarchical column names)
    # (each set of columns for one variable in block_sub_info contains the mean,
    # standard error and estimated error in the standard error for that
    # variable).
    try:
        if 'optimal block' in block_sub_info.name:
            iterator = [('optimal block', block_sub_info)]
        else:
            raise ValueError('No optimal block data')
    except AttributeError:
        # Have DataFrame.
        # 'optimal block' is in the innermost level.
        level = block_sub_info.columns.nlevels - 1
        opt_cols = [col == 'optimal block'
                     for col in block_sub_info.columns.get_level_values(level)]
        if not any(opt_cols):
            raise ValueError('No optimal block data')
        iterator = block_sub_info.loc[:,opt_cols].iteritems()

    opt = -1
    for (name, col) in iterator:
        col_opt = col[col != ''].index
        if len(col_opt) == 0:
            opt = float('inf')
        elif len(col_opt) == 1:
            opt = max(col_opt[0], opt)
        else:
            raise ValueError('Multiple entries listed as optimal.')

    return opt

def reblock_summary(block_sub_info):
    '''Get the data corresponding to the optimal block from the reblocking data.

Parameters
----------
block_sub_info : :class:`pandas.DataFrame` or :class:`pandas.Series`
    Reblocking data (i.e. the first item of the tuple returned by ``reblock``),
    or a subset thereof containing the statistics columns for one or more data
    items.

Returns
-------
summary : :class:`pandas.DataFrame`
    Mean, standard error and estimate of the error in the standard error
    corresponding to the optimal block size in the reblocking data (or largest
    optimal size if multiple data sets are given.  The index is labelled with
    the data name, if known.  An empty DataFrame is returned if no optimal block
    size was found.
'''
    opt = optimal_block(block_sub_info)
    if opt < float('inf'):
        summary = block_sub_info.loc[opt]
        # Convert to DataFrame, with statistics in columns.
        if summary.index.nlevels == 1:
            # Sadly don't know the data name; leave to user.
            summary = pd.DataFrame(summary).T
        else:
            # Have hierarchical index; can pivot into a DataFrame.
            # Each row will be labelled by the data name.
            summary = summary.unstack()
        summary.drop('optimal block', axis=1, inplace=True)
    else:
        summary = pd.DataFrame()
    return summary

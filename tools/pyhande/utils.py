'''Utility procedures'''

import numpy as np

def groupby_beta_loops(data):
    '''Group a HANDE DMQMC data table by beta loop.

Parameters
----------
data : :class:`pandas.DataFrame`
    DMQMC data table (e.g. obtained by :func:`pyhande.extract.extract_data`.

Returns
-------
grouped : :class:`pandas.DataFrameGroupBy`
    GroupBy object with data table grouped by beta loop.
'''

    # Exploit the fact that (except for possibly the last beta loop due to wall
    # time) each beta loop contains the same set of iterations.
    indx = np.arange(len(data)) // len(data['iterations']).unique()
    return data.groupby(indx)

def groupby_iterations(data):
    '''Group a HANDE QMC data table by blocks of iterations.

Parameters
----------
data : :class:`pandas.DataFrame`
    QMC data table (e.g. obtained by :func:`pyhande.extract.extract_data`.

Returns
-------
grouped : :class:`pandas.DataFrameGroupBy`
    GroupBy object with data table grouped into blocks within which the
    iteration count increases monotonically.
'''
    indx = np.zeros(len(data))
    prev_iteration = -1
    curr_indx = 0
    for i in range(len(data)):
        if data['iterations'].iloc[i] < prev_iteration:
            # new block of iterations
            curr_indx += 1
        indx[i] = curr_indx
        prev_iteration = data['iterations'].iloc[i]

    return data.groupby(indx)

#!/usr/bin/env python
'''Analysis of data from state histogram collection'''

import warnings
import numpy as np
import pandas as pd

def read_state_histogram_file(histogram_file, trunc_ex1, trunc_ex2):
    ''' Read in a state histogram file and return the data in
        a pandas DataFrame, as well as the maximum excitation levels.

Parameters
----------
histogram_file : string
    A file which contains state histogram data..
trunc_ex1: integer
    A value to truncate the first excitation level data read in to less than or
    equal to this value. If None is provided no truncation is performed.
trunc_ex2: integer
    A value to truncate the second excitation level data read in to less than or
    equal to this value. If None is provided no truncation is performed.

Returns
-------
histogram : :class:`pandas.DataFrame`
    A pandas data frame of the histogram data from the histogram_file.
mex1 : integer
    The maximum excitation level of the first excitation level.
mex2 : integer
    The maximum excitation level of the second excitation level.
'''
    histogram, allowed_keys, mex1, mex2 = {}, [], 0, 0

    with open(histogram_file, 'r') as of:
        for ln, l in enumerate(of):
            if ln == 0:
                histogram['bin_edges'] = []
                for column in l.split('Ex.Lvl')[1:]:
                    ex1, ex2 = [int(ex) for ex in column.split()]
                    if trunc_ex1 is not None and ex1 > trunc_ex1: continue
                    if trunc_ex2 is not None and ex2 > trunc_ex2: continue
                    if ex1 > mex1: mex1 = ex1
                    if ex2 > mex2: mex2 = ex2
                    histogram[f'Ex.Lvl {ex1} {ex2}'] = []
                    allowed_keys.append(f'Ex.Lvl {ex1} {ex2}')
            else:
                ld = np.array(l.split()).astype(float)
                for column, ndets in zip(list(histogram.keys()), ld):
                    if column not in allowed_keys: continue
                    histogram[column].append(ndets)

    return pd.DataFrame(histogram), mex1, mex2

def average_histograms(state_histograms):
    ''' Generate the averaged profiles of each state.

Parameters
----------
state_histogram : :class:`pandas.DataFrame`
    A pandas data frame with all the histograms concatenated together

Returns
-------
average_state_histogram : :class:`pandas.DataFrame`
    The averaged profiles of each state
'''
    # TODO - WZV
    # Need to determine if covariance is a factor or not based on the RNG seeds.
    # For now we assume each histogram comes from a seperate trajectory
    warnings.warn('\n Assuming the State Histograms are uncorrelated \n'+\
                    ' and so there is no covariance. Be warned!')

    temp_average = state_histograms.groupby('bin_edges').mean()
    temp_std_err = state_histograms.groupby('bin_edges').sem()

    std_err_cols = {col : col + ' error' for col in temp_std_err.columns}
    temp_std_err.rename(columns = std_err_cols, inplace=True)

    average_state_histogram = pd.concat([temp_average, temp_std_err])
    average_state_histogram.reset_index(inplace=True)

    return average_state_histogram

def analyse_state_histograms(state_histogram_outputs, trunc_ex1, trunc_ex2):
    ''' Perform analysis on the state histograms from a DMQMC or FCIQMC
        calculation.

Parameters
----------
state_histogram_outputs : list
    Contains all the state histogram output files which are to be analysed.
trunc_ex1: integer
    A value to truncate the first excitation level data read in to less than or
    equal to this value. If None is provided no truncation is performed.
trunc_ex2: integer
    A value to truncate the second excitation level data read in to less than or
    equal to this value. If None is provided no truncation is performed.

Returns
-------
average_state_histogram : :class:`pandas.DataFrame`
    The averaged profile of all the provided state histograms.
'''
    state_histograms, ex1, ex2 = [], [], []

    for histogram_file in state_histogram_outputs:
        histogram, mex1, mex2 = read_state_histogram_file(histogram_file,
                                                            trunc_ex1,
                                                            trunc_ex2)
        state_histograms.append(histogram)
        ex1.append(mex1)
        ex2.append(mex2)

    # Do a simple sanity check on the data we are processing.
    if any(ex != ex1[0] for ex in ex1) or any(ex != ex2[0] for ex in ex2):
        raise ValueError(' Excitation levels do not match for the \n'+\
                         ' provided histogram files!')

    state_histograms = pd.concat(state_histograms)
    average_state_histogram = average_histograms(state_histograms)

    return average_state_histogram

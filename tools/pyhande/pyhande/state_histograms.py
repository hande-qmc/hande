#!/usr/bin/env python
'''Analysis of data from state histogram collection'''

import numpy as np
import pandas as pd
from collections import OrderedDict


def read_state_histogram_file(histogram_file):
    ''' Read in a state histogram file and return the data in
    a pandas DataFrame, as well as the maximum excitation levels.

    Parameters
    ----------
    histogram_file : string
        A file which contains state histogram data.

    Returns
    -------
    histogram : :class:`pandas.DataFrame`
        A pandas data frame of the histogram data from the histogram_file.
    mex1 : integer
        The maximum excitation level of the first excitation index.
    mex2 : integer
        The maximum excitation level of the second excitation index.
    '''

    histogram, mex1, mex2 = {}, 0, 0

    with open(histogram_file, 'r') as of:
        for ln, l in enumerate(of):
            if ln == 0:
                histogram['bin_edges'] = []
                for column in l.split('Ex.Lvl')[1:]:
                    ex1, ex2 = [int(ex) for ex in column.split()]
                    histogram[f'Ex.Lvl {ex1} {ex2}'] = []
                    if ex1 > mex1:
                        mex1 = ex1
                    if ex2 > mex2:
                        mex2 = ex2
            else:
                ld = np.array(l.split()).astype(float)
                for column, ndets in zip(list(histogram.keys()), ld):
                    histogram[column].append(ndets)

    return pd.DataFrame(histogram), mex1, mex2


def collect_state_histogram_data(state_histogram_outputs):
    ''' Collect all the state histograms data and store in an ordered
    dictionary based on the report (time/temperature) the data was
    generated at. Additionally perform some basic sanity checks on the data.

    Parameters
    ----------
    state_histogram_outputs : list of strings
        The state histogram output files which are to be analysed.

    Returns
    -------
    grouped : tuple of OrderedDicts
        Stores OrderdDicts for the files, data set data-frames,
        and maximum excitation indices indexed by the report the data
        was generated at.

    Raises
    ------
    RuntimeError
        If the maximum excitations for the excitation indices 1/2
        are not consistent across data sets or the data set shapes are
        not consistent across all data sets.
    '''

    # First group the files by the report they are generated from.
    # A file name should follow the format:
    #   EXLEVELPOPS-RNGyyy-IREPORTzzz
    # Where yyy is an rng seed and zzz is a report cycle.
    unique_reports = []
    for hist_file in state_histogram_outputs:
        ireport = hist_file.split('-')[-1].replace('IREPORT', '')
        if ireport not in unique_reports:
            unique_reports.append(ireport)

    unique_reports = np.array(unique_reports)
    unique_reports = unique_reports[np.argsort(unique_reports.astype(int))]

    grouped_files = OrderedDict()
    grouped_data = OrderedDict()
    grouped_mex1 = OrderedDict()
    grouped_mex2 = OrderedDict()
    for ireport in unique_reports:
        grouped_files[ireport] = []
        grouped_data[ireport] = []
        grouped_mex1[ireport] = []
        grouped_mex2[ireport] = []

        updated_state_histogram_outputs = state_histogram_outputs.copy()
        for hist_file in state_histogram_outputs:
            if 'IREPORT' + ireport in hist_file:
                grouped_files[ireport].append(hist_file)

                df, mex1, mex2 = read_state_histogram_file(hist_file)
                grouped_data[ireport].append(df)
                grouped_mex1[ireport].append(mex1)
                grouped_mex2[ireport].append(mex2)

                updated_state_histogram_outputs.remove(hist_file)

        # A minor performance gain: We don't need to loop through
        # files we have already stored.
        state_histogram_outputs = updated_state_histogram_outputs

    # Do a couple simple sanity checks on our data
    report_shapes, report_mex1s, report_mex2s = [], [], []
    for ireport in unique_reports:
        dfs = grouped_data[ireport]
        mex1s = grouped_mex1[ireport]
        mex2s = grouped_mex2[ireport]

        if any(k.shape != dfs[0].shape for k in dfs):
            raise RuntimeError('The shapes of the data sets for files in '
                               f'IREPORT{ireport} do not agree!')
        elif any(k != mex1s[0] for k in mex1s):
            raise RuntimeError('The maximum excitation 1 for files in '
                               f'IREPORT{ireport} do not agree!')
        elif any(k != mex2s[0] for k in mex2s):
            raise RuntimeError('The maximum excitation 2 for files in '
                               f'IREPORT{ireport} do not agree!')

        report_shapes.append(dfs[0].shape)
        report_mex1s.append(mex1s[0])
        report_mex2s.append(mex2s[0])

    # Now do the sanity check across ireports
    if any(k != report_shapes[0] for k in report_shapes):
        raise RuntimeError('The shapes of the data sets for files across '
                           'all IREPORT do not agree!')
    elif any(k != report_mex1s[0] for k in report_mex1s):
        raise RuntimeError('The maximum excitation 1 for files across '
                           'all IREPORT do not agree!')
    elif any(k != report_mex2s[0] for k in report_mex2s):
        raise RuntimeError('The maximum excitation 2 for files across '
                           'all IREPORT do not agree!')

    grouped = (grouped_files, grouped_data, grouped_mex1, grouped_mex2)
    return grouped


def average_histograms(grouped):
    ''' Collapse the state histograms down into an averaged profile
    for each bin of walker populations.

    Parameters
    ----------
    grouped : tuple of OrderedDicts
        Stores OrderdDicts for the files, data set data-frames,
        and maximum excitation indices indexed by the report the data
        was generated at.

    Returns
    -------
    averaged : :class:`pandas.DataFrame`
        Contains the average number of determinants in bin each
        bin for a given temperature/imaginary time.

    Raises
    ------
    RuntimeError
        If the bin edges are not consistent across all data sets
    '''

    grouped_files, grouped_data, grouped_mex1, grouped_mex2 = grouped

    kbe = 'bin_edges'
    averaged = {}
    for ireport, dfs in grouped_data.items():
        if kbe not in averaged.keys():
            averaged[kbe] = dfs[0][kbe]

        # Sanity check our data frames bins agree with one another
        if any(np.count_nonzero(k[kbe] != dfs[0][kbe]) > 0 for k in dfs):
            raise RuntimeError('The bin edges for files in '
                               f'IREPORT{ireport} do not agree!')
        elif any(np.count_nonzero(k[kbe] != averaged[kbe]) > 0 for k in dfs):
            raise RuntimeError('The bin edges for files across '
                               'all IREPORT do not agree!')

        mex1 = grouped_mex1[ireport][0]
        mex2 = grouped_mex2[ireport][0]

        # Averaged data for a given bin, i.e. for: Ex.Lvl x y
        # calculate the average across simulations
        igroup = pd.concat(dfs).groupby('bin_edges', as_index=False)
        iave = igroup.mean()
        isem = igroup.sem()

        # Accumulate the total average across x, i.e. for: Ex.Lvl x y
        # collapse the data set to: Sum Col y = \sum_x Ex.Lvl x y
        kx = [f'Ex.Lvl {x}' for x in range(mex1+1)]
        ky = [[k+f' {y}' for k in kx] for y in range(mex2+1)]
        ex2_sum = [np.sum(iave[k], axis=1) for k in ky]
        ex2_sem = [np.sum(isem[k]**2.0, axis=1)**0.5 for k in ky]

        # Now accumulate data across y, i.e.
        # bin = \sum_y Sum Col y
        bin_sum = np.sum(ex2_sum, axis=0)
        bin_sem = np.sum(np.power(ex2_sem, 2.0), axis=0)**0.5
        averaged[f'{ireport} \\sum'] = bin_sum
        averaged[f'{ireport} sem'] = bin_sem

    averaged = pd.DataFrame(averaged)
    return averaged


def analyse_state_histograms(state_histogram_outputs):
    ''' Perform analysis on the state histograms from a DMQMC or FCIQMC
    calculation.

    Parameters
    ----------
    state_histogram_outputs : list
        The state histogram output files which are to be analysed.

    Returns
    -------
    averaged : :class:`pandas.DataFrame`
        Contains the average number of determinants in bin each
        bin for a given temperature/imaginary time.
    '''

    grouped = collect_state_histogram_data(state_histogram_outputs)

    averaged = average_histograms(grouped)

    return averaged

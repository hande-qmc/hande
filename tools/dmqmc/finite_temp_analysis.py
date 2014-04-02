#!/usr/bin/python
'''finite_temp_analysis.py [options] file_1 file_2 ... file_N

Analyse the output of a HANDE DMQMC calculation by averaging
temperature-dependent data across beta loops.'''

import pandas as pd
import sys
import pyhande
import pyblock
import optparse
import numpy as np
import scipy.interpolate

def perform_averaging(beta_values, estimates, nsamples, options):
    '''Perform analysis of DMQMC data from a HANDE calculation.

Parameters
----------
beta_values : :class:`pandas.Series`
    The beta values at which data was accumulated (for at least one beta loop, but
    not necessarily all loops).
estimates : :class:`pandas.DataFrame`
    All of the estimates from all beta loops, which are to be combined and analysed.
nsamples : :class:`pandas.Series`
    The number of samples at the various values in beta_values.
options : :class:`OptionParser`
    Options read in from command line.

Returns
-------
final_estimates : :class:`pandas.Series`
    All of the final mean and error estimates as a function of beta, ready for output.
None.
'''

    columns = list(estimates.columns.values)

    #DataFrame for all of the final estimates to be printed.
    final_estimates = pd.DataFrame(index=beta_values)
    # DataFrame for the numerator.
    num = pd.DataFrame(columns=['mean','standard error'], index=beta_values)
    # DataFrame for the trace from the first replica.
    tr1 = pd.DataFrame(columns=['mean','standard error'], index=beta_values)
    # DataFrame for the trace from the second replica.
    tr2 = pd.DataFrame(columns=['mean','standard error'], index=beta_values)

    means = estimates.groupby(level=2).mean()
    covariances = estimates.groupby(level=2).cov()

    # The keys hold the names to be output and the values hold the possible
    # columns names in the estimates DataFrame.
    observables = dict([
        ('Tr[Hp]/Tr[p]','\\sum\\rho_{ij}H_{ji}'), 
        ('Tr[Sp]/Tr[p]','\\sum\\rho_{ij}S_{ji}'),
        ('Tr[Cp]/Tr[p]','\\sum\\rho_{ij}C_{ji}'),
        ('Tr[Mp]/Tr[p]','\\sum\\rho_{ij}M2{ji}'),
        ('Tr[H2p]/Tr[p]','\\sum\\rho_{ij}H2{ji}')
    ])

    tr1['mean'] = means['Trace']
    tr1['standard error'] = np.sqrt(covariances.xs('Trace',level=1)['Trace']/nsamples)

    # If requested, add the shift and traces to the final_estimates DataFrame, to be output later.
    if options.with_shift:
        final_estimates['Shift'] = means['Shift']
        final_estimates['Shift s.d.'] = np.sqrt(covariances.xs('Shift',level=1)['Shift'])
    if options.with_trace:
        final_estimates['Trace'] = tr1['mean']
        final_estimates['Trace s.d.'] = np.sqrt(covariances.xs('Trace',level=1)['Trace'])
        if 'Trace_2' in columns:
            final_estimates['Trace 2'] = means['Trace_2']
            final_estimates['Trace 2 s.d.'] = np.sqrt(covariances.xs('Trace_2',level=1)['Trace_2'])

    # Compute the mean and standard error estimates for all quantities of the form Tr(\rho O)/Tr(\rho).
    for (k,v) in observables.items():
        if v in columns:
            num['mean'] = means[v]
            num['standard error'] = np.sqrt(covariances.xs(v,level=1)[v]/nsamples)
            cov_AB = covariances.xs('Trace',level=1)[v]
            stats = pyblock.error.ratio(num, tr1, cov_AB, nsamples)

            final_estimates[k] = stats['mean']
            final_estimates[k+' error'] = stats['standard error']

    # Compute the mean and standard error estimates for the Renyi entropy (S2), for all subsystems.
    nrdms = 0
    for column in columns:
        have_s2 = False
        if 'Renyi_2_numerator' in column: # RDM S2
            nrdms += 1
            have_s2 = True
            num_col = column
            tr1_col = 'RDM'+str(nrdms)+'_trace_1'
            tr2_col = 'RDM'+str(nrdms)+'_trace_2'
            out_str = 'RDM'+str(nrdms)+' S2'
        elif 'Full_R2_numerator' in column: # Full S2
            have_s2 = True
            num_col = column
            tr1_col = 'Trace'
            tr2_col = 'Trace_2'
            out_str = 'Full S2'

        if have_s2:
            num['mean'] = means[num_col]
            tr1['mean'] = means[tr1_col]
            tr2['mean'] = means[tr2_col]
            num['standard error'] = np.sqrt(covariances.xs(num_col,level=1)[num_col]/nsamples)
            tr1['standard error'] = np.sqrt(covariances.xs(tr1_col,level=1)[tr1_col]/nsamples)
            tr2['standard error'] = np.sqrt(covariances.xs(tr2_col,level=1)[tr2_col]/nsamples)

            # A denotes the numerator, B denotes the first trace and C denotes the second trace.
            cov_AB = covariances.xs(num_col,level=1)[tr1_col]
            cov_AC = covariances.xs(num_col,level=1)[tr2_col]
            cov_BC = covariances.xs(tr1_col,level=1)[tr2_col]

            final_estimates[out_str], final_estimates[out_str+' error'] = \
                calc_S2(num, tr1, tr2, cov_AB, cov_AC, cov_BC, nsamples)

    return final_estimates


def calc_S2(stats_A, stats_B, stats_C, cov_AB, cov_AC, cov_BC, data_len):
    '''Calculate the mean and standard error of :math:`f = -log(A/BC)`.

Parameters
----------
stats_A : :class:`pandas.DataFrame`
    Statistics (containing at least the 'mean' and 'standard error' fields) for
    variable :math:`A`.  The rows contain different values of these statistics
    for different Beta values.
stats_B : :class:`pandas.DataFrame`
    Similarly for variable :math:`B`.
stats_C : :class:`pandas.DataFrame`
    Similarly for variable :math:`C`.
cov_AB : :class:`pandas.Series`
    Covariance between variables A and B.
cov_AC : :class:`pandas.Series`
    Covariance between variables A and C.
cov_BC : :class:`pandas.Series`
    Covariance between variables B and C.
data_len : :class:`pandas.Series`
    Number of Beta loops used to obtain the statistics given in ``stats_A``, ``stats_B``
    and ``stats_C``.

Returns
-------
mean : :class:`pandas.Series`
    Mean for :math:`f = -log(A/BC)`.
std_err: :class:`pandas.Series`
    Standard error for :math:`f = -log(A/BC)
'''

    mean = -np.log(stats_A['mean']/(stats_B['mean']*stats_C['mean']))/np.log(2)
    std_err = 1/np.log(2) * np.sqrt(
        (stats_A['standard error']/stats_A['mean'])**2 +
        (stats_B['standard error']/stats_B['mean'])**2 +
        (stats_C['standard error']/stats_C['mean'])**2 -
        2*cov_AB/(data_len*stats_A['mean']*stats_B['mean']) -
        2*cov_AC/(data_len*stats_A['mean']*stats_C['mean']) +
        2*cov_BC/(data_len*stats_B['mean']*stats_C['mean']) )

    return mean, std_err

def calc_spline_fit(column, estimates):
    '''Find a cubic B-spline fit for the requested column in estimates.

Parameters
----------
column: string
    Column name of estimate for which to find a spline fit.
estimates : :class:`pandas.DataFrame`
    Must contain a column with the above name, and another with the
    name column+' error'. The indices are used for the beta values.

Returns
-------
spline_fit : :class:`pandas.Series`
    Cubic B-spline fit.
'''
    
    beta_values = list(estimates.index.values)
    values = list(estimates[column].values)
    weights = list(1/estimates[column+' error'].values)
    
    tck = scipy.interpolate.splrep(beta_values, values, weights, k=3)
    spline_fit = scipy.interpolate.splev(beta_values, tck)
    
    return pd.Series(spline_fit)

def parse_args(args):
    '''Parse command-line arguments.

Parameters
----------
args : list of strings
    command-line arguments.

Returns
-------
filenames : list of strings
    list of QMC output files
options : :class:`OptionParser`
    Options read in from command line.
'''

    parser = optparse.OptionParser(usage = __doc__)
    parser.add_option('-s', '--with-shift', action='store_true', dest='with_shift', 
                      default=False, help='Output the averaged shift profile and '
                      'the standard deviation of these profiles across beta loops.')
    parser.add_option('-t', '--with-trace', action='store_true', dest='with_trace',
                      default=False, help='Output the averaged traces and the '
                      'standard deviation of these traces, for all replicas present.')
    parser.add_option('-b', '--with-spline', action='store_true', dest='with_spline',
                      default=False, help='Output a B-spline fit for each of '
                      ' estimates calculated')

    (options, filenames) = parser.parse_args(args)

    if not filenames:
        parser.print_help()
        sys.exit(1)

    return (filenames, options)

def main(args):
    '''Run data analysis on finite-temperature HANDE output.

Parameters
----------
args : list of strings
    command-line arguments.

Returns
-------
None.
'''

    (files, options) = parse_args(args)
    (metadata, data) = pyhande.extract.extract_data_sets(files)

    # Convert the iteration number to the beta value.
    tau = metadata[0]['tau']
    data.rename(columns={'iterations' : 'Beta'}, inplace=True)
    data['Beta'] = data['Beta']*tau

    data.set_index('Beta', inplace=True, append=True)
    # The number of beta loops contributing to each beta value.
    nsamples = data['Trace'].groupby(level=2).count()
    beta_values = pd.Series(nsamples.index)
    estimates = data.loc[:,'Shift':'# H psips']

    final_estimates = perform_averaging(beta_values, estimates, nsamples, options)

    # If requested, calculate a spline fit for all mean estimates.
    if options.with_spline:
        columns = list(final_estimates.columns.values)
        for column in columns:
            if ('Tr[p]' in column or 'S2' in column) and (not 'error' in column):
                final_estimates[column+' spline'] = calc_spline_fit(column, final_estimates)

    print final_estimates.to_string()

if __name__ == '__main__':

    main(sys.argv[1:])


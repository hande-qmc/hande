#!/usr/bin/python
'''finite_temp_analysis.py [options] file_1 file_2 ... file_N

Analyse the output of a HANDE DMQMC calculation by averaging
temperature-dependent data across beta loops.'''

import pandas as pd
import sys
import pyhande
import pyblock
import optparse
import numpy

def run_dmqmc_analysis(beta_values, estimates, nsamples):
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

Returns
-------
None.
'''

    columns = list(estimates.columns.values)
    final_estimates = pd.DataFrame(index=beta_values)
    trace_1 = pd.DataFrame(columns=['mean','standard error'], index=beta_values)
    trace_2 = pd.DataFrame(columns=['mean','standard error'], index=beta_values)
    numerator = pd.DataFrame(columns=['mean','standard error'], index=beta_values)

    # The keys hold the names to be output and the values hold the possible
    # columns names in the estimates DataFrame.
    observables = dict([
        ('Tr[Hp]/Tr[p]','\\sum\\rho_{ij}H_{ji}'), 
        ('Tr[Sp]/Tr[p]','\\sum\\rho_{ij}S_{ji}'),
        ('Tr[Cp]/Tr[p]','\\sum\\rho_{ij}C_{ji}'),
        ('Tr[Mp]/Tr[p]','\\sum\\rho_{ij}M2{ji}'),
        ('Tr[H2p]/Tr[p]','\\sum\\rho_{ij}H2{ji}')
    ])

    means = estimates.groupby(level=2).mean()
    covariances = estimates.groupby(level=2).cov()

    trace_1['mean'] = means['Trace']
    trace_1['standard error'] = numpy.sqrt(covariances.xs('Trace',level=1)['Trace']/nsamples)

    for (k,v) in observables.items():
        if v in columns:
            numerator['mean'] = means[v]
            numerator['standard error'] = numpy.sqrt(covariances.xs(v,level=1)[v]/nsamples)
            cov_AB = covariances.xs('Trace',level=1)[v]

            stats = pyblock.error.ratio(numerator, trace_1, cov_AB, nsamples)
            final_estimates[k] = stats['mean']
            final_estimates[k+' error'] = stats['standard error']

    v = 'Full_R2_numerator'
    k = 'Full S2'
    if v in columns:
        numerator['mean'] = means[v]
        numerator['standard error'] = numpy.sqrt(covariances.xs(v,level=1)[v]/nsamples)
        trace_2['mean'] = means['Trace_2']
        trace_2['standard error'] = numpy.sqrt(covariances.xs('Trace_2',level=1)['Trace_2']/nsamples)

        # A denotes the numerator, B denotes the first trace and C denotes the second trace.
        cov_AB = covariances.xs(v,level=1)['Trace']
        cov_AC = covariances.xs(v,level=1)['Trace_2']
        cov_BC = covariances.xs('Trace',level=1)['Trace_2']

        final_estimates[k], final_estimates[k+' error'] = calc_S2(numerator, trace_1, trace_2, cov_AB, cov_AC, cov_BC, nsamples)

    print final_estimates.to_string()


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

    mean = -numpy.log(stats_A['mean']/(stats_B['mean']*stats_C['mean']))/numpy.log(2)
    std_err = 1/numpy.log(2) * numpy.sqrt(
        (stats_A['standard error']/stats_A['mean'])**2 +
        (stats_B['standard error']/stats_B['mean'])**2 +
        (stats_C['standard error']/stats_C['mean'])**2 -
        2*cov_AB/(data_len*stats_A['mean']*stats_B['mean']) -
        2*cov_AC/(data_len*stats_A['mean']*stats_C['mean']) +
        2*cov_BC/(data_len*stats_B['mean']*stats_C['mean']) )

    return mean, std_err

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
'''

    parser = optparse.OptionParser(usage = __doc__)

    (options, filenames) = parser.parse_args(args)

    if not filenames:
        parser.print_help()
        sys.exit(1)

    return (filenames)

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

    (files) = parse_args(args)
    (metadata, data) = pyhande.extract.extract_data_sets(files)
    tau = metadata[0]['tau']

    # Convert the iteration number to the beta value.
    data.rename(columns={'iterations' : 'Beta'}, inplace=True)
    data['Beta'] = data['Beta']*tau

    data.set_index('Beta', inplace=True, append=True)
    nsamples = data['Trace'].groupby(level=2).count()
    beta_values = pd.Series(nsamples.index)
    estimates = data.loc[:,'Shift':'# H psips']
    run_dmqmc_analysis(beta_values, estimates, nsamples)

if __name__ == '__main__':

    main(sys.argv[1:])


#!/usr/bin/python
'''finite_temp_analysis.py [options] file_1 file_2 ... file_N

Analyse the output of a HANDE DMQMC calculation by averaging
temperature-dependent data across beta loops.'''

import pandas as pd
import sys
import pyhande
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
    All of the estimates from all beta loops, whcich are to be combined and analysed.
nsamples : :class:`pandas.Series`
    The number of samples at the various values in beta_values.

Returns
-------
None.
'''

    columns = list(estimates.columns.values)
    final_estimates = pd.DataFrame({'Beta' : beta_values})
    trace = pd.DataFrame(columns=['mean','standard_error'], index=beta_values)
    numerator = pd.DataFrame(columns=['mean','standard_error'], index=beta_values)

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

    trace['mean'] = means['Trace']
    trace['standard_error'] = numpy.sqrt(covariances.xs('Trace',level=1)['Trace']/nsamples)

    for (k,v) in observables.items():
        if v in columns:
            numerator['mean'] = means[v]
            numerator['standard_error'] = numpy.sqrt(covariances.xs(v,level=1)[v]/nsamples)
            cov = covariances.xs('Trace',level=1)[v]/nsamples

            final_estimates[k] = numerator['mean']/trace['mean']
            final_estimates[k+' error'] = numpy.sqrt(
                (numerator['standard_error']/numerator['mean'])**2 +
                (trace['standard_error']/trace['mean'])**2 -
                2*cov/(numerator['mean']*trace['mean']))*\
                abs(final_estimates[k])

    print final_estimates.to_string(index=False)

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


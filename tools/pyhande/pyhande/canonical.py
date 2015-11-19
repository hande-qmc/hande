'''Analysis of data from  canonical thermodynamic calculations.'''

import pandas as pd
import pyblock
import pyhande.legacy
import numpy as np


def analyse_hf_observables(means, covariances, nsamples):
    ''' Perform Error analysis for Hartree-Fock estimates which
are the ratio of two quantities.

Parameters
----------
means : :class:`pandas.DataFrame`
    Data frame containing means of verious observables.
covariances : :class:`pandas.DataFrame`
    Data frame containing covariances between various observables.
nsamples : int
    Number of samples contributing to estimates and standard errors

Returns
-------
results : :class:`pandas.DataFrame`
    Averaged Hartree-Fock estimates along with error estimates.
'''

    observables = dict([
        ('<T>_HF', r'Tr(T\rho_HF)'),
        ('<V>_HF', r'Tr(V\rho_HF)'),
        ('<H>_HF', r'Tr(H\rho_HF)'),
    ])

    num = pd.DataFrame()
    trace = pd.DataFrame()
    results = pd.DataFrame()
    trace['mean'] = [means[r'Tr(\rho_HF)']]
    trace['standard error'] = (
            [np.sqrt(covariances[r'Tr(\rho_HF)'][r'Tr(\rho_HF)']/nsamples)])

    for (k, v) in observables.items():
        num['mean'] = [means[v]]
        num['standard error'] = [np.sqrt(covariances[v][v]/nsamples)]
        cov_ab = covariances[v][r'Tr(\rho_HF)']

        stats = pyblock.error.ratio(num, trace, cov_ab, nsamples)

        results[k] = stats['mean']
        results[k+'_error'] = stats['standard error']

    return results


def estimates(metadata, data):
    '''Perform error analysis for canonical thermodynamic estimates.

Parameters
----------
metadata : dict
    metadata (i.e. calculation information, parameters and settings) extracted
    from output files.
data : :class:`pandas.DataFrame`
    HANDE QMC data.

Returns
-------
results : :class:`pandas.DataFrame`
    Averaged estimates.
'''

    data['<H>_0'] = data['<T>_0'] + data['<V>_0']
    data[r'Tr(H\rho_HF)'] = data[r'Tr(T\rho_HF)'] + data[r'Tr(V\rho_HF)']

    means = data.mean()
    covariances = data.cov()
    nsamples = len(data['<T>_0'])

    results = pd.DataFrame()
    if 'beta' in metadata:
        # New, richer JSON-based metadata.
        results['Beta'] = [metadata['beta']]
    else:
        # Hope to find it in the input file...
        results['Beta'] = pyhande.legacy.extract_input(metadata, 'beta')
    # Free estimates contain no denominator so the error is
    # just the standard error.
    results['<H>_0'] = [means['<H>_0']]
    results['<H>_0_error'] = [np.sqrt(covariances['<H>_0']['<H>_0']/nsamples)]
    results['<T>_0'] = [means['<T>_0']]
    results['<T>_0_error'] = [np.sqrt(covariances['<T>_0']['<T>_0']/nsamples)]
    results['<V>_0'] = [means['<V>_0']]
    results['<V>_0_error'] = [np.sqrt(covariances['<V>_0']['<V>_0']/nsamples)]
    results['N_acc/N_att'] = [means['N_ACC/N_ATT']]
    results['N_acc/N_att_error'] = (
                  [np.sqrt(covariances['N_ACC/N_ATT']['N_ACC/N_ATT'])/nsamples])

    # Take care of the correlation between numerator and denominator
    # in Hartree-Fock estimates.
    results = (
            results.join(analyse_hf_observables(means, covariances, nsamples)))

    return results

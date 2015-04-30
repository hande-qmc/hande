'''Analyse output of canonical free-particle thermodynamic calculation.'''

import pandas as pd
import pyblock
import pyhande
import numpy as np

def estimates(metadata, data):
    '''Perform error analysis for canonical thermodynamic estimates.

Parameters
----------
metadata : :class:`pandas.DataFrame`
    metadata (i.e. calculation information, parameters and settings) extracted
    from output files.
data : :class:`pandas.DataFrame`
    HANDE QMC data.

Returns
-------
results : :class:`pandas.DataFrame`
    Averaged estimates.
'''

    num = r'\sum\rho_HF_{ii}H_{ii}'
    denom = r'\sum\rho_HF_{ii}'

    means = data.mean()
    covariances = data.cov()
    nsamples = len(data[num])

    numerator = pd.DataFrame()
    numerator['mean'] = [means[num]]
    numerator['standard error'] = [np.sqrt(covariances[num][num]/nsamples)]
    denominator = pd.DataFrame()
    denominator['mean'] = [means[denom]]
    denominator['standard error'] = [np.sqrt(covariances[denom][denom]/nsamples)]
    cov_thf = covariances[num][denom]
    # The numerator and denominator are correlated for the
    # HF estimate for the total energy.
    e_thf = pyblock.error.ratio(numerator, denominator, cov_thf, nsamples)
    e_thf.reset_index(inplace=True)

    results = pd.DataFrame()
    results['Beta'] = pyhande.extract.extract_input(metadata, 'beta')
    # E_0 and E_HF0 contain no denominator so the error is
    # just the standard error.
    results['E_0'] = [means['E_0']]
    results['E_0-Error'] = [np.sqrt(covariances['E_0']['E_0']/nsamples)]
    results['E_HF0'] = [means['E_HF0']]
    results['E_HF0-Error'] = [np.sqrt(covariances['E_HF0']['E_HF0']/nsamples)]
    results['E_THF'] = list(e_thf['mean'])
    results['E_THF-Error'] = list(e_thf['standard error'])

    return (results)

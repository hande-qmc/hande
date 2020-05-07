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
results : :class:`pandas.Series`
    Averaged Hartree-Fock estimates along with error estimates.
'''

    observables = dict([
        ('T_HF', r'Tr(T\rho_HF)'),
        ('V_HF', r'Tr(V\rho_HF)'),
        ('U_HF', r'Tr(H\rho_HF)'),
    ])

    num = pd.Series({'mean': 0.0, 'standard error': 0.0})
    trace = pd.Series({
        'mean': means[r'Tr(\rho_HF)'],
        'standard error': np.sqrt(covariances[r'Tr(\rho_HF)'][r'Tr(\rho_HF)']/nsamples)
    })
    results = {}

    for (k, v) in observables.items():
        num['mean'] = means[v]
        num['standard error'] = np.sqrt(covariances[v][v]/nsamples)
        cov_ab = covariances[v][r'Tr(\rho_HF)']

        stats = pyblock.error.ratio(num, trace, cov_ab, nsamples)

        results[k] = stats['mean']
        results[k+'_error'] = stats['standard error']

    return pd.Series(results)


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
results : :class:`pandas.Series`
    Averaged estimates.
'''

    data['U_0'] = data['<T>_0'] + data['<V>_0']
    data[r'Tr(H\rho_HF)'] = data[r'Tr(T\rho_HF)'] + data[r'Tr(V\rho_HF)']

    ncycles = len(data) # For standard error.
    # Ensure some backwards compatability.
    if 'N_ACC/N_ATT' in data.columns:
        # Work out weights given that number of configurations generated differs
        # between cycles.
        w = data['N_ACC/N_ATT'] * metadata['nattempts']
    else:
        w = np.ones(ncycles)
    # Number of configurations.
    nsamples = sum(w)
    # Normalise the weights.
    w = w / nsamples
    # Weighted estimate for the means.
    means = (
        pd.Series(np.average(data, axis=0, weights=w), index=list(data.keys()))
    )
    # Weighted estimate for the covariance. See
    # https://en.wikipedia.org/wiki/Weighted_arithmetic_mean and
    # http://stats.stackexchange.com/questions/61225/correct-equation-for-weighted-unbiased-sample-covariance
    # for more details.
    w2 = sum(w**2.0)
    xm = data.sub(means, axis=1)
    covariances = 1.0/(1.0-w2) * xm.mul(w, axis=0).T.dot(xm)

    if 'beta' in metadata:
        # New, richer JSON-based metadata.
        beta = metadata['beta']
    else:
        # Hope to find it in the input file...
        # WARNING: This is not tested for.
        beta = pyhande.legacy.extract_input(metadata, 'beta')
    results = {
        'Beta': beta,
        # Free estimates contain no denominator so the error is
        # just the standard error.
        'U_0': means['U_0'],
        'U_0_error': np.sqrt(covariances['U_0']['U_0']/ncycles),
        'T_0': means['<T>_0'],
        'T_0_error': np.sqrt(covariances['<T>_0']['<T>_0']/ncycles),
        'V_0': means['<V>_0'],
        'V_0_error': np.sqrt(covariances['<V>_0']['<V>_0']/ncycles),
    }
    if 'N_ACC/N_ATT' in data.columns:
        if metadata['fermi_temperature']:
            beta = results['Beta'] / metadata['system']['ueg']['E_fermi']
        else:
            beta = results['Beta']
        correction = metadata['free_energy_corr']
        results.update({
            'N_ACC/N_ATT': means['N_ACC/N_ATT'],
            'N_ACC/N_ATT_error': (
                np.sqrt(covariances['N_ACC/N_ATT']['N_ACC/N_ATT']/ncycles)
            ),
        })
        results['F_0'] = (
                (-1.0/beta)*np.log(results['N_ACC/N_ATT']) + correction
        )
        # Normal error estimate for natural logarithm.
        results['F_0_error'] = (
                    results['N_ACC/N_ATT_error']/(beta*results['N_ACC/N_ATT'])
        )
        results['S_0'] = beta * (results['T_0']-results['F_0'])
        # T_0 and F are correlated so DS = ((DT_0^2 + DF^2 + 2 COV(T_0,F))/T)^.5
        # COV(x,f(y)) = 1/(N*N-1) sum_{i<j}(x_i-X)(y_i-Y)df/dy|y=Y
        # df/dy = -kT/(N_ACC/N_ATT)
        results['S_0_error'] = (beta*np.sqrt(results['T_0_error']**2.0 +
                                results['F_0_error']**2.0 -
                                2.0*covariances['N_ACC/N_ATT']['<T>_0'] /
                                (ncycles*results['N_ACC/N_ATT']*beta)))
    results = pd.Series(results)

    # Take care of the correlation between numerator and denominator
    # in Hartree-Fock estimates.
    results = results.append(analyse_hf_observables(means, covariances, ncycles))

    return results

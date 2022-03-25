'''Analysis of data from dmqmc calculations.'''

import pandas as pd
import numpy as np
import warnings
import scipy.interpolate
import scipy.integrate
import pyblock
import pyhande

def analyse_observables(means, covariances, nsamples):
    '''Calculate all mean and error estimates of the form Tr(\rho O)/Tr(\rho).

Parameters
----------
means : :class:`pandas.DatFrame`
    The mean estimates, obtained by averaging over beta loops, as a function of
    beta, which is used as the index.
covariances : :class:`pandas.DataFrame`
    Estimates of the covariance matrix as a function of the beta, which is used
    as the index.
nsamples : :class:`pandas.Series`
    The number of samples contributing to the various beta values.

Returns
-------
results : :class:`pandas.DataFrame`
    Contains all of the final mean and error estimates.
None.
'''
    columns = list(means.columns.values)
    beta_values = means.index.values

    # The keys hold the names to be output and the values hold the possible
    # columns names in the means and covariances DataFrames.
    observables = dict([
        ('Tr[Hp]/Tr[p]','\\sum\\rho_{ij}H_{ji}'),
        ('Tr[H2p]/Tr[p]','\\sum\\rho_{ij}H2{ji}'),
        ('Tr[Sp]/Tr[p]','\\sum\\rho_{ij}S_{ji}'),
        ('Tr[Mp]/Tr[p]','\\sum\\rho_{ij}M2{ji}'),
        ('Tr[Tp]/Tr[p]','\\sum\\rho_{ij}T_{ji}'),
        ('Tr[Up]/Tr[p]','\\sum\\rho_{ij}U_{ji}'),
        ('Tr[H0p]/Tr[p]','\\sum\\rho_{ij}H0{ji}'),
        ('Tr[HIp]/Tr[p]','\\sum\\rho_{ij}HI{ji}'),
        ('VI', '\\sum\\rho_{ij}VI_{ji}'),
        ('Re{Tr[Hp]/Tr[p]}', r'Re{\sum \rho H}'),
        ('Im{Tr[Hp]/Tr[p]}', r'Re{\sum \rho H}'),
        ('Tr[p0H0]/Tr[p00]', r'\sum\rho_{0j}H_{j0}'),
        ('Re{Tr[p0H0]/Tr[p00]}', r'Re{Sum\rho_0j H_j0}'),
        ('Im{Tr[p0H0]/Tr[p00]}', r'Im{Sum\rho_0j H_j0}'),
    ])
    # A dictionary of all trace keys, so alternative keys are accessible
    denominators = dict([
        ('Tr[Hp]/Tr[p]','Trace'),
        ('Tr[H2p]/Tr[p]','Trace'),
        ('Tr[Sp]/Tr[p]','Trace'),
        ('Tr[Mp]/Tr[p]','Trace'),
        ('Tr[Tp]/Tr[p]','Trace'),
        ('Tr[Up]/Tr[p]','Trace'),
        ('Tr[H0p]/Tr[p]','Trace'),
        ('Tr[HIp]/Tr[p]','Trace'),
        ('VI', 'Trace'),
        ('Re{Tr[Hp]/Tr[p]}', r'Re{Trace}'),
        ('Im{Tr[Hp]/Tr[p]}', r'Im{Trace}'),
        ('Tr[p0H0]/Tr[p00]', r'\rho_00'),
        ('Re{Tr[p0H0]/Tr[p00]}', r'Re{\rho_00}'),
        ('Im{Tr[p0H0]/Tr[p00]}', r'Im{\rho_00}'),
    ])

    # Add momentum distribution to dictionary of observables to be analaysed.
    add_observable_to_dict(observables, columns, 'n_')
    # Add spin-averaged static structure factor to dict of observables to be analaysed.
    add_observable_to_dict(observables, columns, 'S_')
    # Add spin up static structure factor to dict of observables to be analaysed.
    add_observable_to_dict(observables, columns, 'Suu_')
    # Add spin down static structure factor to dict of observables to be analaysed.
    add_observable_to_dict(observables, columns, 'Sud_')
    # DataFrame to hold the final mean and error estimates.
    results = pd.DataFrame(index=beta_values)
    # DataFrame for the numerator.
    num = pd.DataFrame(columns=['mean','standard error'], index=beta_values)
    # DataFrame for the trace from the first replica.
    tr1 = pd.DataFrame(columns=['mean','standard error'], index=beta_values)

    for (k,v) in observables.items():
        if v in columns:
            den_key = denominators[k]

            num['mean'] = means[v]
            num['standard error'] = np.sqrt(covariances.xs(v,level=1)[v]/nsamples)
            tr1['mean'] = means[den_key]
            tr1['standard error'] = np.sqrt(covariances.xs(den_key,level=1)[den_key]/nsamples)
            cov_AB = covariances.xs(den_key,level=1)[v]

            stats = pyblock.error.ratio(num, tr1, cov_AB, nsamples)

            results[k] = stats['mean']
            results[k+'_error'] = stats['standard error']

    return results


def add_observable_to_dict(observables, columns, label):
    '''Add observable to dictionary of analysed data.

    This sets the value to be the same as the key specified by label variable.

Parameters
----------
observables : dict
    Dictionary of observables to be analysed.
columns : list
    Columns in hande output.
label : string
    regex to search for in list of columns.
'''

    new_obs = [c for c in columns if label in c]
    observables.update(dict(zip(new_obs, new_obs)))


def free_energy_error_analysis(data, results, dtau):
    '''Calculate the mean and error estimates for the exchange correlation
       free energy (appropriately defined)

Parameters
----------
data : :class:`pandas.DataFrame`
    Raw dmqmc data.
results : :class:`pandas.DataFrame`
    The mean estimates, obtained by averaging over beta loops, as a function of
    beta, which is used as the index. On output estimates for f_xc and its error
    are appended.
dtau : float
    Time step used in simulation.
'''
    # Integral is evaulated as cumulative integral of integral for each
    # temperature/imaginary time value. Naturally I(tau=0) = 0.
    I = scipy.integrate.cumtrapz(results['VI'], dx=dtau, initial=0)

    # An estimate for the error can be found by considering the variance of each
    # simulations estimate for the Integral.
    grouped = pyhande.utils.groupby_beta_loops(data.reset_index(), name='Beta')
    I_single = [scipy.integrate.cumtrapz(d[r'\sum\rho_{ij}VI_{ji}']/d['Trace'],
                dx=dtau, initial=0) for (i, d) in grouped]

    # Place cumulative integrals in frame so variance can be easily calculated.
    frame = pd.DataFrame(I_single)
    I_single_mean = frame.mean()
    I_error = ((frame.var() / len(grouped)))**0.5

    # Some sort of check for poor estimation of mean/error.
    # Check if last point's 'ratio' error (not quite the normal definition) is
    # significant (so comparible to the stochastic error).
    err = abs(I[-1]-I_single_mean.iloc[-1])
    if (abs(err) > I_error.iloc[-1]):
        warnings.warn("Ratio error of %f is not insignificant - check results."
                      %err)

    # Need to multiply integral by kT. Any scaling factor drops out due to
    # presence in timestep as well.
    results.reset_index(drop=True, inplace=True)
    results['f_xc'] = I / results['Beta'].iloc[-1]
    results['f_xc_error'] = I_error/results['Beta'].iloc[-1]


def analyse_renyi_entropy(means, covariances, nsamples):
    '''Calculate the mean and error estimates for the Renyi entropy (S2), for
       all subsystems, including the entire system if present.

Parameters
----------
means : :class:`pandas.DatFrame`
    The mean estimates, obtained by averaging over beta loops, as a function of
    beta, which is used as the index.
covariances : :class:`pandas.DataFrame`
    Estimates of the covariance matrix as a function of the beta, which is used
    as the index.
nsamples : :class:`pandas.Series`
    The number of samples contributing to the various beta values.

Returns
-------
results : :class:`pandas.DataFrame`
    Contains all of the final mean and error estimates.
None.
'''

    columns = list(means.columns.values)
    beta_values = means.index.values

    # DataFrame to hold the final mean and error estimates.
    results = pd.DataFrame(index=beta_values)
    # DataFrame for the numerator.
    num = pd.DataFrame(columns=['mean','standard error'], index=beta_values)
    # DataFrame for the trace from the first replica.
    tr1 = pd.DataFrame(columns=['mean','standard error'], index=beta_values)
    # DataFrame for the trace from the second replica.
    tr2 = pd.DataFrame(columns=['mean','standard error'], index=beta_values)

    # Compute the mean and standard error estimates for the Renyi entropy (S2),
    # for all subsystems.
    nrdms = 0
    for column in columns:
        have_s2 = False
        if 'RDM' in column and 'S2' in column: # RDM S2
            nrdms += 1
            have_s2 = True
            num_col = column
            tr1_col = 'RDM'+str(nrdms)+' trace 1'
            tr2_col = 'RDM'+str(nrdms)+' trace 2'
            out_str = 'RDM'+str(nrdms)+' S2'
        elif 'Full S2' in column: # Full S2
            have_s2 = True
            num_col = column
            tr1_col = 'Trace'
            tr2_col = 'Trace 2'
            out_str = 'Full S2'

        if have_s2:
            num['mean'] = means[num_col]
            tr1['mean'] = means[tr1_col]
            tr2['mean'] = means[tr2_col]
            num['standard error'] = np.sqrt(covariances.xs(num_col,level=1)[num_col]/nsamples)
            tr1['standard error'] = np.sqrt(covariances.xs(tr1_col,level=1)[tr1_col]/nsamples)
            tr2['standard error'] = np.sqrt(covariances.xs(tr2_col,level=1)[tr2_col]/nsamples)

            # A denotes the numerator, B denotes the first trace and C denotes
            # the second trace.
            cov_AB = covariances.xs(num_col,level=1)[tr1_col]
            cov_AC = covariances.xs(num_col,level=1)[tr2_col]
            cov_BC = covariances.xs(tr1_col,level=1)[tr2_col]

            results[out_str], results[out_str+' error'] = \
                calc_S2(num, tr1, tr2, cov_AB, cov_AC, cov_BC, nsamples)

    return results


def calc_S2(stats_A, stats_B, stats_C, cov_AB, cov_AC, cov_BC, data_len):
    '''Calculate the mean and standard error of :math:`f = -log(A/BC)`.

Parameters
----------
stats_A : :class:`pandas.DataFrame`
    Statistics (containing at least the 'mean' and 'standard error' fields) for
    variable :math:`A`.  The rows contain different values of these statistics
    for different beta values.
stats_B : :class:`pandas.DataFrame`
    Similarly for variable :math:`B`.
stats_C : :class:`pandas.DataFrame`
    Similarly for variable :math:`C`.
cov_AB : :class:`pandas.Series`
    Covariance between variables :math:`A` and :math:`B`.
cov_AC : :class:`pandas.Series`
    Covariance between variables :math:`A` and :math:`C`.
cov_BC : :class:`pandas.Series`
    Covariance between variables :math:`B` and :math:`C`.
data_len : :class:`pandas.Series`
    Number of beta loops used to obtain the statistics given in ``stats_A``,
    ``stats_B`` and ``stats_C``.

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
    Must contain a column with the above name, and another with the name
    column+' error'. The indices are used for the beta values.

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

    return pd.Series(spline_fit, index=beta_values)


def sort_momentum(columns):
    ''' Naturally sort results columns based off of their numeric values.

    e.g. give columns = [n_1, n_10, n_5], return [n_1, n_5, n_10]

Parameters
----------
columns : list
    list of column names

Returns
-------
sort : list
    Sorted list of column names.
'''

    values = [float(c.split('_')[1]) for c in columns]

    sort = [c for (v, c) in sorted(zip(values, columns))]

    return sort


def extract_momentum_correlation_function(data, beta):
    '''Extract correlation function from analysed DMQMC data.

Parameters
----------
data : :class:`pandas.DataFrame`
    Analysed DMQMC output at a specific temperature / imaginary time.
beta : float
    Temperature / imaginary time value to extrat correlation function at.

Returns
-------
full : :class:`pandas.DataFrame`
    Momentum space correlation functions as a function of k at given
    temperature.
'''

    dbeta = data.groupby('Beta').get_group(beta)
    base_headers = ['n_', 'S_', 'Suu_', 'Sud_']

    frames = []

    for h in base_headers:

        mom = [c for c in dbeta.columns.values if h in c and '_error' not in c]
        mome = [c for c in dbeta.columns.values if h in c and '_error' in c]

        kvals = [float(c.split('_')[1]) for c in mom]

        # The columns in the output do not contain the factors of 1
        # necessary for the upup+downdown contributions or the spin-averaged
        # contribution.
        if h == 'S_' or h == 'Suu_':
            kspace_fn = 1.0 + dbeta[mom].values
        else:
            kspace_fn = dbeta[mom].values
        kspace_fn_error = dbeta[mome].values
        kspace_fn_error[np.isnan(kspace_fn_error)] = 0
        frames.append(pd.DataFrame({'Beta': beta, 'k': kvals,
                                    h+'k': kspace_fn.ravel(),
                                    h+'k_error': kspace_fn_error.ravel()}))


    # Want to perform merge so that the correlation function with the most k points
    # defines the total number of k points in the final frame.
    frames.sort(key=lambda x: len(x.k))
    full = frames[0]
    for f in frames[1:]:
        full = full.merge(f, on=['Beta', 'k'], how='right')

    return full


def analyse_data(hande_out, shift=False, free_energy=False, spline=False,
                 trace=False, calc_number=None, energy_numerator=False,
                 population=False, beta_loop_count=False, very_very_loud=False):
    '''Clean up Hande output so that analysis can be performed.

Parameters
----------
hande_out : list of (dict, :class:`pandas.DataFrame` or :class:`pandas.Series`)
    Calculation output represented by a tuple for each calculation, consisting
    of metadata (dict) and a :class:`pandas.DataFrame` (MC calculations) or
    :class:`pandas.Series` (other calculations) containing the calculation
    output/results.
shift : bool
    Perform analysis on the shift.
free_energy : bool
    Perform free energy error analysis, requires appropriate column in HANDE
    output.
spline : bool
    Fit splines to analysed data.
trace : bool
    Perform analysis on the trace of the density matrix.
calc_number : int or None
    If not None then only perform analysis on the calc_number calculation.

Returns
-------
metadata : dict
    Calculation metadata for first simulation in hande_out (or user specified
    calculation number in options).
results : :class:`pandas.DataFrame`
    Analysed DMQMC data.
'''
    if spline:
        raise ValueError("spline option is broken! Please set to default "
                         "(False) or fix code.")
    if calc_number:
        raise ValueError("calc_number option is broken! Please set to default "
                         "(None) or fix code.")

    if very_very_loud:
        shift, trace, energy_numerator, population, beta_loop_count = [True]*5

    (metadata, data) = ([], [])
    for (md, df) in hande_out:
        if 'DMQMC' in md['calc_type'] and not md['dmqmc']['find_weights']:
            metadata.append(md)
            # Convert the iteration number to the beta value.
            cycles = md['qmc']['ncycles']
            tau = md['qmc']['tau']
            df.rename(columns={'iterations' : 'Beta'}, inplace=True)
            df['Beta'] = df['Beta']*tau
            data.append(df)
    if data:
        if calc_number is not None:
            data = data[calc_number]
            metadata = metadata[calc_number]
        else:
            data = pd.concat(data)

    # Sanity check: Same time step used in all calculations?
    # Only check for new metadata format...
    if 'qmc' in metadata[0]:
        tau = metadata[0]['qmc']['tau']
        for md in metadata[1:]:
            if 'tau' in md['qmc'] and md['qmc']['tau'] != tau:
                warnings.warn('Tau values in input files not consistent.')

    # Discard beta loops which didn't reach final iteration.
    grouped = pyhande.utils.groupby_beta_loops(data, name='Beta')
    niterations = len(grouped.get_group(0))
    last_group = len(grouped) - 1
    if len(grouped.get_group(last_group)) < niterations:
        data.drop(grouped.get_group(last_group).index, inplace=True)
    # Make the Beta column a MultiIndex.
    data.set_index('Beta', inplace=True, append=True)
    # The number of beta loops contributing to each beta value.
    if 'Trace' in data.columns:
        nsamples = data['Trace'].groupby(level=1).count()
    else:
        nsamples = data['Re{Trace}'].groupby(level=1).count()
    beta_values = nsamples.index.values
    # The data that we are going to use.
    colnames = [colname for colname in list(data.columns) 
                if colname not in ['# spawn_events','# states','R_spawn','time']]
    estimates = data.loc[:,colnames]

    if free_energy:
        # Deal with legacy input.
        if 'propagate_to_beta' in metadata[0]['ipdmqmc']:
            ipdmqmc = metadata[0]['ipdmqmc']['propagate_to_beta']
        else:
            ipdmqmc = metadata[0]['ipdmqmc']['ipdmqmc']
        # Set up estimator we need to integrate wrt time/temperature to evaluate
        # free energy difference.
        if ipdmqmc:
            if metadata[0]['ipdmqmc']['symmetric']:
                estimates[r'\sum\rho_{ij}VI_{ji}'] = (
                                                 data[r'\sum\rho_{ij}HI{ji}'] -
                                                 data[r'\sum\rho_{ij}H0{ji}']
                                                 )
            else:
                estimates[r'\sum\rho_{ij}VI_{ji}'] = (
                                                  data[r'\sum\rho_{ij}H_{ji}'] -
                                                  data[r'\sum\rho_{ij}H0{ji}']
                                                  )
        else:
            estimates[r'\sum\rho_{ij}VI_{ji}'] = data[r'\sum\rho_{ij}H_{ji}']


    # Compute the mean of all estimates across all beta loops.
    means = estimates.groupby(level=1).mean()
    # Compute the covariances between all pairs of columns.
    covariances = estimates.groupby(level=1).cov()

    columns = list(estimates.columns.values)

    # results will hold all of the final values to be printed.
    results = pd.DataFrame({'Beta' : pd.Series(beta_values, index=beta_values)})
    results = results.join(analyse_observables(means, covariances, nsamples))
    results = results.join(analyse_renyi_entropy(means, covariances, nsamples))

    # If requested, add the averaged shift profile to results.
    if shift:
        results['Shift'] = means['Shift']
        results['Shift s.d.'] = np.sqrt(covariances.xs('Shift',level=1)['Shift'])

    # If requested, calculate a spline fit for all mean estimates.
    if spline:
        results_columns = list(results.columns.values)
        for column in results_columns:
            if ('Tr[p]' in column or 'S2' in column) and (not 'error' in column):
                results[column+' spline'] = calc_spline_fit(column, results)

    # If requested, calculate excess free-energy.
    if free_energy:
        free_energy_error_analysis(estimates, results, cycles*tau)

    # If requested, add the averaged trace profiles to results.
    if trace:
        trace_keys = [
                'Trace',
                'Trace 2',
                'Re{Trace}',
                'Im{Trace}',
                r'\rho_00',
                r'Re{\rho_00}',
                r'Im{\rho_00}',
            ]
        for key in trace_keys:
            if key in columns:
                results[key] = means[key]
                results[key+' s.d.'] = np.sqrt(covariances.xs(key,level=1)[key])

    # If requested, return the averaged energy numerator profile to the results.
    if energy_numerator:
        numerator_key_val = {
                r'\sum\rho_{ij}H_{ji}'  : 'Tr[Hp]',
                r'Re{\sum \rho H}'      : 'Re{Tr[Hp]}',
                r'Im{\sum \rho H}'      : 'Im{Tr[Hp]}',
                r'\sum\rho_{0j}H_{j0}'  : r'Tr[p0H0]',
                r'Re{Sum\rho_0j H_j0}'  : r'Re{Tr[p0H0]}',
                r'Im{Sum\rho_0j H_j0}'  : r'Im{Tr[p0H0]}',
            }
        for (k,v) in numerator_key_val.items():
            if k in columns:
                results[v] = means[k]
                results[v+' s.d.'] = np.sqrt(covariances.xs(k,level=1)[k])

    # If requested, return the averaged walker population profile to the results.
    if population:
        results['# particles'] = means['# H psips']
        results['# particles s.d.'] = np.sqrt(covariances.xs('# H psips',level=1)['# H psips'])
        if r'# \rho_{0j} psips' in columns:
            results[r'# \rho_{0j} psips'] = means[r'# \rho_{0j} psips']
            results[r'# \rho_{0j} psips s.d.'] = \
                np.sqrt(covariances.xs(r'# \rho_{0j} psips',level=1)[r'# \rho_{0j} psips'])

    # If requested, return the number of beta loops in the averaging profile.
    if beta_loop_count:
        results['N_beta'] = nsamples

    return (metadata, results)

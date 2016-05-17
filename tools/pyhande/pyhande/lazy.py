'''Tools for the lazy amongst us: automation of common HANDE analysis tasks.'''

import collections
import math
from os import path
import pkgutil
import sys
import warnings

import matplotlib.pyplot as plt
import pandas as pd

if pkgutil.find_loader('pyblock'):
    sys.path.append(path.join(path.abspath(path.dirname(__file__)), '../../pyblock'))
import pyblock
import pyhande.extract
import pyhande.analysis
import pyhande.weight

def std_analysis(datafiles, start=0, select_function=None, extract_psips=False,
                reweight_history=0, mean_shift=0.0, arith_mean=False,
                calc_inefficiency=False ):
    '''Perform a 'standard' analysis of HANDE output files.

Parameters
----------
datafiles : list of strings
    names of files containing HANDE QMC calculation output.
start : int
    iteration from which the blocking analysis is performed.
select_function : function 
    function which returns a boolean mask for the iterations to include in the
    analysis.  Not used if set to None (default).  Overrides ``start``.  See
    below for examples.
extract_psips : bool
    also extract the mean number of psips from the calculation.
reweight_history : integer
    reweight in an attempt to remove population control bias. According to
    [Umrigar93]_ this should be set to be a few correlation times.
mean_shift : float
    prevent the weights from beoming to large.

Returns
-------
info : list of :func:`collections.namedtuple`
    raw and analysed data, consisting of:

        metadata, data
            from :func:`pyhande.extract.extract_data_sets`.  If ``data``
            consists of several concatenated calculations, then the only
            ``metadata`` object is from the first calculation.
        data_len, reblock, covariance
            from :func:`pyblock.pd_utils.reblock`.  The projected energy
            estimator (evaluated by :func:`pyhande.analysis.projected_energy`)
            is included in ``reblock``.
        opt_block, no_opt_block
            from :func:`pyhande.analysis.qmc_summary`.  A 'pretty-printed'
            estimate string is included in ``opt_block``.

Examples
--------

The following are equivalent and will extract the data from the file called
hande.fciqmc.out, perform a blocking analysis from the 10000th iteration
onwards, calculated the projected energy estimator and find the optimal block
size from the blocking analysis:

>>> std_analysis(['hande.fciqmc.out'], 10000)
>>> std_analysis(['hande.fciqmc.out'],
...              select_function=lambda d: d['iterations'] > 10000)

References
----------
Umrigar93
    Umrigar et al., J. Chem. Phys. 99, 2865 (1993).
'''
    (calcs, calcs_md) = zeroT_qmc(datafiles, reweight_history, mean_shift,
                                  arith_mean)
    return [lazy_block(calc, md, start, select_function, extract_psips,
                      calc_inefficiency) for (calc, md) in zip(calcs, calcs_md)]

def zeroT_qmc(datafiles, reweight_history=0, mean_shift=0.0, arith_mean=False):
    '''Extract zero-temperature QMC (i.e. FCIQMC and CCMC) calculations.

Reweighting information is added to the calculation data if requested.

.. note::

    :func:`std_analysis` is recommended unless custom processing is required
    before blocking analysis is performed.

Parameters
----------
datafiles, reweight_history, mean_shift, arith_mean :
    See :func:`std_analysis`.

Returns
-------
calcs : list of :class:`pandas.DataFrame`
    Calculation outputs for just the zero-temperature/ground-state QMC
    calculations contained in `datafiles`.
metadata : list of dict
    Metadata corresponding to each calculation in `calcs`.
'''

    hande_out = pyhande.extract.extract_data_sets(datafiles)

    # Concat all QMC data (We did say 'lazy', so assumptions are being made...)
    data = []
    metadata = []
    for (md, df) in filter_calcs(hande_out, ('FCIQMC', 'CCMC')):
        if reweight_history > 0:
            df = pyhande.weight.reweight(df, md['mc_cycles'],
                md['qmc']['tau'], reweight_history, mean_shift,
                arith_mean=arith_mean)
            df['W * \sum H_0j N_j'] = df['\sum H_0j N_j'] * df['Weight']
            df['W * N_0'] = df['N_0'] * df['Weight']
        data.append(df)
        metadata.append(md)
    if data:
        calcs_metadata, calcs = concat_calcs(metadata, data)
    else:
        raise ValueError('No data found in '+' '.join(datafiles))
    return (calcs, calcs_metadata)

def lazy_block(calc, md, start=0, select_function=None, extract_psips=False,
               calc_inefficiency=False):
    '''Standard blocking analysis on zero-temperature QMC calcaulations.

.. note::

    :func:`std_analysis` is recommended unless custom processing is required
    before blocking analysis is performed.

Parameters
----------
calc : :class:`pandas.DataFrame`
    Zero-temperature QMC calculation output.
md : dict
    Metadata for the calculation in `calc`.
start, select_function, extract_psips, calc_inefficiency:
    See :func:`std_analysis`.

Returns
--------
info : :func:`collections.namedtuple`
    See :func:`std_analysis`.
'''

    tuple_fields = ('metadata data data_len reblock covariance opt_block '
                   'no_opt_block'.split())
    info_tuple = collections.namedtuple('HandeInfo', tuple_fields)
    # Reblock Monte Carlo data over desired window.
    reweight_calc = 'W * N_0' in calc
    if select_function is None:
        indx = calc['iterations'] > start
    else:
        indx = select_function(calc)
    to_block = []
    if extract_psips:
        to_block.append('# H psips')
    to_block.extend(['\sum H_0j N_j', 'N_0', 'Shift'])
    if reweight_calc:
        to_block.extend(['W * \sum H_0j N_j', 'W * N_0'])

    mc_data = calc.ix[indx, to_block]

    if mc_data['Shift'].iloc[0] == mc_data['Shift'].iloc[1]:
        warnings.warn('The blocking analysis starts from before the shift '
                     'begins to vary.')

    (data_len, reblock, covariance) = pyblock.pd_utils.reblock(mc_data)
    
    proje = pyhande.analysis.projected_energy(reblock, covariance, data_len)
    reblock = pd.concat([reblock, proje], axis=1)
    to_block.append('Proj. Energy')

    if reweight_calc:
        proje = pyhande.analysis.projected_energy(reblock, covariance, 
                    data_len, sum_key='W * \sum H_0j N_j', ref_key='W * N_0',
                    col_name='Weighted Proj. E.')
        reblock = pd.concat([reblock, proje], axis=1)
        to_block.append('Weighted Proj. E.')
    
    # Summary (including pretty printing of estimates).
    (opt_block, no_opt_block) = pyhande.analysis.qmc_summary(reblock, to_block)

    if calc_inefficiency:
        # Calculate quantities needed for the inefficiency.
        dtau = md['qmc']['tau']
        reblocked_iters = calc.ix[indx, 'iterations']
        N = reblocked_iters.iloc[-1] - reblocked_iters.iloc[0]
        
        # This returns a data frame with inefficiency data from the
        # projected energy estimators if available.
        ineff = pyhande.analysis.inefficiency(opt_block, dtau, N)
        if ineff is not None:
            opt_block = opt_block.append(ineff)

    estimates = []
    for (name, row) in opt_block.iterrows():
        estimates.append(
                pyblock.error.pretty_fmt_err(row['mean'], row['standard error'])
                       )
    opt_block['estimate'] = estimates
    info = info_tuple(md, calc, data_len, reblock, covariance, opt_block,
                      no_opt_block)

    return info

def filter_calcs(outputs, calc_types):
    '''Select calculations corresponding to a given list of calculation types.

Parameters
----------
outputs : list of (dict, :class:`pandas.DataFrame` or :class:`pandas.Series`)
    List of (metadata, data) tuples for each calculation, as created in
    :func:`pyhande.extract.extract_data_sets`.
calc_types : iterable of strings
    Calculation types (e.g. 'FCIQMC', 'CCMC', etc.) to select.

Returns
-------
filtered : list of (dict, :class:`pandas.DataFrame` or :class:`pandas.Series`)
    As in :func:`pyhande.extract.extract_data_sets` but containing only the
    desired calculations.
'''

    calc_filter = lambda md: any(md['calc_type'] for calc in calc_types)
    filtered = [(md, calc) for (md, calc) in outputs if calc_filter(md)]
    return filtered

def concat_calcs(metadata, data):
    '''Concatenate data from restarted calculations to analyse together.

Parameters
----------
metadata : list of dicts
    Extracted metadata for each calculation.
data : list of :class:`pandas.DataFrame`
    Output of each QMC calculation.

Returns
-------
calcs_metadata : list of dicts
    Metadata for each calculation, with duplicates from restarting dropped.
calcs : list of :class:`pandas.DataFrame`
    Output of each QMC calculation, with parts of a restarted calculation combined.
'''

    restart_uuids = [md['restart'].get('uuid_restart','') for md in metadata]
    uuids = [md['UUID'] for md in metadata]
    if any(restart_uuids) and all(uuids):
        data = list(data)
        metadata = list(metadata)
        calcs = []
        calcs_metadata = []
        while uuids:
            for indx in range(len(uuids)):
                if uuids[indx] not in restart_uuids:
                    # Found the end of a chain.
                    break
            uuid = uuids.pop(indx)
            restart = restart_uuids.pop(indx)
            calc = [data.pop(indx)]
            calcs_metadata.append(metadata.pop(indx))
            while restart and restart in uuids:
                indx = uuids.index(restart)
                uuid = uuids.pop(indx)
                restart = restart_uuids.pop(indx)
                calc.append(data.pop(indx))
                metadata.pop(indx)
            calcs.append(pd.concat(calc[::-1]))

    else:
        # Don't have UUID information.
        # Assume any restarted calculations are in the right order and contiguous
        # Check concatenating data is at least possibly sane.
        step = data[0]['iterations'].iloc[-1] - data[0]['iterations'].iloc[-2]
        prev_iteration = data[0]['iterations'].iloc[-1]
        calc_type = metadata[0]['calc_type']
        calcs = []
        calcs_metadata = [metadata[0]]
        xcalc = [data[0]]
        for i in range(1, len(data)):
            if metadata[i]['calc_type'] != calc_type or \
                    data[i]['iterations'].iloc[0] - step != prev_iteration or \
                    data[i]['iterations'].iloc[-1] - data[i]['iterations'].iloc[-2] != step:
                # Different calculation
                step = data[i]['iterations'].iloc[-1] - data[i]['iterations'].iloc[-2]
                calc_type = metadata[i]['calc_type']
                calcs.append(pd.concat(xcalc))
                xcalc = [data[i]]
                calcs_metadata.append(metadata[i])
            else:
                # Continuation of same calculation (probably)
                xcalc.append(data[i])
            prev_iteration = data[i]['iterations'].iloc[-1]
        calcs.append(pd.concat(xcalc, ignore_index=True))
    return calcs_metadata, calcs

def find_starting_iteration(data, md, frac_screen_interval=500,
    number_of_reblockings=50, number_of_reblocks_to_cut_off=1, pos_min_frac=0.5,
    verbose=False, show_graph=False):
    '''Find the best iteration to start analysing CCMC/FCIQMC data.

.. warning::

    Use with caution, check whether output is sensible and adjust parameters if
    necessary.

First, consider only data from when the shift begins to vary.  The error in the
error on the shift is artificially high due to the equilibration period.  The
error in the error hence initially decreases as the start of the blocking
analysis is increased before beginning to increase once data from after
equilibration begins to be discarded.  We are thus interested in finding the
minimum in the error in the error of the shift.

The best estimate of the iteration to start the blocking analysis is found by:

1. discard data during the constant shift phase.
2. estimate the error in the error of the shift by blocking the remaining data
   :math:`n` times, where the blocking analysis considers the last :math:`1-i/f`
   fraction of the data and where :math:`i` is the number of blocking analyses
   already performed, :math:`n`  is `number_of_reblockings`  and :math:`f` is
   `frac_screen_interval`.
3. find the iteration which gives the minimum estimate of the error in the error
   of the shift.  If this is in the first `pos_min_frac` fraction of the
   blocking attempts, go to 4, otherwise repeat 2 and perform an additional
   `number_of_reblockings` attempts.
4. To be extra conservative, discard the first `number_of_reblocks_to_cut_off`
   blocks from the start iteration, where each block corresponds to roughly the
   autocorrelation time, and return the resultant iteration number as the
   estimate of the best place to start blocking from.

Parameters
----------
data : :class:`pandas.DataFrame`
    Calculation output for a FCIQMC or CCMC calculation.
md : dict
    Metadata corresponding to the calculation in `data`.
frac_screen_interval : int
    Number of intervals the iterations from where the shift started to vary to
    the end are divided up into. Has to be greater than zero.
number_of_reblockings : int
    Number of reblocking analyses done in steps set by the width of an interval
    before it is checked whether suitable minimum error in the error has been
    found. Has to be greater than zero.
number_of_reblocks_to_cut_off : integer
    Number of reblocking analysis blocks to cut off additionally to the data
    before the best iteration with the lowest error in the error. Has to be non
    negative. It is highly recommended to not set this to zero.
pos_min_frac : float
    The minimum has to be in the first pos_min_frac part of the tested data to
    be taken as the true minimum. Has be to greater than a small number (here
    0.00001) and can at most be equal to one.
verbose : bool
    Determines whether extra information is printed out.
show_graph : bool
    Determines whether a window showing the shift vs iteration graph pops up
    highlighting where the minimum was found and - after also excluding some
    reblocking blocks - which iteration was found as the best starting iteration
    to use in reblocking analyses.

Returns
-------
starting_iteration: integer
    Iteration from which to start reblocking analysis for this calculation.
'''

    if frac_screen_interval <= 0:
        raise RuntimeError("frac_screen_interval <= 0")

    if number_of_reblocks_to_cut_off < 0:
        raise RuntimeError("number_of_reblocks_to_cut_off < 0")

    if pos_min_frac < 0.00001 or pos_min_frac > 1.0:
        raise RuntimeError("0.00001 < pos_min_frac < 1 not satisfied")

    if number_of_reblockings <= 0:
        raise RuntimeError("number_of_reblockings <= 0")

    if number_of_reblockings > frac_screen_interval:
        raise RuntimeError("number_of_reblockings > frac_screen_interval")

    # Find the point the shift began to vary.
    variable_shift = data['Shift'] != data['Shift'].iloc[0]
    if variable_shift.any():
        shift_variation_indx = data[variable_shift]['iterations'].index[0]
    else:
        raise RuntimeError("Shift has not started to vary in dataset!")

    # Check we have enough data to screen:
    if data['Shift'].size - shift_variation_indx < frac_screen_interval:
        # i.e. data where shift is not equal to initial value is less than
        # frac_screen_interval, i.e. we cannot screen adequately.
        warnings.warn("Calculation contains less data than "
            "frac_screen_interval. Will continue but frac_screen_interval is "
            "less than one data point.")

    # Find the MC iteration at which shift starts to vary.
    iteration_shift_variation_start = \
            data['iterations'].iloc[shift_variation_indx]

    step = int((data['iterations'].iloc[-1] - \
            iteration_shift_variation_start )/frac_screen_interval)

    min_index = -1
    min_error_error = float('inf')
    starting_iteration_found = False

    for k in range(int(frac_screen_interval/number_of_reblockings)):

        for j in range(k*number_of_reblockings, (k+1)*number_of_reblockings):
            start = iteration_shift_variation_start + j*step
            info = lazy_block(data, md, start, extract_psips=True)

            if info.no_opt_block:
                # Not enough data to get a best estimate for some values.
                # Don't include, even if the shift is estimated.
                s_err_err = float('inf')
            else:
                s_err_err = info.opt_block["standard error error"]["Shift"]
                if s_err_err <= min_error_error:
                    min_index = j
                    min_error_error = s_err_err
                    opt_ind = pyblock.pd_utils.optimal_block(info.reblock['Shift'])

            if verbose:
                print("Blocking attempt: %i. Error in the shift error: %f"
                        % (j, s_err_err))

        if min_index == -1:
            raise ValueError("Failed to find minimum error in the shift error!")

        if int(min_index * pos_min_frac) < j:
            # Also discard the frst n=number_of_reblocks_to_cut_off of data to
            # be conservative.  This amounts to removing n autocorrelation
            # lengths.
            discard_indx = 2**opt_ind * number_of_reblocks_to_cut_off
            starting_iteration = data['iterations'].iloc[shift_variation_indx+discard_indx] + min_index*step
            if data['iterations'].iloc[-1] <= starting_iteration:
                raise ValueError("Too much cut off! Data is not converged or "
                                 "use a smaller number_of_reblocks_to_cut_off.")
            starting_iteration_found = True
            break

    if not starting_iteration_found:
        raise RuntimeError("Failed to find starting iteration. The calculation "
                           "might not be converged.")

    if show_graph:
        plt.xlabel('Iterations')
        plt.ylabel('Shift')
        plt.plot(data['iterations'], data['Shift'], 'b-', label='data')
        plt.axvline(starting_iteration, color='r',
                    label='Suggested starting iteration')
        plt.axvline(iteration_shift_variation_start + min_index*step, color='g',
                    label='Starting iteration without discarding reblocks')
        plt.legend(loc='best')
        plt.show()

    return starting_iteration

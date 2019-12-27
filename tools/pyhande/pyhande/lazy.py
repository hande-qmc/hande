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
import numpy
import statsmodels.tsa.ar_model as ar_model
import statsmodels.tsa.stattools as tsastats


def find_starting_iteration_mser_min(data, md, start_max_frac=0.9, n_blocks=100, verbose=None, end=None):
    '''Find the best iteration to start analysing CCMC/FCIQMC data based on MSER minimization scheme.

.. warning::
    
    Use with caution, check whether output is sensible and adjust parameters 
    if necessary.

This function gives an optimal estimation of the starting interations 
based on MSER minimization heuristics. 
This methods decides the starting iterations :math:`d` as minimizing an evalualtion 
function 
MSER(:math:`d`) = :math:`\Sigma_{i=1}^{n-d} ( X_{i+d} - X_{mean}(d) ) / (n-d)^2`.
Here, :math:`n` is length of time-series, :math:`X_i` is '\sum H_0j N_j' / 'N_0' of :math:`i`-th step, and
:math:`X_{mean}` is the average of :math:`X_i` after the :math:`d`-th step.

Parameters 
----------
data : :class:`pandas.DataFrame`
    Calculation output for a FCIQMC or CCMC calculation.
md : dict
    Metadata corresponding to the calculation in `data`.
n_blocks : int
    This analysis takes long time when :math:`n` is large.
    Thus, we pick up :math:`d` for every 'n_blocks'
    samples, calculate MSER(:math:`d`), and decide the 
    optimal estimation of the starting iterations only 
    from these `d`.
start_max_frac : float
    MSER(d) may oscillate when become unreanably small 
    when :math:`n-d` is large. Thus, we calculate MSER(:math:`d`) 
    for :math:`d` < (:math:`n` * start_max_frac) and 
    give the optimal estimation of the starting iterations
    only in this range of :math:`d`.
verbose : int
    Inactive. This valuable does not change anything.
end : int or None
    Last iteration included in analysis. If None, the last iteration included
    is the last iteration of the data set.

Returns
-------
starting_iteration: integer
    Iteration from which to start reblocking analysis for this calculation.
'''

    if end is None:
        end = data['iterations'].iloc[-1]
    before_end_indx = data['iterations'] <= end
    data_before_end = data.ix[before_end_indx]
    
    list = data_before_end['Proj. Energy']
    n_data=len(list)

    mser_min = sys.float_info.max
    for i in range(n_blocks): 
        start_line = int(i*(n_data*start_max_frac)/n_blocks)
        mser = numpy.var(list[start_line:n_data]) / (n_data-start_line)
        if ( mser < mser_min):
            mser_min = mser
            starting_iteration = start_line * md['qmc']['ncycles']
            final_start_line = start_line
        
    if ( final_start_line > n_data*start_max_frac*0.8):
        warnings.warn('Proj. energy may not be converged. MSER min. may underestimate the starting iteration. One should check 1:$3/$4 plot.')

    return starting_iteration

def lazy_hybrid(calc, md, start=0, end=None, batch_size=1):
    '''New post-analysis on zero-temperature QMC calcaulations.

.. note::

    :func:`std_analysis` is recommended unless custom processing is required
    before blocking analysis is performed.

This scheme is made by hybridizing two different post-analysis methods,
AR model and Straatsma. The former (the latter) is comparatively 
good at estimating the statistic error for smaller (larger) length
of time-series, respectively. This method just picks up the larger
statistic error from the ones given by both methods. The mathematical 
details of both methods are explained in an upcoming paper.


Parameters
----------
calc : :class:`pandas.DataFrame`
    Zero-temperature QMC calculation output.
md : dict
    Metadata for the calculation in `calc`.
start, end : 
    See :func:`std_analysis`.
batch_size : int
    The energy time-series is coarse-grained by 
    averaging several sequential samples 
    into just one sample and the statistic error
    is calculated for the coarse-grained time-series.
    This variable designates how many sequential 
    samples are averaged together.

Returns
--------
info : :func:`collections.namedtuple`
    See :func:`std_analysis`.
'''

    if end is None:
        end = calc['iterations'].iloc[-1]
    before_end_indx = calc['iterations'] <= end
    data_before_end = calc.ix[before_end_indx]
    
    after_start_indx = data_before_end['iterations'] >= start
    calc_tr          = data_before_end[after_start_indx]
        
    #list = calc_tr['Proj. Energy'].as_matrix()
    #n_data = len(list)
    
    list_org = calc_tr['Proj. Energy'].as_matrix()
    n_data = len(list_org) / batch_size
    list = [0]*n_data
    for i in range(n_data):
        list[i] = numpy.mean(list_org[i*batch_size:(i+1)*batch_size])
    #print list

    #print list    
    mean = numpy.mean(list)
    var  = numpy.var(list)
    acf = tsastats.acf(x=list, unbiased=True, nlags=n_data-1, fft=True)
        
    # ar model
    ar = ar_model.AR(list)
    model_ar = ar.fit(ic='aic', trend='c', method='cmle')
    params = model_ar.params
    denom = nom = 1
    for j in range( len(params)-1 ):
        denom -= params[j+1]
        nom -= params[j+1] * acf[j+1]
    tau =  nom / denom**2
    error_ar = numpy.sqrt(var/n_data*tau)

    # autocorr
    tau = 1.0    
    for i in range(1, n_data-1):
        if(acf[i]<0):
            break
        tau += 2.0*acf[i]
    error_ac = numpy.sqrt(var/n_data*tau)    

    # return value
    error = max(error_ar, error_ac)
    opt_block = pd.DataFrame(
        {'mean': mean,
         'standard error': error,
         'standard error error': None,
         'estimate': pyblock.error.pretty_fmt_err(mean, error)},
        columns=['mean', 'standard error', 'standard error error', 'estimate'],
        index=['Proj. Energy'])
    no_opt_block = ['N_0','Shift','# H psips','\sum H_0j N_j']
    tuple_fields = ('metadata data data_len reblock covariance opt_block '
                    'no_opt_block'.split())
    info_tuple = collections.namedtuple('HandeInfo', tuple_fields)
    info = info_tuple(md, calc, None, None, None, opt_block, no_opt_block)
    return info


def std_analysis(datafiles, start=None, end=None, select_function=None,
        extract_psips=False, reweight_history=0, mean_shift=0.0,
        arith_mean=False, calc_inefficiency=False, verbosity = 1, 
        starts_reweighting=None, extract_rep_loop_time=False,
        analysis_method=None, warmup_detection=None):
    '''Perform a 'standard' analysis of HANDE output files.

Parameters
----------
datafiles : list of strings
    names of files containing HANDE QMC calculation output.
start, end : int or None
    iteration after which/until which the blocking analysis is performed. The
    end iteration is included in analysis, the start iteration is not.
    If start is None, then attempt to automatically determine a good iteration
    using :func:`find_starting_iteration`.  If end is None, the last iteration
    included is the last iteration of the data set.
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
    prevent the weights from becoming to large.
arith_mean : bool
calc_inefficiency : bool
    determines whether inefficiency should be calculated.
verbosity : int
    values greater than 1 print out blocking information when automatically
    finding the starting iteration. 0 and 1 print out the starting iteration if
    automatically found. Negative values print out nothing from the automatic
    starting point search.
starts_reweighting : list of floats
    used by the reweighting_graph function to pass more than one starting 
    iteration
extract_rep_loop_time : bool
    also extract the mean time taken per report loop from the calculation.
analysis_method : string
    determines which post-analysis method is used to estimate the statistic
    error. Currently 'reblocking' and 'hybrid' are prepared.
warmup_detection : string
    determines which method is used to decide the starting iterations 
    to be discarded before calculation the statistic error. Currently
    'hande_org' and 'mser_min' are prepared.

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
    infos = []
    for (calc, md) in zip(calcs, calcs_md):        
        calc_start = start
        calc_end = end
        if calc_start is None:
            if starts_reweighting is None:
                
                if (warmup_detection == 'hande_org'): # added_by_ichibha
                    calc_start = find_starting_iteration(calc, md, verbose=verbosity,
                                                         end=calc_end)
                elif (warmup_detection == 'mser_min'): # added_by_ichibha
                    calc_start = find_starting_iteration_mser_min(calc, md, verbose=verbosity,
                                                                  end=calc_end) # added_by_ichibha                             
            else:
                calc_start = starts_reweighting[len(infos)]
        md['pyhande'] = {'reblock_start': calc_start}
        if (verbosity > -1) :
            print('Block from: %i' % calc_start)
            
        if (analysis_method == 'reblocking'):  # added_by_ichibha
            infos.append(lazy_block(calc, md, calc_start, calc_end,          
                                    select_function, extract_psips, calc_inefficiency,
                                    extract_rep_loop_time))
        elif (analysis_method == 'hybrid'):  # added_by_ichibha
            infos.append(lazy_hybrid(calc, md, calc_start, calc_end))      # added_by_ichibha
    return infos

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
    for (md, df) in filter_calcs(hande_out, ('FCIQMC', 'CCMC', 'Simple FCIQMC','UCCMC')):
        if reweight_history > 0:
            df = pyhande.weight.reweight(df, md['qmc']['ncycles'],
                md['qmc']['tau'], reweight_history, mean_shift,
                arith_mean=arith_mean)
            df['W * \sum H_0j N_j'] = df['\sum H_0j N_j'] * df['Weight']
            df['W * N_0'] = df['N_0'] * df['Weight']
        df['Proj. Energy'] = df['\sum H_0j N_j'] / df['N_0'] 
        data.append(df)
        metadata.append(md)
    if data:
        calcs_metadata, calcs = concat_calcs(metadata, data)
    else:
        raise ValueError('No data found in '+' '.join(datafiles))
    return (calcs, calcs_metadata)

def lazy_block(calc, md, start=0, end=None, select_function=None,
            extract_psips=False, calc_inefficiency=False,
            extract_rep_loop_time=False):
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
start, end, select_function, extract_psips, calc_inefficiency,
    extract_rep_loop_time:
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
    if end is None:
        # Default end is the last iteration.
        end = calc['iterations'].iloc[-1]
    if select_function is None:
        # start+1 due to backwards compatibility
        # (used to be indx = calc['iterations'] > start
        # but pd.Series.between(start,end) gives True at
        # calc['iterations'] == start).
        indx = calc['iterations'].between(start+1, end)
    else:
        indx = select_function(calc)
    to_block = []
    if extract_psips:
        to_block.append('# H psips')
    to_block.extend(['\sum H_0j N_j', 'N_0', 'Shift'])
    if reweight_calc:
        to_block.extend(['W * \sum H_0j N_j', 'W * N_0'])
    if extract_rep_loop_time:
        to_block.append('time')
    if 'uccmc' in md.keys():
        to_block.append('N_0 UCCMC')
    mc_data = calc.ix[indx, to_block]
    if mc_data['Shift'].iloc[0] == mc_data['Shift'].iloc[1]:
        if calc['Shift'][~indx].iloc[-1] == mc_data['Shift'].iloc[0]:
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

    calc_filter = lambda md: md['calc_type'] in calc_types
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
        data = calcs
        metadata = calcs_metadata

    # Don't have UUID information in all calculations.
    # Assume any restarted calculations/set of calculations if sorted by uuids
    # above are in the right order from here and contiguous.
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
            # Different (set of) calculation(s)
            step = data[i]['iterations'].iloc[-1] - data[i]['iterations'].iloc[-2]
            calc_type = metadata[i]['calc_type']
            calcs.append(pd.concat(xcalc))
            xcalc = [data[i]]
            calcs_metadata.append(metadata[i])
        else:
            # Continuation of same (set of) calculation(s) (probably)
            xcalc.append(data[i])
        prev_iteration = data[i]['iterations'].iloc[-1]
    calcs.append(pd.concat(xcalc, ignore_index=True))
    calcs = [ca.drop_duplicates(subset='iterations', keep='last').reset_index(drop=True) for ca in calcs]
    return calcs_metadata, calcs

def find_starting_iteration(data, md, frac_screen_interval=300,
    number_of_reblockings=30, number_of_reblocks_to_cut_off=1, pos_min_frac=0.8,
    verbose=0, show_graph=False, end=None):
    '''Find the best iteration to start analysing CCMC/FCIQMC data.

.. warning::

    Use with caution, check whether output is sensible and adjust parameters if
    necessary.

First, consider only data from when the shift begins to vary. We are interested
in finding the minimum in the fractional error in the error of the shift
weighted by 1/sqrt(number of data points left). The error in the error of the
shift and the error in the shift vary as 1/sqrt(number of data points to
analyse) with the number of data points to analyse. If we were looking for the
minimum in either of these quantities, the minimum would therefore be biased to
the lower iterations as then more data points are included in the analysis.
However, we have noticed that the error in the shift and its error fluctuate as
we have less iterations to analyse which means that our search for the minimum
could get trapped easily in a local minimum. We therefore consider their
fraction. As they are divided by each other in the fractional error, the
1/sqrt(number of data points to analyse) gets removed. It is therefore
artificially included as a weight. To be more conservative, we also find the
minimum in the weighted fractional error in the error of # H psips, N_0,
\sum H_0j N_j. We then consider the minimum out of these four minima which is
at the highest number of iterations.

The best estimate of the iteration to start the blocking analysis is found by:

1. discard data during the constant shift phase.
2. estimate the weighted fractional error in the error of the shift, # H psips,
   N_0, \sum H_0j N_j, by blocking the remaining data :math:`n` times, where
   the blocking analysis considers the last :math:`1-i/f` fraction of the data
   and where :math:`i` is the number of blocking analyses already performed,
   :math:`n`  is `number_of_reblockings`  and :math:`f` is
   `frac_screen_interval`.
3. find the iteration which gives the minimum estimate of the weighted
   fractional error in the error of the shift, numerator of projected energy,
   reference and total population. We then focus on the minimum out of these
   four minima which is at the highest number of iterations. If this is in the
   first `pos_min_frac` fraction of the blocking attempts, go to 4, otherwise
   repeat 2 and perform an additional `number_of_reblockings` attempts.
4. To be conservative, discard the first `number_of_reblocks_to_cut_off` blocks
   from the start iteration, where each block corresponds to roughly the
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
verbose : int
    If greater than 1, prints out which blocking attempt is currently being
    performed.
show_graph : bool
    Determines whether a window showing the shift vs iteration graph pops up
    highlighting where the minimum was found and - after also excluding some
    reblocking blocks - which iteration was found as the best starting iteration
    to use in reblocking analyses.
end : int or None
    Last iteration included in analysis. If None, the last iteration included
    is the last iteration of the data set.

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

    if end is None:
        # Default end is the last iteration.
        end = data['iterations'].iloc[-1]
    before_end_indx = data['iterations'] <= end
    data_before_end = data.ix[before_end_indx]

    # Find the point the shift began to vary.
    variable_shift = data_before_end['Shift'] != data_before_end['Shift'].iloc[0]
    if variable_shift.any():
        shift_variation_indx = data_before_end[variable_shift]['iterations'].index[0]
    else:
        raise RuntimeError("Shift has not started to vary in dataset!")

    # Check we have enough data to screen:
    if data_before_end['Shift'].size - shift_variation_indx < frac_screen_interval:
        # i.e. data where shift is not equal to initial value is less than
        # frac_screen_interval, i.e. we cannot screen adequately.
        warnings.warn("Calculation contains less data than "
            "frac_screen_interval. Will continue but frac_screen_interval is "
            "less than one data point.")

    # Find the MC iteration at which shift starts to vary.
    iteration_shift_variation_start = \
            data_before_end['iterations'].iloc[shift_variation_indx]

    step = int((data_before_end['iterations'].iloc[-1] - \
            iteration_shift_variation_start )/frac_screen_interval)

    step_indx = int((data_before_end['iterations'].index[-1]-\
            shift_variation_indx)/frac_screen_interval)

    min_index = -1
    err_keys = ['Shift',  'N_0', '\sum H_0j N_j', '# H psips']
    min_error_frac_weighted = pd.Series([float('inf')]*len(err_keys), index=err_keys)
    starting_iteration_found = False

    for k in range(int(frac_screen_interval/number_of_reblockings)):

        for j in range(k*number_of_reblockings, (k+1)*number_of_reblockings):
            start = iteration_shift_variation_start + j*step
            info = lazy_block(data_before_end, md, start, end=end, extract_psips=True)

            if info.no_opt_block:
                # Not enough data to get a best estimate for some values.
                # Don't include, even if the shift is estimated.
                s_err_frac_weighted = float('inf')
            else:
                number_of_data_left = data_before_end['Shift'].index[-1] - shift_variation_indx - j*step_indx + 1
                err_err = info.opt_block.loc[err_keys, 'standard error error']
                err = info.opt_block.loc[err_keys, 'standard error']
                err_frac = err_err.divide(err)
                err_frac_weighted = err_frac.divide(math.sqrt(float(number_of_data_left)))
                s_err_frac_weighted = err_frac_weighted['Shift']
                if (err_frac_weighted <= min_error_frac_weighted).any():
                    min_index = j
                    min_error_frac_weighted = err_frac_weighted.copy()
                    opt_ind = pyblock.pd_utils.optimal_block(info.reblock[err_keys])

            if (verbose > 1):
                print("Blocking attempt: %i. Blocking from: %i. "
                      "Upper bound on shift fractional weighted error: %f" % (j, start, s_err_frac_weighted))

        if min_index < int(pos_min_frac * j) and min_index > -1:
            # Also discard the frst n=number_of_reblocks_to_cut_off of data to
            # be conservative.  This amounts to removing n autocorrelation
            # lengths.
            discard_indx = 2**opt_ind * number_of_reblocks_to_cut_off
            starting_iteration = data_before_end['iterations'].iloc[shift_variation_indx+discard_indx] + min_index*step
            if data_before_end['iterations'].iloc[-1] <= starting_iteration:
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
        plt.plot(data_before_end['iterations'], data_before_end['Shift'], 'b-', label='data')
        plt.axvline(starting_iteration, color='r',
                    label='Suggested starting iteration')
        plt.axvline(iteration_shift_variation_start + min_index*step, color='g',
                    label='Starting iteration without discarding reblocks')
        plt.legend(loc='best')
        plt.show()

    return starting_iteration

def reweighting_graph(datafiles, start=None, verbosity=1, mean_shift=0.0,
                      arith_mean=False):
    '''Plot a graph of reweighted projected energy vs. reweighted factor W.
    
Detecting biases by reweighting is described in [Umrigar93]_ and [Vigor15]_ , 
see pyhande.weight for details. The graph produced by this function is similar
to figure 4 in [Vigor15]_.

A similar function has been published in
Neufeld, V., & Thom, A. J. Research data and further information supporting
"A study of the dense uniform electron gas with high orders of coupled cluster" [Dataset].
https://doi.org/10.17863/CAM.14336 under Attribution 4.0 International (CC BY 4.0).
   
Parameters
----------
datafiles : list of strings
    names of files containing HANDE QMC calculation output.
start : int or None
    iteration from which the blocking analysis is performed.  If None, then
    attempt to automatically determine a good iteration using
    :func:`find_starting_iteration`.
verbosity : int
    values greater than 1 print out blocking information when automatically
    finding the starting iteration. 0 and 1 print out the starting iteration if
    automatically found. Negative values print out nothing from the automatic 
    starting point search.
mean_shift : float
    prevent the weights from becoming to large.
arith_mean : bool
        

References
----------
Umrigar93
    C.J. Umirigar et al., J. Chem. Phys. 99, 2865 (1993)
Vigor15
    W.A. Vigor, et al., J. Chem. Phys. 142, 104101 (2015).

Thanks to Will Vigor for original implementation.
'''
    weights = [2,4,8,16,32,64,128,256]
    infos_collection = []
    infos = std_analysis(datafiles=datafiles, start=start, extract_psips=True, 
            verbosity = verbosity)
    infos_collection.append(infos)
    ndiff_calcs = len(infos)
    starts = [info.metadata['pyhande']['reblock_start'] for info in infos]
    
    for weight in weights:
        infos = std_analysis(datafiles=datafiles, start=start, 
                extract_psips=True, reweight_history=weight, 
                mean_shift=mean_shift, arith_mean=arith_mean, 
                verbosity=verbosity, starts_reweighting=starts)
        infos_collection.append(infos)
    
    for k in range(0, ndiff_calcs):
        reweights = [0]
        weighted_proj_energies = []
        weighted_proj_energy_errs = []
        weighted_proj_energies.append(infos_collection[0][k].opt_block['mean']['Proj. Energy'])
        weighted_proj_energy_errs.append(infos_collection[0][k].opt_block['standard error']['Proj. Energy'])
        for j in range(0,len(weights)):
            try:
                weighted_proj_energies.append(infos_collection[j+1][k].opt_block['mean']['Weighted Proj. E.'])
                weighted_proj_energy_errs.append(infos_collection[j+1][k].opt_block['standard error']['Weighted Proj. E.'])
                reweights.append(weights[j])
            except:
                print("Reweighting with " + str(weights[j]) + " and further failed.") 
                break
        
        # Plotting the graph
        plt.xlabel('Reweights')
        plt.ylabel('(Reweighted) Projected Energy')
        plt.errorbar(reweights, weighted_proj_energies, yerr=weighted_proj_energy_errs, fmt='o')
        plt.title('Calculation #' + str(k+1))
        plt.tight_layout()
        plt.show()

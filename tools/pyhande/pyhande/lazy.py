'''Tools for the lazy amongst us: automation of common HANDE analysis tasks.'''

import collections
import pandas as pd
import warnings

import pyblock
import pyhande.extract
import pyhande.analysis
import pyhande.weight
import math

def std_analysis(datafiles, start=0, select_function=None, extract_psips=False,
                reweight_history=0, mean_shift=0.0, arith_mean=False):
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

    hande_out = pyhande.extract.extract_data_sets(datafiles)

    # Concat all QMC data (We did say 'lazy', so assumptions are being made...)
    data = []
    metadata = []
    for (md, df) in hande_out:
        if any(calc in md['calc_type'] for calc in ('FCIQMC', 'CCMC')):
            if reweight_history > 0:
                if 'mc_cycles' in md:
                    df = pyhande.weight.reweight(df, md['mc_cycles'],
                        md['tau'], reweight_history, mean_shift,
                        arith_mean=arith_mean)
                else:
                    df = pyhande.weight.reweight(df, md['qmc']['ncycles'],
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

    tuple_fields = ('metadata data data_len reblock covariance opt_block '
                   'no_opt_block'.split())
    info_tuple = collections.namedtuple('HandeInfo', tuple_fields)
    return_vals = []
    for calc,md in zip(calcs,calcs_metadata):
        # Reblock Monte Carlo data over desired window.
        if select_function is None:
            indx = calc['iterations'] > start
        else:
            indx = select_function(calc)
        to_block = []
        if extract_psips:
            to_block.append('# H psips')
        to_block.extend(['\sum H_0j N_j', 'N_0', 'Shift'])
        if reweight_history > 0:
            to_block.extend(['W * \sum H_0j N_j', 'W * N_0'])

        mc_data = calc.ix[indx, to_block]

        if mc_data['Shift'].iloc[0] == mc_data['Shift'].iloc[1]:
            warnings.warn('The blocking analysis starts from before the shift '
                         'begins to vary.')

        (data_len, reblock, covariance) = pyblock.pd_utils.reblock(mc_data)
        
        proje = pyhande.analysis.projected_energy(reblock, covariance, data_len)
        reblock = pd.concat([reblock, proje], axis=1)
        to_block.append('Proj. Energy')

        if reweight_history > 0:
            proje = pyhande.analysis.projected_energy(reblock, covariance, 
                        data_len, sum_key='W * \sum H_0j N_j', ref_key='W * N_0',
                        col_name='Weighted Proj. E.')
            reblock = pd.concat([reblock, proje], axis=1)
            to_block.append('Weighted Proj. E.')
        
        # Summary (including pretty printing of estimates).
        (opt_block, no_opt_block) = pyhande.analysis.qmc_summary(reblock, to_block)

        # Now try to calculate an inefficiency according to
        # Understanding and improving the efficiency of full configuration interaction quantum Monte Carlo
        # W. A. Vigor, J. S. Spencer, M. J. Bearpark, and A. J. W. Thom
        # J. Chem. Phys. 144, 094110 (2016); doi: 10.1063/1.4943113 

        #  We also try to provide some sort of error estimate for the inefficiency on the basis of the relative standard error error of the projected energy numerator..
        try:
            try:
                dtau = md['tau']
            except:
                dtau = md['qmc']['tau']
            sigmaE = opt_block['standard error']['Proj. Energy']
            Np = opt_block['mean']['# H psips']
            N = data_len[int(pyblock.pd_utils.reblock_summary(reblock.ix[:, "Proj. Energy"]).index)]
            inefficiency = sigmaE * math.sqrt(Np*N*dtau)
            stderr = inefficiency*(0.5*opt_block['standard error']['# H psips']/Np+(opt_block['standard error error']['\sum H_0j N_j']/opt_block['standard error']['\sum H_0j N_j']))
            d = pd.DataFrame(data={'mean':inefficiency, 'standard error':stderr}, index = ['Inefficiency'])
            opt_block = opt_block.append(d)
        except:
            pass

        estimates = []
        for (name, row) in opt_block.iterrows():
            estimates.append(
                    pyblock.error.pretty_fmt_err(row['mean'], row['standard error'])
                           )
        opt_block['estimate'] = estimates

        return_vals.append(info_tuple(md, calc, data_len, reblock,
                                      covariance, opt_block, no_opt_block))

    return return_vals

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

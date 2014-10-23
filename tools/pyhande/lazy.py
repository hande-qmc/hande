'''Tools for the lazy amongst us: automate common analysis tasks.'''

import collections
import pandas as pd
import warnings

import pyblock
import pyhande.extract
import pyhande.analysis
import pyhande.weight

def std_analysis(datafiles, start=0, select_function=None, extract_psips=False,
                reweight_itts=0, mean_shift=0.0, arith_mean=False):
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
reweight_itts : integer
    reweight in an attempt to remove population control bias. According to
    C. J. Umirigar et. al. J. Chem. Phys. 99, 2865 (1993) this should be set
    to be a few correlation times.
mean_shift: float
    prevent the weights from beoming to large.
Returns
-------
info : :func:`collections.namedtuple`
    raw and analysed data, consisting of:

        metadata, data
            from :func:`pyhande.extract.extract_data_sets`.
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
'''

    (metadata, data) = pyhande.extract.extract_data_sets(datafiles)

    # Reblock Monte Carlo data over desired window.
    if select_function is None:
        indx = data['iterations'] > start
    else:
        indx = select_function(data)
    to_block = ['Shift', '\sum H_0j N_j', 'N_0']
    if extract_psips:
        to_block.append('# H psips')

    # Compute and define new weighted columns to reblock.
    if reweight_itts > 0:
        data = pyhande.weight.reweight(data, metadata[0]['mc_cycles'],
            metadata[0]['tau'], reweight_itts, mean_shift, 
            arith_mean=arith_mean)
        data['W * \sum H_0j N_j'] = data['\sum H_0j N_j'] * data['Weight']
        data['W * N_0'] = data['N_0'] * data['Weight']
        to_block.append('W * \sum H_0j N_j')
        to_block.append('W * N_0')

    mc_data = data.ix[indx, to_block]

#    if mc_data['Shift'][1] == mc_data['Shift'][2]:
#        warnings.warn('The blocking analysis starts from before the shift '
#                     'begins to vary.')

    (data_len, reblock, covariance) = pyblock.pd_utils.reblock(mc_data)
    
    proje = pyhande.analysis.projected_energy(reblock, covariance, data_len)
    reblock = pd.concat([reblock, proje], axis=1)

    if reweight_itts > 0:
        proje = pyhande.analysis.projected_energy(reblock, covariance, 
                    data_len, sum_key='W * \sum H_0j N_j', ref_key='W * N_0'
                    ,col_name='Weighted Proj. E.')
        reblock = pd.concat([reblock, proje], axis=1)

    # Summary (including pretty printing of estimates).
    (opt_block, no_opt_block) = pyhande.analysis.qmc_summary(reblock)
    if extract_psips:
        (opt_block, no_opt_block) = pyhande.analysis.qmc_summary(reblock,
                keys=('# H psips',), summary_tuple=(opt_block, no_opt_block))
    if reweight_itts > 0:
        (opt_block, no_opt_block) = pyhande.analysis.qmc_summary(reblock,
                keys=('W * \sum H_0j N_j', 'W * N_0', 'Weighted Proj. E.'),
                summary_tuple=(opt_block, no_opt_block))

    estimates = []
    for (name, row) in opt_block.iterrows():
        estimates.append(
                pyblock.error.pretty_fmt_err(row['mean'], row['standard error'])
                       )
    opt_block['estimate'] = estimates

    tuple_fields = ('metadata data data_len reblock covariance opt_block '
                   'no_opt_block'.split())
    info_tuple = collections.namedtuple('HandeInfo', tuple_fields)

    return info_tuple(metadata, data, data_len, reblock, covariance, opt_block,
                      no_opt_block)

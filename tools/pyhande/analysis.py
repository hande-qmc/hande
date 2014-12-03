'''Analysis helper routines for HANDE data.'''

from os import path
import pandas as pd
import sys

sys.path.append(path.join(path.abspath(path.dirname(__file__)), '../pyblock'))
import pyblock.error
import pyblock.pd_utils

def projected_energy(reblock_data, covariance, data_length,
                     sum_key='\sum H_0j N_j', ref_key='N_0'):
    '''Calculate the projected energy estimator and associated error.

The projected energy estimator is given by

.. math::

    E = \\frac{\\sum_H_0j N_j}{N_0}

The numerator and denominator are correlated and so their covariance must be
taken into account.

Parameters
----------
reblock_data : :class:`pandas.DataFrame`
    reblock data for (at least) the numerator and denominator in the
    projected energy estimator.
covariance: :class:`pandas.DataFrame`
    covariance at each reblock iteration between (at least) the numerator
    and denominator in the projected energy estimator.
data_length: :class:`pandas.DataFrame`
    number of data points in each reblock iteration.
sum_key : string
    column name in reblock_data containing :math:`\\sum_H_0j N_j``, i.e. the sum
    of the population weighted by the Hamiltonian matrix element with the trial
    wavefunction.
ref_key : string
    column name in reblock_data containing :math:`N_0``, i.e. the population of
    the trial wavefunction (often/originally just a single determinant).

Returns
-------
proje : :class:`pandas.DataFrame`
    The projected energy estimator at each reblock iteration.

See also
--------
:func:`pyblock.pd_utils.reblock` for producing the input parameters.
'''

    proje_sum = reblock_data.ix[:, sum_key]
    ref_pop = reblock_data.ix[:, ref_key]
    proje_ref_cov = covariance.xs(ref_key, level=1)[sum_key]
    proje = pyblock.error.ratio(proje_sum, ref_pop, proje_ref_cov, data_length)
    proje.columns = pd.MultiIndex.from_tuples(
            [('Proj. Energy', col) for col in proje.columns])
    return proje

def qmc_summary(data, keys=('Shift', '\sum H_0j N_j', 'N_0',
                            'Proj. Energy'), summary_tuple=None):
    '''Summarise a reblocked data set by the optimal block.

Parameters
----------
data : :class:`pandas.DataFrame`
    reblocked data (i.e. data with the reblock iteration as the index).
keys : list of strings
    columns (by top-level index) of the data table to inspect.  Each top-level
    column must contain an optimal block column.
summary_tuple :  (:class:`pandas.DataFrame`, list of strings)
    Optionally append summary data to this tuple. Allows repeated calling of 
    this function.

Returns
-------
opt_data : :class:`pandas.DataFrame`
    Data for each column from the optimal block size of that column.
no_opt : list of strings
    list of columns for which no optimal block size was found.
'''

    if summary_tuple:
        (opt_data, no_opt) = ([summary_tuple[0]], summary_tuple[1])
    else:
        (opt_data, no_opt) = ([], [])
    for col in keys:
        if col in data:
            summary = pyblock.pd_utils.reblock_summary(data.ix[:, col])
            if summary.empty:
                no_opt.append(col)
            else:
                summary.index = [col]
            opt_data.append(summary)
    opt_data = pd.concat(opt_data)
    return (opt_data, no_opt)


# [review] - JSS: no need for \ to mark a continuation inside a pair of ().
# [review] - JSS: the plateau is a more general feature than the shoulder--bad name?
def shoulder_estimator(data, total_key='# H psips', ref_key='N_0',\
# [review] - JSS: indent this to be in line with the first argument.
 shift_key='Shift', min_pop=10):
    # [review] - JSS: note additional blank line for PEP-8 compliant docstring.
    # [review] - JSS: check rewording of the description.
    ''' Estimate the shoulder (plateau) from a (CCMC) FCIQMC calculation.

The population on the reference starts to grow exponentially during the plateau,
whilst the total population grows exponentially from the start of the
calculation before stabilising (perhaps only briefly) during the plateau phase.
As a result, the ratio of the total population to the population on the
reference is at a maximum at the start of the plateau.

The shoulder estimator is defined to be mean of the ten points with the smallest
proportion of the population on the reference (excluding points when the 
population drops below min_pop excips (psips). The shoulder height is the total
population at this point.

Credit to AJWT for original implementation.

Parameters
----------
data : :class:`pandas.DataFrame`
    HANDE QMC data. The function pyhande.extract.extract_data_sets can 
    be used to extract this from a HANDE output file.
total_key : string
    column name in reblock_data containing the total number of psips.
ref_key : string
    column name in reblock_data containing the number of psips on the
    reference determinant.
shift_key : string
    column name in reblock_data containing the shift.
min_pop : int
    exclude points with less than min_pop on the reference.

Returns
-------
plateau_data : :class:`pandas.DataFrame`
    An estimate of the shoulder (plateau) from a FCIQMC (CCMC) calculation,
    along with the associated standard error.
'''
    plateau_data = []
    shift_first = data[shift_key][0][0]
    # [review] - JSS: select pre-variable-shift data before doing the division.
    # [review] - JSS: what if the population on the reference is negative?
    data['Shoulder'] = data[total_key]/data[ref_key]

    # Select data which both the reference is greater than the min_pop and 
    # shift updating hasn't started. Annoyingly data[data[ref_key] < min_pop &&
    # data[shift_key] == shift_first] doesn't work as the && operation on a 
    # series of booleans is undefined.
    # [review] - JSS: use & rather than &&.
    # [review] - JSS: comments say the population must be < min_pop to be discarded but you discard it if it's <=.
    sorted_data =  data[data[ref_key] > min_pop]\
                   [data[shift_key][data[ref_key]>min_pop] == shift_first ]\
                   .sort('Shoulder')
    # [review] - JSS: no need for \.
    plateau_data.append([sorted_data[-10:]['Shoulder'].mean(),\
                       sorted_data[-10:]['Shoulder'].sem()])
    plateau_data.append([sorted_data[-10:][total_key].mean(),\
                       sorted_data[-10:][total_key].sem()])
    plateau_data = pd.DataFrame(data=plateau_data, \
                               columns=['mean', 'standard error'],\
                               index=['shoulder estimator', 'shoulder height'])
    return plateau_data

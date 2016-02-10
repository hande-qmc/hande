'''Analysis of data from FCIQMC and CCMC calculations.'''

from os import path
import numpy
import pandas as pd
import sys

try:
    import pyblock
except ImportError:
    sys.path.append(path.join(path.abspath(path.dirname(__file__)), '../../pyblock'))
    import pyblock

def projected_energy(reblock_data, covariance, data_length,
                     sum_key='\sum H_0j N_j', ref_key='N_0',
                     col_name='Proj. Energy'):
    '''Calculate the projected energy estimator and associated error.

The projected energy estimator is given by

.. math::

    E = \\frac{\\sum H_0j N_j}{N_0}

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
    column name in reblock_data containing :math:`\\sum H_0j N_j`, i.e. the sum
    of the population weighted by the Hamiltonian matrix element with the trial
    wavefunction.
ref_key : string
    column name in reblock_data containing :math:`N_0`, i.e. the population of
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
            [(col_name, col) for col in proje.columns])
    return proje

def qmc_summary(data, keys=('\sum H_0j N_j', 'N_0', 'Shift', 'Proj. Energy'),
                            summary_tuple=None):
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

def extract_pop_growth(data, ref_key='N_0', shift_key='Shift', min_ref_pop=10):
    '''Select QMC data during which the population was allowed to grow.

We define the region of population growth as the period in which the shift is
held constant.

Parameters
----------
data : :class:`pandas.DataFrame`
    HANDE QMC data. :func:`pyhande.extract.extract_data_sets` can be used to
    extract this from a HANDE output file.
ref_key : string
    column name in reblock_data containing the number of psips on the reference
    determinant.
shift_key : string
    column name in reblock_data containing the shift.
min_pop : int
    discard data entries with fewer than min_pop on the reference.

Returns
-------
pop_data : :class:`pandas.DataFrame`
    The subset of data prior to the shift being varied.
'''
    pop_data = data[data[shift_key] == data[shift_key].iloc[0]]
    # Now discard any data with a smaller absolute population on reference.
    # Don't combine this with the above as the total number of iterations might
    # be much larger than the number before the shift started to vary.
    pop_data = pop_data[pop_data[ref_key] >= min_ref_pop]

    return pop_data

def plateau_estimator(data, total_key='# H psips', ref_key='N_0',
                      shift_key='Shift', min_ref_pop=10, pop_data=None):
    ''' Estimate the (plateau) shoulder from a FCIQMC/CCMC calculation.

The population on the reference starts to grow exponentially during the plateau,
whilst the total population grows exponentially from the start of the
calculation before stabilising (perhaps only briefly) during the plateau phase.
As a result, the ratio of the total population to the population on the
reference is at a maximum at the start of the plateau.

The shoulder estimator is defined to be mean of the ten points with the smallest
proportion of the population on the reference (excluding points when the 
population drops below min_pop excips (psips). The shoulder height is the total
population at this point.

Credit to Alex Thom for original implementation.

Parameters
----------
data : :class:`pandas.DataFrame`
    HANDE QMC data. :func:`pyhande.extract.extract_data_sets` can be used to
    extract this from a HANDE output file.
total_key : string
    column name in reblock_data containing the total number of psips.
ref_key : string
    column name in reblock_data containing the number of psips on the
    reference determinant.
shift_key : string
    column name in reblock_data containing the shift.
min_ref_pop : int
    exclude points with less than min_ref_pop on the reference.
pop_data: :class:`pandas.DataFrame`
    The subset of data prior to the shift being varied.  Calculated if not
    supplied from extract_pop_growth.

Returns
-------
plateau_data : :class:`pandas.DataFrame`
    An estimate of the shoulder (plateau) from a FCIQMC (CCMC) calculation,
    along with the associated standard error.
'''
    plateau_data = []
    if pop_data is None:
        pop_data = extract_pop_growth(data, ref_key, shift_key, min_ref_pop)

    pop_data['Shoulder'] = pop_data[total_key]/abs(pop_data[ref_key])
    sorted_data = pop_data.sort('Shoulder')
    plateau_data.append([sorted_data[-10:]['Shoulder'].mean(),
                         sorted_data[-10:]['Shoulder'].sem()])
    plateau_data.append([sorted_data[-10:][total_key].mean(),
                         sorted_data[-10:][total_key].sem()])
    plateau_data = pd.DataFrame(data=plateau_data, 
                               columns=['mean', 'standard error'],\
                               index=['shoulder estimator', 'shoulder height'])
    return plateau_data

def plateau_estimator_hist(data, total_key='# H psips', shift_key='Shift',
                          pop_data=None, bin_width_fn=None):
    '''Estimate the plateau height via a histogram of the population.

The population (approximately) stabilises during the plateau phase.  By taking
a histogram of the population, the plateau can be estimated from the histogram
bin with greatest frequency.  Due to the exponential population growth outside
of the plateau, we histogram the logarithm of the population.

This tends to give similar numbers to shoulder_estimator, though may be less
useful for shoulder-like plateaus.  Detecting a plateau automatically is tricky
so having multiple approaches for comparison helps with corner cases.

Used in [Shepherd14]_.

Credit to James Shepherd for the idea and original (perl) implementation.

Parameters
----------
data : :class:`pandas.DataFrame`
    HANDE QMC data. func:`pyhande.extract.extract_data_sets` can be used to
    extract this from a HANDE output file.
total_key : string
    column name in reblock_data containing the total number of psips.
shift_key : string
    column name in reblock_data containing the shift.
pop_data: :class:`pandas.DataFrame`
    The subset of data prior to the shift being varied.  Calculated if not
    supplied from extract_pop_growth.
bin_width_fn : function
    A function which calculates the bin width in the histogram based upon
    pop_data.  12500/len(data)^2 (obtained empirically) is used if not supplied.

Returns
-------
plateau : float
    An estimate of the population at the plateau.

References
----------
Shepherd14
    J.J. Shepherd et al., Phys. Rev. B 90, 155130 (2014).
'''
    if pop_data is None:
        pop_data = extract_pop_growth(data, shift_key=shift_key, min_ref_pop=0)

    log_pop = numpy.log10(pop_data[total_key])

    if bin_width_fn:
        bin_width = bin_width_fn(pop_data)
    else:
        # bin size determined empirically for a range of 10^1 to 10^7.
        bin_width = 12500.0 / len(pop_data)**2
    nbins = int(max(log_pop) / bin_width)

    (hist, bin_edges) = numpy.histogram(log_pop, bins=nbins)
    hist_max = hist.argmax()
    if hist_max == nbins-1:
        # If the histogram peaks in the last bin, then most likely the
        # simulation has not reached a plateau.
        return numpy.nan
    else:
        return numpy.mean(10**bin_edges[hist_max:hist_max+2])

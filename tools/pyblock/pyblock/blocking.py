'''Tools for reblocking of data to remove serial correlation from data sets.'''

# copyright: (c) 2014 James Spencer
# license: modified BSD license; see LICENSE for further details.

import numpy
import collections

def reblock(data, rowvar=1, ddof=None, weights=None):
    '''Blocking analysis of correlated data.

Repeatedly average neighbouring data points in order to remove the effect of
serial correlation on the estimate of the standard error of a data set, as
described by Flyvbjerg and Petersen [Flyvbjerg]_.  The standard error is constant
(within error bars) once the correlation has been removed.

If a weighting is provided then the weighted variance and standard error of
each variable is calculated, as described in [Pozzi]_. Bessel correction is
obtained using the "effective sample size" from [Madansky]_.

.. default-role:: math

Parameters
----------
data : :class:`numpy.ndarray`
    1D or 2D array containing multiple variables and data points.  See ``rowvar``.
rowvar : int
    If ``rowvar`` is non-zero (default) then each row represents a variable and
    each column a data point per variable.  Otherwise the relationship is
    swapped.  Only used if data is a 2D array.
ddof : int
    If not ``None``, then the standard error and covariance are normalised by
    `(N - \\text{ddof})`, where `N` is the number of data points per variable.
    Otherwise, the numpy default is used (i.e. `(N - 1)`).
weights : :class:`numpy.array`
    A 1D weighting of the data to be reblocked. For multidimensional data an
    identical weighting is applied to the data for each variable.

Returns
-------
block_info : :class:`list` of :func:`collections.namedtuple`
    Statistics from each reblocking iteration.  Each tuple contains:

        block : int
            blocking iteration.  Each iteration successively averages neighbouring
            pairs of data points.  The final data point is discarded if the number
            of data points is odd.
        ndata: int
            number of data points in the blocking iteration.
        mean : :class:`numpy.ndarray`
            mean of each variable in the data set.
        cov : :class:`numpy.ndarray`
            covariance matrix.
        std_err : :class:`numpy.ndarray`
            standard error of each variable.
        std_err_err : :class:`numpy.ndarray`
            an estimate of the error in the standard error, assuming a Gaussian
            distribution.

References
----------
.. [Flyvbjerg]  "Error estimates on averages of correlated data", H. Flyvbjerg,
   H.G. Petersen, J. Chem. Phys. 91, 461 (1989).
.. [Pozzi]  "Exponential smoothing weighted correlations", F. Pozzi, T. Matteo,
   T. Aste, Eur. Phys. J. B. 85, 175 (2012).
.. [Madansky]  "An Analysis of WinCross, SPSS, and Mentor Procedures for
   Estimating the Variance of a Weighted Mean", A. Madansky, H. G. B. Alexander,
   www.analyticalgroup.com/download/weighted_variance.pdf
'''

    if ddof is not None and ddof != int(ddof):
        raise ValueError("ddof must be integer")
    if ddof is None:
        ddof = 1

    if data.ndim > 2:
        raise RuntimeError("do not understand how to reblock in more than two dimensions")

    if data.ndim == 1 or data.shape[0] == 1:
        rowvar = 1
        axis = 0
        nvar = 1
    elif rowvar:
        nvar = data.shape[0]
        axis = 1
    else:
        nvar = data.shape[1]
        axis = 0

    if weights is not None:
        if weights.ndim > 1:
            raise RuntimeError("cannot handle multidimensional weights")
        if weights.shape[0] != data.shape[axis]:
            raise RuntimeError("incompatible numbers of weights and samples")
        if numpy.any(weights < 0):
            raise RuntimeError("cannot handle negative weights")

    iblock = 0
    stats = []
    block_tuple_fields = 'block ndata mean cov std_err std_err_err'.split()
    block_tuple = collections.namedtuple('BlockTuple', block_tuple_fields)
    while data.shape[axis] >= 2:

        data_len = data.shape[axis]

        if weights is None:
            nsamp = data_len
            mean = numpy.array(numpy.mean(data, axis=axis))
            cov = numpy.cov(data, rowvar=rowvar, ddof=ddof)
        else:
            mean, tot_weight = numpy.average(data,
                                             axis=axis,
                                             weights=weights,
                                             returned=True)
            mean = numpy.array(mean)
            if nvar > 1:
                tot_weight = tot_weight[0]
            norm_wts = weights/tot_weight
            nsamp = 1.0/numpy.sum(norm_wts*norm_wts)
            bessel = nsamp/(nsamp - ddof)
            if data.ndim == 1:
                ds = data - mean
                cov = bessel*numpy.sum(norm_wts*ds*ds)
            else:
                ds = numpy.empty((nvar, data_len))
                for i in range(nvar):
                    d = data[i, :] if rowvar else data[:, i]
                    ds[i] = d - mean[i]
                cov = numpy.zeros((nvar, nvar))
                for i in range(nvar):
                    for j in range(i, nvar):
                        cov[i, j] = bessel*numpy.sum(norm_wts*ds[i]*ds[j])
            if nvar > 1:
                cov = cov + cov.T - numpy.diag(cov.diagonal())

        if cov.ndim < 2:
            std_err = numpy.array(numpy.sqrt(cov/nsamp))
        else:
            std_err = numpy.sqrt(cov.diagonal()/nsamp)
        std_err_err = numpy.array(std_err/(numpy.sqrt(2*(nsamp - ddof))))
        stats.append(
            block_tuple(iblock, data_len, mean, cov, std_err, std_err_err)
        )

        # last even-indexed value (ignore the odd one, if relevant)
        half = int(data.shape[axis]/2)
        last = 2*half
        if weights is None:
            if data.ndim == 1 or not rowvar:
                data = (data[:last:2] + data[1:last:2])/2
            else:
                data = (data[:,:last:2] + data[:,1:last:2])/2
        else:
            weights = norm_wts[:last:2] + norm_wts[1:last:2]
            if data.ndim == 1:
                wt_data = data[:last]*norm_wts[:last]
                data = (wt_data[::2] + wt_data[1::2])/weights
            elif rowvar:
                wt_data = data[:,:last]*norm_wts[:last]
                data = (wt_data[:,::2] + wt_data[:,1::2])/weights
            else:
                idxs = 'ij,i->ij'
                wt_data = numpy.einsum(idxs, data[:last], norm_wts[:last])
                summed_wt_data = wt_data[::2] + wt_data[1::2]
                data = numpy.einsum(idxs, summed_wt_data, 1.0/weights)

        iblock += 1

    return stats

def find_optimal_block(ndata, stats):
    '''Find the optimal block length from a reblocking calculation.

Inspect a reblocking calculation and find the block length which minimises the
stochastic error and removes the effect of correlation from the data set.  This
follows the procedures detailed by [Wolff]_ and [Lee]_ et al.

.. default-role:: math

Parameters
----------
ndata : int
    number of data points ('observations') in the data set.
stats : list of tuples
    statistics in the format as returned by :func:`pyblock.blocking.reblock`.

Returns
-------
list of int
    the optimal block index for each variable (i.e. the first block index in
    which the correlation has been removed).  If NaN, then the statistics
    provided were not sufficient to estimate the correlation length and more
    data should be collected.

Notes
-----
[Wolff]_ (Eq 47) and [Lee]_ et al. (Eq 14) give the optimal block size to be

.. math::

    B^3 = 2 n n_{\\text{corr}}^2

where `n` is the number of data points in the data set, `B` is the number of
data points in each 'block' (ie the data set has been divided into `n/B`
contiguous blocks) and `n_{\\text{corr}}`.
[todo] - describe n_corr.
Following the scheme proposed by [Lee]_ et al., we hence look for the largest
block size which satisfies 

.. math::

    B^3 >= 2 n n_{\\text{corr}}^2.

From Eq 13 in [Lee]_ et al. (which they cast in terms of the variance):

.. math::

    n_{\\text{err}} SE = SE_{\\text{true}}

where the 'error factor', `n_{\\text{err}}`, is the square root of the
estimated correlation length,  `SE` is the standard error of the data set and
`SE_{\\text{true}}` is the true standard error once the correlation length has
been taken into account.  Hence the condition becomes:

.. math::

    B^3 >= 2 n (SE(B) / SE(0))^4

where `SE(B)` is the estimate of the standard error of the data divided in
blocks of size `B`.

I am grateful to Will Vigor for discussions and the initial implementation.

References
----------
.. [Wolff] "Monte Carlo errors with less errors", U. Wolff, Comput. Phys. Commun.
       156, 143 (2004) and arXiv:hep-lat/0306017.
.. [Lee] "Strategies for improving the efficiency of quantum Monte Carlo
       calculations", R.M. Lee, G.J. Conduit, N. Nemec, P. Lopez Rios,
       N.D.  Drummond, Phys. Rev. E. 83, 066706 (2011).
'''

    # Get the number of variables by looking at the number of means calculated
    # in the first stats entry.
    nvariables = stats[0][2].size
    optimal_block = [float('NaN')]*nvariables
    # If the data was just of a single variable, then the numpy arrays returned
    # by blocking are all 0-dimensions.  Make sure they're 1-D so we can use
    # enumerate safely (just to keep the code short).
    std_err_first = numpy.array(stats[0][4], ndmin=1)
    for (iblock, data_len, mean, cov, std_err, std_err_err) in reversed(stats):
        # 2**iblock data points per block.
        B3 = 2**(3*iblock)
        std_err = numpy.array(std_err, ndmin=1)
        for (i, var_std_err) in enumerate(std_err):
            if B3 > 2*ndata*(var_std_err/std_err_first[i])**4:
                optimal_block[i] = iblock

    return optimal_block

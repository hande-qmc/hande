"""Creating mock dataframes/series for testing."""
from typing import List, Union
import numpy as np
import pandas as pd


def create_qmc_frame(rng, cols: List[str], means: List[float],
                     sine_periods: List[float], noise_facs: List[float],
                     frac_not_convergeds: Union[int, List[float]] = 0,
                     num_mc_its: int = 30) -> pd.DataFrame:
    """Create correlated time series mock qmc data table.

    Parameters
    ----------
    rng : np.random._generator.Generator
        Random number generator used.  Passed in seeded.
    cols : list of strings
        Elements are the names of the columns.  Number of elements
        is equal to the number of columns in qmc_df.
    means : list of floats
        The values are the rough means of the these columns in cols
        when converged.
    sine_periods : list of floats
        sine_period in sin(counter * np.pi/sine_period), modelling
        correlated time series with length of number of columns in
        qmc_df.
    noise_facs : list of floats
        Prefactor of Gaussian random noise.
    frac_not_convergeds : Union[int, List[float]]
        What fraction of the correlated time series should still
        increase/decrease.  If =0, then a list with zeros is created
        with lengths of number of columns in qmc_df.
    num_mc_its : integer
        Number of rows/Monte Carlo iterations.

    Returns
    -------
    qmc_df : :class:`pandas.DataFrame`
        Correlated time series mock qmc data table with passed in
        cols as columns.
        Format:
        cols[0] cols[1]     ...     cols[-1]
        ---------------     ...     -----------------
        all elements of type <class 'numpy.float64'>
    """
    if frac_not_convergeds == 0:
        frac_not_convergeds = [0.0 for _ in range(len(means))]
    qmc_dict = {}
    for (col, mean, sine_period, noise_fac, frac_not_converged) in\
            zip(cols, means, sine_periods, noise_facs, frac_not_convergeds):
        if frac_not_converged > 1e-9:
            # frac_not_converged of the data plus/minus some random
            # number of data points are not converged.
            x_stdn = rng.standard_normal()
            not_converged_its = min(
                int(frac_not_converged*num_mc_its + 0.1*num_mc_its*x_stdn),
                num_mc_its
            )
            log_mean = np.log(abs(mean))
            ncits = not_converged_its
            # Model an exponential increase in magnitude.
            qmc_dict[col] = [
                np.sign(mean) *
                np.exp(log_mean * (i+0.2*rng.standard_normal())/ncits)
                for i in range(ncits)
            ]
        else:
            not_converged_its = 0
            qmc_dict[col] = []
        # Create mock correlated noisy data with a sine curve plus
        # random, Gaussian noise, shifted by a constant.
        qmc_dict[col].extend([
            mean
            + 0.1*abs(int(mean))*np.sin(i*np.pi/sine_period)
            + noise_fac*rng.standard_normal()
            for i in range(num_mc_its - not_converged_its)
        ])
    return pd.DataFrame.from_dict(qmc_dict)


def create_reblock_frame(rng, cols: List[str], means: List[float],
                         it_optbls: List[int],
                         num_reblock_its: int = 10) -> pd.DataFrame:
    """Create mock covariance dataframe.

    The reblocked means are generated by a given mean plus Gaussian
    noise which decays with reblock iteration.
    Their standard errors at first increase, until it_optbl-1, then
    they stay roughly constant as a fraction of the mean plus
    Gaussian noise.
    The standard error errors are a fraction of the standard errors
    plus Gaussian noise.
    The optimal block has empty strings except for '<---' at the
    it_optbl.

    Parameters
    ----------
    rng : np.random._generator.Generator
        Random number generator used.  Passed in seeded.
    cols : list of strings
        Elements are the names of the columns/rows.  Number of
        elements is equal to the number of columns in cov_df.
    means : list of floats
        The values are the rough means of the these columns in cols
        when converged.
    it_optbls : list of integers
        iterations where optimal block is.  Note that if it_optbls =
        num_reblock_its, there is no optimal block for that column.
    num_reblock_its : integer
        Number of reblock iterations.

    Returns
    -------
    reblock_df : :class:`pandas.DataFrame`
        DataFrame containg the mean and its error at every reblock
        iteration forthe cols passed.
        Format:
        (entries are of type <class 'numpy.float64'> for columns
         'mean', 'standard error', 'standard error error' and of
         type  <class 'str'> for 'optimal block')

                            cols[0]                             ...
        reblock    [mean standard error standard error error
                   optimal block]
        0
        1
        ...
    """
    col_dfs = []
    n_reb_its = num_reblock_its
    for (col, mean, it_optbl) in zip(cols, means, it_optbls):
        reb_means = [
            mean + (n_reb_its+1-i)*0.01*mean*rng.standard_normal()
            for i in range(n_reb_its)
        ]
        reb_stes = (
            [abs((i+1)*0.005*mean + 0.001*mean*rng.standard_normal())
             for i in range(it_optbl-1)]
            + [abs((it_optbl+1)*0.005*mean
                   + 0.001*mean*rng.standard_normal())
               for _ in range(min(n_reb_its-it_optbl+1, n_reb_its))]
        )
        # where the min is there for the case optbl=0
        if col != 'Proj. Energy':
            reb_stees = [
                abs(0.1*reb_ste*(1.0 + rng.standard_normal()))
                for reb_ste in reb_stes
            ]
        else:
            reb_stees = None
        reb_optb = [
            '' if i != it_optbl else '<---' for i in range(n_reb_its)
        ]
        col_dfs.append(pd.DataFrame.from_dict({
            'mean': reb_means,
            'standard error': reb_stes,
            'standard error error': reb_stees,
            'optimal block': reb_optb
        }))
    reblock_df = pd.concat(col_dfs, axis=1, keys=cols)
    reblock_df.index.name = 'reblock'
    return reblock_df


def create_cov_frame(rng, cols: List[str], means: List[float],
                     num_reblock_its: int = 10,
                     index_name: str = 'reblock') -> pd.DataFrame:
    """Create mock covariance dataframe.

    Parameters
    ----------
    rng : np.random._generator.Generator
        Random number generator used.  Passed in seeded.
    cols : list of strings
        Elements are the names of the columns/rows.  Number of
        elements is equal to the number of columns in cov_df.
    means : list of floats
        The values are the rough means of the these columns in cols
        when  converged.
    num_reblock_its : integer
        Number of reblock iterations.
    index_name : string
        Name of index (outer level). By default 'reblock'.

    Returns
    -------
    cov_df : :class:`pandas.DataFrame`
        DataFrame containg the covariance matrix between the cols
        for each reblock iteration.
        Format:
        (with example numbers, entries are of
         <class 'numpy.float64'>)
                            cols[0]          cols[1]   ...  cols[-1]
        reblock
        0       cols[0]    3495.988001    131.070144   ... 30.573263
                cols[1]     159.113260    365.103745   ...  2.866295
                ...           ...           ...        ...   ...
                cols[-1]     57.778278     19.470779   ...  0.925906
        1       cols[0]    3495.988001    131.070144   ... 30.573263
                cols[1]     159.113260    365.103745   ...  2.866295
                ...           ...           ...        ...    ...
                cols[-1]     57.778278     19.470779   ...  0.925906
        ...
    """
    covs = []
    keys = []
    for i in range(num_reblock_its):
        # Start with a uniform array.
        cov = rng.uniform(size=(len(cols), len(cols)))
        # Make values a bit more realistic by scaling the columns
        # and rows by the abs. values of the means of the involved
        # reblock columns.  Ignore signs here, just use magnitudes.
        for j, mean in enumerate(means):
            cov[j] *= abs(mean)
            cov[:, j] *= abs(mean)
        covs.append(pd.DataFrame(cov, columns=cols, index=cols))
        keys.append(i)
    return pd.concat(covs, keys=keys, names=[index_name, None])


def create_date_length_series(rng, num_reblock_its: int = 10) -> pd.DataFrame:
    """Create mock data_length dataseries.

    Parameters
    ----------
    rng : np.random._generator.Generator
        Random number generator used.  Passed in seeded.
    num_reblock_its : integer
        Number of reblock iterations.

    Returns
    -------
    data_length : :class:`pandas.Series`
        Series containg the number of blocks for reblocking, where
        each row represents another reblock iteration.
        Format:
        - index named "reblock"
        - all elements of type <class 'numpy.int64'>
    """
    # Numbers get divided by 2 and then rounded down as we move
    # down the dataseries.
    data_l = []
    prev = 1
    for _ in range(num_reblock_its):
        new_element = 2*prev
        if rng.uniform() > 0.5:
            new_element += 1
        prev = new_element
        data_l.append(new_element)
    # reverse as we have created the series the wrong way round
    data_length = pd.Series(data_l[::-1])
    data_length.index.name = 'reblock'
    return data_length

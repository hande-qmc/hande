"""Functions to find starting iteration for analysis."""
from typing import Dict, List
import warnings
import math
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pyblock
from pyhande.helpers.simple_callables import RaiseValueError


def _show_starting_iterations_graph(
        data: pd.DataFrame, it_key: str, col_to_show: str,
        starting_it: int) -> None:
    """
    Show a plot of data of specified column with starting iteration.

    Parameters
    ----------
    data : pd.DataFrame
        Data where starting iteration was found.
    it_key : str
        Key of iteration column.
    col_to_show : str
        Key of column to plot here.
    starting_it : int
        Suggested/found starting iteration.
    """
    plt.xlabel(it_key)
    plt.ylabel(col_to_show)
    plt.plot(data[it_key], data[col_to_show], 'b-', label='data')
    plt.axvline(starting_it, color='r',
                label='Suggested starting iteration')
    plt.legend(loc='best')
    plt.show()


def _get_blocking_loss(start_ind: int, data: pd.DataFrame) -> pd.Series:
    """
    Get blocking loss which is to be minimised.

    Loss here is defined as the fractional standard error over the
    square root of the number of data points used for the blocking.

    Parameters
    ----------
    start_ind : int
        Index where blocking will start.  The blocking loss of that
        analysis will be found.
    data : pd.DataFrame
        Data to be blocked, only including columns that are blocked.

    Returns
    -------
    pd.Series
        Index consists of the columns in `data` considered when finding
        starting iteration.  Values are the losses for each column
        respectively.  If loss calculation fails, respective value is
        NaN.
    """
    (_, reblock, _) = pyblock.pd_utils.reblock(data.iloc[start_ind:])
    opt_block = pyblock.pd_utils.reblock_summary(reblock)
    if all([col in opt_block.index for col in data.columns]):
        loss = ((opt_block['standard error error'] /
                 opt_block['standard error']) /
                math.sqrt(float(len(data.iloc[start_ind:]))))
        if not loss.isna().any():
            return loss.astype('float64')
    # 'nan' values can easily be ignored by _grid_search.
    return pd.Series([float('nan')]*len(data.columns), index=data.columns)


def _grid_search(data: pd.DataFrame, grid_size: int, min_ind: int,
                 max_ind: int) -> int:
    """
    Do log adaptive grid search between `min_ind` and `max_ind`.

    Parameters
    ----------
    data : pd.DataFrame
         Data to be blocked, only including columns that are blocked.
    grid_size : int
        Number of logarithmically spaced grid points per run.
    min_ind : int
        Minimum value of to be found `start_ind`.
    max_ind : int
        Maximum value of to be found `start_ind`.

    Returns
    -------
    int
        Found `start_ind`, index where blocking will start.
    """
    while grid_size > 2:
        grid_pts = np.logspace(
            np.log10(min_ind), np.log10(int(max_ind)), grid_size)
        # The grid points correspond to possible - discrete - start
        # indices so resolution is 1 at best.
        while int(grid_pts[0]) == int(grid_pts[1]) and grid_size > 2:
            grid_pts = np.delete(grid_pts, 0)
            grid_size -= 1
        if int(grid_pts[0]) == int(grid_pts[1]) and grid_size == 2:
            # Only two grid points remain and they are identical.
            # Can stop looping and return their value.
            return int(grid_pts[0])
        losses = pd.concat(
            [_get_blocking_loss(int(grid_pt), data) for grid_pt in grid_pts],
            keys=list(map(int, grid_pts)), axis=1)
        # idxmin(axis=1) will compare values in each row, giving the
        # losses column name where the minimum is.  The losses column
        # names correspond to the grid points.  Of those, we want the
        # lowest grid point (.min()).
        poss_min = losses.idxmin(axis=1).min()
        # Find next minimum if the one above is ignored.  Of the grid
        # point found, select highest to be most conservative.
        poss_max = losses.drop(columns=poss_min).idxmin(axis=1).max()
        # Sort.
        poss_min, poss_max = ((poss_min, poss_max) if poss_min < poss_max
                              else (poss_max, poss_min))
        min_ind = max(min_ind, poss_min)
        max_ind = min(max_ind, poss_max)
        if min_ind == int(grid_pts[0]) and max_ind == int(grid_pts[-1]):
            break
    return max_ind


def find_starting_iteration_blocking(
        data: pd.DataFrame, end_it: int, it_key: str, cols: List[str],
        hybrid_col: str, start_max_frac: float = 0.8,
        grid_size: int = 10, number_of_reblocks_to_cut_off: int = 1,
        show_graph: bool = False) -> int:
    """
    Find the best iteration to start analysing CCMC/FCIQMC data.

    It first excludes data before not all data in all columns specified
    in `cols` are varying and after `end_it`.  Then it searches for the
    starting iteration using an adaptive grid search on a log scale
    since we assume that the starting iteration is closer to the
    beginning than the end of the available data.  During the search, a
    loss function is minimised.  The loss is the fractional error over
    number of data involved in the blocking for each data column in
    `cols`.

    This implementation is based on an older version in pyhande/lazy.py.

    .. warning::

        Use with caution, check whether output is sensible and adjust
        parameters if necessary.

    Parameters
    ----------
    data : pd.DataFrame
        QMC data, e.g. as extracted by extract.py.  Has to contain
        columns with key `it_key` and columns in `cols`, used for
        blocking.
    end_it : int
        Last iteration to be considered in blocking.
    it_key : str
        Key of column containing MC iterations.
    cols : List[str]
        List of keys of columns involved in blocking.
    hybrid_col : str
        Ignored here, for common interface.
    start_max_frac : float, optional
        The start iterations found has to be in the first
        `start_max_frac` fraction of the data between
        the point where all columns in `cols` have started varying and
        `end_it`.  This prevents finding a starting iteration too close
        to the end.  Has to be between 0.00001 and 1.0.
        The default is 0.8.
    grid_size : int, optional
        Number of logarithmically spaced grid points per run.
        The default is 10.
    number_of_reblocks_to_cut_off : int, optional
        To be extra sure, cut off a few reblocks to make sure data after
        starting iteration is truly in equilibrium. Cannot be negative.
        The default is 1.
    show_graph : bool, optional
        If True, show a graph showing the columns with key `cols[0]` as
        a function of iterations.  The suggested starting iteration is
        highlighted.  The default is False.

    Raises
    ------
    ValueError
        If `start_max_frac` or
        `number_of_reblocks_to_cut_off` are out of range.
    RuntimeError
        If not all columns with keys in `cols` have started varying in
        `data` or if suitable starting iteration was not found.

    Returns
    -------
    int
        Suggestion iteration in columns `it_key` where analysis should
        start.
    """
    # Check some inputs.
    if (start_max_frac < 0.00001 or
            start_max_frac > 1.0):
        raise ValueError("0.00001 < start_max_frac < 1 not "
                         "satisfied!")
    if number_of_reblocks_to_cut_off < 0:
        raise ValueError("'number_of_reblocks_to_cut_off' can't be negative!")

    # Data cleaning.
    # Remove iterations passed the specified end iteration and make
    # sure all cols have started varying in dataset excluding data
    # before.
    data = data[data[it_key] <= end_it]
    max_varying_it = 0
    for col in cols:
        if not any(data[data[col] != data[col].iloc[0]]):
            raise RuntimeError(f"{col} has not started varying in considered "
                               "dataset.")
        max_varying_it = max(
            max_varying_it,
            data[data[col] != data[col].iloc[0]][it_key].iloc[0])
    data = data[data[it_key] >= max_varying_it]
    # Finding starting iteration.
    # Do grid search to find the index of the starting iteration.
    start_ind = _grid_search(
        data[cols], grid_size, 1, int(start_max_frac*len(data)) + 1)
    # Search has failed if index is too close to the end.
    if start_ind > int(start_max_frac*len(data)):
        raise RuntimeError(f"Failed to find starting iteration. ")
    # Discarding number_of_reblocks_to_cut_off reblocks.
    (_, reblock, _) = pyblock.pd_utils.reblock(data[cols].iloc[start_ind:])
    opt_ind = pyblock.pd_utils.optimal_block(reblock)
    discard_indx = 2**opt_ind * number_of_reblocks_to_cut_off
    # Converting to iteration.
    starting_it = data['iterations'].iloc[start_ind + discard_indx]

    # Show plot if desired, aiding judgment whether to trust estimate.
    if show_graph:
        # Note that the data has non varying phase cut off!
        _show_starting_iterations_graph(data, it_key, cols[0], starting_it)
    return starting_it


def find_starting_iteration_mser_min(
        data: pd.DataFrame, end_it: int, it_key: str, cols: List[str],
        hybrid_col: str, start_max_frac: float = 0.84,
        n_blocks: int = 100) -> int:
    r'''Estimate starting iteration with MSER minimization scheme.

    .. warning::

        Use with caution, check whether output is sensible and adjust
        parameters if necessary.

    This function gives an optimal estimation of the starting
    interations based on MSER minimization heuristics.
    This methods decides the starting iterations :math:`d` as minimizing
    an evaluation function
    MSER(:math:`d`) =
    :math:`\Sigma_{i=1}^{n-d} ( X_{i+d} - X_{mean}(d) ) / (n-d)^2`.
    Here, :math:`n` is length of time-series, :math:`X_i` is
    `eval_ratio['num']` / `eval_ratio['denom']` of :math:`i`-th step,
    and :math:`X_{mean}` is the average of :math:`X_i` after the
    :math:`d`-th step.

    This is a reformatted and altered version of a previous
    implementation in lazy.py by Tom Ichibha.

    Parameters
    ----------
    data : :class:`pandas.DataFrame`
        Calculation output of a FCIQMC or CCMC calculation.
    end_it : int
        Last iteration to be considered in blocking.
    it_key : str
        Key of column containing MC iterations.
    cols : List[str]
        Ignored here.  Keep for common interface.
    hybrid_col: str
        Column in data to be analysed here, e.g. 'Inst. Proj. Energy'.
    start_max_frac : float
        MSER(d) may oscillate when become unreanably small
        when :math:`n-d` is large. Thus, we calculate MSER(:math:`d`)
        for :math:`d` < (:math:`n` * start_max_frac) and
        give the optimal estimation of the starting iterations
        only in this range of :math:`d`.
        The default is 0.84.
    n_blocks : int
        This analysis takes long time when :math:`n` is large.
        Thus, we pick up :math:`d` for every 'n_blocks' samples,
        calculate MSER(:math:`d`), and decide the optimal estimation of
        the starting iterations only from these `d`.
        The default is 100.

    Returns
    -------
    starting_it: int
        Iteration from which to start reblocking analysis for this
        calculation.
    '''
    data = data[data[it_key] <= end_it]
    inst_ratio = data[hybrid_col]

    mser_min = float('inf')
    for i in range(n_blocks):
        start_ind = int(i*(len(inst_ratio)*start_max_frac)/n_blocks)
        mser = (np.var(inst_ratio[start_ind:len(inst_ratio)]) /
                (len(inst_ratio)-start_ind))
        if mser < mser_min:
            mser_min = mser
            starting_it = data[it_key].loc[start_ind]
            final_start_ind = start_ind

    if final_start_ind > len(inst_ratio)*(start_max_frac**2):
        warnings.warn(
            f"Instantaneous ratio '{hybrid_col}' may not be "
            "converged.  MSER min. may underestimate the starting iteration.  "
            "Check!")
    return starting_it


def select_find_start(key: str):
    """Select find_starting_iteration function to use.

    Parameters
    ----------
    key : str
        Key linked to find_starting_iteration.

    Returns
    -------
    Find_starting_iteration function.
    """
    return {'blocking': find_starting_iteration_blocking,
            'mser': find_starting_iteration_mser_min}.get(
                key, RaiseValueError("The find start iteration selected in "
                                     f"'start_its', '{key}', is not "
                                     "available!"))

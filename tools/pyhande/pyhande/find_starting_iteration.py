"""Functions to find starting iteration for analysis."""
from typing import List
import copy
import warnings
import math
import pandas as pd
import matplotlib.pyplot as plt
import pyblock
import pyhande.analysis as analysis

def _check_find_starting_iteration_blocking_inputs(
        number_of_intervals, min_number_of_blockings,
        number_of_reblocks_to_cut_off, start_position_in_data_max_frac):
    """Check parameters of find_starting_iteration_blocking."""
    if number_of_intervals <= 0:
        raise ValueError("'number_of_intervals' has to be greater than zero!")

    if number_of_reblocks_to_cut_off < 0:
        raise ValueError("'number_of_reblocks_to_cut_off' can't be negative!")

    if (start_position_in_data_max_frac < 0.00001 or
            start_position_in_data_max_frac > 1.0):
        raise ValueError("0.00001 < start_position_in_data_max_frac < 1 not "
                         "satisfied!")

    if min_number_of_blockings <= 0:
        raise ValueError("'min_number_of_blockings' has to be greater than "
                         "zero!")

    if min_number_of_blockings > number_of_intervals:
        raise ValueError("'min_number_of_blockings' can't be greater than "
                         "'number_of_intervals'!")

def find_starting_iteration_blocking(
        data: pd.DataFrame, end_it: int, it_key: str, cols: List[str],
        number_of_intervals: int = 300, min_number_of_blockings: int = 30,
        number_of_reblocks_to_cut_off: int = 1,
        start_position_in_data_max_frac: float = 0.8, verbose: int = 0,
        show_graph: bool = False) -> int:
    '''Find the best iteration to start analysing CCMC/FCIQMC data.

    This is a modification of a previous function in pyhande.lazy.py
    'find_starting_iteration'.

    .. warning::

        Use with caution, check whether output is sensible and adjust
        parameters if necessary.
    '''
    # Check inputs
    _check_find_starting_iteration_blocking_inputs(
        number_of_intervals, min_number_of_blockings,
        number_of_reblocks_to_cut_off, start_position_in_data_max_frac)

    data = data[data[it_key] <= end_it]

    # Make sure all cols have started varying in dataset and exclude data
    # before.
    max_varying_it = end_it
    for col in cols:
        if not any(data[data[col] != data[col].iloc[0]]):
            raise RuntimeError(f"{col} has not started varying in considered "
                               "dataset.")
        max_varying_it = min(
            max_varying_it,
            data[data[col] != data[col].iloc[0]][it_key].iloc[0])
    data = data[data[it_key] >= max_varying_it]

    # Check we have enough data to screen:
    if len(data) < number_of_intervals:
        warnings.warn(f"Length of data to be analysed, {len(data)}, is less "
                      f"than 'number_of_intervals', {number_of_intervals}. "
                      "Setting 'number_of_intervals' equal to length of data.")
        number_of_intervals = len(data)

    interval_step = int(len(data)/number_of_intervals)
    min_index = -1
    min_error_frac_weighted = pd.Series([float('inf')]*len(cols), index=cols)
    starting_it_found = False
    dat_c = copy.copy(data[cols])
    for k in range(int(number_of_intervals/min_number_of_blockings)):
        for j in range(
                k*min_number_of_blockings, (k+1)*min_number_of_blockings):
            (_, reblock, _) = pyblock.pd_utils.reblock(dat_c)
            (opt_block, no_opt_block) = analysis.qmc_summary(reblock, cols)

            if not no_opt_block:
                err_frac_weighted = ((
                    opt_block.loc['standard error error']/
                    opt_block.loc['standard error']
                    )/math.sqrt(float(len(dat_c))))
                if (err_frac_weighted < min_error_frac_weighted).any():
                    min_index = j
                    min_error_frac_weighted = err_frac_weighted.copy()
                    opt_ind = pyblock.pd_utils.optimal_block(reblock[cols])

            if verbose > 1:
                print(f"Blocking attempt: {j}. Blocking from: "
                      f"{dat_c[it_key].iloc[0]}.")
            dat_c = dat_c.iloc[interval_step:]

        if -1 < min_index < (start_position_in_data_max_frac * j):
            # Also discard the frst n=number_of_reblocks_to_cut_off of data to
            # be conservative.  This amounts to removing n autocorrelation
            # lengths.
            discard_indx = 2**opt_ind * number_of_reblocks_to_cut_off
            starting_it = data[it_key].iloc[
                discard_indx + min_index*interval_step
                ]
            if data['iterations'].iloc[-1] <= starting_it:
                raise ValueError("Too much cut off! Data is not converged or "
                                 "use a smaller "
                                 "'number_of_reblocks_to_cut_off'.")
            starting_it_found = True
            break

    if not starting_it_found:
        raise RuntimeError("Failed to find starting iteration. The "
                           "calculation might not be converged.")

    if show_graph:
        plt.xlabel(it_key)
        plt.ylabel(cols[0])
        plt.plot(data[it_key], data[cols[0]], 'b-', label='data')
        plt.axvline(starting_it, color='r',
                    label='Suggested starting iteration')
        plt.legend(loc='best')
        plt.show()

    return starting_it

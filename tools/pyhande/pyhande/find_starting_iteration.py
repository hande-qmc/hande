from os import path
import warnings
import sys

try:
    import pyblock
except ImportError:
    sys.path.append(path.join(path.abspath(path.dirname(__file__)), '../../pyblock'))
    import pyblock

import math
import matplotlib.pyplot as plt
import pyhande.lazy

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
            info = pyhande.lazy.lazy_block(data, md, start, extract_psips=True)

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

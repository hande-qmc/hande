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
import pyhande.extract

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

    # Check whether the input is valid. If not, set to default values.
    if (frac_screen_interval <= 0):
        warnings.warn("frac_screen_interval > 0 is not satisfied, \
            frac_screen_interval set to default: 500")
        frac_screen_interval = 500

    if (number_of_reblocks_to_cut_off < 0):
        warnings.warn("number_of_reblocks_to_cut_off >= 0 is not satisfied, \
            number_of_reblocks_to_cut_off set to default: 1")
        number_of_reblocks_to_cut_off = 1

    if (pos_min_frac < 0.00001) or (pos_min_frac > 1.0):
        warnings.warn("0.00001 < pos_min_frac < 1 not satisfied, \
            pos_min_frac set to default: 0.5")
        pos_min_frac  = 0.5

    if (number_of_reblockings <= 0):
        warnings.warn("number_of_reblockings > 0 not satisfied, \
            number_of_reblockings set to default: 50")
        number_of_reblockings = 50

    if (number_of_reblockings > frac_screen_interval):
        warnings.warn("number_of_reblockings has to be less than or equal to \
            frac_screen_interval, it is now set to be equal to \
            frac_screen_interval.")
        number_of_reblockings = frac_screen_interval

    # Find the point the shift began to vary.
    variable_shift = data['Shift'] != data['Shift'].iloc[0]
    if variable_shift.any():
        shift_variation_indx = data[variable_shift]['iterations'].index[0]
    else:
        raise RuntimeError("Shift has not started to vary in dataset!")

    # Check we have enough data to screen:
    if(data['Shift'].size - shift_variation_indx) < \
            (frac_screen_interval):
        # i.e. data where shift is not equal to initial value is less than
        # frac_screen_interval, i.e. we cannot screen adequately.
        warnings.warn("Files %s contain less data than \
            we wanted to screen. Will continues but frac_screen_interval \
            is less than one data point." % (','.join(outputfiles)))

    # Find the MC iteration at which shift starts to vary.
    iteration_shift_variation_start = \
            data['iterations'].iloc[shift_variation_indx]

    step = int((data['iterations'].iloc[-1] - \
            iteration_shift_variation_start )/frac_screen_interval)

    shift_error_errors = []
    starting_iteration_found = False

    for k in range(int(frac_screen_interval/number_of_reblockings)):

        for j in range(k*number_of_reblockings, (k+1)*number_of_reblockings):
            start = iteration_shift_variation_start + j*step
            info = pyhande.lazy.lazy_block(data, md, start, extract_psips=True)

            if info.no_opt_block:
                # Add a large enough (imaginary) shift error error
                # so this does not win in this search.
                shift_error_errors.append(float('inf'))
            else:
                s_err_err = info.opt_block["standard error error"]["Shift"]
                shift_error_errors.append(s_err_err)

            if verbose:
                print("Blocking attempt: %i. Error in the shift error: %f"
                        % (j, shift_error_errors[-1]))

        min_error_error = float('inf')
        min_index = len(shift_error_errors)
        for j_min_ind in range(len(shift_error_errors) - 1, 0, -1):
            if shift_error_errors[j_min_ind] < min_error_error:
                min_index = j_min_ind
                min_error_error = shift_error_errors[j_min_ind]
        if min_index == len(shift_error_errors):
            raise ValueError("Failed to find minimum error in the shift error!")

        if int(min_index * pos_min_frac) < j:
            start = iteration_shift_variation_start + (min_index*step)
            info = pyhande.lazy.lazy_block(data, md, start, extract_psips=True)
            opt_ind = pyblock.pd_utils.optimal_block(info.reblock['Shift'])
            data_points_in_block = 2**opt_ind
            # [review] - JSS: if statement over 5 lines with various indent levels - no chance of easy comprehension.
            if ((iteration_shift_variation_start + min_index*step + \
                number_of_reblocks_to_cut_off*data_points_in_block* \
                (data['iterations'].iloc[1] - \
                data['iterations'].iloc[0])) \
                >= data['iterations'].iloc[-1]) :
                raise ValueError("Too much cut off! Try again with \
                    a smaller (positive) number_of_reblocks_to_cut_off, \
                    or data is not fit for analysis.")
            starting_iteration = (min_index*step + \
                data['iterations'].iloc[shift_variation_indx + \
                number_of_reblocks_to_cut_off*data_points_in_block])
            starting_iteration_found = True
            break

    if not starting_iteration_found:
        raise RuntimeError("Failed to find starting iteration. The data \
        might not be converged.")

    if show_graph:
        starting_iteration_list = []
        starting_no_correlation = []
        for iteration in range(0,data['Shift'].size):
            starting_iteration_list.append(starting_iteration)
            starting_no_correlation.append(iteration_shift_variation_start \
                + min_index*step)
        plt.xlabel('Iterations')
        plt.ylabel('Shift')
        plt.plot(data['iterations'], data['Shift'], 'xb', label='data')
        plt.plot(starting_iteration_list, data['Shift'], 'r',
                 label='suggested starting iteration')
        plt.plot(starting_no_correlation, data['Shift'], 'g', \
                 label='starting iteration without subtracting reblocks')
        plt.legend(loc='best')
        plt.show()

    return starting_iteration

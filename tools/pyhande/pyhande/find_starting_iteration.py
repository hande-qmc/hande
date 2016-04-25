import warnings

try:
    import pyblock
except ImportError:
    sys.path.append(path.join(path.abspath(path.dirname(__file__)), '../../pyblock'))
    import pyblock

import reblock_hande
import math
import matplotlib.pyplot as plt
import pyhande.lazy
import pyhande.extract

def find_starting_iteration(outputfiles, frac_screen_interval=500,  
    number_of_reblockings = 50, number_of_reblocks_to_cut_off = 1, 
    pos_min_frac = 0.5, verbose = False, show_graph = False) :

    '''
!!!Use with caution, check whether output is sensible and adjust 
parameters if necessary!!!

Function to find the best iteration to start analysing CCMC/FCIQMC 
data. The output of this function can be passed to reblock_hande as the
starting iteration. 
The best iteration is found by first finding the iteration that if used as the 
starting iteration gives the lowest error in the error of the shift. To be
conservative, an extra few blocks from the reblocking analysis can be cut off as
well, the number can be set by the user as "number_of_reblocks_to_cut_off". To
be efficient, first the point is found where the shift started to vary. Any
iteration before that point is discarded. The remaining iterations are then
divided up into "frac_screen_interval" intervals, set by the user. A certain
number ("number_of_reblockings", set by user) of reblocking analyses are done,
starting at the point where the shift started to vary and then progressing in
steps set by the lengths of an interval. If, after having done this certain
number of reblocking analyses with different starting points, the minimum error
in the error is in the first "pos_min_frac" fraction of the points used as
starting points so far, the search for the minimum is ended. If this is not the
case, the search continues with another "number_of_reblocking" trial starting
points. If the interval has been reached and the minimum is not in the first
"pos_min_frac" fraction of the data (excluding the iterations before the shift),
the search has failed.

Parameters
----------
outputfiles: list of strings
    List of CCMC/FCIQMC output file names. If more than one, they should be
    restarted versions from each other for this to make sense.
frac_screen_interval: integer
    Number of intervals the iterations from where the shift started to vary to
    the end are divided up into. Has to be greater than zero.
number_of_reblockings: integer
    Number of reblocking analyses done in steps set by the width of an interval
    before it is checked whether suitable minimum error in the error has been
    found. Has to be greater than zero.
number_of_reblocks_to_cut_off: integer
    Number of reblocking analysis blocks to cut off additionally to the data
    before the best iteration with the lowest error in the error. Has to be non
    negative. It is highly recommended to not set this to zero.
pos_min_frac: float/double 
    The minimum has to be in the first pos_min_frac part of the tested data to
    be taken as the true minimum. Has be to greater than a small number (here
    0.00001) and can at most be
    equal to one.
verbose: boolean
    Determines whether extra information is printed out.
show_graph: boolean
    Determines whether a window showing the shift vs iteration graph pops up
    highlighting where the minimum was found and - after also excluding some
    reblocking blocks - which iteration was found as the best starting iteration
    to use in reblocking analyses.

Returns
-------
starting_iterations: list of float/double
    Iterations to start reblocking analyses for this output from. Each iteration
    corresponds to an outputfile/chain of restarted outputfiles given as
    input in outputfiles.
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

    # [review] - JSS: this assumes all data sets are from CCMC or FCIQMC
    # [review] - JSS: calculations.  Fix to only get CCMC/FCIQMC data.
    # [review] - JSS: actually, I think the best way is to assume the user has
    # [review] - JSS: already done these three lines and require just the data
    # [review] - JSS: list from concat_calcs to be passed in.  This also removes
    # [review] - JSS: the need for the data to be extracted multiple times...
    hande_outs = pyhande.extract.extract_data_sets(outputfiles)
    metadata_sets, data_sets = zip(*hande_outs)
    metadata_list, data = \
        pyhande.lazy.concat_calcs(metadata_sets, data_sets)		

    starting_iterations = []

    for it in range(0,len(data)):
        #finding the point the shift began to vary, at the 
        #"shift_variation_start" iteration.
        # [review] - JSS: simpler way to do this -- see std_analysis.
        start_not_found = True
        shift_variation_start = 0
        i = 0
        # [review] - JSS: don't include debug output (and not compatible with python 2 and 3)
        print data[it]['Shift'].size
        # [review] - JSS: bold move not indenting the while block...
        while ((start_not_found == True) and (data[it]['Shift'].size > i)):  
	    if(math.fabs(data[it]['Shift'].iloc[i] != data[it]['Shift'].iloc[0])):  
                start_not_found = False
	        shift_variation_start = i
	    i += 1
        # [review] - JSS: also bold here.
        if (start_not_found == True) :
            # [review] - JSS: use RuntimeError for this kind of issue.
	    raise BaseException("ERROR: Shift has not started to vary in \
                dataset!")
        # [review] - JSS: not necessary to delete variables explicitly (at least usually)
        del start_not_found
        del i

        #Check we have enough data to screen:
        if(data[it]['Shift'].size - shift_variation_start) < \
            (frac_screen_interval):
            # i.e. data where shift is not equal to initial value is less than
            # frac_screen_interval, i.e. we cannot screen adequately.
            warnings.warn("Files " + outputfiles + " contain less data than \
	        we wanted to screen. Will continues but frac_screen_interval \
	        is less than one data point.")

        # Find the MC iteration at which shift starts to vary.
        # [review] - JSS: indent continuation lines for clarity (your editor should do this automatically...)
        iteration_shift_variation_start = \
	    data[it]['iterations'].iloc[shift_variation_start]

        # [review] - JSS: indent
        step = int((data[it]['iterations'].iloc[-1] - \
	    iteration_shift_variation_start )/frac_screen_interval) 
		
        shift_error_errors = []
        starting_iteration_found = False

        for k in range(0, int(frac_screen_interval/number_of_reblockings)):
		
	    for j in range(k*number_of_reblockings, (k+1)*number_of_reblockings):
	        if verbose:
            	    print "Looping through point #" + str(j) + \
		        " to find a starting iteration."
                info = pyhande.lazy.std_analysis(outputfiles, start = \
		    (iteration_shift_variation_start + j*step), \
                    extract_psips=True)

                if info[it].no_opt_block:
                    # Add a large enough (imaginary) shift error error
                    # so this does not win in this search.
            	    shift_error_errors.append(float('inf'))
	        else:
            	    shift_error_error = \
		        info[it].opt_block["standard error error"]["Shift"]
                    shift_error_errors.append(shift_error_error)
            
            min_error_error = float('inf')
            min_index = len(shift_error_errors)
            for j_min_ind in range(len(shift_error_errors) - 1, 0, -1):
                if(shift_error_errors[j_min_ind] < min_error_error):
                    min_index = j_min_ind
                    min_error_error = shift_error_errors[j_min_ind]
            if min_index == len(shift_error_errors):
                raise ValueError("Calculation to find minimum error in error \
                    failed!")
		    
	    if int(min_index * pos_min_frac) < j:
                info = pyhande.lazy.std_analysis(outputfiles, start =  
		    (iteration_shift_variation_start + (min_index*step)), \
		    extract_psips=True)
	        opt_ind = \
                    pyblock.pd_utils.optimal_block(info[it].reblock['Shift'])
                data_points_in_block = int(info[it].data_len[0] / \
		    info[it].data_len[opt_ind])
                # [review] - JSS: if statement over 5 lines with various indent levels - no chance of easy comprehension.
                if ((iteration_shift_variation_start + min_index*step + \
		    number_of_reblocks_to_cut_off*data_points_in_block* \
		    (data[it]['iterations'].iloc[1] - \
                    data[it]['iterations'].iloc[0])) \
		    >= data[it]['iterations'].iloc[-1]) :
                    raise ValueError("Too much cut off! Try again with \
		        a smaller (positive) number_of_reblocks_to_cut_off, \
	                or data is not fit for analysis.")
	        starting_iteration = (min_index*step + \
		    data[it]['iterations'].iloc[shift_variation_start + \
		    number_of_reblocks_to_cut_off*data_points_in_block])
                starting_iteration_found = True
	        break


        if not starting_iteration_found :
            raise BaseException("Failed to find starting iteration. The data \
            might not be converged.")

        if show_graph == True:
            starting_iteration_list = []
            starting_no_correlation = []
            for iteration in range(0,data[it]['Shift'].size):
                starting_iteration_list.append(starting_iteration)
                starting_no_correlation.append(iteration_shift_variation_start \
                    + min_index*step)
            plt.xlabel('Iterations')
            plt.ylabel('Shift')
            plt.plot(data[it]['iterations'], data[it]['Shift'], 'xb', \
                label='data') 
            plt.plot(starting_iteration_list, data[it]['Shift'], 'r', \
                label = 'suggested starting iteration')
            plt.plot(starting_no_correlation, data[it]['Shift'], 'g', \
                label = 'starting iteration without subtracting reblocks')
            plt.legend(loc = 'best')
            plt.show()

        starting_iterations.append(starting_iteration)
		
    return starting_iterations


#!/usr/bin/python

import math
import re
import sys
import scipy
import scipy.stats

def extract_data(data_file):
    # Holds the estimate of the trace from the 'first simulation'.
    trace1 = []
    # Holds the estimate of the trace from the 'replica simulation'.
    trace2 = []
    # Element-wise multiplication of the above two lists.
    trace_prod = []

    # These lists hold several estimates of only one RDM element at a time.
    # The unnormalised RDM elements from the 'first simulation'.
    elements1 = []
    # The unnormalised RDM elements from the 'replica simulation'.
    elements2 = []
    # Element-wise multiplication of the above two lists.
    elem_prod = []

    # Holds \sum_{ij} q_{ij} w{ij}, where q_{ij} are the unnormalised RDM elements
    # from the first simulation, w_{ij} from the replica simulation.
    # Each element of this list holds this value from one 'pair' of simulations.
    num_estimates = []

    trace_regex = '^Trace'
    # The number of elements above and including the diagonal.
    nelem = 0

    a = 1
    b = 1

    for data_file in data_file:
    
        f = open(data_file)

        # Count the number of lines.
        for line in f:
            nelem = nelem + 1

        # The first line counted was the trace, so uncount it.
        nelem = nelem - 1
        # Calculate number of elements in each row from this value.
        rdm_row_size = int(-0.5 + 0.5*math.sqrt(1 + 8*nelem))

        # Close and reopen so we can loop over all lines again.
        f.close()
        f = open(data_file)

        for line in f:
            if re.match(trace_regex, line):
                i = 0
                values = line.split()
                for value in values[1:]:
                    i += 1
                    # Values assigned to the first simulation.
                    if i%2 == 1:
                        trace1.append(float(value)) 
                    # Values assigned to the replica simulation.
                    if i%2 == 0:
                        trace2.append(float(value)) 

                # If don't have the same number of estimates in both lists, throw the extra one away.
                if len(trace1) != len(trace2):
                    del trace2[-1]

                # For later, build the list and initialise to zero.
                for elem in trace1:
                    num_estimates.append(0.0)

                trace_prod = [x*y for x,y in zip(trace1,trace2)] 

                denom_mean = scipy.mean(trace_prod)
                denom_se = scipy.stats.sem(trace_prod)
            else:
                i = 0
                elements1 = []
                elements2 = []
                values = line.split()
                for value in values[1:]:
                    i += 1
                    # Values assigned to the first simulation.
                    if i%2 == 1:
                        elements1.append(float(value)) 
                    # Values assigned to the replica simulation.
                    if i%2 == 0:
                        elements2.append(float(value)) 

                elem_prod = [x*y for x,y in zip(elements1,elements2)] 

                # Each of the estimates in elem_prod would have come from a different beta loop or iteration
                # when done properly, so add them in as they would be in DMQMC.
                iteration = 1
                for elem in elem_prod:
                    # If on the diagonal of the RDM.
                    if a == b:
                        num_estimates[iteration-1] += elem
                    # If on the off-diagonal, count the element twice.
                    else:
                        num_estimates[iteration-1] += 2*elem
                    iteration += 1

                b += 1
                # If true, the next element will be on a new line of the RDM, on the diagonal.
                if b%rdm_row_size == 1:
                    a += 1
                    b = a

            num_mean = scipy.mean(num_estimates)
            num_se = scipy.stats.sem(num_estimates)
            num_cov = calculate_covariance(num_estimates, num_mean, trace_prod, denom_mean)

    f.close()

    return denom_mean, denom_se, num_mean, num_se, num_cov

def calculate_covariance(num, num_mean, denom, denom_mean):
    cov = 0
    for i in range(len(denom)):
        cov += (num[i] - num_mean)*(denom[i] - denom_mean)
    return cov/(len(denom)*(len(denom)-1))

if __name__ == '__main__':

    data_file = sys.argv[1:]

    denom_mean, denom_se, num_mean, num_se, num_cov = extract_data(data_file)

    # The estimate of Renyi-2 itself.
    r2_estimate = -math.log(num_mean/denom_mean, 2)

    # Error, using the propgation of errors formula for the above expression for Renyi-2
    error = (1/math.log(2)) * scipy.sqrt( (num_se/num_mean)**2 + (denom_se/denom_mean)**2 + (2*num_cov/(num_mean*denom_mean)) )

    print "%.18g" % r2_estimate, error

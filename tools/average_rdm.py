#!/usr/bin/python

import math
import re
import sys
import scipy
import scipy.stats

def extract_data(data_file):
    trace = []
    elements = []
    elem_mean = []
    elem_se = []
    covariances = []

    trace_regex = '^Trace'
    n = 0

    for data_file in data_file:
    
        f = open(data_file)

        for line in f:
            if re.match(trace_regex, line):
                values = line.split()
                for value in values[1:]:
                   trace.append(float(value)) 
                trace_mean = scipy.mean(trace)
                trace_se = scipy.stats.sem(trace)
            else:
                elements = []
                values = line.split()
                for value in values[1:]:
                   elements.append(float(value)) 
                elem_mean.append(scipy.mean(elements))
                elem_se.append(scipy.stats.sem(elements))
                covariances.append(calculate_covariance(elements, elem_mean[n-1], trace, trace_mean))
            n += 1

    # We want the first line of the output to hold the number of elements in each row of the RDM.
    # The number of elements above and including the diagonal of this d-by-d RDM is n-1.
    # It is also true that n-1 = d*(d-1)/2. Putting these together and solving for d (taking
    # the positive solution) gives the following expression.
    rdm_row_size = int(-0.5 + 0.5*math.sqrt(1 + 8*(n-1)))

    # For the first line of the output, print the row size.
    print rdm_row_size

    f.close()

    return trace_mean, trace_se, elem_mean, elem_se, covariances

def calculate_covariance(numerator, numerator_mean, trace, trace_mean):
    cov = 0
    for i in range(len(trace)):
        cov += (numerator[i] - numerator_mean)*(trace[i] - trace_mean)
    return cov/(len(trace)*(len(trace)-1))

if __name__ == '__main__':

    data_file = sys.argv[1:]

    trace_mean, trace_se, elem_mean, elem_se, covariances = extract_data(data_file)

    # Print the averaged value of each element, and the corresponding standard errors.
    for i in range(len(elem_mean)):
        normalised_elem = elem_mean[i]/trace_mean
        if elem_se[i] == 0:
            error = 0.0
        else:
            error = scipy.sqrt((elem_se[i]/elem_mean[i])**2 + (trace_se/trace_mean)**2 - (2*covariances[i]/(elem_mean[i]*trace_mean)) )*abs(normalised_elem)
        print "%.18g" % normalised_elem, error

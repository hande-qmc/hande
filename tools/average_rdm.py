#!/usr/bin/python

import math
import re
import sys
import scipy
import scipy.stats
import linecache

def extract_data(data_file):
    trace = []
    elements = []
    elem_mean = []
    elem_se = []
    covariances = []
    
    first_line_regex = '^# Number of elements in each RDM row:'
    trace_regex = '^# RDM trace'
    
    for data_file in data_file:

        f = open(data_file)

        # Get the system size from the first line of the file.
        line = f.readline()
        if not re.match(first_line_regex, line):
            print 'First line of file is not as expected.'
        else:
            words = line.split()
            rdm_row_size = int(words[8])
            # Number of elements above and including the diagonal.
            rdm_unique_elems = rdm_row_size*(rdm_row_size+1)/2
            num_lines_per_loop = rdm_unique_elems + 3

        # First, get the trace mean and se.
        have_data = True
        num_beta_loops = 0
        # The line of the first RDM trace.
        line_num = 4
        while have_data:
            line = linecache.getline(data_file,line_num)
            linecache.clearcache()
            if not re.match(trace_regex, line):
                have_data = False
            if have_data:
                words = line.split()    
                trace.append(float(words[3]))
                num_beta_loops = num_beta_loops + 1
            # The line of the next RDM trace.
            line_num = line_num + num_lines_per_loop
        trace_mean = scipy.mean(trace)
        trace_se = scipy.stats.sem(trace)

        # Next get the mean and se for each RDM element.

        # Over all RDM elements in a particular loop.
        for i in range(1, rdm_unique_elems+1):
            elements = []
            # Over all beta loops for a given element.
            for j in range(1, num_beta_loops+1):
                line_num = num_lines_per_loop*(j-1) + i + 4
                line = linecache.getline(data_file,line_num)
                linecache.clearcache()
                if len(line) == 0:
                    print "Unexpected empty line found."
                words = line.split()
                elements.append(float(words[0]))
            elem_mean.append(scipy.mean(elements))
            elem_se.append(scipy.stats.sem(elements))
            covariances.append(calculate_covariance(elements, elem_mean[i-1], trace, trace_mean))

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

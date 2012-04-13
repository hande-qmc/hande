#!/usr/bin/python

import sys
import re
import scipy
import scipy.stats

def extract_data(data_files):
    unormalised_concurrence = []
    entropy = []
    trace = []
    entropy_regex = '^ # Von-Neumann Entropy='
    concurrence_regex = '^ # Unnormalised concurrence='
    trace_regex = '^ # RDM trace=' 
    
    for data_file in data_files:

        f = open(data_file)

        have_data = False
        for line in f:
            if re.match(entropy_regex, line):
                words = line.split()
                entropy.append(float(words[3]))
            elif re.match(concurrence_regex, line):
                words = line.split()
                unormalised_concurrence.append(float(words[3]))
            elif re.match(trace_regex, line):
                words = line.split()
                trace.append(float(words[3]))
        f.close()
    return entropy, unormalised_concurrence, trace

def calculate_covariance(numerator, numerator_mean, trace, trace_mean):
    cov = 0
    for i in range(len(trace)):
        cov += (numerator[i] - numerator_mean)*(trace[i] - trace_mean)
    return cov/(len(trace)**2)

def calculate_stats_ratio(numerator, trace):
    numerator_mean = scipy.mean(numerator)
    trace_mean = scipy.mean(trace)
    numerator_se = scipy.stats.sem(numerator)
    trace_se = scipy.stats.sem(trace)
    mean = numerator_mean/trace_mean
    cov = calculate_covariance(numerator, numerator_mean, trace, trace_mean)
    error = scipy.sqrt((numerator_se/numerator_mean)**2 + (trace_se/trace_mean)**2 - (2*cov/(numerator_mean*trace_mean)) )*abs(mean)
    return mean, error

if __name__ == '__main__':

    data_files = sys.argv[1:]
    entropy, concurrence, trace = extract_data(data_files)

    concurrence_mean, concurrence_se = calculate_stats_ratio(concurrence, trace)
 
    entropy_mean = scipy.mean(entropy)
    entropy_se = scipy.stats.sem(entropy)
     
    if len(entropy) > 0: print "Average Von Neumann Entropy = ", entropy_mean, "    s.e. = ", entropy_se, "   beta loops = ", len(entropy)
    if len(concurrence) > 0: print "Average concurrence = ", concurrence_mean, "   s.e. = ", concurrence_se, "   beta loops = ", len(concurrence)

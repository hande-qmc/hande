#!/usr/bin/python

import sys
import re
import math
import scipy
import scipy.stats

def extract_data(data_files):
    unormalised_concurrence = []
    entropy = []
    trace = []
    entropy_regex = '^ # Unnormalised von Neumann entropy='
    concurrence_regex = '^ # Unnormalised concurrence='
    trace_regex = '^ # RDM trace=' 
    
    for data_file in data_files:

        f = open(data_file)

        for line in f:
            if re.match(entropy_regex, line):
                words = line.split()
                entropy.append(float(words[5]))
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
    return cov/(len(trace)*(len(trace)-1))

def calculate_stats_ratio(numerator, trace):
    numerator_mean = scipy.mean(numerator)
    trace_mean = scipy.mean(trace)
    numerator_se = scipy.stats.sem(numerator)
    trace_se = scipy.stats.sem(trace)
    mean = numerator_mean/trace_mean
    cov = calculate_covariance(numerator, numerator_mean, trace, trace_mean)
    error = scipy.sqrt((numerator_se/numerator_mean)**2 + (trace_se/trace_mean)**2 - (2*cov/(numerator_mean*trace_mean)) )*abs(mean)
    return mean, error

def calculate_entropy_stats(numerator, trace):
    numerator_mean = scipy.mean(numerator)
    trace_mean = scipy.mean(trace)
    numerator_se = scipy.stats.sem(numerator)
    trace_se = scipy.stats.sem(trace)
    mean = numerator_mean/trace_mean
    cov = calculate_covariance(numerator, numerator_mean, trace, trace_mean)
    error = scipy.sqrt((numerator_se*mean/numerator_mean)**2 + ((trace_se/trace_mean)**2)*((math.log(2)-mean)**2) - (2*cov/(numerator_mean*trace_mean))*(mean-math.log(2))*mean )
    mean = mean + math.log(trace_mean, 2)
    return mean, error

if __name__ == '__main__':

    data_files = sys.argv[1:]
    entropy, concurrence, trace = extract_data(data_files)

    if len(entropy) > 0: entropy_mean, entropy_se = calculate_entropy_stats(entropy, trace)
    if len(concurrence) > 0: concurrence_mean, concurrence_se = calculate_stats_ratio(concurrence, trace)
 
    if len(entropy) > 0: print "Average Von Neumann Entropy = ", entropy_mean, "    s.e. = ", entropy_se, "   beta loops = ", len(entropy)
    if len(concurrence) > 0: print "Average concurrence = ", concurrence_mean, "   s.e. = ", concurrence_se, "   beta loops = ", len(concurrence)

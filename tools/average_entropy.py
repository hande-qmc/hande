#!/usr/bin/python

import sys
import re
import scipy
import scipy.stats

def extract_data(data_files):
    concurrence = []
    entropy = []
    entropy_regex = '^ # Von-Neumann Entropy='
    concurrence_regex = '^ # Concurrence='
    for data_file in data_files:

        f = open(data_file)

        have_data = False
        for line in f:
            if re.match(entropy_regex, line):
                words = line.split()
                entropy_i = float(words[3])
                entropy.append(entropy_i)
            if re.match(concurrence_regex, line):
                words = line.split()
                concurrence_i = float(words[2])
                concurrence.append(concurrence_i)
        f.close()
    return entropy, concurrence



if __name__ == '__main__':

    data_files = sys.argv[1:]
    entropy, concurrence = extract_data(data_files)
    entropy_mean = scipy.mean(entropy)
    entropy_se = scipy.stats.sem(entropy)
    concurrence_mean = scipy.mean(concurrence)
    concurrence_se = scipy.stats.sem(concurrence)
    if len(entropy) > 0: print "Average Von Neumann Entropy = ", entropy_mean, "    s.e. = ", entropy_se, "   beta loops = ", len(entropy)
    if len(concurrence) > 0: print "Average concurrence = ", concurrence_mean, "   s.e. = ", concurrence_se, "   beta loops = ", len(concurrence) 

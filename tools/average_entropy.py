#!/usr/bin/python

import sys
import re
import scipy
import scipy.stats

def extract_data(data_files):
    data = []
    # timestep value from echoed input data
    entropy_regex = '^ # Von-Neumann Entropy='

    for data_file in data_files:

        f = open(data_file)

        have_data = False
        for line in f:
            if re.match(entropy_regex, line):
                words = line.split()
                entropy = float(words[3])
                data.append(entropy)
        f.close()
    return data



if __name__ == '__main__':

    data_files = sys.argv[1:]
    data = extract_data(data_files)
    mean = scipy.mean(data)
    se = scipy.stats.sem(data)
    print "Average Von Neumann Entropy = ", mean, "    s.e. = ", se, "   beta loops = ", len(data)

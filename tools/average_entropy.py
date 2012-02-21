#!/usr/bin/python

import sys
import re
import scipy
import scipy.stats

def extract_data(data_files):
    concurrence = []
    entropy = []
    entanglement = []
    entropy_regex = '^ # Von-Neumann Entropy='
    concurrence_regex = '^ # Concurrence='
    entanglement_regex = '^ # Entanglement of formation='
    for data_file in data_files:

        f = open(data_file)

        have_data = False
        for line in f:
            if re.match(entropy_regex, line):
                words = line.split()
                entropy.append(float(words[3]))
            elif re.match(concurrence_regex, line):
                words = line.split()
                concurrence.append(float(words[2]))
            elif re.match(entanglement_regex, line):
                words = line.split()
                entanglement.append(float(words[4]))
        f.close()
    return entropy, concurrence, entanglement



if __name__ == '__main__':

    data_files = sys.argv[1:]
    entropy, concurrence, entanglement = extract_data(data_files)
    entropy_mean = scipy.mean(entropy)
    entropy_se = scipy.stats.sem(entropy)
    concurrence_mean = scipy.mean(concurrence)
    concurrence_se = scipy.stats.sem(concurrence)
    entanglement_mean = scipy.mean(entanglement)
    entanglement_se = scipy.stats.sem(entanglement)
    if len(entropy) > 0: print "Average Von Neumann Entropy = ", entropy_mean, "    s.e. = ", entropy_se, "   beta loops = ", len(entropy)
    if len(concurrence) > 0: print "Average concurrence = ", concurrence_mean, "   s.e. = ", concurrence_se, "   beta loops = ", len(concurrence)
    if len(entanglement) > 0: print "Average Entanglement of Formation = ", entanglement_mean, "   s.e. = ", entanglement_se, "   beta loops = ", len(entanglement)  

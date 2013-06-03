#!/usr/bin/python

import math
import re
import sys

def extract_data(data_file):
    eigv = []
    regex = '^ Reduced density matrix eigenvalues:'
    
    for data_file in data_file:

        f = open(data_file)

        have_data = False
        for line in f:
            if not line.strip():
                have_data = False
            if have_data:
                if '----' not in line:
                    words = line.split()
                    eigv.append(float(words[1]))
            if re.match(regex, line):
                have_data = True
            
        f.close()
    return eigv

def calculate_vn_entropy(eigvs):
    vn_entropy = 0

    for eigv in eigvs:
        vn_entropy = vn_entropy - eigv*math.log(eigv,2)

    return vn_entropy

def calculate_renyi_2(eigvs):
    renyi_2 = 0

    for eigv in eigvs:
        renyi_2 = renyi_2 + eigv*eigv

    renyi_2 = -math.log(renyi_2,2)

    return renyi_2

if __name__ == '__main__':

    data_file = sys.argv[1:]

    eigv = extract_data(data_file)

    if len(eigv) > 0: 
        vn_entropy = calculate_vn_entropy(eigv)
        print "Von Neumann entropy = ", vn_entropy
        renyi_2 = calculate_renyi_2(eigv)
        print "Renyi 2 entropy = ", renyi_2
    else:
        print "No RDM eigenvalues were found."

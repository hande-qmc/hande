#!/usr/bin/python
'''Analysis of rdm_entropy_estimates.x output.

The estimate of the von Neumann and Renyi entropy measure comes directly from
the mean provided by rdm_entropy_estimates.x, the mean and standard error are
calculated by evaluating the entropy measure over the noisy RDMs produced by
rdm_entropy_estimates.x.'''
import scipy
import sys

def extract_data(data_file):

    vn_data = []
    renyi_data = []
    vn_estimate = 0
    renyi_estimate = 0

    f = open(data_file)

    vn_entropy = True
    for line in f:
        if 'Renyi' in line:
            vn_entropy = False
        if vn_entropy:
            if '#' in line:
                words = line.split()
                vn_estimate = float(words[-1])
            else:
                words = line.split()
                vn_data.append(float(words[0]))
        else:
            if '#' in line:
                words = line.split()
                renyi_estimate = float(words[-1])
            else:
                words = line.split()
                renyi_data.append(float(words[0]))

    return vn_data, vn_estimate, renyi_data, renyi_estimate

if __name__ == '__main__':

    if len(sys.argv) != 2:
        print('''%s: %s

No arguments supplied.

Usage:

%s FILE

where FILE contains the output from rdm_entropy_estimates.x.''' % (sys.argv[0], __doc__, sys.argv[0]))
    else:
        data_file = sys.argv[1]

        vn_data, vn_estimate, renyi_data, renyi_estimate = extract_data(data_file)

        if (len(vn_data)) > 0:
            vn_mean = scipy.mean(vn_data)
            vn_se = numpy.std(vn_data,ddof=1)
            print('vn_estimate: %s  mean: %s  se: %s' % (vn_estimate, vn_mean, vn_se))

        if (len(renyi_data)) > 0:
            renyi_mean = scipy.mean(renyi_data)
            renyi_se = numpy.std(renyi_data,ddof=1)
            print('renyi_estimate: %s  mean: %s  se: %s' % renyi_estimate, renyi_mean, renyi_se)

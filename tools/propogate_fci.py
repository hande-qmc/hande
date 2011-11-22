#!/usr/bin/python

import numpy
try:
    import matplotlib.pyplot as plt
    USE_MATPLOTLIB = True
except ImportError:
    USE_MATPLOTLIB = False
import sys

def extract_spectrum(fci_file):

    f = open(fci_file)

    have_data = False
    spectrum = []
    for line in f:
        if not line.strip():
            # empty line signals end of spectrum output
            have_data = False
        elif have_data:
            spectrum.append(float(line.split()[-1]))
        elif 'Total energy' in line:
            # start of spectrum output
            have_data = True

    f.close()

    return spectrum

def finite_temp_energy(beta, spectrum):

    s1 = 0.0
    s2 = 0.0
    for eigv in spectrum:
        e = numpy.exp(-beta*eigv)
        s1 += e*eigv
        s2 += e
    return s1/s2

def propogate_spectrum(beta_min, beta_max, nbeta, spectrum):

    beta = numpy.arange(beta_min, beta_max, float(beta_max-beta_min)/(nbeta-1))

    energies = [finite_temp_energy(b, spectrum) for b in beta]

    print '#     beta             E(beta)'
    for i in range(len(beta)):
        print '%16.8f %16.8f' % (beta[i], energies[i])

    if USE_MATPLOTLIB:
        plt.plot(beta, energies)
        plt.show()

if __name__ == '__main__':

    if len(sys.argv) != 5:
        print 'Usage:', sys.argv[0], 'fci_file beta_min beta_max nbeta'
        print r'Evaluate E(\beta) from the output of an FCI calculation contained in fci_file produced by hubbard.x, between beta_min and beta_max in steps of (beta_max-beta_min)/(nbeta-1).'
        sys.exit(1)

    (fci_file, beta_min, beta_max, nbeta) = sys.argv[1:]
    beta_min = float(beta_min)
    beta_max = float(beta_max)
    nbeta = float(nbeta)

    spectrum = extract_spectrum(fci_file)
    propogate_spectrum(beta_min, beta_max, nbeta, spectrum)

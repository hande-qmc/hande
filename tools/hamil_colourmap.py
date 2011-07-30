#!/usr/bin/env python2

import numpy
import pylab
import sys

def colormap(filename, N):

    # parallel output or serial?
    f = open(filename)
    parallel_output = 'hamil' in f.readline()
    f.seek(0)

    HMat = pylab.zeros((N,N))

    if parallel_output:
        for line in f:
            i, j, hamil = line.split()[1:]
            i = numpy.float(i[:-1]) - 1
            j = numpy.float(j[:-2]) - 1
            if i != j:
                hamil = numpy.float(hamil)
                HMat[i][j] = hamil
                HMat[j][i] = hamil
    else:
        for line in f:
            i, j, hamil = line.split()
            i = numpy.float(i) - 1
            j = numpy.float(j) - 1
            if i != j:
                hamil = numpy.float(hamil)
                HMat[i][j] = hamil
                HMat[j][i] = hamil

    f.close()

    pylab.imsave('%s.png' % filename, HMat, cmap=pylab.get_cmap('RdBu'))

if __name__ == '__main__':

    if len(sys.argv) == 3: 
        colormap(sys.argv[1], int(sys.argv[2]))
    else:
        print 'Usage: hamil_colormap.py hamil_file N'
        sys.exit(1)

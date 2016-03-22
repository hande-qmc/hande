#!/usr/bin/env python
'''Plot excitation distribution from dmqmc output.'''

import os
import sys
import matplotlib.pyplot as pl
import argparse
try:
    import pyhande as ph
except ImportError:
    _script_dir = os.path.dirname(os.path.abspath(__file__))
    sys.path.append(os.path.join(_script_dir, '../pyhande'))
    import pyhande as ph


def plot_excit_dist(filename, plotfile, calc, max_excit):
    ''' Plot excitation distribution.

Paramters
---------
filename : string
    file to plot from.
plotfile : string or None
    name of file to plot to.
calc : int
   calculation number to plot.
max_excit : int or None
    maximum excitation level to plot to.

'''

    data = ph.extract.extract_data(filename)

    (m, d) = data[calc]
    if not max_excit:
        max_excit = int(m['system']['max_number_excitations'])
    for e in range(0, max_excit):
        pl.plot(d['iterations']*m['qmc']['tau'], d['Excit. level %s'%e],
                label=r'$n_{\mathrm{ex}} = %s$'%e)
    pl.legend(numpoints=1, loc='best')
    if m['ipdmqmc']['propagate_to_beta']:
        pl.xlabel(r'$\tau$')
    else:
        pl.xlabel(r'$\beta$')
    pl.ylabel('Weight')
    if plotfile:
        pl.savefig(plotfile+'.png', fmt='png')
    else:
        pl.show()


def parse_args(args):

    parser = argparse.ArgumentParser(description='Plot the excitation'
                                     ' distribution of DMQMC calulation.')
    parser.add_argument('-p', '--plotfile', default=None, help='File to save '
                        'the graphs to. The graphs are shown interactively by '
                        'default.')
    parser.add_argument('calc', type=int, help='Calculation number to plot. '
                        'C indexed.', default=0)
    parser.add_argument('-m', '--max-excit', action='store', dest='max_excit',
                        type=int, help='Plot up to maximum excitation '
                        'distribution. Plot all by default.', default=None)
    parser.add_argument('file', help='File to plot.')

    opts = parser.parse_args(args)

    return (opts.file, opts.plotfile, opts.calc, opts.max_excit)


if __name__ == '__main__':

    (datafile, plotfile, calc, max_excit) = parse_args(sys.argv[1:])
    plot_excit_dist(datafile, plotfile, calc, max_excit)

#!/usr/bin/env python
'''Extract the momentum distribution from an analysed DMQMC simulation.'''

import pandas as pd
import numpy as np
import sys
import argparse


def parse_args(args):
    '''Parse command-line arguments.

Parameters
----------
args : list of strings
    command-line arguments.

Returns
-------
options : :class:`ArgumentParser`
    Options read in from command line.
'''
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument('-b', '--beta-val', action='store', required=True,
                        type=float, dest='beta', help='Inverse temperature '
                        'to extract the mometnum distribution at.')
    parser.add_argument('filename', nargs='+', help='Analysed DMQMC data. '
                        'i.e. the output of running finite_temp_analysis.py '
                        'on a HANDE calculation.')

    options = parser.parse_args(args)

    if not options.filename:
        parser.print_help()
        sys.exit(1)

    return options


def main(args):

    options = parse_args(args)
    data = pd.read_csv(options.filename[0],
                       sep=r'\s+').groupby('Beta').get_group(options.beta)

    mom = [c for c in data.columns.values if 'n_' in c and '_error' not in c]
    mome = [c for c in data.columns.values if 'n_' in c and '_error' in c]

    vals = [float(c.split('_')[1]) for c in mom]

    n_k = (data[mom].transpose()).values
    n_k_error = (data[mome].transpose()).values
    n_k_error[np.isnan(n_k_error)] = 0
    frame = pd.DataFrame({'Beta': options.beta, 'k': vals, 'n_k': n_k.ravel(),
                          'n_k_error': n_k_error.ravel()})

    print (frame.to_string(index=False))


if __name__ == '__main__':

    main(sys.argv[1:])

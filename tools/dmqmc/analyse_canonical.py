#!/usr/bin/env python

import pandas as pd
import os
import pkgutil
import sys
import warnings
import argparse

if not pkgutil.find_loader('pyhande'):
    _script_dir = os.path.dirname(os.path.abspath(__file__))
    sys.path.append(os.path.join(_script_dir, '../pyhande'))
import pyhande


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
    parser.add_argument('-s', '--sim', action='store_true', default=False,
                        dest='multi_sim', help='Do not average over multiple '
                        'simulations in the same or from multiple data files.')
    parser.add_argument('filename', nargs='+', help='HANDE output.')

    options = parser.parse_args(args)

    if not options.filename:
        parser.print_help()
        sys.exit(1)

    return options


def main(args):
    ''' Analyse the output from a canonical estimates calculation.

Parameters
----------
filename : list of strings
    files to be analysed.
'''

    args = parse_args(args)

    hande_out = pyhande.extract.extract_data_sets(args.filename)

    (metadata, data) = ([], [])
    for (md, df) in hande_out:
        # Handle old output with incorrect title...
        if md['calc_type'] == 'Canonical energy' or md['calc_type'] == 'RNG':
            metadata.append(md)
            data.append(df)
    if data and not args.multi_sim:
        data = pd.concat(data)

    # Sanity check: are all the calculations from the same calculation?
    # Only check for new metadata format...
    if 'beta' in metadata[0]:
        beta = metadata[0]['beta']
        for md in metadata[1:]:
            if 'beta' in md and md['beta'] != beta:
                warnings.warn('Beta values in input files not consistent.')

    if args.multi_sim:
        results = pd.DataFrame([pyhande.canonical.estimates(m, d) for (m, d)
                                in zip(metadata, data)])
    else:
        results = pd.DataFrame(pyhande.canonical.estimates(metadata[0], data)).T

    try:
        float_fmt = '{0:-#.8e}'.format
        float_fmt(1.0)
    except ValueError:
        # GAH.  Alternate formatting only added to format function after
        # python 2.6..
        float_fmt = '{0:-.8e}'.format

    # Work around bug in to_string alignment with index=False
    lines = results.to_string(float_format=float_fmt).split('\n')
    index_wdith = max([len(str(s)) for s in results.index]) + 1
    print('\n'.join(l[index_wdith:] for l in lines))

if __name__ == '__main__':

    main(sys.argv[1:])

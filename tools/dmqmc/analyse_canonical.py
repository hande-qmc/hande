#!/usr/bin/env python

import pandas as pd
import os
import sys
import warnings
script_dir = os.path.dirname(__file__)
sys.path.extend([os.path.join(script_dir, '../')])
import pyhande


def main(filename):
    ''' Analyse the output from a canonical kinetic energy calculation.

Parameters
----------
filename : list of strings
    files to be analysed.
'''

    if len(filename) < 1:
        print("Usage: ./analyse_canonical.py files")
        sys.exit()

    hande_out = pyhande.extract.extract_data_sets(filename)

    (metadata, data) = ([], [])
    for (md, df) in hande_out:
        # Handle old output with incorrect title...
        if md['calc_type'] == 'Canonical energy' or md['calc_type'] == 'RNG':
            metadata.append(md)
            data.append(df)
    if data:
        data = pd.concat(data)

    # Sanity check: are all the calculations from the same calculation?
    # Only check for new metadata format...
    if 'beta' in metadata[0]:
        beta = metadata[0]['beta']
        for md in metadata[1:]:
            if 'beta' in md and md['beta'] != beta:
                warnings.warn('Beta values in input files not consistent.')

    results = pyhande.canonical.estimates(metadata[0], data)

    try:
        float_fmt = '{0:-#.8e}'.format
        float_fmt(1.0)
    except ValueError:
        # GAH.  Alternate formatting only added to format function after
        # python 2.6..
        float_fmt = '{0:-.8e}'.format

    print(results.to_string(index=False, float_format=float_fmt))


if __name__ == '__main__':

    main(sys.argv[1:])

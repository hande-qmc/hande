#!/usr/bin/env python

import pandas as pd
import os
import sys
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

    (metadata, data) = pyhande.extract.extract_data_sets(filename)

    results = pyhande.canonical.estimates(metadata, data)

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

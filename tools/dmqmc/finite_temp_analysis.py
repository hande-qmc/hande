#!/usr/bin/env python
'''finite_temp_analysis.py [options] file_1 file_2 ... file_N

Analyse the output of a HANDE DMQMC calculation by averaging
temperature-dependent data across beta loops.'''

import pandas as pd
import os
import pkgutil
import sys

_script_dir = os.path.dirname(os.path.abspath(__file__))
if not pkgutil.find_loader('pyhande'):
    sys.path.append(os.path.join(_script_dir, '../pyhande'))

import pyhande


def main(args):
    '''Run data analysis on finite-temperature HANDE output.

Parameters
----------
args : list of strings
    command-line arguments.

Returns
-------
None.
'''

    (files, options) = pyhande.dmqmc.parse_args(args)
    hande_out = pyhande.extract.extract_data_sets(files)

    # Finally, output the results!
    results = pyhande.dmqmc.analyse_data(hande_out, options)[1]

    # For anal-retentiveness, print the energy first after beta and then all
    # columns in alphabetical order.
    columns = sorted(results.columns.values)
    columns.insert(1, columns.pop(columns.index('Tr[Hp]/Tr[p]')))
    columns.insert(2, columns.pop(columns.index('Tr[Hp]/Tr[p]_error')))

    print(results.to_string(index=False, columns=columns))


if __name__ == '__main__':

    main(sys.argv[1:])

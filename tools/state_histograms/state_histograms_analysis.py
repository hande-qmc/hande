#!/usr/bin/env python
'''state_histograms_analysis.py [options] file_1 file_2 ... file_N

Analyse the output of a HANDE DMQMC or FCIQMC state histogram calculation by
averaging the provided state histograms by bin/excitation indexes.'''

import os
import sys
import pkgutil
import argparse
import pandas as pd

_script_dir = os.path.dirname(os.path.abspath(__file__))
if not pkgutil.find_loader('pyhande'):
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
filenames : list of strings
    list of HANDE DMQMC or FCIQMC state histogram output files.
options : :class:`ArgumentParser`
    Options read in from command line.
'''

    parser = argparse.ArgumentParser(usage = __doc__)
    parser.add_argument('-o', '--output', default='csv', choices=['csv', 'txt'],
                        help='Format for data table.  Default: %(default)s.')
    parser.add_argument('-fltfmt', '--float-format', action='store', default=None,
                        type=str, dest='fltfmt', help='Format the values from '
                        'the resulting analysis before reporting, only applies '
                        'to "txt" format. I.E., %%6.4f, %%12.8E, ..., etc.')
    parser.add_argument('-trex1', '--trunc-exlevel1', action='store', default=None,
                        type=int, dest='trunc_ex1', help='Truncate the analysis to '
                        'the first excitation levels less than or equal to '
                        'the provided excitation level.')
    parser.add_argument('-trex2', '--trunc-exlevel2', action='store', default=None,
                        type=int, dest='trunc_ex2', help='Truncate the analysis to '
                        'the second excitation levels less than or equal to '
                        'the provided excitation level.')
    parser.add_argument('-tol', '--tolerance-number', action='store', default=None,
                        type=int, dest='csv_tol', help='Set the tolerance of the '
                        'reported data from analysis. For example, 12 would round '
                        'the final data report (after analysis) to the 12th decimal '
                        'place before reporting. Currently only for the "csv" format.')
    parser.add_argument('filenames', nargs='+', help='HANDE files to analyse.')

    options = parser.parse_args(args)

    if not options.filenames:
        parser.print_help()
        sys.exit(1)

    return (options.filenames, options)

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

    (files, options) = parse_args(args)

    results = pyhande.state_histograms.analyse_state_histograms(files,
                                                            options.trunc_ex1,
                                                            options.trunc_ex2)

    if options.output is 'csv':
        if options.csv_tol: results = results.round(options.csv_tol)
        print(results.to_csv(index=False))
    elif options.output is 'txt' and options.fltfmt is not None:
        print(results.to_string(index=False, float_format=options.fltfmt))
    else:
        print(results.to_string(index=False))

if __name__ == '__main__':
    main(sys.argv[1:])

#!/usr/bin/env python
'''finite_temp_analysis.py [options] file_1 file_2 ... file_N

Analyse the output of a HANDE DMQMC calculation by averaging
temperature-dependent data across beta loops.'''

import pandas as pd
import os
import pkgutil
import sys
import argparse

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
    list of HANDE DMQMC output files.
options : :class:`ArgumentParser`
    Options read in from command line.
'''

    parser = argparse.ArgumentParser(usage = __doc__)
    parser.add_argument('-s', '--with-shift', action='store_true', dest='with_shift',
                      default=False, help='Output the averaged shift profile and '
                      'the standard deviation of these profiles across beta loops.')
    parser.add_argument('-t', '--with-trace', action='store_true', dest='with_trace',
                      default=False, help='Output the averaged traces and the '
                      'standard deviation of these traces, for all replicas present.')
    parser.add_argument('-b', '--with-spline', action='store_true', dest='with_spline',
                      default=False, help='Output a B-spline fit for each of '
                      ' estimates calculated')
    parser.add_argument('-c', '--calc-number', action='store', default=None, type=int,
                      dest='calc_number', help='Calculation number to analyse. '
                      'Note any simulation using find_weights option should not '
                      'be included.')
    parser.add_argument('-f', '--with-free-energy', action='store_true',
                      dest='with_free_energy', default=False,
                      help='Calculate Free energy')
    parser.add_argument('filenames', nargs='+', help='HANDE files to analyse.')

    options = parser.parse_args(args)

    if not options.filenames:
        parser.print_help()
        sys.exit(1)

    return (options.filenames, options)


def remove_observable(columns, label):
    '''Remove observable from list.

Parameters
----------
columns : list
    List of column names.
label : string
    regex to search for removal.
'''

    return [c for c in columns if label not in c]


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
    hande_out = pyhande.extract.extract_data_sets(files)

    # Finally, output the results!
    results = pyhande.dmqmc.analyse_data(hande_out, options.with_shift,
                                         options.with_free_energy,
                                         options.with_spline,
                                         options.with_trace,
                                         options.calc_number)[1]

    # We want to sort momentum distribution array naturally so extract them here
    # before sorting the rest of the column names by their label.
    columns = sorted(remove_observable(results.columns.values, 'n_'))
    # For anal-retentiveness, print the energy first after beta and then all
    # columns in alphabetical order.
    columns = sorted([c for c in results.columns.values if ('n_' not in c) and ('S_'
                      not in c) and ('Suu_' not in c) and ('Sud_' not in c)])

    columns.insert(1, columns.pop(columns.index('Tr[Hp]/Tr[p]')))
    columns.insert(2, columns.pop(columns.index('Tr[Hp]/Tr[p]_error')))
    momentum_dist = pyhande.dmqmc.sort_momentum([c for c in
                                                 results.columns.values
                                                 if 'n_' in c])
    structure_factor = pyhande.dmqmc.sort_momentum([c for c in
                                                    results.columns.values
                                                    if ('S_' in c) or
                                                    ('Suu_' in c) or
                                                    ('Sud_' in c)])
    print(results.to_string(index=False, columns=columns+momentum_dist+structure_factor))


if __name__ == '__main__':

    main(sys.argv[1:])

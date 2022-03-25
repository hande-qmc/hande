#!/usr/bin/env python
''' state_histograms_analysis.py [options] file_1 file_2 ... file_N

Analyse the output of a HANDE DMQMC or FCIQMC state histogram calculation by
averaging the provided state histograms by bin/excitation indexes.
'''

import os
import sys
import pkgutil
import argparse

try:
    import pyhande
except ModuleNotFoundError:
    _script_dir = os.path.dirname(os.path.abspath(__file__))
    if not pkgutil.find_loader('pyhande'):
        sys.path.append(os.path.join(_script_dir, '../pyhande'))
    import pyhande


def parse_arguments(arguments):
    ''' Parse command-line arguments.

    Parameters
    ----------
    arguments : list of strings
        command-line arguments.

    Returns
    -------
    filenames : list of strings
        list of HANDE DMQMC or FCIQMC state histogram output files.
    options : :class:`ArgumentParser`
        Options read in from command line.
    '''

    parser = argparse.ArgumentParser(usage=__doc__)

    parser.add_argument('-o', '--output', default='csv',
                        choices=['csv', 'txt'], help='Format for data table. '
                        '  Default: %(default)s.')
    parser.add_argument('-ff', '--float-format', action='store',
                        default=None, type=str, dest='float_format',
                        help='Format the values from the resulting analysis '
                        'before reporting, only used for the "txt" format. '
                        'I.E., %%6.4f, %%12.8E, ..., etc.')
    parser.add_argument('-tol', '--tolerance-number', action='store',
                        default=None, type=int, dest='csv_tol', help='Set the '
                        'tolerance of the reported data from analysis. For '
                        'example, 12 would round the final data report '
                        '(after analysis) to the 12th decimal place before '
                        'reporting. Only used for the "csv" format.')
    parser.add_argument('filenames', nargs='+', help='HANDE files to analyse.')
    parser.parse_args(args=None if arguments else ['--help'])

    options = parser.parse_args(arguments)

    return (options.filenames, options)


def main(arguments):
    ''' Run data analysis on state histogram files from DMQMC or FCIQMC.

    Parameters
    ----------
    arguments : list of strings
        command-line arguments.

    Returns
    -------
    None.
    '''

    (files, options) = parse_arguments(arguments)

    results = pyhande.state_histograms.analyse_state_histograms(files)

    if options.output == 'csv':
        if options.csv_tol:
            results = results.round(options.csv_tol)
        print(results.to_csv(index=False))
    elif options.output == 'txt':
        if options.float_format is not None:
            print(results.to_string(index=False,
                                    float_format=options.float_format))
        else:
            print(results.to_string(index=False))


if __name__ == '__main__':
    main(sys.argv[1:])

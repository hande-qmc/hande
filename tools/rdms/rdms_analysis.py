#!/usr/bin/env python
''' rdms_analysis.py [options] file_1 file_2 ... file_N

Perform reblocking analysis and/or average the RDM file written from
an FCIQMC or DMQMC simulation accumulating statistics for the 1- and 2-RDM.
'''

import os
import sys
import pkgutil
import argparse

try:
    from pyhande.extract import extract_data
    from pyhande.lazy import find_starting_iteration
    from pyhande.analysis import projected_energy
    from pyblock.pd_utils import reblock
except ModuleNotFoundError:
    _script_dir = os.path.dirname(os.path.abspath(__file__))
    if not pkgutil.find_loader('pyhande'):
        sys.path.append(os.path.join(_script_dir, '../pyhande'))
    from pyhande.extract import extract_data
    from pyhande.lazy import find_starting_iteration
    from pyhande.analysis import projected_energy
    if not pkgutil.find_loader('pyblock'):
        sys.path.append(os.path.join(_script_dir, '../pyblock'))
    from pyblock.pd_utils import reblock

def parse_arguments(arguments):
    ''' Parse command-line arguments.

    Parameters
    ----------
    arguments : list of strings
        command-line arguments.

    Returns
    -------
    list_of_filenames : list of list of strings
        list of lists with HANDE DMQMC or FCIQMC state histogram output files.
    options : :class:`ArgumentParser`
        Options read in from command line.
    '''

    parser = argparse.ArgumentParser(usage=__doc__)

    parser.add_argument('-o', '--output', default='csv',
                        choices=['csv', 'txt'], help='Format for data table. '
                        '  Default: %(default)s.')
    parser.add_argument('-of', '--output_file', action='store', default=None,
                        type=str, dest='output_file', help='Define a file '
                        'name to store the results of analysis to.')
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

    return options.filenames, options


def reblock_rdm(filename):
    ''' Perform reblocking analysis on the RDM numerator and denominator.
    Behaves very similar to the traditional analysis of the projected
    energy from FCIQMC, see related routines for more information.

    Parameters
    ----------
    filename : str
        A filename containing the RDM quantities to perform reblocking
        analysis on.

    Returns
    -------
    '''
    meta_data, data = extract_data(filename)[0]

    calc_end = data['iterations'].iloc[-1]

    tmp_data = data.copy()

    tmp_data[['sum H_0j N_j', 'N_0']] = data[['RDM energy num.', 'RDM trace']]

    calc_start = find_starting_iteration(tmp_data, meta_data, end=calc_end)

    del tmp_data

    data = data.set_index('iterations').loc[calc_start+1:calc_end]
    data = data.reset_index(drop=False)

    data_length, reblock_data, covariance = reblock(data)

    erdm = projected_energy(reblock_data, covariance, data_length,
                            sum_key='RDM energy num.', ref_key='RDM trace',
                            col_name='RDM energy')

    print(reblock_data[['RDM energy num.', 'RDM trace']])
    print(erdm)


def main(arguments):
    ''' Run data analysis to produce an estimate of the energy from
    the RDM numerator and RDM denominator.

    Parameters
    ----------
    arguments : list of strings
        command-line arguments.

    Returns
    -------
    None.
    '''

    filenames, options = parse_arguments(arguments)

    for filename in filenames:
        reblock_rdm(filename)


if __name__ == '__main__':
    main(sys.argv[1:])

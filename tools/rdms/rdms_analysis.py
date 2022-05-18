#!/usr/bin/env python
''' rdms_analysis.py [options] file_1 file_2 ... file_N

Perform reblocking analysis on the RDM energy num./trace from an
FCIQMC calculation to generate an estimate of the correlation energy.
'''

import os
import sys
import pkgutil
import argparse

try:
    from pyhande.extract import extract_data
except ModuleNotFoundError:
    _script_dir = os.path.dirname(os.path.abspath(__file__))
    if not pkgutil.find_loader('pyhande'):
        sys.path.append(os.path.join(_script_dir, '../pyhande'))
    if not pkgutil.find_loader('pyblock'):
        sys.path.append(os.path.join(_script_dir, '../pyblock'))
    from pyhande.extract import extract_data

from numpy import nan
from pandas import DataFrame, concat
from pyhande.lazy import find_starting_iteration
from pyhande.analysis import projected_energy, qmc_summary
from pyblock.pd_utils import reblock
from pyblock.plot import plot_reblocking


def parse_arguments(arguments):
    ''' Parse command-line arguments.

    Parameters
    ----------
    arguments : list of strings
        command-line arguments.

    Returns
    -------
    filenames : list of list of strings
        list of lists with HANDE FCIQMC RDM data.
    options : :class:`ArgumentParser`
        Options read in from command line.
    '''

    parser = argparse.ArgumentParser(usage=__doc__)

    parser.add_argument('-o', '--output', default='txt',
                        choices=['csv', 'txt'], help='Format for data table. '
                        '  Default: %(default)s.')
    parser.add_argument('-ff', '--float-format', action='store',
                        default=None, type=str, dest='float_format',
                        help='Format the values from the resulting analysis '
                        'before reporting, only used for the "txt" format. '
                        'I.E., %%6.4f, %%12.8E, ..., etc.')
    parser.add_argument('-plt', '--plot', default=None, type=str,
                        action='store', help='Provide a file name for the '
                        'reblocking analysis to be plotted and saved to.')
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
    opt : :class:`pandas.DataFrame`
        The estimates for the RDM energy numerator, RDM trace and RDM Energy
        based on the reblocking analysis.
    rdm : :class:`pandas.DataFrame`
        Reblockign data for RDM energy num. and RDM trace in a csv
        friendlier format.
    block_info : :class:`pandas.DataFrame`
        Reblocking data for RDM energy num. and RDM trace, used when
        generating a plot of the reblocking analysis.
    '''
    rdm_keys = ['RDM energy num.', 'RDM trace']
    prj_keys = [r'\sum H_0j N_j', 'N_0']
    rep_keys = [
            'iterations', '# H psips', 'Shift',
            '# states', '# spawn_events', 'R_spawn', 'time',
        ]

    meta_data, data = extract_data(filename)[0]

    calc_start = 0
    calc_end = data['iterations'].iloc[-1]

    # Find the starting iteration by considering both replicas
    # and taking the conservative estimate of the two.
    # Note we replace the typical projected energy terms with the RDM
    # terms so they carry the same influence in the starting iteration
    # as the projected energy terms normally do.
    for ireplica in [1, 2]:
        tmp_keys = rep_keys.copy()
        tmp_keys[1] = tmp_keys[1] + f'_{ireplica}'
        tmp_keys[2] = tmp_keys[2] + f'_{ireplica}'

        tmp_data = DataFrame()
        tmp_data[rep_keys] = data.loc[:, tmp_keys].copy()
        tmp_data[prj_keys] = data.loc[:, rdm_keys].copy()
        tmp_data = tmp_data.loc[data['RDM trace'] != 0.0]

        tmp_start = find_starting_iteration(tmp_data, meta_data, end=calc_end)

        if tmp_start > calc_start:
            calc_start = tmp_start

        del tmp_data, tmp_keys, tmp_start

    data = data.set_index('iterations').loc[calc_start+1:calc_end]
    data = data.reset_index(drop=False)

    data_length, reblock_data, covariance = reblock(data)

    erdm = projected_energy(reblock_data, covariance, data_length,
                            sum_key='RDM energy num.', ref_key='RDM trace',
                            col_name='RDM energy')

    opt = []
    rdm = reblock_data.loc[:, rdm_keys].copy()
    block_info = reblock_data.loc[:, rdm_keys].copy()
    rdm_opt, rdm_no_opt = qmc_summary(rdm, keys=rdm_keys)

    opt.append(rdm_opt)

    erdm_opt, erdm_no_opt = qmc_summary(erdm, keys=['RDM energy'])
    erdm_opt['standard error error'] = [nan]

    opt.append(erdm_opt)
    opt = concat(opt)

    # Convert the "<---   "/"" notation to 1/0 respectively for
    # a csv read in friendly format.
    nopt = (rdm.loc[:, ('RDM energy num.', 'optimal block')] != '')
    topt = (rdm.loc[:, ('RDM trace', 'optimal block')] != '')

    rdm.loc[:, ('RDM energy num.', 'optimal block')] = nopt.astype(int)
    rdm.loc[:, ('RDM trace', 'optimal block')] = topt.astype(int)

    return opt, rdm, block_info


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
        opt, rdm, block_info = reblock_rdm(filename)

        if options.output == 'csv':
            print('Reblocked RDM estimates:')
            print(opt.to_csv(), '\n')
            print('Raw RDM reblock data:')
            print(rdm.reset_index(level=[0]).to_csv(index=False))
        else:
            print('Reblocked RDM estimates:')
            print(opt.to_string(float_format=options.float_format), '\n')
            print('Raw RDM reblock data:')
            print(rdm.to_string(float_format=options.float_format))

        if options.plot is not None:
            if options.plot[-4:] != '.pdf':
                options.plot += '.pdf'

            plot_reblocking(block_info, plotfile=options.plot, plotshow=False)


if __name__ == '__main__':
    main(sys.argv[1:])

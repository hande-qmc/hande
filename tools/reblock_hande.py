#!/usr/bin/env python
'''reblock_hande.py [options] file_1 file_2 ... file_N

Run a reblocking analysis on HANDE QMC output files.  CCMC and FCIQMC
calculations only are supported.'''

import pandas as pd
from os import path
import sys
sys.path.append(path.join(path.abspath(path.dirname(sys.argv[0])), 'pyblock'))
import pyblock
import pyhande

# Still supporting 2.6.  *sigh*
import optparse

def run_hande_blocking(files, start_iteration, reblock_plot=None, verbose=1):
    '''Run a reblocking analysis on HANDE output and print to STDOUT.

See :func:`pyblock.pd_utils.reblock` and :func:`pyblock.blocking.reblock` for
details on the reblocking procedure.

Parameters
----------
files : list of strings
    names of files containing HANDE QMC calculation output.
start_iteration : int
    QMC iteration from which statistics are gathered.
reblock_plot : string
    Filename to which the reblocking convergence plot (standard error vs reblock
    iteration) is saved.  The plot is not created if None and shown
    interactively if '-'.
verbose : int
    Level of verbosity.

    <0: print nothing
    0: print only the recommended statistics from the optimal block length.
    1: print blocking analysis and recommended statistics.
    2: print calculation metadata, blocking analysis and recommended statistics.

Returns
-------
metadata : :class:`pandas.DataFrame`
    Metadata extracted from the calculation output file.
data : :class:`pandas.DataFrame`
    QMC data extracted from the calculation output file.
opt_block: :class:`pandas.DataFrame`
    Recommended statistics based upon the estimated 'optimal' block size
    as suggested by Wolff and Lee et al. (see
    :func:`pyblock.blocking.find_optimal_block`).
'''

    try:
        float_fmt = '{0:-#.8e}'.format
        float_fmt(1.0)
    except ValueError:
        # GAH.  Alternate formatting only added to format function after
        # python 2.6..
        float_fmt = '{0:-.8e}'.format

    info = pyhande.lazy.std_analysis(files, start_iteration)

    if verbose >= 2:

        col_name = info.metadata.columns.name
        for (calc_name, calc) in info.metadata.iteritems():
            calc_local = calc.copy()
            # problems with pop on pandas 0.13?  It seems to return a list...
            calc_input = calc_local['input']
            calc_local = calc_local.drop('input')
            # Add the calc index to the series and make it come first.
            calc_local[col_name] = calc_name
            indx = calc_local.index.copy()
            indx = indx.delete(indx.get_loc(col_name)).insert(0, col_name)
            calc_local = calc_local.reindex(indx)
            print(calc_local.to_string(na_rep='n/a'))
            if verbose >= 3:
                print('\nFull input options:\n\n%s' % ('\n'.join(calc_input)))
            print('')
    if verbose >= 1:
        print(info.reblock.to_string(float_format=float_fmt, line_width=80))
        print('')
    if not info.opt_block.empty and verbose >= 0:
        print('Recommended statistics from optimal block size:')
        print('')
        print(info.opt_block.to_string(float_format=float_fmt, na_rep='n/a'))
    if info.no_opt_block and verbose >= 0:
        print('WARNING: could not find optimal block size.')
        print('Insufficient statistics collected for the following variables: '
              '%s.' % (', '.join(info.no_opt_block)))
    if reblock_plot:
        pyblock.pd_utils.plot_reblocking(reblock, reblock_plot)

    return info

def parse_args(args):
    '''Parse command-line arguments.

Parameters
----------
args : list of strings
    command-line arguments.

Returns
-------
(filenames, start_iteration, reblock_plot)

where

filenames : list of strings
    list of QMC output files
start_iteration : int
    iteration number from which statistics should be gathered.
reblock_plot : string
    filename for the reblock convergence plot output.
'''

    parser = optparse.OptionParser(usage = __doc__)
    parser.add_option('-p', '--plot', default=None, dest='plotfile',
                      help='Filename to which the reblocking convergence plot '
                      'is saved.  Use \'-\' to show plot interactively.  '
                      'Default: off.')
    parser.add_option('-q', '--quiet', dest='verbose', action='store_const',
                      const=0, default=1,
                      help='Output only the final summary table.  '
                      'Overrides --verbose.')
    parser.add_option('-s', '--start', type='int', dest='start_iteration',
                      default=0, help='Iteration number from which to gather '
                           'statistics.  Default: %default.')
    parser.add_option('-v', '--verbose', dest='verbose', action='count',
                      default=1, help='Increase verbosity of the output.  Can '
                      'be specified multiple times.')

    (options, filenames) = parser.parse_args(args)

    if not filenames:
        parser.print_help()
        sys.exit(1)

    return (filenames, options.start_iteration, options.plotfile,
            options.verbose)

def main(args):
    '''Run reblocking and data analysis on HANDE output.

Parameters
----------
args : list of strings
    command-line arguments.

Returns
-------
None.
'''

    (files, start_iteration, reblock_plot, verbose) = parse_args(args)
    run_hande_blocking(files, start_iteration, reblock_plot, verbose)

if __name__ == '__main__':

    main(sys.argv[1:])

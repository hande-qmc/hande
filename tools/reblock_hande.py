#!/usr/bin/python
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
opt_data: :class:`pandas.DataFrame`
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
    (metadata, data) = pyhande.extract.extract_data_sets(files)
    if verbose >= 2:
        print(metadata.to_string())
        print('')

    # Reblock over desired window.
    indx = data['iterations'] >= start_iteration
    mc_data =  data.ix[indx, ['Instant shift', '\sum H_0j Nj', '# D0']]
    (data_length, reblock, covariance) = pyblock.pd_utils.reblock(mc_data)

    # Calculate projected energy.
    proje_sum = reblock.ix[:, '\sum H_0j Nj']
    ref_pop = reblock.ix[:, '# D0']
    proje_ref_cov = covariance.xs('# D0', level=1)['\sum H_0j Nj']
    proje = pyblock.error.ratio(proje_sum, ref_pop, proje_ref_cov, data_length)

    if verbose >= 1:
        print(reblock.to_string(float_format=float_fmt, line_width=80))
        print('')

    # Data summary: suggested data to use from reblocking analysis.
    opt_data = []
    no_opt = []
    for col in ('Instant shift', '\sum H_0j Nj', '# D0'):
        summary = pyblock.pd_utils.reblock_summary(reblock.ix[:, col])
        if summary.empty:
            no_opt.append(col)
        else:
            summary.index = [col]
        opt_data.append(summary)
    summary = pyblock.pd_utils.reblock_summary(proje)
    if summary.empty:
        no_opt.append('Proj. Energy')
    else:
        summary.index = ['Proj. Energy']
    opt_data.append(summary)
    opt_data = pd.concat(opt_data)
    if not opt_data.empty and verbose >= 0:
        print('Recommended statistics from optimal block size:')
        print('')
        print(opt_data.to_string(float_format=float_fmt))
    if no_opt and verbose >= 0:
        print('WARNING: could not find optimal block size.')
        print('Insufficient statistics collected for the following variables: '
              '%s.' % (', '.join(no_opt)))

    if reblock_plot:
        pyblock.pd_utils.plot_reblocking(reblock, reblock_plot)

    return (metadata, data, opt_data)

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
                      const=0, default=1, help='')
    parser.add_option('-s', '--start', type='int', dest='start_iteration',
                      default=0, help='Iteration number from which to gather '
                           'statistics.  Default: %default.')
    parser.add_option('-v', '--verbose', dest='verbose', action='store_const',
                      const=2, help='')

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

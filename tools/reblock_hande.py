#!/usr/bin/env python
'''Run a reblocking analysis on HANDE QMC output files.  Files may be compressed
with either gzip, bzip2 or xz (python 3 only).  CCMC and FCIQMC calculations
only are supported; other calculations have specific analysis scripts and/or
should be analysed directly with pyhande.'''

import argparse
import os
import pkgutil
import pprint
import sys

import pandas as pd

_script_dir = os.path.abspath(os.path.dirname(__file__))
if not pkgutil.find_loader('pyblock'):
    sys.path.append(os.path.join(_script_dir, 'pyblock'))
if not pkgutil.find_loader('pyhande'):
    sys.path.append(os.path.join(_script_dir, 'pyhande'))

import pyblock
import pyhande

def run_hande_blocking(files, start_iteration=None, end_iteration=None,
                        reblock_plot=None, verbose=1, width=0,
                        out_method='to_string', inefficiency=False,
                        reweight_plot=False, extract_rl_time=False,
                        analysis_method=None, warmup_detection=None):
    '''Run a reblocking analysis on HANDE output and print to STDOUT.

See :func:`pyblock.pd_utils.reblock` and :func:`pyblock.blocking.reblock` for
details on the reblocking procedure.

Parameters
----------
files : list of list of strings
    names of files containing HANDE QMC calculation output.  Each list contains
    the a set of files which are analysed together (ie a series of calculations
    restarted from the previous calculation).
start_iteration : int or None (Default)
    QMC iteration from which statistics are gathered. While the end_iteration
    is included in analysis, the start_iteration is not.
end_iteration : int or None (Default)
    QMC iteration until which statistics are gathered. If None, the last QMC
    iteration included is the last iteration of the data set.
reblock_plot : string
    Filename to which the reblocking convergence plot (standard error vs reblock
    iteration) is saved.  The plot is not created if None and shown
    interactively if '-'.
verbose : int
    Level of verbosity.

    <0: print nothing
    0: print only the estimate from the optimal block length.
    1: print only the recommended statistics from the optimal block length.
    2: print search for automatic starting iteration (if required), blocking
    analysis and recommended statistics.
    3: print calculation metadata, search for automatic starting iteration
    (if required), blocking analysis and recommended statistics.

width : int
    Maximum width (in characters) of lines to print out for
    :class:`pandas.DataFrame` objects; exceeding this results in line wrapping.
    A non-positive value corresponds to disabling line wrapping.
out_method : string
    Output method for printing out tables.  Either 'to_string' to print a
    space-separate table or 'to_csv' to print a CSV table.
inefficiency : bool
    Attempt to calculate the inefficiency factor for the calculations, and
    include it in the output.
reweight_plot: do reweighting the projected energy and show plot to determine
    population bias.
extract_rl_time: extract times taken for a report loop and find mean and errors.

Returns
-------
info :
    Output from :func:`pyhande.lazy.std_analysis`.
opt_block: :class:`pandas.DataFrame`
    Recommended statistics based upon the estimated 'optimal' block size
    as suggested by Wolff and Lee et al. (see
    :func:`pyblock.blocking.find_optimal_block`).
'''

    try:
        float_str = '%-#.8e'
        float_fmt = '{0:-#.8e}'.format
        float_fmt(1.0)
    except ValueError:
        # GAH.  Alternate formatting only added to format function after
        # python 2.6..
        float_str = '%-.8e'
        float_fmt = '{0:-.8e}'.format

    def df_to_x(df, out_method, float_fmt, float_str, width):
        tbl_fn = getattr(df, out_method)
        try:
            return tbl_fn(float_format=float_fmt, na_rep='n/a',
                          line_width=width)
        except TypeError:
            # Annoyingly to_csv and to_string take different types for
            # their float_fmt arguments.
            # See: https://github.com/pydata/pandas/issues/9448
            return tbl_fn(float_format=float_str, na_rep='n/a',
                          line_width=width)

    if width <= 0:
        width = None

    # verbosity levels
    v_silent = -1
    (v_estimate, v_rec_stats, v_analysis, v_meta, v_input) = (0, 1, 2, 3, 4)

    infos = []
    indices = []
    for calc in files:
        try:
            info = pyhande.lazy.std_analysis(calc, start_iteration,
                                             end=end_iteration,
                                             extract_psips=True,
                                             calc_inefficiency=inefficiency,
                                             verbosity = verbose,
                                             extract_rep_loop_time=extract_rl_time,
                                             analysis_method=analysis_method,
                                             warmup_detection=warmup_detection)
            for (i, i_info) in enumerate(info):
                if verbose >= v_analysis:
                    msg = 'Analysing file(s): %s.' % (' '.join(calc))
                    if len(info) > 1:
                        msg += '\nCalculation: %i.' % (i,)
                    msg += ('\nReblocking from iteration: %i.' %
                                (i_info.metadata['pyhande']['reblock_start'],))
                    print(msg)
                if verbose >= v_meta:
                    md = i_info.metadata
                    calc_type = md.pop('calc_type')
                    calc_input = md.pop('input')
                    print('calc_type: %s.\n' % (calc_type))
                    pprint.pprint(md)
                    if verbose >= v_input:
                        print('\nFull input options:\n%s' % '\n'.join(calc_input))
                    print('')
                if verbose >= v_analysis:
                    print(df_to_x(i_info.reblock, out_method, float_fmt, float_str,
                                  width))
                    print('')
            infos.extend(info)
            if len(info) == 1:
                indices.append(','.join(calc))
            else:
                indices.extend((','.join(calc),i) for i in range(len(info)))

            if reweight_plot:
                pyhande.lazy.reweighting_graph(calc, start=start_iteration,
                        verbosity=verbose)

        except ValueError:
            print('WARNING: No data found in file '+' '.join(calc)+'.')
        except RuntimeError as err:
            print('WARNING: Reblocking failed for file '+' '.join(calc)+' with error '+format(err)+'.')
    opt_blocks = [pd.DataFrame(data=
        {'iteration':info.metadata['pyhande']['reblock_start']},
        index = ['Block from']).append(info.opt_block) for info in infos]
    if verbose < v_rec_stats:
        for opt_block in opt_blocks:
            if not opt_block.empty:
                levels = ['mean', 'standard error', 'standard error error']
                for level in levels:
                    opt_block.drop(level, axis=1, inplace=True)
    opt_blocks = [opt_block.stack() for opt_block in opt_blocks]
    opt_block = pd.DataFrame(dict(zip(indices, opt_blocks))).T
    if verbose < v_rec_stats and not opt_block.empty:
        opt_block.columns = opt_block.columns.droplevel(1)

    if not opt_block.empty and verbose > v_silent:
        print('Recommended statistics from optimal block size:')
        print('')
        print(df_to_x(opt_block, out_method, float_fmt, float_str, width))

    for (calc, info) in zip(indices, infos):
        if info.no_opt_block and verbose > v_silent:
            fnames = ''
            if (len(indices) > 1):
                try:
                    fnames = ' in ' + calc.replace(',',' ')
                except AttributeError:
                    # if there is more than one calculation in the file calc is a tuple
                    fnames = ' in ' + calc[0] + ' ' + str(calc[1])
            print('WARNING: could not find optimal block size%s.' % (fnames))
            print('Insufficient statistics collected for the following '
                  'variables: %s.' % (', '.join(info.no_opt_block)))

    if reblock_plot:
        for info in infos:
            pyblock.plot.plot_reblocking(info.reblock, reblock_plot)

    return infos

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

    try:
        cols = pd.util.terminal.get_terminal_size()[0]
    except AttributeError:
        # terminal module moved in pandas 0.20
        cols = pd.io.formats.terminal.get_terminal_size()[0]
    if not sys.stdout.isatty():
        cols = -1

    parser = argparse.ArgumentParser(description = __doc__)
    parser.add_argument('-m', '--merge', default=False, action='store_true',
                        help='Combine data from each file before analysing. '
                        'Separate calculations can be denoted by placing \'--\''
                        ' between groups of files.  Default: treat each file as'
                        ' an independent calculation.')
    parser.add_argument('-o', '--output', default='txt', choices=['txt', 'csv'],
                        help='Format for data table.  Default: %(default)s.')
    parser.add_argument('-p', '--plot', default=None, dest='plotfile',
                        help='Filename to which the reblocking convergence plot '
                        'is saved.  Use \'-\' to show plot interactively.  '
                        'Default: off.')
    parser.add_argument('-r', '--reweight', default=False, dest='reweight_plot', 
                        action='store_true', help='For each independent passed '
                        'calculation show a reweighting plot')
    parser.add_argument('-q', '--quiet', dest='verbose', action='store_const',
                        const=0, default=1,
                        help='Output only the final summary table.  '
                        'Overrides --verbose.')
    parser.add_argument('-s', '--start', type=int, dest='start_iteration',
                        default=None, help='Iteration number from which to '
                        'gather statistics. The start iteration itself is not '
                        'included in the analysis. Default: Try finding '
                        'starting iteration automatically. ')
    parser.add_argument('-e', '--end', type=int, dest='end_iteration',
                        default=None, help='Iteration number until which to '
                        'gather statistics.  Default: Last iteration in data '
                        'set. ')
    parser.add_argument('-v', '--verbose', dest='verbose', action='count',
                        default=1, help='Increase verbosity of the output.  Can '
                        'be specified multiple times.')
    parser.add_argument('-w', '--width', type=int, default=cols,
                        help='Width (in characters) of data to print out '
                        'before wrapping them.  A non-positive value disables '
                        'wrapping.  Default: current terminal width if printing '
                        'to a terminal, -1 if redirecting.')
    parser.add_argument('-i','--inefficiency', default=False, action='store_true',
                        help='Calculate the inefficiency factor for the calculation '
                        'if possible.')
    parser.add_argument('-t','--extract_rl_time', default=False, action='store_true',
                        help='Find the mean time taken for a report loop.')
    parser.add_argument('-a','--analysis_method',  dest='analysis_method', #added_by_ichibha
                        default='reblocking', help='Designate post-analysis method '
                        'to estimate statistic error. [reblocking or hybrid] '
                        'Default: reblocking')
    parser.add_argument('-b','--warmup_detection', dest='warmup_detection', #added_by_ichibha 
                        default='hande_org', help='Designate method to determine '
                        'starting iteration. [hande_org or mser_min] Default: hande_org')
    parser.add_argument('filenames', nargs=argparse.REMAINDER,
                        help='Space-separated list of files to analyse.')


    options = parser.parse_args(args)

    if not options.filenames:
        parser.print_help()
        sys.exit(1)

    if options.merge:
        merged = [[]]
        for fname in options.filenames:
            if fname == '--':
                if merged[-1]:
                    merged.append([])
            else:
                merged[-1].append(fname)
        options.filenames = merged
    else:
        options.filenames = [[fname] for fname in options.filenames]

    out_methods = {'txt': 'to_string', 'csv': 'to_csv'}
    options.output = out_methods[options.output]

    return options

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

    options = parse_args(args)
    run_hande_blocking(options.filenames, options.start_iteration,
                       options.end_iteration, options.plotfile,
                       options.verbose, options.width, options.output,
                       options.inefficiency, options.reweight_plot,
                       options.extract_rl_time, 
                       options.analysis_method, options.warmup_detection) # added_by_ichibha

if __name__ == '__main__':

    main(sys.argv[1:])

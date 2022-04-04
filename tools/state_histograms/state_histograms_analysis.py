#!/usr/bin/env python
''' state_histograms_analysis.py [options] -f file_1 file_2 ... file_N

Analyse the output of a HANDE DMQMC or FCIQMC state histogram calculation by
averaging the provided state histograms by bin/excitation indexes.
'''

import os
import sys
import pkgutil
import argparse
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

try:
    from pyhande.state_histograms import analyse_state_histograms
except ModuleNotFoundError:
    _script_dir = os.path.dirname(os.path.abspath(__file__))
    if not pkgutil.find_loader('pyhande'):
        sys.path.append(os.path.join(_script_dir, '../pyhande'))
    from pyhande.state_histograms import analyse_state_histograms


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
    parser.add_argument('-pdf', '--pdf-plot', action='store', default=None,
                        type=str, dest='pdf_plot', help='Define a file name '
                        'for the state histograms to plotted and saved to.')
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
    parser.add_argument('-fci', '--fciqmc-data', action='store_true',
                        dest='fciqmc', default=False, help='Average '
                        'calculations across iterations which can be more '
                        'appropriate for data from an FCIQMC simulation which '
                        'has reached the ground state.')
    parser.add_argument('-f', '--files', type=str, nargs='+', action='append',
                        dest='list_of_filenames', help='EXLEVEL files to '
                        'analyse. For multiple analysis separate sets of '
                        'files with the -f flag.')
    parser.parse_args(args=None if arguments else ['--help'])

    options = parser.parse_args(arguments)

    return (options.list_of_filenames, options)


def plot_state_histograms(list_of_results, pdf_name):
    ''' Plot the user provided state histograms if they requested it.

    Parameters
    ----------
    list_of_results : list of :class:`pandas.DataFrame`
        A list containing the average number of determinants in bin each
        bin for a given temperature/imaginary time.
    pdf_name : str
        A filename for the plot(s) to be saved to.

    Results
    -------
    None.

    Raises
    ------
    RuntimeError
        When there are multiple data sets, if the data sets have a different
        number of columns.
    '''
    font = {'family': 'serif', 'sans-serif': 'Computer Modern Roman'}
    mpl.rc('font', **font)
    mpl.rc('savefig', dpi=300)
    mpl.rc('lines', lw=2, markersize=5)
    mpl.rc('legend', fontsize=8, numpoints=1)
    mpl.rc(('axes', 'xtick', 'ytick'), labelsize=8)
    mpl.rc('figure', dpi=300, figsize=(3.37, 3.37*((5)**0.5-1)/2))

    if '.pdf' not in pdf_name:
        pdf_name += '.pdf'

    bkey = 'bin_edges'
    sum_keys, sem_keys, xmins, xmaxs, ymins = [], [], [], [], []
    for results in list_of_results:
        sum_columns = [c for c in results.columns if '\\sum' in c]
        sem_columns = [c.replace('\\sum', 'sem') for c in sum_columns]
        sum_keys.append(sum_columns)
        sem_keys.append(sem_columns)

        min_nonzero_sums = results[sum_columns][results[sum_columns] != 0.0]
        xmin = 10**round(np.log10(np.nanmin(min_nonzero_sums)) - 1)
        xmax = 10**round(np.log10(np.max(results[sum_columns].max())) + 1)

        imaxs = [np.max(np.nonzero(np.array(results[c]))) for c in sum_columns]
        ymin = results[bkey].min()/results[bkey].iloc[np.max(imaxs)]
        ymin = 10**round(np.log10(ymin) - 1)

        xmins.append(xmin)
        xmaxs.append(xmax)
        ymins.append(ymin)

    if any(len(sum_keys[0]) != len(keys) for keys in sum_keys):
        raise RuntimeError('Data sets must have the same shape to be '
                           'plotted together!')

    xmin, xmax, ymin = min(xmins), max(xmaxs), min(ymins)
    max_data_set_length = max([len(sum_key) for sum_key in sum_keys])
    ndata_sets = len(list_of_results)

    with PdfPages(pdf_name) as pdf:

        for iplot in range(max_data_set_length):
            plt.clf()

            for idata in range(ndata_sets):

                iresults = list_of_results[idata]
                sum_key = sum_keys[idata][iplot]
                sem_key = sem_keys[idata][iplot]

                these_data = iresults[[bkey, sum_key, sem_key]]
                these_data = these_data[these_data[sum_key] != 0.0]
                these_data[bkey] /= these_data[bkey].max()

                ireport = int(sem_key.strip("sem"))

                plt.plot(
                        these_data[sum_key],
                        these_data[bkey],
                        color=f'C{idata}',
                        lw=1,
                        zorder=2,
                    )
                plt.fill_betweenx(
                        these_data[bkey],
                        these_data[sum_key] - these_data[sem_key],
                        these_data[sum_key] + these_data[sem_key],
                        label=f'Data set: {idata+1}, Report: {ireport}',
                        color=f'C{idata}',
                        alpha=0.5,
                    )

                plt.axhline(
                        y=1E0,
                        color='k',
                        lw=1,
                        ls=':',
                        alpha=0.6,
                    )

            plt.xlim(left=xmin, right=xmax)
            plt.ylim(bottom=ymin, top=1.5E0)
            plt.xscale('log')
            plt.yscale('log')
            plt.xlabel(r'Density matrix elements')
            plt.ylabel(r'Normalized population')

            plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.25), ncol=1)
            pdf.savefig(bbox_inches='tight')


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

    (list_of_filenames, options) = parse_arguments(arguments)

    list_of_results = []
    for ifiles, files in enumerate(list_of_filenames):
        ifiles += 1
        results = analyse_state_histograms(files, options.fciqmc)
        list_of_results.append(results)

        if options.output_file is not None:
            file_stream_name = str(ifiles).zfill(5) + options.output_file
            file_stream = open(file_stream_name + '.' + options.output, 'w')
        else:
            file_stream = sys.stdout

        if options.output == 'csv':
            if options.csv_tol:
                results = results.round(options.csv_tol)
            print(results.to_csv(index=False), file=file_stream)
        elif options.output == 'txt':
            results_str = results.to_string(index=False,
                                            float_format=options.float_format)
            print(results_str, file=file_stream)

        if options.output_file is not None:
            file_stream.close()

    if options.pdf_plot is not None:
        plot_state_histograms(list_of_results, options.pdf_plot)


if __name__ == '__main__':
    main(sys.argv[1:])

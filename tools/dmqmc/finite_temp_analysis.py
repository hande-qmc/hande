#!/usr/bin/python
'''finite_temp_analysis.py [options] file_1 file_2 ... file_N

Analyse the output of a HANDE DMQMC calculation by averaging
temperature-dependent data across beta loops.'''

import pandas as pd
from os import path
import sys
sys.path.append(path.join(path.abspath(path.dirname(sys.argv[0])), 'pyblock'))
import pyblock
import pyhande
import optparse

def run_dmqmc_analysis(estimates):
    '''Perform analysis of DMQMC data from a HANDE calculation.

Parameters
----------
estimates : :class:`pandas.DataFrame`
    All of the estimates from all beta loops, whcich are to be combined and analysed.

Returns
-------
None.
'''

    means = estimates.mean(level=2)

def parse_args(args):
    '''Parse command-line arguments.

Parameters
----------
args : list of strings
    command-line arguments.

Returns
-------
filenames : list of strings
    list of QMC output files
'''

    parser = optparse.OptionParser(usage = __doc__)

    (options, filenames) = parser.parse_args(args)

    if not filenames:
        parser.print_help()
        sys.exit(1)

    return (filenames)

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

    (files) = parse_args(args)
    (metadata, data) = pyhande.extract.extract_data_sets(files)
    data.set_index('iterations', inplace=True, append=True)
    estimates = data.loc[:,'Shift':'# H psips']
    run_dmqmc_analysis(estimates)

if __name__ == '__main__':

    main(sys.argv[1:])


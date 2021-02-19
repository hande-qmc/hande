#!/usr/bin/env python
"""Analyse CCMC/FCIQMC calculations run with HANDE QMC.

Modified from reblock_hande.py.
Files may be compressed with either gzip, bzip2 or xz.  CCMC and FCIQMC
calculations only are supported; other calculations have specific
analysis scripts and/or should be analysed directly with pyhande.
"""

from typing import Any, List
import argparse
import sys
import pandas as pd
import pyhande


def parse_args(args: List[str]) -> Any:
    """Parse command-line arguments.

    Parameters
    ----------
    args : List[str]
        command-line arguments.

    Returns
    -------
    (merge, start, inefficiency, shoulder, analysis_method,
    warmup_detection, filenames)

    where

    See individual descriptions for details.
    """

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-m", "--merge", default="uuid",
                        choices=["uuid", "legacy", "no"],
                        help="Combine data from each file before analysing. "
                        "Default: Merge using 'uuid' strategy, ie. using "
                        "UUIDs.")
    parser.add_argument("-s", "--start", type=int, dest="start_iteration",
                        default=None, help="Iteration number from which to "
                        "gather statistics. Note that this applies to all "
                        "calculations passed.  Default: Try finding "
                        "starting iteration automatically. ")
    parser.add_argument("-i", "--inefficiency", default=False,
                        action="store_true", help="Calculate the inefficiency "
                        "factor for the calculation(s) if possible.")
    parser.add_argument("-sh", "--shoulder", default=False,
                        action="store_true", help="Determine the shoulder for "
                        "the calculation(s) if possible.")
    parser.add_argument("-a", "--analysis_method",  dest="analysis_method",
                        default="blocking", choices=["blocking", "hybrid"],
                        help="Designate the post-analysis method "
                        "to estimate the statistic error. Default: "
                        "%(default)s")
    parser.add_argument("-wd", "--warmup_detection", dest="warmup_detection",
                        default="blocking", choices=["blocking", "mser"],
                        help="Designate the method to determine "
                        "the starting iterations to be discarded before "
                        "calculating the statistic error.  Ignored if start "
                        "is given.  Default: %(default)s")
    parser.add_argument("filenames", nargs=argparse.REMAINDER,
                        help="Space-separated list of files to analyse.")

    options = parser.parse_args(args)

    if not options.filenames:
        parser.print_help()
        sys.exit(1)

    return options


def main(args: List[str]) -> None:
    """Run reblocking and data analysis on HANDE output.

    Parameters
    ----------
    args : list of strings
        command-line arguments.

    Returns
    -------
    None.
    """

    options = parse_args(args)
    start_its = ([options.start_iteration]*len(options.filenames)
                 if options.start_iteration or options.start_iteration == 0
                 else options.warmup_detection)
    results = pyhande.results_viewer.get_results.get_results(
        options.filenames, merge_type=options.merge,
        analyser=options.analysis_method, start_its=start_its)
    if options.inefficiency:
        results.add_inefficiency()
    if options.shoulder:
        results.add_shoulder()
    # See https://pandas.pydata.org/pandas-docs/stable/user_guide/options.html
    # Set to a high number.
    pd.options.display.max_columns = 300
    print(results.summary_pretty)


if __name__ == "__main__":

    main(sys.argv[1:])

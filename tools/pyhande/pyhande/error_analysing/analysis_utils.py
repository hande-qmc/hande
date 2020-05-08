"""Shared helper functions for analysers."""
from typing import Dict, List, Union
import pandas as pd


def check_data_input(data: List[pd.DataFrame], cols: List[str],
                     eval_ratio: Dict[str, str], hybrid_col: str,
                     start_its: Union[List[int], str], end_its: List[int]):
    """Check data input against other, previous, input.

    Parameters
    ----------
    data : List[pd.DataFrame]
        List of QMC data.
    cols : List[str]
        Columns to be analysed when blocking/ finding starting iteration
        with 'blocking'.
    eval_ratio : Dict[str, str]
        Contains information to evaluate ratio (e.g. projected energy)
        when doing blocking analysis.
    hybrid_col : str
        Column name when doing hybrid analysis/ finding starting
        iteration with 'mser'.
    start_its : Union[List[int], str]
        Starting iterations for analysis or information on type of
        find_starting_it function.
    end_its : List[int]
        Last iterations for analysis.

    Raises
    ------
    ValueError
        If cols, eval_ratio, hybrid_cols are specified but not in data
        respectively.
        If start_its/end_its are lists of iterations but the list has a
        different length than data.
    """
    if cols and not all([col in dat for dat in data for col in cols]):
        raise ValueError("'cols' parameter must only contain columns names "
                         "present in all dataframes in 'data'.")
    if eval_ratio and not all([col in dat for dat in data for col in
                               [eval_ratio['num'], eval_ratio['denom']]]):
        raise ValueError("When 'eval_ratio' is defined, its values must be "
                         "column names present in all dataframes in 'data'.")
    if hybrid_col and not all(hybrid_col in dat for dat in data):
        raise ValueError("When 'hybrid_col' is defined (necessary when doing "
                         "hybrid analysis or using 'mser' starting iteration "
                         "finder), it has to be a column in all elements of "
                         "'data'.")
    if (isinstance(start_its, list) and len(start_its) != len(data)):
        raise ValueError(f"If 'start_its' (here of length {len(start_its)}) "
                         "is specified as list of start iterations, it has to "
                         "be a list of the same length as 'data' (here of "
                         f"length {len(data)}).")
    if end_its and len(end_its) != len(data):
        raise ValueError(f"If 'end_its' (here of length {len(end_its)}) is "
                         "specified, it has to be a list of the same length "
                         f"as 'data' (here of length {len(data)}).")


def set_start_and_end_its(
        data: List[pd.DataFrame], it_key: str, cols: List[str],
        hybrid_col: str, find_starting_it_func,
        find_start_kw_args: Dict[str, Union[bool, float, int]],
        start_its: List[int], end_its: List[int]) -> (List[int], List[int]):
    """Find end and start iteration if required.

    If `end_its` is None, the last iteration will be chosen respectively
    and `end_its` updated.
    If `start_its` is None, the start iteration is found using
    `find_starting_it_func`.
    If not None, they are returned as passed in.

    Parameters
    ----------
    data : List[pd.DataFrame]
        List of data DataFrames with QMC data.
    it_key : str
        Name of iteration columns, e.g. 'iterations'.
    cols : List[str]
        Columns to be analysed if doing 'blocking' starting it finding.
    hybrid_col : str
        Column to be analysed if doing 'mser' starting it finding.
    find_starting_it_func : Function from find_starting_iteration.py
        Function to find starting iteration.  Adheres to certain common
        interface.
    find_start_kw_args : Dict[str, Union[bool, float, int]
        Extra input parameters for finding starting iteration.
    start_its : List[int]
        List of starting iterations if already present.
    end_its : List[int]
        List of final iterations if already present.

    Returns
    -------
    (List[int], List[int])
        Lists of start and final iterations.
    """
    if not end_its:
        end_its = [dat[it_key].iloc[-1] for dat in data]
    if not start_its:
        start_its = [
            find_starting_it_func(dat, end_it, it_key, cols, hybrid_col,
                                  **find_start_kw_args)
            for dat, end_it in zip(data, end_its)
        ]
    return start_its, end_its

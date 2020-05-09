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

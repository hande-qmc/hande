"""Shared helper functions for analysers."""
from typing import Dict, List, Optional, Tuple, Union
import pandas as pd


def check_data_input(data: List[pd.DataFrame],
                     cols: Union[List[str], Optional[List[str]]],
                     eval_ratio: Optional[Dict[str, str]],
                     hybrid_col: Union[Optional[str], str],
                     start_its: Union[List[int], str],
                     end_its: Optional[List[int]]) -> None:
    """Check data input against other, previous, input.

    Parameters
    ----------
    data : List[pd.DataFrame]
        List of QMC data.
    cols : Union[List[str], Optional[List[str]]]
        Columns to be analysed when blocking/ finding starting iteration
        with 'blocking'.
    eval_ratio : Optional[Dict[str, str]]
        Contains information to evaluate ratio (e.g. projected energy)
        when doing blocking analysis.
    hybrid_col : Union[Optional[str], str]
        Column name when doing hybrid analysis/ finding starting
        iteration with 'mser'.
    start_its : Union[List[int], str]
        Starting iterations for analysis or information on type of
        find_starting_it function.
    end_its : Optional[List[int]]
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


def _set_value(observables: Dict[str, str], col_item: Union[Optional[str], str]
               ) -> Union[Optional[str], str]:
    """Helper to set value of column names. e.g. if col_item='N_0', return 'N_0'.
       If col_item='obs:ref_key', return observables['ref_key']"""
    return (
        observables[col_item[4:]] if (col_item and col_item.startswith('obs:'))
        else col_item
    )


def set_cols(observables: Dict[str, str], it_key: str,
             cols: Union[Optional[List[str]], List[str]], replica_col: str,
             eval_ratio: Optional[Dict[str, str]],
             hybrid_col: Union[Optional[str], str]
             ) -> Tuple[str, Union[Optional[List[str]], List[str]], str,
                        Optional[Dict[str, str]], str]:
    """Set various columns and observable names.

    Either the input is simply returned or set to observables[input] if
    input starts with 'obs:'.

    Parameters
    ----------
    observables : Dict[str, str]
        Map of key to column/observable name, e.g. {'ref_key': 'N_0'}
    it_key : str
        Key or actual name for iterations.
    cols : Union[Optional[List[str]], List[str]]
        Keys or actual names of columns/observables to be analysed in
        blocking.
    replica_col : str
        Key or actual name for replica column.
    eval_ratio : optional[Dict[str, str]]
        Keys or actual names of elements in observable ratio to be
        evaluated.
    hybrid_col : Union[Optional[str], str]
        Key or actual name of column/observable to be analysed in hybrid
        analysis.

    Returns
    -------
    (Set) values from above (except `observables`).
    """
    it_key = _set_value(observables, it_key)
    cols = [_set_value(observables, col) for col in cols] if cols else None
    replica_col = _set_value(observables, replica_col)
    eval_ratio = {ev_key: _set_value(observables, eval_ratio[ev_key])
                  for ev_key in eval_ratio.keys()} if eval_ratio else None
    hybrid_col = _set_value(observables, hybrid_col)
    return it_key, cols, replica_col, eval_ratio, hybrid_col

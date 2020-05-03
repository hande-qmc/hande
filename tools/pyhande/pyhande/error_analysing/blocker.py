"""Analyse Monte Carlo correlated output using reblocking."""
from typing import Dict, List, Union
import copy
import pandas as pd
import pyblock
import pyhande.analysis as analysis
from pyhande.error_analysing.abs_error_analyser import AbsErrorAnalyser
from pyhande.error_analysing.find_starting_iteration import select_find_start


class Blocker(AbsErrorAnalyser):
    """Reblock specified columns from HANDE output using pyblock.

    [todo] - cite Flyberg?
    Can be used instead of Hybrid (which is not implemented yet in this
    form).
    """

    def __init__(
            self, cols: List[str] = None, eval_ratio: Dict[str, str] = None,
            it_key: str = 'iterations',
            start_its: Union[List[int], str] = 'blocking',
            end_its: List[int] = None,
            find_start_kw_args: Dict[str, Union[bool, float, int]]
            = None) -> None:
        r"""
        Initialise a Blocking instance.

        Parameters
        ----------
        cols : List[str], optional
            Columns in QMC data to be reblocked.
            The default is
            ['Shift', '\sum H_0j N_j', 'N_0', '# H psips'].
        eval_ratio: Dict, optional
            Evaluate mean and standard error of ratio of column 'num'
            with column 'denom'. Add this to the blocking results as
            column with name 'name'.
            The default is {'name': 'Proj. Energy', 'num': 'N_0',
                            'denom': '# H psips'}.
        it_key: str, optional
            MC iterations columns key.
            The default is 'iterations'.
        start_its : Union[List[int], str], optional
            Either list of start iterations which has to be of same
            length as data.
            Or string specifying find starting iteration function to use
            to estimate starting iterations automatically, either
            'blocking' (`find_starting_iteration_blocking()`) or
            'mser', (`find_starting_iteration_mser_min()`).  Note that
            this choice then applies to all calculations in passed in
            `data` to `.exe()`.
            The default is 'blocking'.
        end_its : List[int], optional
            List of end iterations. Has to be of same length as
            data. If None, end iteration is the last iteration for each
            calculation.
            The default is None.
        find_start_kw_args : Dict[str, Union[bool, float, int]],
                optional
            Possible extra arguments that can be passed to
            find_starting_iterations functions.  See their definitions
            for details.  E. g. {'show_graph' : True} shows a graph
            highlighting the starting iteration found when
            start_its = 'blocking'.
            The default is None, which defaults to an empty dictionary.
        """
        self._cols: List[str] = cols if cols else [
            'Shift', r'\sum H_0j N_j', 'N_0', '# H psips'
        ]
        self._eval_ratio: Dict[str, str] = eval_ratio if eval_ratio else {
            'name': 'Proj. Energy', 'num': r'\sum H_0j N_j', 'denom': 'N_0'
        }
        if any([v not in self._cols for v in [self._eval_ratio['num'],
                                              self._eval_ratio['denom']]]):
            raise ValueError(f"'eval_ratio' requires "
                             f"'{self._eval_ratio['num']}' and "
                             f"'{self._eval_ratio['denom']}' columns to be "
                             "specified in 'cols'.")
        self._it_key = it_key
        if isinstance(start_its, list):
            self._start_its = start_its
        else:
            self._start_its = None
            self._find_starting_it = select_find_start(start_its)
        self._end_its: List[int] = end_its
        self._find_start_kw_args: Dict[str, Union[bool, float, int]] = (
            find_start_kw_args if find_start_kw_args else {})
        # These attributes are set later:
        self._reblock: List[pd.DataFrame]
        self._covariance: List[pd.DataFrame]
        self._data_len: List[pd.Series]
        self._opt_block: pd.DataFrame
        self._no_opt_block: List[List[str]]

    @property
    def cols(self) -> List[str]:
        """Access _cols attribute, the columns (to be) blocked."""
        return self._cols

    @property
    def reblock(self) -> List[pd.DataFrame]:
        """Access _reblock attribute if available. Else raise error."""
        try:
            return self._reblock
        except AttributeError:
            print("First do analysis by running 'exe' instance method.")
            raise

    @property
    def data_len(self) -> List[pd.DataFrame]:
        """Access _data_len attribute if available. Else raise error."""
        try:
            return self._data_len
        except AttributeError:
            print("First do analysis by running 'exe' instance method.")
            raise

    @property
    def covariance(self) -> List[pd.DataFrame]:
        """Access _covariance attribute if available. Else error."""
        try:
            return self._covariance
        except AttributeError:
            print("First do analysis by running 'exe' instance method.")
            raise

    def exe(self, data: List[pd.DataFrame]):
        """
        Do reblocking (first finding starting iteration if required).

        Parameters
        ----------
        data : List[pd.DataFrame]
            HANDE calculation Monte Carlo output data.

        Raises
        ------
        ValueError
            If not all columns to be blocked appear in 'data'
            or if the length of 'data' is different to length of
            'start_its' or 'end_its' if they are defined.

        """
        self._check_data_input(data)
        self._set_start_and_end_its(data)

        # Do blocking.
        self._data_len = []
        self._reblock = []
        self._covariance = []
        self._no_opt_block = []
        self._opt_block = []
        for dat, start_it, end_it in zip(data, self.start_its,
                                         self.end_its):
            # Select subset of data to block:
            dat_c = dat.copy()
            dat_c = dat_c[dat_c['iterations'].between(start_it, end_it)]
            dat_c = dat_c[self.cols]
            # Block:
            (data_len, reblock, covariance) = pyblock.pd_utils.reblock(dat_c)
            cols_in_opt = copy.copy(self.cols)
            # Add ratio if required:
            if self.eval_ratio:
                ratio = analysis.projected_energy(
                    reblock, covariance, data_len,
                    sum_key=self.eval_ratio['num'],
                    ref_key=self.eval_ratio['denom'],
                    col_name=self.eval_ratio['name'])
                reblock = pd.concat([reblock, ratio], axis=1)
                cols_in_opt.append(self.eval_ratio['name'])
            # Get summary of blocking:
            (opt_block, no_opt_block) = analysis.qmc_summary(reblock,
                                                             cols_in_opt)
            self._data_len.append(data_len)
            self._reblock.append(reblock)
            self._covariance.append(covariance)
            self._no_opt_block.append(no_opt_block)
            self.opt_block.append(opt_block)

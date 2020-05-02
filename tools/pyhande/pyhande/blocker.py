"""Analyse Monte Carlo correlated output using reblocking."""
from typing import Dict, List
import copy
import pandas as pd
import pyblock
import pyhande.analysis as analysis
import pyhande.lazy as lazy


class Blocking:
    """Reblock specified columns from HANDE output using pyblock.

    [todo] - cite Flyberg?
    Can be used instead of Hybrid (which is not implemented yet in this
    form).
    """

    def __init__(
            self,
            cols: List[str] = None, eval_ratio: Dict[str, str] = None,
            start_its: List[int] = None, end_its: List[int] = None) -> None:
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
        start_its : List[int]
            List of start iterations. Has to be of same length as
            data. If None, estimate starting iterations automatically.
            The default is None.
        end_its : List[int], optional
            List of end iterations. Has to be of same length as
            data. If None, end iteration is the last iteration for each
            calculation.
            The default is None.

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
        self._start_its: List[int] = start_its
        self._end_its: List[int] = end_its
        # These attributes are set later:
        self._data: List[pd.DataFrame]
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
    def eval_ratio(self) -> Dict[str, str]:
        """Access _eval_ratio attribute, ratio to evaluate."""
        return self._eval_ratio

    @property
    def start_its(self) -> List[int]:
        """Access _start_its attribute, blocking start iterations."""
        return self._start_its

    @property
    def end_its(self) -> List[int]:
        """Access _end_its attribute, blocking end iterations."""
        return self._end_its

    @property
    def data(self) -> List[pd.DataFrame]:
        """Access _data attribute if available. Else raise error."""
        try:
            return self._data
        except AttributeError:
            print("First pass data to 'exe' instance method.")
            raise

    @property
    def reblock(self) -> List[pd.DataFrame]:
        """Access _reblock attribute if available. Else raise error."""
        try:
            return self._reblock
        except AttributeError:
            print("First do reblocking by running 'exe' instance method.")
            raise

    @property
    def data_len(self) -> List[pd.DataFrame]:
        """Access _data_len attribute if available. Else raise error."""
        try:
            return self._data_len
        except AttributeError:
            print("First do reblocking by running 'exe' instance method.")
            raise

    @property
    def covariance(self) -> List[pd.DataFrame]:
        """Access _covariance attribute if available. Else error."""
        try:
            return self._covariance
        except AttributeError:
            print("First do reblocking by running 'exe' instance method.")
            raise

    @property
    def opt_block(self) -> pd.DataFrame:
        """Access _opt_block attribute if available. Else error."""
        try:
            return self._opt_block
        except AttributeError:
            print("First do reblocking by running 'exe' instance method.")
            raise

    @property
    def no_opt_block(self) -> List[List[str]]:
        """Access _no_opt_block attribute if available. Else error."""
        try:
            return self._no_opt_block
        except AttributeError:
            print("First do reblocking by running 'exe' instance method.")
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
        # Check and set parameter.
        if not all([col in dat for dat in data for col in self.cols]):
            raise ValueError("'cols' parameter must only contain columns "
                             "names present in all dataframes in 'data'. If "
                             "you have not specified 'cols', specify it "
                             "explicity in another instance of Blocking to "
                             "only contain columns in all of your 'data'.")
        if self.start_its and len(self.start_its) != len(data):
            raise ValueError("If 'start_its' (here of length "
                             f"{len(self.start_its)}) is specified, it has to "
                             "be a list of the same length as 'data' (here of "
                             f"length {len(data)}).")
        if self.end_its and len(self.end_its) != len(data):
            raise ValueError("If 'end_its' (here of length "
                             f"{len(self.end_its)}) is specified, it has to "
                             "be a list of the same length as 'data' (here of "
                             f"length {len(data)}).")
        self._data = data

        # Find end and start iteration if required.
        if not self.end_its:
            self._end_its = [dat['iterations'].iloc[-1] for dat in self.data]
        if not self.start_its:
            # Modify! Don't use lazy.
            self._start_its = [
                lazy.find_starting_iteration(dat, {}) for dat in self.data
            ]

        # Do blocking.
        self._data_len = []
        self._reblock = []
        self._covariance = []
        self._no_opt_block = []
        self._opt_block = []
        for dat, start_it, end_it in zip(self._data, self.start_its,
                                         self.end_its):
            # Select subset of data to block:
            dat_c = dat.copy()
            dat_c = dat_c[dat_c['iterations'].between(start_it + 1, end_it)]
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

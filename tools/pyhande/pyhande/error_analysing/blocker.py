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

    This uses pyblock [1] to do reblocking, see Ref. [2] for more
    details on the reblocking algorithm.
    Can be used instead of Hybrid.

    [1] - pyblock, James Spencer, http://github.com/jsspencer/pyblock
    [2] - Flyvbjerg, H., Petersen, H. G.. J. Chem. Phys. 91, 461, 1989.
    """

    def __init__(
            self, it_key: str, cols: List[str],
            eval_ratio: Dict[str, str] = None, hybrid_col: str = None,
            start_its: Union[List[int], str] = 'blocking',
            end_its: List[int] = None,
            find_start_kw_args: Dict[str, Union[bool, float, int]]
            = None) -> None:
        r"""
        Initialise a Blocking instance.

        Parameters
        ----------
        it_key : str
            Column name of MC iterations, e.g. 'iterations'.
        cols : List[str]
            Columns in QMC data to be reblocked, e.g.
            ['Shift', r'\sum H_0j N_j', 'N_0', '# H psips'].
        eval_ratio : Dict[str, str], optional
            After blocking, evaluate ratio (e.g. projected energy).
            Contains:
            `num` : Numerator of ratio, e.g. r'\sum H_0j N_j'.
            `denom` : Denominator of ratio, e.g. 'N_0'.
            `name` : Name of ratio, e.g. 'Proj. Energy'.
            The default is None, in which case no ratio is evaluated.
        hybrid_col : str, optional
            Column to be analysed in a hybrid way, e.g.
            'Inst. Proj. Energy'.  Here, only used if
            `start_its` = 'mser', ignored otherwise.
            The default is None.  Has to be specified if `start_its` =
            'mser'.
        start_its : Union[List[int], str], optional
            Either list of start iterations which has to be of same
            length as data.
            Or string specifying find starting iteration function to use
            to estimate starting iterations automatically, either
            'blocking' (`find_starting_iteration_blocking()`) or
            'mser', (`find_starting_iteration_mser_min()`).  In case of
            'mser', need to specify `hybrid_col`.  Note that this choice
            then applies to all calculations in passed in `data` to
            `.exe()`.
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
        # Set input.
        Blocker._check_input(it_key, cols, hybrid_col, start_its)
        self._it_key: str = it_key
        self._cols: List[str] = cols
        self._eval_ratio: Dict[str, str] = eval_ratio
        self._hybrid_col: str = hybrid_col
        if isinstance(start_its, list):
            self._start_its = start_its
        else:
            self._start_its = None
            self._find_starting_it = select_find_start(start_its)
        self._end_its: List[int] = end_its
        self._find_start_kw_args: Dict[str, Union[bool, float, int]] = (
            find_start_kw_args if find_start_kw_args else {})

        # Set helper attribute.
        self._pre_exe_error_message = (
            "First do analysis by running 'exe' instance method."
        )

        # These attributes are set later:
        self._reblock: List[pd.DataFrame]
        self._covariance: List[pd.DataFrame]
        self._data_len: List[pd.Series]
        self._opt_block: pd.DataFrame
        self._no_opt_block: List[List[str]]

    @staticmethod
    def _check_input(it_key: str, cols: List[str], hybrid_col: str,
                     start_its: Union[List[int], str]):
        """Check some input parameters."""
        if not it_key:
            raise ValueError("'it_key' must be specified!")
        if not cols:
            raise ValueError("'cols' has to be specified.")
        if start_its == 'mser' and not hybrid_col:
            raise ValueError("When starting iterations should be found with "
                             "'mser', i.e. 'start_its' == 'mser', 'hybrid_col' "
                             "has to be specified.")

    @property
    def start_its(self) -> List[int]:
        """Access _start_its attribute, analysis start iterations."""
        return self._start_its

    @property
    def end_its(self) -> List[int]:
        """Access _end_its attribute, analysis end iterations."""
        return self._end_its

    @property
    def opt_block(self) -> pd.DataFrame:
        """Access _opt_block attribute if available. Else error."""
        try:
            return self._opt_block
        except AttributeError:
            print(self._pre_exe_error_message)
            raise

    @property
    def no_opt_block(self) -> List[List[str]]:
        """Access _no_opt_block attribute if available. Else error."""
        try:
            return self._no_opt_block
        except AttributeError:
            print(self._pre_exe_error_message)
            raise

    @property
    def reblock(self) -> List[pd.DataFrame]:
        """Access _reblock attribute if available. Else raise error."""
        try:
            return self._reblock
        except AttributeError:
            print(self._pre_exe_error_message)
            raise

    @property
    def data_len(self) -> List[pd.DataFrame]:
        """Access _data_len attribute if available. Else raise error."""
        try:
            return self._data_len
        except AttributeError:
            print(self._pre_exe_error_message)
            raise

    @property
    def covariance(self) -> List[pd.DataFrame]:
        """Access _covariance attribute if available. Else error."""
        try:
            return self._covariance
        except AttributeError:
            print(self._pre_exe_error_message)
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
            dat_c = dat_c[dat_c[self._it_key].between(start_it, end_it)]
            dat_c = dat_c[self._cols]
            # Block:
            (data_len, reblock, covariance) = pyblock.pd_utils.reblock(dat_c)
            cols_in_opt = copy.copy(self._cols)
            # Add ratio if required:
            if self.eval_ratio:
                ratio = analysis.projected_energy(
                    reblock, covariance, data_len,
                    sum_key=self._eval_ratio['num'],
                    ref_key=self._eval_ratio['denom'],
                    col_name=self._eval_ratio['name'])
                reblock = pd.concat([reblock, ratio], axis=1)
                cols_in_opt.append(self._eval_ratio['name'])
            # Get summary of blocking:
            (opt_block, no_opt_block) = analysis.qmc_summary(reblock,
                                                             cols_in_opt)
            self._data_len.append(data_len)
            self._reblock.append(reblock)
            self._covariance.append(covariance)
            self._no_opt_block.append(no_opt_block)
            self._opt_block.append(opt_block)

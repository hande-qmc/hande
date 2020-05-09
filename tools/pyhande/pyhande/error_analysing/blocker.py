"""Analyse Monte Carlo correlated output using reblocking."""
from typing import Dict, List, Tuple, Union
import copy
import pandas as pd
import pyblock
import pyhande.analysis as analysis
from pyhande.error_analysing.abs_error_analyser import AbsErrorAnalyser
from pyhande.error_analysing.find_starting_iteration import select_find_start
from pyhande.error_analysing.analysis_utils import check_data_input


class Blocker(AbsErrorAnalyser):
    """Reblock specified columns from HANDE output using pyblock.

    Can be used instead of HybridAnalyser.

    This uses pyblock [1] to do reblocking, see Ref. [2] for more
    details on the reblocking algorithm.

    [1] - pyblock, James Spencer, http://github.com/jsspencer/pyblock
    [2] - Flyvbjerg, H., Petersen, H. G., 1989, J. Chem. Phys. 91, 461.
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
            it is up to the user to specify correct keys here.
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

    @classmethod
    def inst_hande_ccmc_fciqmc(
            cls, observables: Dict[str, str],
            start_its: Union[List[int], str] = 'blocking',
            end_its: List[int] = None,
            find_start_kw_args: Dict[str, Union[bool, float, int]] = None
    ):
        """Return Blocker instance for a HANDE CCMC/FCIQMC calculation.

        Parameters
        ----------
        observables : Dict[str, str]
            Maps generic column names, 'it_key', 'shift_key', etc, to
            their HANDE CCMC/FCIQMC column names, e.g.
            PrepHandeCcmcFciqmc.observables in data_preparing.
            'it_key', 'shift_key', 'sum_key', 'ref_key', 'total_key',
            'proje_key' and 'inst_proje_key' are required here.
        For the other arguments, see __init__().

        Returns
        -------
        Blocker
            Instance of the Blocker class, customised for a HANDE CCMC/
            FCIQMC calculation.
        """
        return Blocker(
            observables['it_key'],
            [observables['shift_key'], observables['sum_key'],
             observables['ref_key'], observables['total_key']],
            eval_ratio={'name': observables['proje_key'],
                        'num': observables['sum_key'],
                        'denom': observables['ref_key']},
            hybrid_col=observables['inst_proje_key'], start_its=start_its,
            end_its=end_its, find_start_kw_args=find_start_kw_args)

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
                             "'mser', i.e. 'start_its' == 'mser', 'hybrid_col'"
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

    def _do_blocking_dat(self, dat: pd.DataFrame, dat_ind: int
                         ) -> Tuple(pd.Series, pd.DataFrame, pd.DataFrame,
                                    List[str], pd.DataFrame):
        """Block one QMC calculation, one replica if replica tricks.

        Parameters
        ----------
        dat : pd.DataFrame
            QMC data for one calculation (of one replica if doing that).
        dat_in : int
            Index of `dat` in all `data` passed to .exe().

        Returns
        -------
        pd.Series, pd.DataFrame, pd.DataFrame, List[str], pd.DataFrame
            data_len, reblock, covariance, no_opt_block, opt_block of
            analysis.
        """
        # Find start and end iteration:
        end_it = (self._end_its[dat_ind] if self._end_its
                  else dat[self._it_key].iloc[-1])
        start_it = (self._start_its[dat_ind] if self._start_its
                    else self._find_starting_it(dat, end_it, self._it_key,
                                                self._cols,
                                                self._hybrid_col,
                                                **self._find_start_kw_args)
                    )
        # Select subset of data to block:
        dat = dat[dat[self._it_key].between(start_it, end_it)]
        dat = dat[self._cols]
        # Block:
        (data_len, reblock, covariance) = pyblock.pd_utils.reblock(dat)
        cols_in_opt = copy.copy(self._cols)
        # Add ratio if required:
        if self._eval_ratio:
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
        return (data_len, reblock, covariance, no_opt_block, opt_block)

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
        check_data_input(
            data, self._cols, self._eval_ratio, self._hybrid_col,
            self._start_its, self._end_its)

        # Do blocking.
        self._data_len = []
        self._reblock = []
        self._covariance = []
        self._no_opt_block = []
        self._opt_block = []
        for i, dat in enumerate(data):
            dat_c = dat.copy()
            (data_len, reblock, covariance, no_opt_block, opt_block) = (
                self._do_blocking_dat(dat_c, i)
            )
            self._data_len.append(data_len)
            self._reblock.append(reblock)
            self._covariance.append(covariance)
            self._no_opt_block.append(no_opt_block)
            self._opt_block.append(opt_block)

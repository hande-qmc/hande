"""Analyse Monte Carlo correlated output with hybrid analyser."""
from typing import Dict, List, Optional, Tuple, Union
import copy
import pandas as pd
import numpy as np
import statsmodels.tsa.ar_model as ar_model
import statsmodels.tsa.stattools as tsastats
from pyhande.error_analysing.abs_error_analyser import AbsErrorAnalyser
from pyhande.error_analysing.find_starting_iteration import select_find_start
from pyhande.error_analysing.analysis_utils import check_data_input, set_cols


class HybridAna(AbsErrorAnalyser):
    """Analyse ratio observable, such as projected energy.

    Can be used instead of Blocker.

    This scheme is made by hybridizing two different post-analysis
    methods, AR model and Straatsma. The former (the latter) is
    comparatively good at estimating the statistic error for smaller
    (larger) length of time-series, respectively. This method just
    picks up the larger statistic error from the ones given by both
    methods. The mathematical details of both methods are explained
    in (please cite if you use this):

    Ichibha, T., Hongo, K., Maezono, R., Thom, A. J. W., 2019
    arXiv:1904.09934 [physics.comp-ph]
    """

    def __init__(
            self, it_key: str, hybrid_col: str, replica_col: str,
            cols: List[str] = None,
            start_its: Union[List[int], str] = 'mser',
            end_its: List[int] = None, batch_size: int = 1,
            find_start_kw_args: Dict[str, Union[bool, float, int]]
            = None) -> None:
        r"""
        Initialise a HybridAnalyser instance.

        Parameters
        ----------
        it_key : str
            Column name of MC iterations, e.g. 'iterations'.
            Can be set to 'obs:'+key, so it will be set later in
            .exe using passed in observables as
            it_key = observables[key].
        hybrid_col : str
            Column name to be analysed here, e.g. 'Inst. Proj. Energy'.
            Values can be set to 'obs:'+key, so it will be set later in
            .exe using passed in observables as
            hybrid_col = observables[key].
        replica_col : str
            Name of replica columns, e.g. 'replica id'.
        cols : List[str], optional
            Columns in QMC data potentially used when finding start
            iteration (if 'blocking' start_its selected), e.g.
            ['Shift', '\sum H_0j N_j', 'N_0', '# H psips'].
            The default is None, compulsory when
            `start_its` = 'blocking' though.
            All elements can be set to 'obs:'+key, so it will be set
            later in .exe using passed in observables as
            obs_key = observables[key].
        start_its : Union[List[int], str], optional
            Either list of start iterations which has to be of same
            length as data.
            Or string specifying find starting iteration function to use
            to estimate starting iterations automatically, either
            'blocking' (`find_starting_iteration_blocking()`) or
            'mser', (`find_starting_iteration_mser_min()`).  Note that
            this choice then applies to all calculations in passed in
            `data` to `.exe()`.  'blocking' requires `cols` to be set.
            The default is 'mser'.
        end_its : List[int], optional
            List of end iterations. Has to be of same length as
            data. If None, end iteration is the last iteration for each
            calculation.
            The default is None.
        batch_size : int
            The energy time-series is coarse-grained by averaging
            several sequential samples into just one sample and the
            statistic error is calculated for the coarse-grained
            time-series. This variable designates how many sequential
            samples are averaged together.
            The default is 1.
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
        HybridAnalyser._check_input(it_key, cols, hybrid_col, start_its)
        self._input_it_key = it_key
        self._input_hybrid_col = hybrid_col
        self._input_replica_col = replica_col
        self._input_cols = cols
        self._input_start_its = start_its
        if not isinstance(start_its, list):
            self._find_starting_it = select_find_start(start_its)
        self._input_end_its: Optional[List[int]] = end_its
        self._batch_size: int = batch_size
        self._find_start_kw_args: Dict[str, Union[bool, float, int]] = (
            find_start_kw_args if find_start_kw_args else {})

        # Set helper attribute.
        self._pre_exe_error_message = (
            "First do analysis by running 'exe' instance method."
        )

        # These attributes are set later:
        self._it_key: str
        self._cols: List[str]
        self._replica_col: str
        self._hybrid_col: str
        self._start_its: List[int]
        self._end_its: List[int]
        self._opt_block: pd.DataFrame
        self._no_opt_block: List[List[str]]

    @classmethod
    def inst_hande_ccmc_fciqmc(
            cls, start_its: Union[List[int], str] = 'mser',
            end_its: List[int] = None, batch_size: int = 1,
            find_start_kw_args: Dict[str, Union[bool, float, int]] = None
    ):
        """Return HybridAnalyser instance for a HANDE CCMC/FCIQMC calc.

        Parameters
        ----------
        See __init__().

        Returns
        -------
        HybridAnalyser
            Instance of the HybridAnalyser class, customised for a HANDE
            CCMC/FCIQMC calculation.
        """
        return HybridAna(
            'obs:it_key', 'obs:inst_proje_key', 'obs:replica_key',
            ['obs:shift_key', 'obs:sum_key', 'obs:ref_key', 'obs:total_key'],
            start_its=start_its, end_its=end_its, batch_size=batch_size,
            find_start_kw_args=find_start_kw_args)

    @property
    def start_its(self) -> List[int]:
        """Access _start_its attribute if available, else error."""
        try:
            return self._start_its
        except AttributeError:
            print(self._pre_exe_error_message)
            raise

    @property
    def end_its(self) -> List[int]:
        """Access _end_its attribute if available, else error."""
        try:
            return self._end_its
        except AttributeError:
            print(self._pre_exe_error_message)
            raise

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

    @staticmethod
    def _check_input(it_key: str, cols: Optional[List[str]], hybrid_col: str,
                     start_its: Union[List[int], str]):
        """Check some input parameters."""
        if not it_key:
            raise ValueError("'it_key' must be specified!")
        if not cols and start_its == 'blocking':
            raise ValueError("'cols' has to be specified when 'start_its' == "
                             "'blocking', i.e. 'blocking' find starting point "
                             "function is used.")
        if not hybrid_col:
            raise ValueError("'hybrid_col' has to be specified.")

    def _do_hybrid_analysis(self, dat: pd.DataFrame) -> Tuple[
            pd.DataFrame, List[str]]:
        """Do hybrid analysis.

        This implementation is reformatted from the original
        implementation by Tom Ichibha in pyhande.lazy.py.
        """

        ratio_values = dat[self._hybrid_col].values
        n_data = len(ratio_values) // self._batch_size
        means_batch = [0]*n_data
        for i in range(n_data):
            means_batch[i] = np.mean(
                ratio_values[i*self._batch_size:(i+1)*self._batch_size])
        mean = np.mean(means_batch)
        var = np.var(means_batch)
        acf = tsastats.acf(x=means_batch, unbiased=True,
                           nlags=n_data-1, fft=True)

        # ar model
        ar = ar_model.AR(means_batch)
        model_ar = ar.fit(ic='aic', trend='c', method='cmle')
        params = model_ar.params
        denom = nom = 1
        for j in range(len(params)-1):
            denom -= params[j+1]
            nom -= params[j+1] * acf[j+1]
        tau = nom / denom**2
        error_ar = np.sqrt(var/n_data*tau)

        # autocorr
        tau = 1.0
        for i in range(1, n_data-1):
            if acf[i] < 0:
                break
            tau += 2.0*acf[i]
        error_ac = np.sqrt(var/n_data*tau)

        # return value
        error = max(error_ar, error_ac)
        opt_block = pd.DataFrame(
            {'mean': mean,
             'standard error': error,
             'standard error error': None},
            columns=['mean', 'standard error', 'standard error error'],
            index=[self._hybrid_col])
        if self._cols:
            no_opt_block = copy.copy(self._cols)
            if self._hybrid_col in self._cols:
                # [todo] unnecessary?
                no_opt_block.remove(self._hybrid_col)
        else:
            no_opt_block = []
        return (opt_block, no_opt_block)

    def _do_hybrid_dat(self, dat: pd.DataFrame, dat_ind: int
                       ) -> Tuple[pd.DataFrame, List[str]]:
        """Do hybrid analysis for one calc/one replica if replica.

        Parameters
        ----------
        dat : pd.DataFrame
            QMC data for one calculation (of one replica if doing that).
        dat_in : int
            Index of `dat` in all `data` passed to .exe().

        Returns
        -------
        pd.DataFrame, List[str]
            opt_block and no_opt_block of analysis.
        """
        # Find start and end iteration:
        self._end_its.append(
            self._input_end_its[dat_ind] if self._input_end_its
            else dat[self._it_key].iloc[-1]
        )
        self._start_its.append(
            self._input_start_its[dat_ind]
            if (self._input_start_its and
                not isinstance(self._input_start_its, str))
            else self._find_starting_it(dat, self._end_its[-1],
                                        self._it_key, self._cols,
                                        self._hybrid_col,
                                        **self._find_start_kw_args)
        )

        dat = dat[dat[self._it_key].between(
            self._start_its[-1], self._end_its[-1])]
        (opt_block, no_opt_block) = self._do_hybrid_analysis(dat)
        return (opt_block, no_opt_block)

    def exe(self, data: List[pd.DataFrame], observables: Dict[str, str]):
        """
        Do analysis (first finding starting iteration if required).

        Parameters
        ----------
        data : List[pd.DataFrame]
            HANDE calculation Monte Carlo output data.
        observables : Dict[str, str]
            Mapping of column key to column name, e.g. 'ref_key': 'N_0'.
            The default is None.  If any of `it_key`, `cols`,
            `eval_ratio`, `hybrid_col`, `replica_col` were instantiated
            as 'obs:key' to be overwritten with observables['key'],
            observables can't be None and those keys have to be present.

        Raises
        ------
        ValueError
            If not all columns to be blocked appear in 'data'
            or if the length of 'data' is different to length of
            'start_its' or 'end_its' if they are defined.

        """
        # Set column keys.  They might have been specified as 'obs:col'
        # to be set now as observables['col'].  Else, leave as it.
        (self._it_key, self._cols, self._replica_col, _, self._hybrid_col) = (
            set_cols(observables, self._input_it_key,
                     copy.deepcopy(self._input_cols),
                     self._input_replica_col, None, self._input_hybrid_col)
        )

        check_data_input(
            data, self._cols, None, self._hybrid_col, self._input_start_its,
            self._input_end_its)

        # Do analysis.
        self._start_its = []
        self._end_its = []
        self._no_opt_block = []
        self._opt_block = []
        for i, dat in enumerate(data):
            # Select subset of data to analyse:
            dat_c = dat.copy()
            if self._replica_col in dat_c:
                # Have used replica tricks.  Analyse them one by one.
                dat_c_repl = dat_c.groupby([self._replica_col])
                no_opt_block, opt_block = [], []
                for rep_id in range(1, len(dat_c_repl)+1):
                    (opt_bl, no_opt_bl) = (
                        self._do_hybrid_dat(dat_c_repl.get_group(rep_id), i)
                    )
                    no_opt_block.append(no_opt_bl)
                    opt_block.append(opt_bl)
                opt_block = pd.concat(opt_block, keys=dat_c_repl.groups.keys(),
                                      names=[self._replica_col])
            else:
                (opt_block, no_opt_block) = self._do_hybrid_dat(dat_c, i)
            self._no_opt_block.append(no_opt_block)
            self.opt_block.append(opt_block)

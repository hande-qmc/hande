"""Analyse Monte Carlo correlated output with hybrid analyser."""
from typing import Dict, List, Tuple, Union
import copy
import pandas as pd
import numpy as np
import statsmodels.tsa.ar_model as ar_model
import statsmodels.tsa.stattools as tsastats
from pyhande.error_analysing.abs_error_analyser import AbsErrorAnalyser
from pyhande.error_analysing.find_starting_iteration import select_find_start


class HybridAnalyser(AbsErrorAnalyser):
    """Analyse ratio observable, such as projected energy.

    Can be used instead of Blocker.

    This scheme is made by hybridizing two different post-analysis
    methods, AR model and Straatsma. The former (the latter) is
    comparatively good at estimating the statistic error for smaller
    (larger) length of time-series, respectively. This method just
    picks up the larger statistic error from the ones given by both
    methods. The mathematical details of both methods are explained
    in (please cite if you use this):

    Tom Ichibha, Kenta Hongo, Ryo Maezono, Alex J.W. Thom (2019),
    arXiv:1904.09934 [physics.comp-ph]
    """

    def __init__(
            self, cols: List[str] = None, eval_ratio: Dict[str, str] = None,
            it_key: str = 'iterations',
            start_its: Union[List[int], str] = 'mser',
            end_its: List[int] = None, batch_size: int = 1,
            find_start_kw_args: Dict[str, Union[bool, float, int]]
            = None) -> None:
        r"""
        Initialise a HybridAnalyser instance.

        Parameters
        ----------
        cols : List[str], optional
            Columns in QMC data potentially used when finding start
            iteration (if 'blocking' start_its selected).
            The default is
            ['Shift', '\sum H_0j N_j', 'N_0', '# H psips'].
        eval_ratio: Dict, optional
            Evaluate mean and standard error of ratio of column 'num'
            with column 'denom'. Add this to the analysis results as
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
            start_its = 'analysis'.
            The default is None, which defaults to an empty dictionary.
        """
        self._cols: List[str] = cols if cols else [
            'Shift', r'\sum H_0j N_j', 'N_0', '# H psips'
        ]
        self._eval_ratio: Dict[str, str] = eval_ratio if eval_ratio else {
            'name': 'Proj. Energy', 'num': r'\sum H_0j N_j', 'denom': 'N_0'
        }
        self._it_key = it_key
        if isinstance(start_its, list):
            self._start_its = start_its
        else:
            self._start_its = None
            self._find_starting_it = select_find_start(start_its)
        self._end_its: List[int] = end_its
        self._batch_size: int = batch_size
        self._find_start_kw_args: Dict[str, Union[bool, float, int]] = (
            find_start_kw_args if find_start_kw_args else {})
        # These attributes are set later:
        self._opt_block: pd.DataFrame
        self._no_opt_block: List[List[str]]

    def _do_hybrid_analysis(self, dat: pd.DataFrame) -> Tuple[
            pd.DataFrame, List[str]]:
        """Do hybrid analysis.

        This implementation is reformatted from the original
        implementation by Tom Ichibha in pyhande.lazy.py.
        """

        ratio_values = (dat[self.eval_ratio['num']] /
                        dat[self.eval_ratio['denom']].values)
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
            index=[self.eval_ratio['name']])
        no_opt_block = copy.copy(self._cols)
        if self.eval_ratio['name'] in self._cols:
            # [todo] unnecessary?
            no_opt_block.remove(self.eval_ratio['name'])
        return (opt_block, no_opt_block)

    def exe(self, data: List[pd.DataFrame]):
        """
        Do analysis (first finding starting iteration if required).

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

        # Find end and start iteration if required.
        if not self.end_its:
            self._end_its = [dat['iterations'].iloc[-1] for dat in data]
        if not self.start_its:
            self._start_its = [
                self._find_starting_it(
                    dat, end_it, self._it_key, self._cols, self.eval_ratio,
                    **self._find_start_kw_args)
                for dat, end_it in zip(data, self.end_its)]

        # Do analysis.
        self._no_opt_block = []
        self._opt_block = []
        for dat, start_it, end_it in zip(data, self.start_its,
                                         self.end_its):
            # Select subset of data to analyse:
            dat_c = dat.copy()
            dat_c = dat_c[dat_c['iterations'].between(start_it, end_it)]
            (opt_block, no_opt_block) = self._do_hybrid_analysis(dat_c)
            self._no_opt_block.append(no_opt_block)
            self.opt_block.append(opt_block)

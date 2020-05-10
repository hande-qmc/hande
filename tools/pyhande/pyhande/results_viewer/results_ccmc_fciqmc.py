"""Access and investigate CCMC/FCIQMC results from HANDE QMC."""
from typing import Dict, List, Union
import copy
import pandas as pd
import matplotlib.pyplot as plt
from pyhande.data_preparing.hande_ccmc_fciqmc import PrepHandeCcmcFciqmc
from pyhande.extracting.extractor import Extractor
from pyhande.error_analysing.blocker import Blocker
from pyhande.error_analysing.hybrid_analyser import HybridAnalyser
from pyhande.results_viewer.results import Results
import pyhande.analysis as analysis


class ResultsCcmcFciqmc(Results):
    """Show CCMC and FCIQMC HANDE results and allow further analysis."""

    def __init__(
            self, extractor: Extractor,
            preparator: PrepHandeCcmcFciqmc = None,
            analyser: Union[Blocker, HybridAnalyser] = None) -> None:
        """
        Initialise `ResultsCcmcFciqmc` instance.

        `extractor.exe()` (`preparator.exe()` and `analyser.exe()`) have
        already be executed before creating an instance of this class.

        Parameters
        ----------
        extractor : Extractor
            Extractor instance which has extracted HANDE QMC data.
        preparator : PrepHandeCcmcFciqmc, optional
            If present, contains prepared data (e.g. dealing with
            complex, replica tricks and adding an inst. projected energy
            column) which might have then be passed to the `analyser`.
        analyser : Union[Blocker, HybridAnalyser], optional
            If present, information on Blocker or HybridAnalyser.
        """
        super().__init__(extractor)
        self._preparator = preparator
        self._analyser: Union[Blocker, HybridAnalyser] = analyser
        if analyser.opt_block:
            self.summary = self._opt_block()
        # To be set later:
        self._shoulder: pd.DataFrame
        self._inefficiency: pd.DataFrame

    @property
    def preparator(self) -> PrepHandeCcmcFciqmc:
        """Access preparator used to prepare data for analysis."""
        return self._preparator

    @property
    def analyser(self) -> Union[Blocker, HybridAnalyser]:
        """Access analyser used to supply the analysed results."""
        return self._analyser

    @staticmethod
    def _concat_reset_rename(df: List[pd.DataFrame]) -> pd.DataFrame:
        """Concat df from diff calcs, reset index, rename new cols.

        Parameters
        ----------
        df : List[pd.DataFrame]
            List of dataframes, length of list is the number of QMC
            calculations.

        Returns
        -------
        pd.DataFrame
            Concatenated dataframe from df with reset index.  Will then
            have column 'calc id' corresponding to calculation id and
            columns 'observables' which used to be the index of each
            element of passed in df.
        """
        df = pd.concat(df, keys=list(range(len(df))), names=['calc id'])
        df = df.reset_index()
        return df.rename(
            columns={'level_1': 'observable', 'level_2': 'observable',
                     'mean': 'value/mean'}
        )  # Only one level will be used, depending on whether replica.

    def _opt_block(self) -> pd.DataFrame:
        return ResultsCcmcFciqmc._concat_reset_rename(
            self.analyser.opt_block)

    def _calc_shoulder(self, dat: pd.DataFrame) -> pd.DataFrame:
        """Calculate shoulder using analysis.plateau_estimator."""
        return analysis.plateau_estimator(
            dat, total_key=self.preparator.observables['total_key'],
            ref_key=self.preparator.observables['ref_key'],
            shift_key=self.preparator.observables['shift_key'])

    @property
    def shoulder(self) -> pd.DataFrame:
        """Access shoulder. For now, not hist shoulder [todo]."""
        try:
            return self._shoulder
        except AttributeError:
            self._shoulder = []
            for dat in self.preparator.data:
                if self.preparator.observables['replica_key'] in dat:
                    datg = dat.groupby(
                        [self.preparator.observables['replica_key']])
                    shoulder = [
                        self._calc_shoulder(datg.get_group(replica_ind))
                        for replica_ind in range(1, len(datg)+1)
                    ]
                    shoulder = pd.concat(
                        shoulder, keys=list(range(1, len(datg)+1)),
                        names=[self.preparator.observables['replica_key']])
                else:
                    shoulder = self._calc_shoulder(dat)
                self._shoulder.append(shoulder)
            self._shoulder = ResultsCcmcFciqmc._concat_reset_rename(
                self._shoulder)
            return self._shoulder

    def add_shoulder(self):
        """Add shoulder to summary. [todo]: allow hist shoulder."""
        self.summary = pd.concat([self.summary, self.shoulder])
        self.summary.sort_values(
            by=['calc id', self.preparator.observables['replica_key']],
            inplace=True, ignore_index=True)

    @property
    def inefficiency(self) -> pd.DataFrame:
        """Access inefficiency."""
        try:
            return self._inefficiency
        except AttributeError:
            if isinstance(self._analyser, Blocker):
                ineffs = []
                for dat_ind in range(len(self.extractor.data)):
                    opt_block = self.analyser.opt_block[dat_ind]
                    dtau = self.extractor.metadata[dat_ind][0]['qmc']['tau']
                    n_its = (self.analyser.end_its[dat_ind]
                             - self.analyser.start_its[dat_ind])
                    ineff = analysis.inefficiency(opt_block, dtau, n_its)
                    ineffs.append(ineff)
                # Idea for reshaping from ../../reblock_hande.py.
                self._inefficiency = pd.concat(
                    [pd.DataFrame(ineff.stack()).T for ineff in ineffs],
                    ignore_index=True)
                return self._inefficiency
            raise AttributeError("Inefficiency cannot be evaluated. "
                                 "Blocker info required.")

    def add_inefficiency(self):
        """Add inefficiency to summary."""
        self.summary = pd.concat([self.summary, self.inefficiency], axis=1)
        return self

    def plot_shoulder(
            self, inds: List[int] = None, show_shoulder: bool = True,
            log_scale: bool = True,
            show: bool = True, savefig: str = None) -> None:
        """Plot shoulder.

        Parameters
        ----------
        inds : List[int]
            Indices of calculations to plot. If None, plot all.
            The default is None.
        show_shoulder : bool
            Show positions of shoulder height with vertical lines.
        log_scale : bool
            Set x and y axis on log scale.
        show : bool
            If True, show plot. The default is True.
        savefig : str
            If not None, save fig to specified path 'savefig'.
            The default is None.
        """
        if not inds:
            inds = list(range(len(self.extractor.data)))
        fig = plt.figure()
        # See https://matplotlib.org/3.2.1/gallery/animation/
        # double_pendulum_sgskip.html#sphx-glr-gallery-animation-double-
        # pendulum-sgskip-py
        shou_plot = fig.add_subplot()
        # [todo] Improve
        max_ratio = -1000000
        min_ratio = 1000000
        for ind in inds:
            tot_part_ref_part = (self.extractor.data[ind]['# H psips'] /
                                 self.extractor.data[ind]['N_0'])
            if max_ratio < max(tot_part_ref_part):
                max_ratio = max(tot_part_ref_part)
            if min_ratio > min(tot_part_ref_part):
                min_ratio = min(tot_part_ref_part)
            shou_plot.plot(self.extractor.data[ind]['# H psips'],
                           tot_part_ref_part, color='C'+str(ind),
                           label=str(ind))
        if show_shoulder:
            for ind in inds:
                shou_plot.vlines(
                    self.shoulder['shoulder height']['mean'].iloc[ind],
                    0.9*min_ratio, 1.1*max_ratio, color='C'+str(ind),
                    linestyle='--')
        shou_plot.set_xlabel(r'# Particles')
        shou_plot.set_ylabel(r'# Particles/# Particles on Reference')
        if log_scale:
            shou_plot.set_xscale('log')
            shou_plot.set_yscale('log')
        if show:
            plt.show()
        if savefig:
            # [todo]
            pass

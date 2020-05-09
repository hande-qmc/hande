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
            self.summary = self._add_opt_block()
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

    def _add_opt_block(self) -> pd.DataFrame:
        opt_block = pd.concat(
            self.analyser.opt_block,
            keys=list(range(len(self.analyser.opt_block))), names=['calc id'])
        opt_block.reset_index(inplace=True)
        return opt_block.rename(
            columns={'level_1': 'observables', 'level_2': 'observables',
                     'mean': 'value/mean'}
        )  # only one will be used, depending on if replica.

    @property
    def shoulder(self) -> pd.DataFrame:
        """Access shoulder. For now, not hist shoulder [todo]."""
        try:
            return self._shoulder
        except AttributeError:
            shoulders = [
                analysis.plateau_estimator(dat) for dat in self.extractor.data
            ]
            # Idea for reshaping from ../../reblock_hande.py.
            self._shoulder = pd.concat(
                [pd.DataFrame(shoulder.stack()).T for shoulder in shoulders],
                ignore_index=True)
            return self._shoulder

    def add_shoulder(self):
        """Add shoulder to summary. [todo]: allow hist shoulder."""
        self.summary = pd.concat([self.summary, self.shoulder], axis=1)

    @property
    def inefficiency(self) -> pd.DataFrame:
        """Access inefficiency."""
        try:
            return self._inefficiency
        except AttributeError:
            if isinstance(self._analyser, Blocker):
                ineffs = []
                for ind in range(len(self.extractor.data)):
                    opt_block = self.analyser.opt_block[ind]
                    dtau = self.extractor.metadata[ind]['qmc']['tau']
                    # [todo] - This is one iteration more than in lazy!
                    dat_c = self.extractor.data[ind].copy()
                    dat_c = dat_c[dat_c['iterations'].between(
                        self.analyser.start_its[ind] + 1,
                        self.analyser.end_its[ind])]
                    its = (dat_c['iterations'].iloc[-1] -
                           dat_c['iterations'].iloc[0])
                    ineff = analysis.inefficiency(opt_block, dtau, its)
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

"""Access and investigate CCMC/FCIQMC results from HANDE QMC."""
from typing import Dict, List, Union
import warnings
import copy
import pandas as pd
import matplotlib.pyplot as plt
import pyblock
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

    @property
    def summary_pretty(self) -> pd.DataFrame:
        """Access self._summary but prettify for viewing data.

        Combine value in "value/mean" column with "standard error"
        columns for easy viewing, e.g. '0.123(4)'.  If not possible,
        due to type or not present values, fill in value in
        "value/mean".

        Returns
        -------
        pd.DataFrame
            Prettified summary table for viewing (not further analysis).
        """
        pretty_fmt_errs = []
        for sum_ind in range(len(self.summary)):
            try:
                pretty_fmt_err = pyblock.error.pretty_fmt_err(
                    self._summary.iloc[sum_ind]['value/mean'],
                    self._summary.iloc[sum_ind]['standard error'])
            except (ValueError, TypeError):
                # Not all observables have an error, e.g. some metadata.
                pretty_fmt_err = self._summary.iloc[sum_ind]['value/mean']
            pretty_fmt_errs.append(pretty_fmt_err)
        summary_pretty = pd.concat([self.summary, pd.DataFrame(
            {'pretty': pretty_fmt_errs})], axis=1)
        if self.preparator.observables['replica_key'] in summary_pretty:
            summary_pretty.set_index(
                ['calc id', self.preparator.observables['replica_key']],
                inplace=True)
            n_calcs = len(self.preparator.data)
            return pd.concat(
                [summary_pretty.loc[i_calc].pivot(columns='observable',
                                                  values='pretty')
                 for i_calc in range(n_calcs)], keys=list(range(n_calcs)),
                names=['calc id'])
        else:
            summary_pretty.set_index(['calc id'], inplace=True)
            return summary_pretty.pivot(columns='observable', values='pretty')

    def compare_obs(self, observables: List[str]) -> pd.DataFrame:
        """Compare observables from .summary where obs are columns.

        Parameters
        ----------
        observables : List[str]
            Observables from .summary to compare.

        Returns
        -------
        pd.DataFrame
            DataFrame where easier comparisons are possible.
        """
        comp_df = self.summary.query("observable in @observables")
        if self.preparator.observables['replica_key'] in comp_df:
            # replica tricks used
            comp_df.set_index(
                ['calc id', self.preparator.observables['replica_key']],
                inplace=True)
            n_calcs = len(self.preparator.data)
            return pd.concat(
                [comp_df.loc[i_calc].pivot(columns='observable', values=[
                    'value/mean', 'standard error'])
                 for i_calc in range(n_calcs)], keys=list(range(n_calcs)),
                names=['calc id'])
        else:
            comp_df.set_index(['calc id'], inplace=True)
            comp_df = comp_df.pivot(columns='observable',
                                    values=['value/mean', 'standard error'])
        return comp_df

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
        """Return opt_block in more normalised form."""
        return ResultsCcmcFciqmc._concat_reset_rename(
            self.analyser.opt_block)

    def _add_to_summary(self, df: pd.DataFrame) -> None:
        """Add data to summary."""
        for obs in df['observable']:
            if obs in self.summary['observable'].values:
                df = df.query('observable != @obs')
                warnings.warn("Add attempt failed: summary already "
                              f"contains {obs}.")
        if df.empty:
            # all metadata items already were in summary!
            return
        self.summary = pd.concat([self.summary, df])
        try:
            self.summary.sort_values(
                by=['calc id', self.preparator.observables['replica_key']],
                inplace=True, ignore_index=True)
        except KeyError:  # no replica tricks
            self.summary.sort_values(
                by=['calc id'], inplace=True, ignore_index=True)

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
            # Calculate shoulder.
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
        self._add_to_summary(self.shoulder)

    def _calc_ineff(self, dat_ind: int, opt_block: pd.DataFrame) -> float:
        """Calculate inefficiency using analysis.inefficiency."""
        dtau = self.extractor.metadata[dat_ind][-1]['qmc']['tau']
        n_its = (self.analyser.end_its[dat_ind]
                 - self.analyser.start_its[dat_ind])
        ineff = analysis.inefficiency(
            opt_block, dtau, n_its,
            sum_key=self.preparator.observables['sum_key'],
            ref_key=self.preparator.observables['ref_key'],
            total_key=self.preparator.observables['total_key'],
            proje_key=self.preparator.observables['proje_key'])
        return ineff

    @property
    def inefficiency(self) -> pd.DataFrame:
        """Access inefficiency."""
        try:
            return self._inefficiency
        except AttributeError:
            # Calculate inefficiency.
            self._inefficiency = []
            for dat_ind in range(len(self.preparator.data)):
                if (self.preparator.observables['replica_key'] in
                        self.preparator.data[dat_ind]):
                    opt_block = self.analyser.opt_block[dat_ind]
                    ineff = [
                        self._calc_ineff(dat_ind, opt_block.loc[replica_ind])
                        for replica_ind in range(1, opt_block.index.max()[0]+1)
                    ]
                    ineff = pd.concat(
                        ineff, keys=list(range(1, opt_block.index.max()[0]+1)),
                        names=[self.preparator.observables['replica_key']])
                else:
                    ineff = self._calc_ineff(
                        dat_ind, self.analyser.opt_block[dat_ind])
                self._inefficiency.append(ineff)
            self._inefficiency = ResultsCcmcFciqmc._concat_reset_rename(
                self._inefficiency)
            return self._inefficiency

    def add_inefficiency(self):
        """Add inefficiency to summary."""
        self._add_to_summary(self.inefficiency)

    def add_metadata(self, meta_keys: List[str]):
        """Overwritten version of Results.add_metadata.

        Parameters
        ----------
        meta_keys : List[str]
            List of metadata to add in strings where different level
            keys are separated by colons. E.g.
            ['qmc:tau', 'system:ueg:r_s'] adds
            extractor.metadata[:]['qmc']['tau'] as well as
            extractor.metadata[:]['system']['ueg']['r_s'] to summary
            (if they exist).
        """
        metadata = self.get_metadata(meta_keys)
        replica_key = self.preparator.observables['replica_key']
        if (replica_key in self.preparator.data[0]):
            # replica tricks. Assume if doing replica tricks, doing it for all.
            # Just duplicated metadata for each replica.
            metadata = pd.concat(
                [metadata]*int(self.summary[replica_key].max()),
                keys=list(range(1, int(self.summary[replica_key].max())+1)),
                names=[replica_key])
            metadata.reset_index(inplace=True)
            metadata.drop(columns=['level_1'], inplace=True)
        self._add_to_summary(metadata)

    def plot_shoulder(
            self, inds: List[int] = None, show_shoulder: bool = True,
            log_scale: bool = True) -> None:
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
        """
        if not inds:
            inds = list(range(len(self.preparator.data)))
        fig = plt.figure()
        # See https://matplotlib.org/3.2.1/gallery/animation/
        # double_pendulum_sgskip.html#sphx-glr-gallery-animation-double-
        # pendulum-sgskip-py
        shou_plot = fig.add_subplot()
        # [todo] Improve
        max_ratio = -1000000
        min_ratio = 1000000
        repl_id = self.preparator.observables['replica_key']
        for ind in inds:
            if repl_id in self.preparator.data[ind]:
                # doing replica tricks, only show first replica
                warnings.warn("Only showing shoulder plot of first replica.")
                data = self.preparator.data[ind].groupby(repl_id).get_group(1)
            else:
                data = self.preparator.data[ind]
            ref_key = self.preparator.observables['ref_key']
            total_key = self.preparator.observables['total_key']
            tot_part_ref_part = data[total_key] / data[ref_key]
            if max_ratio < max(tot_part_ref_part):
                max_ratio = max(tot_part_ref_part)
            if min_ratio > min(tot_part_ref_part):
                min_ratio = min(tot_part_ref_part)
            shou_plot.plot(data[total_key], tot_part_ref_part,
                           color='C'+str(ind), label=str(ind))
        if show_shoulder:
            if repl_id in self.shoulder:
                shoulder = self.shoulder.groupby(repl_id).get_group(1)
            else:
                shoulder = self.shoulder
            shoulder_h = (
                shoulder[shoulder['observable'] == 'shoulder height']
            )
            for ind in inds:
                shou_plot.vlines(
                    shoulder_h[shoulder_h['calc id'] == ind]['value/mean'],
                    0.9*min_ratio, 1.1*max_ratio, color='C'+str(ind),
                    linestyle='--')
        shou_plot.set_xlabel(total_key)
        shou_plot.set_ylabel(total_key+"/"+ref_key)
        if log_scale:
            shou_plot.set_xscale('log')
            shou_plot.set_yscale('log')
        plt.legend(loc='best')
        plt.show()

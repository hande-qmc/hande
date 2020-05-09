"""Access and investigate generic results from HANDE QMC."""
from typing import Dict, List
import warnings
import pandas as pd
from pyhande.extracting.extractor import Extractor


class Results:
    """Show and allow investigation of HANDE QMC results.

    Extraction has already happened.
    This is a base class, used for now for all non CCMC and non
    FCIQMC calculations who use a more specific class.
    """

    def __init__(self, extractor: Extractor) -> None:
        """
        Initialise Results instance.

        Parameters
        ----------
        extractor : Extractor
            Extractor instance which has extracted HANDE QMC data.
        """
        self._extractor: Extractor = extractor
        self.summary: pd.DataFrame = pd.DataFrame()

    @property
    def extractor(self) -> Extractor:
        """Access extractor used to supply these results."""
        return self._extractor

    @property
    def summary(self) -> pd.DataFrame:
        """Access summary."""
        return self._summary

    @summary.setter
    def summary(self, summary) -> None:
        """Set summary which has to be a pandas DataFrame."""
        if not isinstance(summary, pd.DataFrame):
            raise TypeError("Cannot set summary. It has to be a pd.DataFrame.")
        self._summary = summary

    @staticmethod
    def _access_meta_and_check(calc_ind: int, metadata: List[Dict],
                               meta_key: str):
        """Access value of metadata, returning None if non-existent.

        Parameters
        ----------
        calc_ind : int
            Index of calculation.
        metadata : List[Dict]
            Metadata for one calculation.  Can contain more than one
            dictionary, if calculation was merged from multiple
            calculations.
        meta_key : str
            'keyOuter:keyInner:...', e.g. ['qmc:tau', 'system:ueg:r_s']
            adds extractor.metadata[:]['qmc']['tau'] as well as
            extractor.metadata[:]['system']['ueg']['r_s'].

        Returns
        -------
        Tuple or single entry
            Tuple if requested metadata entries differ in metadata,
            otherwise collapsed to single value of metadata value.
            None if value not found.
        """
        same_calc_meta_list = []
        for metadat in metadata:
            try:
                for key in meta_key.split(':'):
                    metadat = metadat[key]
                same_calc_meta_list.append(metadat)
            except KeyError:
                warnings.warn(
                    f"Metadata for calc #{calc_ind} has no key {meta_key}.")
                same_calc_meta_list.append(None)
        # Is this a safe comparison for all datatypes?
        if all([same_calc_meta_list[i] == same_calc_meta_list[0]
                for i in range(len(same_calc_meta_list))]):
            return same_calc_meta_list[0]
        return tuple(same_calc_meta_list)

    def get_metadata(self, meta_keys: List[str]) -> pd.DataFrame:
        """Get part(s) of metadata in pandas DataFrame.

        Parameters
        ----------
        meta_keys : List[str]
            List of metadata items to put into DataFrame.  Each item as
            'keyOuter:keyInner:...', e.g. ['qmc:tau', 'system:ueg:r_s']
            adds extractor.metadata[:]['qmc']['tau'] as well as
            extractor.metadata[:]['system']['ueg']['r_s'].

        Returns
        -------
        pd.DataFrame
            Contains metadata requested for all calculations.
        """
        return pd.DataFrame([
            [self._access_meta_and_check(ind, metadata, meta_key)
             for meta_key in meta_keys]
            for ind, metadata in enumerate(self.extractor.metadata)],
            columns=[meta_key.split(':')[-1] for meta_key in meta_keys])

    def add_metadata(self, meta_keys: List[str]):
        """Add metadata to summary.

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
        self.summary = pd.concat(
            [self.summary, self.get_metadata(meta_keys)], axis=1)

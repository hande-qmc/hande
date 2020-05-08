"""Abstract base class for error analyser objects."""
import abc
from typing import Dict, List
import pandas as pd


class AbsErrorAnalyser(metaclass=abc.ABCMeta):
    """Define required/common attributes for error analyser objects."""

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
            print("First do analysis by running 'exe' instance method.")
            raise

    @property
    def no_opt_block(self) -> List[List[str]]:
        """Access _no_opt_block attribute if available. Else error."""
        try:
            return self._no_opt_block
        except AttributeError:
            print("First do analysis by running 'exe' instance method.")
            raise

    @abc.abstractmethod
    def exe(self):
        """Do error analysis."""

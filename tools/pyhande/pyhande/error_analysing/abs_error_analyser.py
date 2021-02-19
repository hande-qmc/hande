"""Abstract base class for error analyser objects."""
import abc
from typing import List, Dict
import pandas as pd


class AbsErrorAnalyser(metaclass=abc.ABCMeta):
    """Define required/common attributes for error analyser objects."""

    @property
    @abc.abstractmethod
    def start_its(self) -> List[int]:
        """Access _start_its attribute, analysis start iterations."""

    @property
    @abc.abstractmethod
    def end_its(self) -> List[int]:
        """Access _end_its attribute, analysis end iterations."""

    @property
    @abc.abstractmethod
    def opt_block(self) -> pd.DataFrame:
        """Access _opt_block attribute if available. Else error."""

    @property
    @abc.abstractmethod
    def no_opt_block(self) -> List[List[str]]:
        """Access _no_opt_block attribute if available. Else error."""

    @abc.abstractmethod
    def exe(self, data: List[pd.DataFrame], observables: Dict[str, str]
            ) -> None:
        """Do error analysis."""

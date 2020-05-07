"""Abstract base class for mapping and preparing data."""
import abc
from typing import Dict, List
import pandas as pd


class AbsDataPreparator(metaclass=abc.ABCMeta):
    """Prepare and map data columns."""

    @property
    @abc.abstractmethod
    def observables(self) -> Dict[str, str]:
        pass

    @property
    @abc.abstractmethod
    def data(self) -> List[pd.DataFrame]:
        pass

    @abc.abstractmethod
    def exe(self):
        pass

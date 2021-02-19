"""Abstract base class for extractor objects."""
import abc
from typing import List


class AbsExtractor(metaclass=abc.ABCMeta):
    """Define required attributes for extractor objects."""

    @property
    @abc.abstractmethod
    def data(self):
        """Should contain data in pandas DataFrame format."""

    @property
    @abc.abstractmethod
    def metadata(self):
        """Should contain a dictionary with metadata."""

    @abc.abstractmethod
    def exe(self, out_files: List[str]):
        """Do the extraction and possibly merging of calcs."""

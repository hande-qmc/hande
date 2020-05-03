"""Abstract base class for extractor objects."""
import abc


class AbsExtractor(metaclass=abc.ABCMeta):
    """Define required attributes for extractor objects."""

    @property
    @abc.abstractmethod
    def data(self):
        pass

    @property
    @abc.abstractmethod
    def metadata(self):
        pass

    @abc.abstractmethod
    def exe(self):
        pass

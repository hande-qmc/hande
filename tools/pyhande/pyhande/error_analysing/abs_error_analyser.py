"""Abstract base class for error analyser objects."""
import abc


class AbsErrorAnalyser(metaclass=abc.ABCMeta):
    """Define required attributes for error analyser objects."""

    @property
    @abc.abstractmethod
    def cols(self):
        """Columns of observables analysed."""
        pass

    @property
    @abc.abstractmethod
    def start_its(self):
        """Start iterations for analysis."""
        pass

    @property
    @abc.abstractmethod
    def end_its(self):
        """End iterations for analysis."""
        pass

    @property
    @abc.abstractmethod
    def opt_block(self):
        """Observables mean and standard error (error)."""
        pass

    @property
    @abc.abstractmethod
    def no_opt_block(self):
        """Observables not (successfully) error analysed."""
        pass

    @abc.abstractmethod
    def exe(self):
        """Do error analysis."""
        pass

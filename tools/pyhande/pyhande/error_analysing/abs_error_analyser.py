"""(Abstract) base class for error analyser objects."""
import abc
from typing import Dict, List
import pandas as pd


class AbsErrorAnalyser(metaclass=abc.ABCMeta):
    """Define required/common attributes for error analyser objects."""

    def __init__(self):
        # This is never called.  Just to keep linter happy.
        self._cols = None
        self._eval_ratio = None
        self._start_its = None
        self._end_its = None
        self._opt_block = None
        self._no_opt_block = None
        self._find_starting_it = lambda x: x
        self._find_start_kw_args = None
        self._it_key = None

    @property
    def eval_ratio(self) -> Dict[str, str]:
        """Access _eval_ratio attribute, ratio to evaluate."""
        return self._eval_ratio

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
        pass

    def _check_data_input(self, data):
        """Check data input."""
        if not all([col in dat for dat in data for col in self._cols]):
            raise ValueError("'cols' parameter must only contain columns "
                             "names present in all dataframes in 'data'.")
        if not all([col in dat for dat in data for col in
                    [self.eval_ratio['num'], self.eval_ratio['denom']]]):
            raise ValueError("Keys in 'eval_ratio' must only contain columns "
                             "names present in all dataframes in 'data'.")
        if (isinstance(self.start_its, list) and
                len(self.start_its) != len(data)):
            raise ValueError("If 'start_its' (here of length "
                             f"{len(self.start_its)}) is specified as list of "
                             "start iterations, it has to be a list of the "
                             "same length as 'data' (here of "
                             f"length {len(data)}).")
        if self.end_its and len(self.end_its) != len(data):
            raise ValueError("If 'end_its' (here of length "
                             f"{len(self.end_its)}) is specified, it has to "
                             "be a list of the same length as 'data' (here of "
                             f"length {len(data)}).")

    def _set_start_and_end_its(self, data):
        """Find end and start iteration if required."""
        if not self.end_its:
            self._end_its = [dat['iterations'].iloc[-1] for dat in data]
        if not self.start_its:
            self._start_its = [
                self._find_starting_it(
                    dat, end_it, self._it_key, self._cols, self.eval_ratio,
                    **self._find_start_kw_args)
                for dat, end_it in zip(data, self.end_its)]

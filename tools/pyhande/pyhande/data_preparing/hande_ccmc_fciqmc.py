"""CCMC and FCIQMC HANDE data preparation for analysis."""
import copy
from typing import Dict, List
import pandas as pd
from pyhande.data_preparing.abs_data_preparator import AbsDataPreparator


class PrepareHandeCcmcFciqmc(AbsDataPreparator):
    """Prepare HANDE CCMC/FCIQMC data for analysis."""

    def __init__(self):
        """Initialise an instance preparing HANDE CCMC/FCIQMC data."""
        self._observables: Dict[str, str] = {
            'it_key': 'iterations',
            'shift_key': 'Shift',
            'sum_key': r'\sum H_0j N_j',
            'ref_key': 'N_0',
            'total_key': '# H psips',
            'proje_key': 'Proj. Energy'
        }
        self._not_extracted_message: str = (
            'Data not prepared yet. Use the ".exe" method for preparing.')

        # Set by .exe()
        self._data: List[pd.DataFrame]
        self._metadata: List[List[Dict]]
        self._complex_data: bool
        self._replica_tricks: bool

    @property
    def observables(self) -> Dict[str, str]:
        """Access observables, key mapping."""
        return self._observables

    @property
    def data(self) -> List[pd.DataFrame]:
        """Access (prepared) data.

        Raises
        ------
        AttributeError
            If data has not been prepared yet.

        Returns
        -------
        List[pd.DataFrame]
            QMC Data.
            Cleaned list over merged calculations.
        """
        try:
            return self._data
        except AttributeError:
            print(self._not_extracted_message)
            raise

    @property
    def metadata(self) -> List[List[Dict]]:
        """Access (prepared) metadata.

        Raises
        ------
        AttributeError
            If metadata has not been prepared yet.

        Returns
        -------
        List[List[Dict]]
            Copy of metadata.  Outer list over (merged) calculations,
            inner list over calculations in different output files when
            merged.
        """
        try:
            return self._metadata
        except AttributeError:
            print(self._not_extracted_message)
            raise

    @property
    def complex_data(self) -> bool:
        """True if data is complex.

        Raises
        ------
        AttributeError
            If preparation has not been done yet.

        Returns
        -------
        bool
        """
        try:
            return self._complex_data
        except AttributeError:
            print(self._not_extracted_message)
            raise

    @property
    def replica_data(self) -> bool:
        """True if replica tricks were used.

        Raises
        ------
        AttributeError
            If preparation has not been done yet.

        Returns
        -------
        bool
        """
        try:
            return self._replica_data
        except AttributeError:
            print(self._not_extracted_message)
            raise

    def _gen_complex_cols(self, col_key: str) -> Dict[str, str]:
        """Parses columns key in, returns real and imag col keys.

        Parameters
        ----------
        col_key : str
            Name of column.

        Returns
        -------
        Dict[str, str]
            Real and imaginary version of `col_key`.
        """
        return {'real': 'Re{'+col_key+'}', 'imag': 'Im{'+col_key+'}'}

    def _add_complex_obs(self):
        """Add (neg) magnitude of ref and sum to data.

        Take magnitude of real and imaginary part of `sum_key` and
        `ref_key` respectively, whose ratio is then an estimator of the
        projected energy and add that to data.  Update list of
        observables to replace previous ref and sum keys with
        magnitudes.

        Credits of first implementation to Charlie Scott.
        Note a comment by Charlie Scott in an ealier version (albeit not
        in the master branch):
        "We just use the magnitude information to obtain all estimates.
        This won't detect if the projected energy isn't real, but this
        is nontrivial to do. If we consider the phase of the reference
        population/projection through the Hamiltonian onto the
        reference, we would expect this to stabilise but nothing in our
        dynamics provides an effective restoring force if it deviates.
        As such we can expect the phase to gradually random walk, so
        extracting any information via reblocking is difficult. We could
        instead attempt to reblock the ratio of the instantaneous phases
        to provide an estimate of the error in the overall phase, as we
        would expect this to be randomly distributed around zero. This
        will however be at best a biased estimator."
        """
        comp_sum_key = self._gen_complex_cols(self.observables['sum_key'])
        comp_ref_key = self._gen_complex_cols(self.observables['ref_key'])
        mag_sum_key = "-|"+self.observables['sum_key']+"|"
        mag_ref_key = "|"+self.observables['ref_key']+"|"
        for i in range(len(self._data)):
            sum_mag_neg = pd.DataFrame(
                - (self._data[i][comp_sum_key['real']]**2 +
                   self._data[i][comp_sum_key['imag']]**2),
                columns=[mag_sum_key])
            ref_mag = pd.DataFrame(
                (self._data[i][comp_ref_key['real']]**2 +
                 self._data[i][comp_ref_key['imag']]**2),
                columns=[mag_ref_key])
            self._data[i] = pd.concat(
                [self._data[i], sum_mag_neg, ref_mag], axis=1)
        self._observables['sum_key'] = mag_sum_key
        self._observables['ref_key'] = mag_ref_key

    def _add_inst_proje(self):
        """Add inst. projected energy, maybe required for analysis.

        Note that its mean is not the true mean of the projected energy
        which is the mean of `sum_key` over `ref_key`.
        """
        for i in range(len(self._data)):
            inst_proje = pd.DataFrame(
                (self.data[i][self.observables['sum_key']] /
                 self.data[i][self.observables['ref_key']]),
                columns=['Inst. '+self.observables['proje_key']])
            self._data[i] = pd.concat([self._data[i], inst_proje], axis=1)

    def exe(self, data: List[pd.DataFrame], metadata: List[List[Dict]],
            make_copy: bool = True):
        """Prepare data; deal with complex, replica and add inst. proje.

        Parameters
        ----------
        data: List[pd.DataFrame]
            List of output data.  Should be all of same type
            (complex/non-complex, replica-tricks/no replica tricks,
            calc_type).
        metadata: List[List[Dict]]
            metadata.  Outer list over (merged) calculations, inner
            list over parts of calculation in different output files
            when merged.
        make_copy: bool, optional
            If true, deepcopy data and metadata so that passed in data
            and metadata are not altered by any changes here.
            The default is True.
        """
        if make_copy:
            self._data = copy.deepcopy(data)
            self._metadata = copy.deepcopy(metadata)
        else:
            self._data = data
            self._metadata = metadata

        # Complex?
        if (self._gen_complex_cols(self.observables['ref_key'])['real']
                in data[0]):
            print("complex!", self._gen_complex_cols(
                self.observables['ref_key'])['real'])
            # Probably complex calculations!
            # Test that all are complex.
            if not all(
                    self._gen_complex_cols(self.observables['ref_key'])['real']
                    in dat for dat in data
            ):
                raise ValueError("Some but not all data is complex! Either "
                                 "pass data where all calculations are either "
                                 "complex or not complex.")
            self._add_complex_obs()
            self._complex_data = True
        else:
            # Not complex.
            self._complex_data = False

        # Calculate instantaneous projected energy, required by some
        # analysis schemes (e.g. hybrid).
        self._add_inst_proje()

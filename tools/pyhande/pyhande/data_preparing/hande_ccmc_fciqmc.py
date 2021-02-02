"""CCMC and FCIQMC HANDE data preparation for analysis."""
import copy
from typing import Dict, List
import pandas as pd
from pyhande.data_preparing.abs_data_preparator import AbsDataPreparator


class PrepHandeCcmcFciqmc(AbsDataPreparator):
    """Prepare HANDE CCMC/FCIQMC data for analysis."""

    def __init__(self):
        """Initialise an instance preparing HANDE CCMC/FCIQMC data."""
        self._observables_init: Dict[str, str] = {
            'it_key': 'iterations',
            'shift_key': 'Shift',
            'sum_key': r'\sum H_0j N_j',
            'ref_key': 'N_0',
            'total_key': '# H psips',
            'proje_key': 'Proj. Energy',
            'inst_proje_key': 'Inst. Proj. Energy',
            'replica_key': 'replica id'
        }
        self._not_extracted_message: str = (
            'Data not prepared yet. Use the ".exe" method for preparing.')

        # Set by .exe()
        self._observables: Dict[str, str]
        self._data: List[pd.DataFrame]
        self._complex_data: bool
        self._replica_tricks: bool

    @property
    def observables(self) -> Dict[str, str]:
        """Access observables, key mapping.

        Raises
        ------
        AttributeError
            If data has not been prepared yet.

        Returns
        -------
        Dict[str,str]
            Map of observables property to column/observables name.
        """
        try:
            return self._observables
        except AttributeError:
            print(self._not_extracted_message)
            raise

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

    @staticmethod
    def _replica_ending(replica_id: int) -> str:
        """Defines replica's ending.

        Parameters
        ----------
        replica_id : int
            ID of replica, i.e. first one has 0, next 1, etc.

        Returns
        -------
        str
            replica ending
        """
        return "_"+str(replica_id)

    @staticmethod
    def _curly_wrapper(string_to_add: str, col_name: str) -> str:
        """Wrap `col_name` in `string_to_add{col_name}`.

        Parameters
        ----------
        string_to_add : str
            String before curly brackets.
        col_name : str
            Name of column.

        Returns
        -------
        str
            Name of wrapped column.
        """
        return string_to_add + "{" + col_name + "}"

    def _re_cols(self, col_key: str) -> str:
        """Wrap `col_key` in `Re{self.observables[col_key]}`.

        Parameters
        ----------
        col_key : str
            Name of column.

        Returns
        -------
        str
            Name of wrapped column.
        """
        return PrepHandeCcmcFciqmc._curly_wrapper(
            'Re', self.observables[col_key])

    def _im_cols(self, col_key: str) -> str:
        """Wrap `col_key` in `Im{self.observables[col_key]}`.

        Parameters
        ----------
        col_key : str
            Name of column.

        Returns
        -------
        str
            Name of wrapped column.
        """
        return PrepHandeCcmcFciqmc._curly_wrapper(
            'Im', self.observables[col_key])

    def _add_replica_grouping(self, max_replica_id: int) -> None:
        """Add extra column for replica id which can be grouped by.

        Parameters
        ----------
        max_replica_id : int
            Highest id of replica.  Assume that replica ids are
            consecutive.  With two replicas, ids are 1 and 2, so
            `max_replica_id` would be 2.
        """
        for i in range(len(self.data)):
            replica_dats = []
            for replica_id in range(1, max_replica_id+1):
                dat = pd.concat([
                    pd.DataFrame({
                        self.observables['replica_key']:
                            len(self.data[i])*[replica_id]
                    }), self.data[i]
                ], axis=1)
                cols_to_drop = [
                    col
                    for rep_id in list(range(1, max_replica_id+1))
                    for col in dat.columns
                    if (rep_id != replica_id and col.endswith(
                        PrepHandeCcmcFciqmc._replica_ending(rep_id)))
                ]
                dat.drop(columns=cols_to_drop, inplace=True)
                dat.rename(
                    columns={
                        col: col[:-len(
                            PrepHandeCcmcFciqmc._replica_ending(replica_id))
                        ] for col in dat.columns if col.endswith(
                            PrepHandeCcmcFciqmc._replica_ending(replica_id))
                    }, inplace=True
                )
                replica_dats.append(dat)
            self._data[i] = pd.concat(replica_dats, ignore_index=True)

    def _check_add_replica(self) -> None:
        """Check if replica exist, if yes, add replica col and group.

        Replica columns contain the index of the replica (see
        `gen_replica_col`).  Here, reduce this column back to original
        form but each row will only contain one's replica's observables.
        That replica is identified in the extra replica column.

        i.e. df with two columns `cols_1` and `cols_2` and one row
        gets changed to df with two columns `replica id` and `cols` with
        two rows, one for each of the two replicas.
        """
        # If replica:
        # Add extra column for i so that replicas can be groupedby.
        if any(col.endswith(PrepHandeCcmcFciqmc._replica_ending(1))
               for dat in self.data for col in dat.columns):
            # Find highest replica_id.
            max_id = 0
            while (any(
                col.endswith(PrepHandeCcmcFciqmc._replica_ending(max_id+1))
                for dat in self.data for col in dat.columns
            )):
                max_id += 1
            # Check that this highest replica_id exist in all data.
            if not (all(
                any(col.endswith(PrepHandeCcmcFciqmc._replica_ending(max_id))
                    for col in dat.columns) for dat in self.data
            )):
                raise ValueError(f"Some but not all data have {max_id} "
                                 "replicas! Make sure that number of replicas "
                                 "is identical!")
            self._add_replica_grouping(max_id)
            self._replica_data = True
        else:
            self._replica_data = False

    def _add_complex_obs(self) -> None:
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
        mag_sum_key = "-|"+self.observables['sum_key']+"|"
        mag_ref_key = "|"+self.observables['ref_key']+"|"
        for i in range(len(self._data)):
            sum_mag_neg = pd.DataFrame(
                -(self._data[i][self._re_cols('sum_key')]**2
                  + self._data[i][self._im_cols('sum_key')]**2)**0.5,
                columns=[mag_sum_key])
            ref_mag = pd.DataFrame(
                (self._data[i][self._re_cols('ref_key')]**2
                 + self._data[i][self._im_cols('ref_key')]**2)**0.5,
                columns=[mag_ref_key])
            self._data[i] = pd.concat(
                [self._data[i], sum_mag_neg, ref_mag], axis=1)
        self._observables['sum_key'] = mag_sum_key
        self._observables['ref_key'] = mag_ref_key

    def _check_add_complex(self) -> None:
        """Check whether complex, if yes, add magnitude.

        Complex calculations have ref and sum columns in Re{key},
        Im{key} format. If complex, take (negative) magnitudes and add
        those as columns to each self._data elements respectively.
        """
        if any(self._re_cols('ref_key') in dat for dat in self.data):
            # Probably complex calculations!
            # Test that all are complex.
            if not all(self._re_cols('ref_key') in dat for dat in self.data):
                raise ValueError("Some but not all data is complex! Either "
                                 "pass data where all calculations are either "
                                 "complex or not complex.")
            self._add_complex_obs()
            self._complex_data = True
        else:
            # Not complex.
            self._complex_data = False

    def _add_inst_proje(self) -> None:
        """Add inst. projected energy, maybe required for analysis.

        Note that its mean is not the true mean of the projected energy
        which is the mean of `sum_key` over `ref_key`.
        """
        for i in range(len(self._data)):
            inst_proje = pd.DataFrame(
                (self.data[i][self.observables['sum_key']] /
                 self.data[i][self.observables['ref_key']]),
                columns=[self.observables['inst_proje_key']])
            self._data[i] = pd.concat([self._data[i], inst_proje], axis=1)

    def exe(self, data: List[pd.DataFrame], make_copy: bool = True):
        """Prepare data; deal with complex, replica and add inst. proje.

        Parameters
        ----------
        data: List[pd.DataFrame]
            List of output data.  Should be all of same type
            (complex/non-complex, replica-tricks/no replica tricks,
            calc_type).
        make_copy: bool, optional
            If true, deepcopy data so that passed in data are not
            altered by any changes here.
            The default is True.
        """
        # Copy observables so that a new run of .exe can start fresh.
        self._observables = copy.deepcopy(self._observables_init)
        self._data = copy.deepcopy(data) if make_copy else data

        # Replicas?
        # i.e. do we have some repeated columns, same name expect for
        # "_i" at the end where i is in list(range(1,#replicas+1))?
        # Then add a column denoting the index of the replica, remove
        # quasi-duplicate columns (only differing by replica id) so that
        # each row only contains results for one replica. This implies
        # that the number of rows in each dataframe will double if there
        # are two replicas for example.
        self._check_add_replica()

        # Complex?
        # i.e. are the ref and sum columns in Re{key}, Im{key}?
        # Then take (negative) magnitudes and add those as columns.
        self._check_add_complex()

        # Calculate instantaneous projected energy, ratio of sum to key
        # columns, required by some analysis schemes (e.g. hybrid).
        self._add_inst_proje()

"""Extract and merge (meta)data from (multiple) HANDE output files."""
from typing import Dict, List
import pandas as pd
import pyhande.extract as extract
import pyhande.lazy as lazy


class Extractor:
    """Extract data/metadata from HANDE output files and merge.

    Merge if desired/sensible, e.g. when calculation was restarted.
    This expands the functionality of extract.py and is more compactly
    represented as a class.
    """

    def __init__(self, out_files: List[str], merge: str = 'uuid') -> None:
        """
        Initialise instance.

        Parameters
        ----------
        out_files : List[str]
            List of HANDE output filenames to be extracted here.
        merge : str, optional
            * 'uuid': Do merge based on UUIDs of calculations.
                Restarted calculation has UUID of calculation restarted
                from in metadata. `out_files` do not need to be ordered
                in that case. Default.
            * 'legacy': Do merge if iterations in consecutive
                calculations in `out_files` are consecutive. Order in
                `out_files` matters therefore. Mainly here for legacy
                reasons when restart UUID was not saved in metadata yet.
            * 'no': Don't merge.
                Treat calculations in different files in `out_files` as
                separate calculations.
            The default is 'uuid'.
            [todo] For now, 'legacy' and 'uuid' do the same thing - Fix
        """
        self._out_files: List[str] = out_files
        if merge not in ['uuid', 'legacy', 'no']:
            raise ValueError('Invalid merge parameter value "{}". Choose from '
                             '"uuid", "legacy" and "no".'.format(merge))
        self._merge: str = merge
        self._not_extracted_message: str = (
            '(Meta)Data has not been extracted yet. Use the ".exe" method to '
            'start extracting.')
        # The following will be set later (by "exe" attribute method):
        self._data: List[pd.DataFrame]
        self._metadata: List[Dict]
        self._all_ccmc_fciqmc: bool

    @property
    def out_files(self) -> List[str]:
        """
        Access (read only) out_files property.

        Returns
        -------
        List[str]
            List of `out_files` names the data is extracted from.
        """
        return self._out_files

    @property
    def data(self) -> List[pd.DataFrame]:
        """
        Access (extracted) data property.

        Raises
        ------
        AttributeError
            If data has not been extracted yet.

        Returns
        -------
        List[pd.DataFrame]
            QMC Data.
            List over merged calculations.
        """
        try:
            return self._data
        except AttributeError:
            print(self._not_extracted_message)
            raise

    @property
    def metadata(self) -> List[Dict]:
        """
        Access (extracted) metadata property.

        Raises
        ------
        AttributeError
            If metadata has not been extracted yet.

        Returns
        -------
        List[Dict]
            Metadata.
            List over merged calculations.
            [todo] - update so that not only one calc's metadata
            [todo] - survives upon merging! Make it a List[List[Dict]]!

        """
        try:
            return self._metadata
        except AttributeError:
            print(self._not_extracted_message)
            raise

    @property
    def all_ccmc_fciqmc(self) -> bool:
        """
        Are all calculations extracted either CCMC or FCIQMC.

        This will affect what postprocessing can be done.

        Raises
        ------
        AttributeError
            If data has not been extracted yet.

        Returns
        -------
        bool
            True if all calculations extracted are either CCMC or
            FCIQMC. False if at least one is of another type, such as
            FCI or Hilbert space estimation.

        """
        try:
            return self._all_ccmc_fciqmc
        except AttributeError:
            print(self._not_extracted_message)
            raise

    def exe(self):
        """Extract and merge.

        [todo] -  Don't use lazy for concatentation! At the moment
        merge = 'uuid' and 'legacy' do the same.
        """
        # Extract.
        raw_data_metadata = extract.extract_data_sets(self._out_files)
        self._metadata = [rawdm[0] for rawdm in raw_data_metadata]
        self._data = [rawdm[1] for rawdm in raw_data_metadata]

        # Merging?
        # Are all calculations either CCMC or FCIQMC?  Only then we merge.
        calc_types = [md['calc_type'] for md in self.metadata]
        self._all_ccmc_fciqmc = all(
            [ct in ['FCIQMC', 'CCMC'] for ct in calc_types])
        if self._merge != 'no' and self._all_ccmc_fciqmc:
            # Calculations are potentially merged.
            # [todo] Don't use lazy.  Keep all metadata when merging.
            # [todo] And do more checks when merging.
            self._metadata, self._data = lazy.concat_calcs(
                self._metadata, self._data)

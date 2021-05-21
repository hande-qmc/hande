"""Test extracting.extractor.Extractor."""
import unittest
import pandas as pd
from pyhande.extracting.extractor import Extractor


class TestExtractingProperties(unittest.TestCase):
    """Test extracting.extractor.Extractor.

    More detailed tests on extracting in ../test_extract.py
    """

    def setUp(self):
        """Setting up."""
        self._extractor = Extractor()

    def test_basic_extracting(self):
        """Testing a basic extraction."""
        self._extractor.exe(["tests/hande_files/ueg.out"])
        # data
        self.assertEqual(len(self._extractor.data), 1)
        # Sample a row.
        pd.testing.assert_series_equal(
            self._extractor.data[0].iloc[1],
            pd.Series([2.00000000e+01, -1.10017479e-01, -8.20941670e-01,
                       3.00000000e+00, 1.28200000e+03,  7.82000000e+02,
                       4.63000000e+02,  7.62000000e-01, 8.00000000e-04],
                      index=['iterations', 'Shift', r'\sum H_0j N_j', 'N_0',
                             '# H psips', '# states', '# spawn_events',
                             'R_spawn', 'time'], name=1))
        # metadata
        self.assertEqual(len(self._extractor.metadata), 1)
        self.assertEqual(len(self._extractor.metadata[0]), 1)
        # Sample a piece.
        self.assertEqual(self._extractor.metadata[0][0]['qmc']['tau'], 0.1)
        # out_files
        self.assertListEqual(self._extractor.out_files,
                             ["tests/hande_files/ueg.out"])
        # calc_to_outfile_ind
        self.assertListEqual(self._extractor.calc_to_outfile_ind, [[0]])
        # all_ccmc_fciqmc
        self.assertTrue(self._extractor.all_ccmc_fciqmc)

    def test_undefined_property_access(self):
        """Testing undefined property access (before .exe())."""
        # These should all raise AttributeError since .exe() has not
        # been executed yet.
        with self.assertRaises(AttributeError):
            self._extractor.data
        with self.assertRaises(AttributeError):
            self._extractor.metadata
        with self.assertRaises(AttributeError):
            self._extractor.out_files
        with self.assertRaises(AttributeError):
            self._extractor.calc_to_outfile_ind
        with self.assertRaises(AttributeError):
            self._extractor.all_ccmc_fciqmc

    def test_mix_type_extraction(self):
        """Testing extraction of FCIQMC and FCI calc."""
        self._extractor.exe(["tests/hande_files/ueg.out",
                             "tests/hande_files/fci_ueg.out"])
        # data
        self.assertEqual(len(self._extractor.data), 2)
        # metadata
        self.assertEqual(len(self._extractor.metadata), 2)
        self.assertEqual(len(self._extractor.metadata[0]), 1)
        # out_files
        self.assertListEqual(self._extractor.out_files,
                             ["tests/hande_files/ueg.out",
                              "tests/hande_files/fci_ueg.out"])
        # calc_to_outfile_ind
        self.assertListEqual(self._extractor.calc_to_outfile_ind, [[0], [1]])
        # all_ccmc_fciqmc
        self.assertFalse(self._extractor.all_ccmc_fciqmc)

    def test_same_type_extraction(self):
        """Testing extraction of two different FCIQMC calcs."""
        self._extractor.exe(["tests/hande_files/ueg.out",
                             "tests/hande_files/ueg2.out"])
        # data
        self.assertEqual(len(self._extractor.data), 2)
        # metadata
        self.assertEqual(len(self._extractor.metadata), 2)
        self.assertEqual(len(self._extractor.metadata[0]), 1)
        # out_files
        self.assertListEqual(self._extractor.out_files,
                             ["tests/hande_files/ueg.out",
                              "tests/hande_files/ueg2.out"])
        # calc_to_outfile_ind
        self.assertListEqual(self._extractor.calc_to_outfile_ind, [[0], [1]])
        # all_ccmc_fciqmc
        self.assertTrue(self._extractor.all_ccmc_fciqmc)


class TestMerging(unittest.TestCase):
    """Test extracting.extractor.Extractor, focus on merging.

    [todo] - test 'it_key' and 'shift_key' merge options.
    """

    def test_uuid_merge(self):
        """Merging by UUID."""
        extractor = Extractor()  # UUID merging is default.
        # First and third passed output file are to be merged
        # (the first being a continuation of the third).
        # This merged data should be same as data from second passed in
        # calculation (except possibly for the 'time' column).
        extractor.exe(['tests/hande_files/long_calc_split2_ueg.out',
                       'tests/hande_files/long_calc_ueg.out',
                       'tests/hande_files/long_calc_split1_ueg.out'])
        pd.testing.assert_frame_equal(extractor.data[0].drop(columns='time'),
                                      extractor.data[1].drop(columns='time'),
                                      check_exact=True)
        self.assertListEqual(extractor.calc_to_outfile_ind, [[1], [2, 0]])

    def test_legacy_merge(self):
        """Merging by continuation check of iter. of neighbour calcs."""
        extractor = Extractor(merge={'type': 'legacy'})
        # First and second passed output file are to be merged.
        # (Note that order matters with legacy merging - unlike UUID.)
        # This merged data should be same as data from third passed in
        # calculation (except possibly for the 'time' column).
        extractor.exe(['tests/hande_files/long_calc_split1_ueg.out',
                       'tests/hande_files/long_calc_split2_ueg.out',
                       'tests/hande_files/long_calc_ueg.out'])
        pd.testing.assert_frame_equal(extractor.data[0].drop(columns='time'),
                                      extractor.data[1].drop(columns='time'),
                                      check_exact=True)
        self.assertListEqual(extractor.calc_to_outfile_ind, [[0, 1], [2]])

    def test_no_merge(self):
        """Don't merge."""
        extractor = Extractor(merge={'type': 'no'})
        extractor.exe(['tests/hande_files/long_calc_split1_ueg.out',
                       'tests/hande_files/long_calc_split2_ueg.out'])
        self.assertEqual(len(extractor.data), 2)
        self.assertListEqual(extractor.calc_to_outfile_ind, [[0], [1]])

    def test_merge_md_shift_restriction_positive(self):
        """Merge since `qmc:quasi_newton` did not change."""
        extractor = Extractor(merge={'md_shift': 'qmc:quasi_newton'})
        # Note that shift has starting varying between split1 and split2.
        extractor.exe(['tests/hande_files/long_calc_split1_ueg.out',
                       'tests/hande_files/long_calc_split2_ueg.out',
                       'tests/hande_files/long_calc_ueg.out'])
        pd.testing.assert_frame_equal(extractor.data[0].drop(columns='time'),
                                      extractor.data[1].drop(columns='time'),
                                      check_exact=True)
        self.assertListEqual(extractor.calc_to_outfile_ind, [[0, 1], [2]])

    def test_merge_md_shift_restriction_negative(self):
        """Don't merge since `cpu_time` did change."""
        extractor = Extractor(merge={'md_shift': 'cpu_time'})
        # Note that shift has starting varying between split1 and split2.
        with self.assertWarns(UserWarning):
            extractor.exe(['tests/hande_files/long_calc_split1_ueg.out',
                           'tests/hande_files/long_calc_split2_ueg.out'])
            self.assertEqual(len(extractor.data), 2)
            self.assertListEqual(extractor.calc_to_outfile_ind, [[0], [1]])

    def test_merge_md_always_restriction_negative(self):
        """Don't merge since `cpu_time` did change."""
        extractor = Extractor(merge={'md_always': 'cpu_time'})
        with self.assertWarns(UserWarning):
            extractor.exe(['tests/hande_files/long_calc_split1_ueg.out',
                           'tests/hande_files/long_calc_split2_ueg.out'])
            self.assertEqual(len(extractor.data), 2)
            self.assertListEqual(extractor.calc_to_outfile_ind, [[0], [1]])

    def test_merge_md_shift_restriction_shift_not_varying(self):
        """Merge - shift hasn't varied, though `UUID` did change."""
        extractor = Extractor(merge={'md_shift': 'UUID'})
        # Note that shift has not starting varying between split1 and split2.
        # Doing UUID merge, so order of output files can be random, "*split2*"
        # being the continuation of "*split1*".
        extractor.exe([
            'tests/hande_files/short_calc_no_shift_varying_split1_ueg.out',
            'tests/hande_files/short_calc_no_shift_varying_ueg.out',
            'tests/hande_files/short_calc_no_shift_varying_split2_ueg.out'
        ])
        pd.testing.assert_frame_equal(extractor.data[0].drop(columns='time'),
                                      extractor.data[1].drop(columns='time'),
                                      check_exact=True)
        self.assertListEqual(extractor.calc_to_outfile_ind, [[0, 2], [1]])

    def test_merge_md_always_restriction_shift_not_varying(self):
        """Don't merge, `UUID` change and always indep. of shift."""
        extractor = Extractor(merge={'md_always': 'UUID'})
        # Note that shift has not starting varying between split1 and split2.
        # Doing UUID merge, so order of output files can be random, "*split2*"
        # being the continuation of "*split1*".
        with self.assertWarns(UserWarning):
            extractor.exe([
                'tests/hande_files/short_calc_no_shift_varying_split1_ueg.out',
                'tests/hande_files/short_calc_no_shift_varying_ueg.out',
                'tests/hande_files/short_calc_no_shift_varying_split2_ueg.out'
            ])
            self.assertEqual(len(extractor.data), 3)
            self.assertListEqual(extractor.calc_to_outfile_ind,
                                 [[0], [1], [2]])

    def test_invalid_merge_type(self):
        """Don't merge since `cpu_time` did change."""
        extractor = Extractor(merge={'type': 'error'})
        with self.assertRaises(ValueError):
            extractor.exe(['tests/hande_files/long_calc_split1_ueg.out',
                           'tests/hande_files/long_calc_split2_ueg.out'])

    def test_invalid_metadata(self):
        """Pass in invalid piece of metadata."""
        extractor = Extractor(merge={'md_always': 'error'})
        with self.assertRaises(KeyError):
            extractor.exe(['tests/hande_files/long_calc_split1_ueg.out',
                           'tests/hande_files/long_calc_split2_ueg.out'])

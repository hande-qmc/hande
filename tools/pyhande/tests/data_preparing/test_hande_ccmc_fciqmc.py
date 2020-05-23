"""Test data_preparing.hande_ccmc_fciqmc.PrepHandeCcmcFciqmc."""
import unittest
import copy
import numpy as np
import pandas as pd
from pyhande.data_preparing.hande_ccmc_fciqmc import PrepHandeCcmcFciqmc
from tests.create_mock_df import create_qmc_frame


class TestAccessPropertiesPreExe(unittest.TestCase):
    """Test data_preparing.hande_ccmc_fciqmc.PrepHandeCcmcFciqmc"""

    def test_undefined_property_access(self):
        """Testing undefined property access (before .exe())."""
        # These should all raise AttributeError since .exe() has not
        # been executed yet.
        self._preparator = PrepHandeCcmcFciqmc()
        with self.assertRaises(AttributeError):
            self._preparator.observables
        with self.assertRaises(AttributeError):
            self._preparator.data
        with self.assertRaises(AttributeError):
            self._preparator.complex_data
        with self.assertRaises(AttributeError):
            self._preparator.replica_data


class TestExe(unittest.TestCase):
    """Test data_preparing.hande_ccmc_fciqmc.PrepHandeCcmcFciqmc.exe"""

    def setUp(self):
        """Setting up."""
        # Standard, non complex non replica.
        self._preparator = PrepHandeCcmcFciqmc()
        rng = np.random.default_rng(3910)
        cols = ['Shift', 'test', 'N_0', r'\sum H_0j N_j', 'alt']
        means = [-0.1, 10.0, 11.0, -9.0, 20.0]
        sine_periods = list(range(1, len(cols)+1))
        noise_facs = [0.1*mean for mean in means]
        self._df_mock = create_qmc_frame(rng, cols, means, sine_periods,
                                         noise_facs)

        # Replica but non complex.
        cols = ['Shift_1', 'Shift_2', 'test', 'N_0_1', 'N_0_2',
                r'\sum H_0j N_j_1', r'\sum H_0j N_j_2', 'alt_1', 'alt_2']
        means = [-0.1, -10.0, 11.0, 9.0, 20.0, -2.0, -100.0, 0.1, -23.0]
        sine_periods = list(range(1, len(cols)+1))
        noise_facs = [0.1*mean for mean in means]
        self._df_mock_replica = create_qmc_frame(rng, cols, means,
                                                 sine_periods, noise_facs)

        # Complex but no replica tricks.
        cols = ['Shift', 'test', 'Re{N_0}', 'testi', 'Im{N_0}',
                r'Re{\sum H_0j N_j}', r'Im{\sum H_0j N_j}', 'alt']
        means = [-0.1, -10.0, 11.0, 9.0, 20.0, -2.0, -100.0, 0.1]
        sine_periods = list(range(1, len(cols)+1))
        noise_facs = [0.1*mean for mean in means]
        self._df_mock_comp = create_qmc_frame(rng, cols, means,
                                              sine_periods, noise_facs)

        # Complex and replica tricks.
        cols = ['Shift', 'test', 'Re{N_0}_2', 'Re{N_0}_1', 'testi',
                'Im{N_0}_1', 'Im{N_0}_2',
                r'Re{\sum H_0j N_j}_1', r'Im{\sum H_0j N_j}_1', 'alt_1',
                r'Re{\sum H_0j N_j}_2', r'Im{\sum H_0j N_j}_2', 'alt_2']
        means = [-0.1, -10.0, 11.0, 9.0, 20.0, 1, 2, -2.0, -100.0, 0.1, -0.9,
                 -32.3, 2.53]
        sine_periods = list(range(1, len(cols)+1))
        noise_facs = [0.1*mean for mean in means]
        self._df_mock_comp_replica = create_qmc_frame(rng, cols, means,
                                                      sine_periods, noise_facs)

    def test_single_calc_non_complex_non_replica(self):
        """Tests calculation that is not complex and has no replica."""
        observables_copy = copy.deepcopy(self._preparator._observables_init)
        self._preparator.exe([self._df_mock])
        inst_proje = pd.DataFrame(
            {'Inst. Proj. Energy':
                self._df_mock[r'\sum H_0j N_j']/self._df_mock['N_0']
             })
        df_exp = pd.concat([self._df_mock, inst_proje], axis=1)
        self.assertDictEqual(observables_copy, self._preparator.observables)
        self.assertDictEqual(observables_copy,
                             self._preparator._observables_init)
        self.assertEqual(len(self._preparator.data), 1)
        self.assertFalse(self._preparator.complex_data)
        self.assertFalse(self._preparator.replica_data)
        pd.testing.assert_frame_equal(self._preparator.data[0], df_exp)
        # and an explicit number test:
        self.assertEqual(
            self._preparator.data[0]['Inst. Proj. Energy'].iloc[17],
            -0.9408831928478553)

    def test_single_calc_non_complex_replica(self):
        """Tests calculation that is not complex and but has replica."""
        observables_copy = copy.deepcopy(self._preparator._observables_init)
        self._preparator.exe([self._df_mock_replica])
        self.assertDictEqual(observables_copy, self._preparator.observables)
        self.assertDictEqual(observables_copy,
                             self._preparator._observables_init)
        self.assertEqual(len(self._preparator.data), 1)
        self.assertFalse(self._preparator.complex_data)
        self.assertTrue(self._preparator.replica_data)
        # Test a row.
        pd.testing.assert_series_equal(
            self._preparator.data[0].iloc[23],
            pd.Series([1.0, -0.10855666764042851, 8.404875958869997,
                       9.709934275877028, -2.021602334081873,
                       0.12271035279340871, -0.20819938391388093],
                      index=['replica id', 'Shift', 'test', 'N_0',
                             r'\sum H_0j N_j', 'alt', 'Inst. Proj. Energy'],
                      name=23))

    def test_single_calc_complex_non_replica(self):
        """Tests complex calculation with no replica tricks."""
        observables_copy = copy.deepcopy(self._preparator._observables_init)
        self._preparator.exe([self._df_mock_comp])
        self.assertDictEqual(observables_copy,
                             self._preparator._observables_init)
        observables_copy['ref_key'] = '|N_0|'
        observables_copy['sum_key'] = r'-|\sum H_0j N_j|'
        self.assertDictEqual(observables_copy, self._preparator.observables)
        self.assertEqual(len(self._preparator.data), 1)
        self.assertTrue(self._preparator.complex_data)
        self.assertFalse(self._preparator.replica_data)
        # Test a row.
        pd.testing.assert_series_equal(
            self._preparator.data[0].iloc[8],
            pd.Series([-0.0965111601859209, -9.18311932648693,
                       13.036315180456127, 10.898009962430494,
                       17.06201092747081, -2.280271892256964,
                       -99.79911877252728, 0.11322312194464978,
                       -99.82516590357174, 21.472254897269735,
                       -4.6490304060364345],
                      index=['Shift', 'test', 'Re{N_0}', 'testi', 'Im{N_0}',
                             r'Re{\sum H_0j N_j}', r'Im{\sum H_0j N_j}', 'alt',
                             r'-|\sum H_0j N_j|', '|N_0|',
                             'Inst. Proj. Energy'],
                      name=8))

    def test_single_calc_complex_replica(self):
        """Tests complex calculation with replica tricks."""
        observables_copy = copy.deepcopy(self._preparator._observables_init)
        self._preparator.exe([self._df_mock_comp_replica])
        self.assertDictEqual(observables_copy,
                             self._preparator._observables_init)
        observables_copy['ref_key'] = '|N_0|'
        observables_copy['sum_key'] = r'-|\sum H_0j N_j|'
        self.assertDictEqual(observables_copy, self._preparator.observables)
        self.assertEqual(len(self._preparator.data), 1)
        self.assertTrue(self._preparator.complex_data)
        self.assertTrue(self._preparator.replica_data)
        # Test a row.
        pd.testing.assert_series_equal(
            self._preparator.data[0].iloc[43],
            pd.Series([2.0, -0.09926766491509308, -8.102674301521086,
                       12.180642349662337, 20.208054064807346,
                       1.8520272792994423, -0.920897999551087,
                       -39.24994243422395, 2.650080260365803,
                       -39.26074418825337, 12.320635255280342,
                       -3.186584406954756],
                      index=['replica id', 'Shift', 'test', 'Re{N_0}', 'testi',
                             'Im{N_0}', r'Re{\sum H_0j N_j}',
                             r'Im{\sum H_0j N_j}', 'alt', r'-|\sum H_0j N_j|',
                             '|N_0|', 'Inst. Proj. Energy'],
                      name=43))

    def test_multi_calcs_non_complex_non_replica(self):
        """Tests 2x calculation, not complex and no replica."""
        observables_copy = copy.deepcopy(self._preparator._observables_init)
        self._preparator.exe([self._df_mock, self._df_mock])
        inst_proje = pd.DataFrame(
            {'Inst. Proj. Energy':
                self._df_mock[r'\sum H_0j N_j']/self._df_mock['N_0']
             })
        df_exp = pd.concat([self._df_mock, inst_proje], axis=1)
        self.assertDictEqual(observables_copy, self._preparator.observables)
        self.assertDictEqual(observables_copy,
                             self._preparator._observables_init)
        self.assertEqual(len(self._preparator.data), 2)
        self.assertFalse(self._preparator.complex_data)
        self.assertFalse(self._preparator.replica_data)
        pd.testing.assert_frame_equal(self._preparator.data[0], df_exp)
        pd.testing.assert_frame_equal(self._preparator.data[1], df_exp)

    def test_multi_calcs_complex_replica(self):
        """Tests 2x complex calculation with replica tricks."""
        observables_copy = copy.deepcopy(self._preparator._observables_init)
        self._preparator.exe(2*[self._df_mock_comp_replica])
        self.assertDictEqual(observables_copy,
                             self._preparator._observables_init)
        observables_copy['ref_key'] = '|N_0|'
        observables_copy['sum_key'] = r'-|\sum H_0j N_j|'
        self.assertDictEqual(observables_copy, self._preparator.observables)
        self.assertEqual(len(self._preparator.data), 2)
        self.assertTrue(self._preparator.complex_data)
        self.assertTrue(self._preparator.replica_data)
        # Test a row.
        row_exp = pd.Series([2.0, -0.09926766491509308, -8.102674301521086,
                             12.180642349662337, 20.208054064807346,
                             1.8520272792994423, -0.920897999551087,
                             -39.24994243422395, 2.650080260365803,
                             -39.26074418825337, 12.320635255280342,
                             -3.186584406954756],
                            index=['replica id', 'Shift', 'test', 'Re{N_0}',
                                   'testi', 'Im{N_0}', r'Re{\sum H_0j N_j}',
                                   r'Im{\sum H_0j N_j}', 'alt',
                                   r'-|\sum H_0j N_j|', '|N_0|',
                                   'Inst. Proj. Energy'],
                            name=43)
        pd.testing.assert_series_equal(self._preparator.data[0].iloc[43],
                                       row_exp)
        pd.testing.assert_series_equal(self._preparator.data[1].iloc[43],
                                       row_exp)

    def test_mix_calc_standard_and_replica(self):
        """Tests standard calc and replic calc passed in."""
        # Don't allow mixing.
        with self.assertRaises(ValueError):
            self._preparator.exe([self._df_mock, self._df_mock_replica])

    def test_mix_calc_standard_and_complex(self):
        """Tests standard calc and complex calc passed in."""
        # Don't allow mixing.
        with self.assertRaises(ValueError):
            self._preparator.exe([self._df_mock, self._df_mock_comp])

    def test_unchanged_mutable(self):
        """Test that data passed doesn't change when make_copy (exe)."""
        df_mock_comp_replica_copy = self._df_mock_comp_replica.copy()
        data_list = [self._df_mock_comp_replica]
        data_list_copy = [df_mock_comp_replica_copy]
        self._preparator.exe(data_list)
        self.assertEqual(len(data_list), len(data_list_copy))
        for dat, dat_copy in zip(data_list, data_list_copy):
            pd.testing.assert_frame_equal(dat, dat_copy, check_exact=True)


class TestStaticMethods(unittest.TestCase):
    """Test PrepHandeCcmcFciqmc's static methods."""

    def test_replica_ending(self):
        """Test PrepHandeCcmcFciqmc._replica_ending."""
        self.assertEqual("_14", PrepHandeCcmcFciqmc._replica_ending(14))

    def test_replica_ending_negative(self):
        """Test PrepHandeCcmcFciqmc._replica_ending, neg. int."""
        self.assertEqual("_-1", PrepHandeCcmcFciqmc._replica_ending(-1))

    def test_curly_wrapper(self):
        """Test PrepHandeCcmcFciqmc._curly_wrapper."""
        self.assertEqual("testi{blubb}",
                         PrepHandeCcmcFciqmc._curly_wrapper("testi", "blubb"))

    def test_curly_wrapper_empty_strs(self):
        """Test PrepHandeCcmcFciqmc._curly_wrapper, empty str inputs."""
        self.assertEqual("{}", PrepHandeCcmcFciqmc._curly_wrapper("", ""))

    def test_curly_wrapper_raw_strings(self):
        """Test PrepHandeCcmcFciqmc._curly_wrapper, pass raw str."""
        self.assertEqual(r"\$%{()}",
                         PrepHandeCcmcFciqmc._curly_wrapper(r"\$%", "()"))

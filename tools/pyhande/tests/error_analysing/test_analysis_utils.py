"""Test pyhande.error_analysing.analysis_utils.py."""
import unittest
import pandas as pd
import numpy as np
import pyhande.error_analysing.analysis_utils as analysis_utils
from tests.create_mock_df import create_qmc_frame


class TestCheckDataInput(unittest.TestCase):
    """Test `check_data_input()`."""

    def setUp(self):
        """Setting up common mock df."""
        rng = np.random.default_rng(7826)
        cols = ['Shift', 'test', 'N_0', r'\sum H_0j N_j', 'alt',
                'Proj. Energy']
        means = [-0.1, 10.0, 11.0, -9.0, 20.0, -0.2]
        sine_periods = list(range(1, len(cols)+1))
        noise_facs = [0.1*mean for mean in means]
        self._df_mock = create_qmc_frame(rng, cols, means, sine_periods,
                                         noise_facs,
                                         frac_not_convergeds=[0.2]*len(cols),
                                         num_mc_its=300)

    def test_cols_no_error(self):
        """All elements of `cols` are in `data`."""
        analysis_utils.check_data_input(
            [self._df_mock], ['Proj. Energy', 'Shift'], None, None, None, None)

    def test_cols_error(self):
        """Not all elements of `cols` are in all elements of `data`."""
        with self.assertRaisesRegex(ValueError, "'cols' parameter must only "
                                    "contain columns names present in all "
                                    "dataframes in 'data'."):
            analysis_utils.check_data_input(
                [self._df_mock,
                 self._df_mock.rename(columns={'Shift': 'Shifti'})],
                ['Proj. Energy', 'Shift'], None, None, None, None)

    def test_eval_ratio_no_error(self):
        """'num' and 'denom' value of `eval_ratio` are in `data`."""
        analysis_utils.check_data_input(
            [self._df_mock, self._df_mock], None,
            {'num': 'Shift', 'denom': 'alt', 'name': 'P'}, None, None, None)

    def test_eval_ratio_error(self):
        """'num' and 'denom' value of `eval_ratio` not in all `data`."""
        with self.assertRaisesRegex(ValueError, "When 'eval_ratio' is "
                                    "defined, its values must be column names "
                                    "present in all dataframes in 'data'."):
            analysis_utils.check_data_input(
                [self._df_mock, self._df_mock.rename(columns={'alt': 'a'})],
                None, {'num': 'Shift', 'denom': 'alt', 'name': 'Shift'}, None,
                None, None)

    def test_hybrid_col_no_error(self):
        """`hybrid_col` is in `data`."""
        analysis_utils.check_data_input(
            [self._df_mock, self._df_mock], None, None, 'alt', None, None)

    def test_hybrid_col_error(self):
        """`hybrid_col` is not in `data`."""
        with self.assertRaises(ValueError):
            # [todo] add check for actual error message.
            analysis_utils.check_data_input(
                [self._df_mock], None, None, 'al', None, None)

    def test_start_its_auto(self):
        """`start_its` is a string, defining starting it function."""
        analysis_utils.check_data_input(
            [self._df_mock, self._df_mock], None, None, None, 'blo', None)

    def test_start_its_same_length(self):
        """`start_its` has the same number of elements as `data`."""
        analysis_utils.check_data_input(
            [self._df_mock, self._df_mock], None, None, None, [0, -22], None)

    def test_start_its_diff_length(self):
        """`start_its` has a different number of elements as `data`."""
        with self.assertRaises(ValueError):
            # [todo] add check for actual error message.
            analysis_utils.check_data_input(
                [self._df_mock, self._df_mock], None, None, None, [6123], None)

    def test_end_its_same_length(self):
        """`end_its` has the same number of elements as `data`."""
        analysis_utils.check_data_input(
            [self._df_mock, self._df_mock], None, None, None, None, [0, -22])

    def test_end_its_diff_length(self):
        """`end_its` has a different number of elements as `data`."""
        with self.assertRaises(ValueError):
            # [todo] add check for actual error message.
            analysis_utils.check_data_input(
                [self._df_mock, self._df_mock], None, None, None, None, [6123])


class TestSetCols(unittest.TestCase):
    """Test `set_cols()`."""

    def setUp(self):
        """Setting up."""
        self._observables = {'afromobs': 'a', 'x_from_o': 'x',
                             'c_key': r'$%\X', 'test_key': 'test'}

    def test_modify_all(self):
        """All input has to be translated."""
        its, cols, rep_cols, eval_ratio, hybrid_col = analysis_utils.set_cols(
            self._observables, 'obs:x_from_o', ['obs:test_key', 'obs:c_key'],
            'obs:afromobs', {'num': 'obs:x_from_o', 'denom': 'obs:c_key',
                             'name': 'obs:test_key'}, 'obs:afromobs')
        self.assertEqual(its, 'x')
        self.assertListEqual(cols, ['test', r'$%\X'])
        self.assertEqual(rep_cols, 'a')
        self.assertDictEqual(eval_ratio, {'num': 'x', 'denom': r'$%\X',
                                          'name': 'test'})
        self.assertEqual(hybrid_col, 'a')

    def test_modify_none(self):
        """No input has to be translated."""
        its, cols, rep_cols, eval_ratio, hybrid_col = analysis_utils.set_cols(
            self._observables, 'xy', ['xobs:test_key', 'c_key'],
            'test', {'num': 'n', 'denom': 'd', 'name': 'z'}, 'Obs:afromobs')
        self.assertEqual(its, 'xy')
        self.assertListEqual(cols, ['xobs:test_key', 'c_key'])
        self.assertEqual(rep_cols, 'test')
        self.assertDictEqual(
            eval_ratio, {'num': 'n', 'denom': 'd', 'name': 'z'})
        self.assertEqual(hybrid_col, 'Obs:afromobs')

    def test_modify_some(self):
        """Some input has to be translated."""
        its, cols, rep_cols, eval_ratio, hybrid_col = analysis_utils.set_cols(
            self._observables, 'xy', ['obstest_key', 'obs:c_key'],
            'test', {'num': 'n', 'denom': 'obs:test_key', 'name': 'z'},
            'Obs:afromobs')
        self.assertEqual(its, 'xy')
        self.assertListEqual(cols, ['obstest_key', r'$%\X'])
        self.assertEqual(rep_cols, 'test')
        self.assertDictEqual(
            eval_ratio, {'num': 'n', 'denom': 'test', 'name': 'z'})
        self.assertEqual(hybrid_col, 'Obs:afromobs')

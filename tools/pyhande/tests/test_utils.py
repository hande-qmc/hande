"""Tests for utils.py."""
import unittest
import pandas as pd
import numpy as np
import pyhande.utils as utils
import create_mock_df


class TestGroupbyBetaLoops(unittest.TestCase):
    """Test utils.groupby_beta_loops"""

    def setUp(self):
        # Create mock partial (dm)qmc dataframe
        rng = np.random.default_rng(375)
        cols = ['N_0', 'alt1']
        means = [100.0, -11.0]
        sine_periods = [5.0, 2.9]
        noise_facs = [0.1, 0.9]
        num_mc_its = 6
        data = create_mock_df.create_qmc_frame(
            rng, cols, means, sine_periods, noise_facs,
            frac_not_convergeds=[0.5 for _ in range(5)], num_mc_its=num_mc_its)
        # Need 'iterations' and alternative column which increase in
        # steps
        mc_cycles1 = 10
        mc_cycles2 = 3
        iterations = pd.DataFrame.from_dict({
            'iterations': 2*[mc_cycles1*i for i in range(1, num_mc_its//2+1)]
            })
        alt2_iterations = pd.DataFrame.from_dict({
            'alt2': 3*[mc_cycles2*i for i in range(1, num_mc_its//3+1)]
            })
        self.data = pd.concat([iterations, data, alt2_iterations], axis=1)

    def test_basic_input(self):
        """Test with basic input."""
        group_by = utils.groupby_beta_loops(self.data)
        group0_mock = self.data.loc[0:2]
        pd.testing.assert_frame_equal(
            group_by.get_group(0.0), group0_mock, check_exact=False)

    def test_custom_col(self):
        """Test with custom column."""
        group_by = utils.groupby_beta_loops(self.data, name='alt2')
        group0_mock = self.data.loc[0:1]
        pd.testing.assert_frame_equal(
            group_by.get_group(0.0), group0_mock, check_exact=False)

    def test_unchanged_mutable(self):
        """Check that mutable objects, such as pd DataFrames, don't
        change when they shouldn't.
        """
        data_copy = self.data.copy()
        _ = utils.groupby_beta_loops(self.data)
        pd.testing.assert_frame_equal(self.data, data_copy, check_exact=True)


class TestGroupbyIterations(unittest.TestCase):
    """Test utils.groupby_iterations"""

    def setUp(self):
        # Create mock partial qmc dataframe
        rng = np.random.default_rng(2871)
        cols = ['N_0', 'alt1']
        means = [100.0, -11.0]
        sine_periods = [5.0, 2.9]
        noise_facs = [0.1, 0.9]
        num_mc_its = 4
        data = create_mock_df.create_qmc_frame(
            rng, cols, means, sine_periods, noise_facs,
            frac_not_convergeds=[0.5 for _ in range(5)], num_mc_its=num_mc_its
            )
        # Need 'iterations' column which increases in steps (mc_cycles)
        # Have two groups of sequences where the iteration increases
        mc_cycles = 10
        iterations = pd.DataFrame.from_dict({
            'iterations': 2*[mc_cycles*i for i in range(1, num_mc_its//2+1)]
            })
        self.data = pd.concat([iterations, data], axis=1)

    def test_basic_input(self):
        """Test with basic input."""
        group_by = utils.groupby_iterations(self.data)
        group0_mock = self.data.loc[0:1]
        pd.testing.assert_frame_equal(
            group_by.get_group(0.0), group0_mock, check_exact=False)

    def test_unchanged_mutable(self):
        """Check that mutable objects, such as pd DataFrames, don't
        change when they shouldn't.
        """
        data_copy = self.data.copy()
        _ = utils.groupby_beta_loops(self.data)
        pd.testing.assert_frame_equal(self.data, data_copy, check_exact=True)

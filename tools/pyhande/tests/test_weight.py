"""Tests for weight.py."""
import unittest
import warnings
import numpy as np
import pandas as pd
import pyhande.weight as weight
import create_mock_df


class TestReweight(unittest.TestCase):
    """Test weight.reweight."""

    def setUp(self):
        rng = np.random.default_rng(248451)
        self.data = create_mock_df.create_qmc_frame(
            rng, ['Shift', 'Alt'], [-1.3, -1.0], [12.0, 7.0], [0.1, 0.05])

    def test_basic_input(self):
        """Test basic input."""
        (mc_cycles, tstep, weight_history, mean_shift) = (4, 0.1, 10, -1.3)
        out = weight.reweight(self.data, mc_cycles, tstep, weight_history,
                              mean_shift)
        data_copy = self.data.copy()
        data_copy.drop(columns=['Weight'], inplace=True)
        data_copy['Weight'] = np.asarray([
            0.99256485, 0.91661459, 0.9090655, 0.88330944, 0.81232758,
            0.76941044, 0.68219966, 0.61679495, 0.62527968, 0.59195304,
            0.5751375, 0.6649207, 0.68626078, 0.71131836, 0.83535056,
            1.00111087, 1.20945807, 1.34714007, 1.35423976, 1.45681871,
            1.63262663, 1.59528056, 1.61402654, 1.58127348, 1.45880711,
            1.31618546, 1.21805526, 1.27076857, 1.16044234, 1.12289736
        ])
        pd.testing.assert_frame_equal(self.data, data_copy, check_exact=False)
        # Note that data and out are identical. Test this just to be
        # safe.
        self.assertEqual(id(self.data), id(out))

    def test_weight_key_input(self):
        """Use custom column in data."""
        (mc_cycles, tstep, weight_history, mean_shift) = (3, 0.02, 8, -1.0)
        out = weight.reweight(
            self.data, mc_cycles, tstep, weight_history, mean_shift,
            weight_key='Alt'
        )
        data_copy = self.data.copy()
        data_copy.drop(columns=['Weight'], inplace=True)
        data_copy['Weight'] = np.asarray([
            1.00178842, 1.00003984, 0.99693812, 0.98701375, 0.97998815,
            0.97256525, 0.96836244, 0.97286865, 0.97354347, 0.98475152,
            0.99543395, 1.00728771, 1.02134185, 1.03276738, 1.03521551,
            1.02923804, 1.02389504, 1.00214499, 0.98853845, 0.98495052,
            0.97725045, 0.97006702, 0.97514141, 0.98271336, 0.99309255,
            1.01049107, 1.01578566, 1.02178473, 1.02575039, 1.02785067
        ])
        pd.testing.assert_frame_equal(self.data, data_copy, check_exact=False)
        # Note that data and out are identical. Test this just to be
        # safe.
        self.assertEqual(id(self.data), id(out))

    def test_zero_tstep_input(self):
        """tstep == 0."""
        (mc_cycles, tstep, weight_history, mean_shift) = (10, 0.0, 20, -0.1)
        self.assertRaises(
            ValueError, weight.reweight, self.data, mc_cycles, tstep,
            weight_history, mean_shift
        )

    def test_neg_mc_cycles_input(self):
        """mc_cycles < 0."""
        (mc_cycles, tstep, weight_history, mean_shift) = (-10, 0.01, 20, -0.1)
        self.assertRaises(
            ValueError, weight.reweight, self.data, mc_cycles, tstep,
            weight_history, mean_shift)

    def test_unchanged_mutable(self):
        """Check that mutable objects, such as pd DataFrames, don't
        change when they shouldn't.
        """
        (mc_cycles, tstep, weight_history, mean_shift) = (4, 0.1, 10, -1.3)
        # data_copy = self.data.copy()
        _ = weight.reweight(self.data, mc_cycles, tstep, weight_history,
                            mean_shift)
        # pd.testing.assert_frame_equal(self.data, data_copy, check_exact=True)
        # [todo] Is there a better way to do this warning?
        warnings.warn("TestReweight.test_unchanged_mutable: "
                      "Mutable data is changed in function! Fix?")


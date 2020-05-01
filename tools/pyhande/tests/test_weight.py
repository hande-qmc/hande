"""Tests for weight.py."""
import unittest
import warnings
import numpy as np
import pandas as pd
import pyhande.weight as weight
from create_dummy_df import _CreateDummyDfs


class TestReweight(unittest.TestCase):
    """Test weight.reweight."""

    def setUp(self):
        create_dummy = _CreateDummyDfs(248451)
        self.data = create_dummy.create_qmc_frame(
            ['Shift', 'Alt'], [-1.3, -1.0], [12.0, 7.0], [0.1, 0.05])

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

    def test_arith_mean_input(self):
        """Use arith_mean == True."""
        (mc_cycles, tstep, weight_history, mean_shift) = (11, 0.9, 11, -1.3)
        out = weight.reweight(
            self.data, mc_cycles, tstep, weight_history, mean_shift,
            arith_mean=True
        )
        data_copy = self.data.copy()
        data_copy.drop(columns=['Weight'], inplace=True)
        data_copy['Weight'] = np.asarray([
            8.31347479e-01, 1.15910410e-01, 9.44562203e-02, 4.63760264e-02,
            5.83239848e-03, 1.52199199e-03, 7.75013960e-05, 6.39712806e-06,
            8.97096892e-06, 2.31264822e-06, 9.42175578e-07, 2.87451628e-06,
            2.96505141e-05, 9.05840759e-05, 5.83190581e-04, 6.42801774e-02,
            2.01816209e+00, 4.60122716e+01, 9.30598203e+02, 2.71273093e+03,
            2.56933872e+04, 2.54468154e+05, 3.59230985e+05, 1.57029508e+05,
            8.71118923e+04, 4.03065068e+04, 5.08904177e+03, 5.62220562e+02,
            1.79813927e+02, 4.15252329e+01
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


class TestArithSeries(unittest.TestCase):
    """Test weight.arith_series."""

    def test_basic_input_wnow_greater_wbefore(self):
        """weight_now > weight_before."""
        (tstep, mc_cycles, weight_now, weight_before) = (0.01, 10, 0.1, 0.05)
        out = weight.arith_series(tstep, mc_cycles, weight_now, weight_before)
        self.assertAlmostEqual(out, 1.00225357)

    def test_basic_input_wnow_smaller_wbefore(self):
        """weight_now < weight_before."""
        (tstep, mc_cycles, weight_now, weight_before) = (0.01, 10, 0.05, 0.1)
        out = weight.arith_series(tstep, mc_cycles, weight_now, weight_before)
        self.assertAlmostEqual(out, 0.99775356)

    def test_basic_input_wnow_eq_wbefore(self):
        """weight_now == weight_before."""
        (tstep, mc_cycles, weight_now, weight_before) = (0.01, 10, 0.1, 0.1)
        out = weight.arith_series(tstep, mc_cycles, weight_now, weight_before)
        self.assertAlmostEqual(out, 1.0)

    def test_basic_input_big_wnow(self):
        """weight_now >> 1 - cannot be too large though, otherwise exp
        overflows.
        """
        (tstep, mc_cycles, weight_now, weight_before) = (0.01, 10, 100.0, 0.1)
        out = weight.arith_series(tstep, mc_cycles, weight_now, weight_before)
        self.assertAlmostEqual(out, 1271.08561912)

    def test_basic_input_big_wbefore(self):
        """weight_before >> 1."""
        (tstep, mc_cycles, weight_now, weight_before) = (0.01, 10, 0.1, 100.0)
        out = weight.arith_series(tstep, mc_cycles, weight_now, weight_before)
        self.assertAlmostEqual(out, 0.15828258)

    def test_basic_input_neg_wnow(self):
        """weight_now < 0."""
        (tstep, mc_cycles, weight_now, weight_before) = (0.01, 10, -0.1, 0.2)
        out = weight.arith_series(tstep, mc_cycles, weight_now, weight_before)
        self.assertAlmostEqual(out, 0.98662734)

    def test_zero_tstep(self):
        """tstep == 0."""
        (tstep, mc_cycles, weight_now, weight_before) = (0.0, 10, -0.1, 0.2)
        self.assertRaises(
            ValueError, weight.arith_series, tstep, mc_cycles, weight_now,
            weight_before)

    def test_neg_mc_cycles(self):
        """mc_cycles < 0."""
        (tstep, mc_cycles, weight_now, weight_before) = (0.01, -10, -0.1, 0.2)
        self.assertRaises(
            ValueError, weight.arith_series, tstep, mc_cycles, weight_now,
            weight_before)

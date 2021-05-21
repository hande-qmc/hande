"""Tests for weight.py."""
import unittest
import warnings
import numpy as np
import pandas as pd
import pyhande.weight as weight
import tests.create_mock_df as create_mock_df


class TestReweight(unittest.TestCase):
    """Test weight.reweight."""

    def setUp(self):
        rng = np.random.default_rng(248451)
        self.data = create_mock_df.create_qmc_frame(
            rng, ['Shift', 'Alt'], [-1.3, -1.0], [12.0, 7.0], [0.1, 0.05])

    def test_basic_input(self):
        """Test basic input."""
        (mc_cycles, tstep, weight_history, mean_shift) = (4, 0.1, 10, -1.3)
        weights = weight.reweight(self.data, mc_cycles, tstep, weight_history,
                                  mean_shift)
        weights_exp = [
            0.9925648521738141, 0.9166145889890867, 0.9090654960602974,
            0.8833094364811825, 0.8123275817820274, 0.7694104408919243,
            0.6821996558442511, 0.616794950851474, 0.6252796826000567,
            0.5919530405088881, 0.5751375025400743, 0.6649207049102872,
            0.686260779187176, 0.7113183626102892, 0.8353505572867969,
            1.001110870216326, 1.2094580656366476, 1.3471400723768445,
            1.3542397576893366, 1.4568187130633774, 1.6326266341297668,
            1.5952805608339418, 1.6140265447518267, 1.5812734786176506,
            1.458807107328307, 1.3161854614506716, 1.2180552596788903,
            1.2707685717179618, 1.1604423403147808, 1.1228973649952922
        ]
        assert np.allclose(weights, weights_exp)

    def test_weight_key_input(self):
        """Use custom column in data."""
        (mc_cycles, tstep, weight_history, mean_shift) = (3, 0.02, 8, -1.0)
        weights = weight.reweight(
            self.data, mc_cycles, tstep, weight_history, mean_shift,
            weight_key='Alt'
        )
        weights_exp = [
            1.0017884205082834, 1.0000398418888259, 0.9969381155705571,
            0.9870137492389269, 0.9799881549626391, 0.972565247591884,
            0.9683624354090272, 0.9728686476670971, 0.9735434671372174,
            0.9847515208233709, 0.9954339536983114, 1.0072877088130072,
            1.0213418528613503, 1.032767375832493, 1.0352155113661938,
            1.0292380393457086, 1.0238950390664328, 1.0021449886133607,
            0.9885384533169285, 0.9849505201506198, 0.977250452926066,
            0.9700670184395724, 0.9751414058334306, 0.9827133605966373,
            0.993092551463694, 1.01049106979456, 1.015785662792946,
            1.0217847276206045, 1.0257503886293498, 1.0278506705546948
        ]
        assert np.allclose(weights, weights_exp)

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
        data_copy = self.data.copy()
        (mc_cycles, tstep, weight_history, mean_shift) = (4, 0.1, 10, -1.3)
        _ = weight.reweight(self.data, mc_cycles, tstep, weight_history,
                            mean_shift)
        pd.testing.assert_frame_equal(self.data, data_copy, check_exact=True)


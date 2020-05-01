"""Tests for canonical.py."""
import warnings
import unittest
import pandas as pd
import numpy as np
import pyhande.canonical as canonical
import create_mock_df


class TestAnalyseHFObservables(unittest.TestCase):
    """Test canonical.analyse_hf_observables."""

    def setUp(self):
        # Create mock series containing means of various DMQMC
        # observables as well as their convariances.
        rng = np.random.default_rng(375)
        means = [-90.0, -99.0, -198.0, 2.8, 10.0]
        self.cols = [
            r'Tr(T\rho_HF)', r'Tr(V\rho_HF)', r'Tr(H\rho_HF)', r'Tr(\rho_HF)',
            'alt1'
            ]
        self.means_series = pd.Series(means, index=self.cols)
        cov_orig = create_mock_df.create_cov_frame(rng, self.cols, means, 1)
        # reblock iterations are not a concept here:
        self.cov_orig = cov_orig.loc[0]
        # Shared mock results
        self.result_mock = pd.Series(
            [-32.142857, 4.004706, -35.357143, 4.122886, -70.714286, 10.74825],
            index=[
                'T_HF', 'T_HF_error', 'V_HF', 'V_HF_error', 'U_HF',
                'U_HF_error'
                ]
            )

    def test_basic_input(self):
        """Basic input."""
        nsamples = 100
        result = canonical.analyse_hf_observables(
            self.means_series.loc[self.cols[:-1]], self.cov_orig.drop(
                columns=['alt1'], index=['alt1']
                ), nsamples
            )
        pd.testing.assert_series_equal(
            result, self.result_mock, check_exact=False
            )

    def test_ignore_col(self):
        """Input with a column ('alt1') which should be ignored."""
        nsamples = 100
        result = canonical.analyse_hf_observables(
            self.means_series, self.cov_orig, nsamples
            )
        pd.testing.assert_series_equal(
            result, self.result_mock, check_exact=False
            )

    def test_unchanged_mutable(self):
        """Check that mutable objects, such as pd DataFrames, don't
        change when they shouldn't.
        """
        means_series_copy = self.means_series.copy()
        cov_orig_copy = self.cov_orig.copy()
        nsamples = 100
        _ = canonical.analyse_hf_observables(
            self.means_series, self.cov_orig, nsamples
            )
        pd.testing.assert_series_equal(
            means_series_copy, self.means_series, check_exact=True
            )
        pd.testing.assert_frame_equal(
            cov_orig_copy, self.cov_orig, check_exact=True
            )


class TestEstimates(unittest.TestCase):
    """Test canonical.estimates."""

    def setUp(self):
        # Create mock DMQMC data as well as metadata.
        rng = np.random.default_rng(375)
        means = [-36.1, -39.3, -90.0, -99.0, 2.8, 0.13]
        cols = [
            '<T>_0', '<V>_0', r'Tr(T\rho_HF)', r'Tr(V\rho_HF)', r'Tr(\rho_HF)',
            'N_ACC/N_ATT'
            ]
        sine_periods = [3.2, 4.2, 7.1, 2.22, 5.73, 1.67]
        noise_facs = [0.1, 0.9, 1.2, 0.2, 0.76, 0.04]
        self.data = create_mock_df.create_qmc_frame(
            rng, cols, means, sine_periods, noise_facs
            )
        # Note that metadata is most likely a dictionary not
        # pd.DataFrame (which is not what the function docstring says)
        # See test_suite/mc_canonical_estimates/np1/H2O_ccpvdz_RHF for
        # similar values (magnitude)
        self.metadata1 = {
            'nattempts': 8000, 'beta': 0.44, 'free_energy_corr': -77.4,
            'fermi_temperature': False, 'alt1': 1
            }
        # See test_suite/mc_canonical_estimates/np1/ueg_n7_ec20_rs1 for
        # similar values (magnitude)
        self.metadata2 = {
            'nattempts': 5000, 'beta': 0.9, 'free_energy_corr': -26.1,
            'fermi_temperature': True, 'alt1': 1, 'system': {
                'ueg': {'E_fermi': 3.2}
                }
            }

    def test_not_ueg(self):
        """Not UEG."""
        result = canonical.estimates(self.metadata1, self.data)
        ind_mock = [
            'U_0', 'T_0', 'V_0', 'N_ACC/N_ATT', 'F_0', 'S_0', 'T_HF', 'V_HF',
            'U_HF'
            ]
        ind_mock_full = []
        for i in ind_mock:
            ind_mock_full.extend([i, i+'_error'])
        ind_mock_full = ['Beta'] + ind_mock_full
        result_mock = pd.Series([
            4.40000000e-01, -7.45375108e+01, 7.95735784e-01, -3.57175099e+01,
            4.73795670e-01, -3.88200009e+01, 5.72059025e-01, 1.41132320e-01,
            5.98715584e-03, -7.29498696e+01, 9.64142895e-02, 1.63822382e+01,
            2.01610646e-01, -3.18572634e+01, 1.43161978e+00, -3.51519450e+01,
            1.61877967e+00, -6.70092084e+01, 2.98339505e+00
            ], index=ind_mock_full)
        pd.testing.assert_series_equal(result, result_mock, check_exact=False)

    def test_ueg(self):
        """UEG."""
        result = canonical.estimates(self.metadata2, self.data)
        ind_mock = ['U_0', 'T_0', 'V_0', 'N_ACC/N_ATT', 'F_0', 'S_0',
                     'T_HF', 'V_HF', 'U_HF']
        ind_mock_full = []
        for i in ind_mock:
            ind_mock_full.extend([i, i+'_error'])
        ind_mock_full = ['Beta'] + ind_mock_full
        result_mock = pd.Series([
            9.00000000e-01, -7.45375108e+01, 7.95735784e-01, -3.57175099e+01,
            4.73795670e-01, -3.88200009e+01, 5.72059025e-01, 1.41132320e-01,
            5.98715584e-03, -19.138018, 0.1508347996, -4.662982060842,
            0.12887270046, -3.18572634e+01, 1.43161978e+00, -3.51519450e+01,
            1.61877967e+00, -6.70092084e+01, 2.98339505e+00
            ], index=ind_mock_full)
        pd.testing.assert_series_equal(result, result_mock, check_exact=False)

    def test_not_ueg_no_naccnatt_col(self):
        """No 'N_ACC/N_ATT' in columns (not UEG)
        [todo] Should the metadata be adjusted as well?
        """
        self.data.drop(columns=(['N_ACC/N_ATT']), inplace=True)
        result = canonical.estimates(self.metadata1, self.data)
        ind_mock = ['U_0', 'T_0', 'V_0', 'T_HF', 'V_HF', 'U_HF']
        ind_mock_full = []
        for i in ind_mock:
            ind_mock_full.extend([i, i+'_error'])
        ind_mock_full = ['Beta'] + ind_mock_full
        result_mock = pd.Series([
            4.40000000e-01, -74.574386427, 0.78720024, -35.86308032,
            0.46069296, -38.7113061, 0.58596805, -32.23777087, 1.44645184,
            -35.44043608, 1.63181849, -67.67820695, 3.01213831
            ], index=ind_mock_full)
        pd.testing.assert_series_equal(result, result_mock, check_exact=False)

    def test_unchanged_mutable(self):
        """Check that mutable objects, such as pd DataFrames, don't
        change when they shouldn't.
        """
        # data_copy = self.data.copy()
        metadata1_copy = self.metadata1.copy()
        _ = canonical.estimates(self.metadata1, self.data)
        # pd.testing.assert_frame_equal(self.data, data_copy, check_exact=True)
        # [todo] Is there a better way to do this warning?
        warnings.warn("TestEstimates.test_unchanged_mutable: "
                      "Mutable data is changed in function! Fix?")
        self.assertDictEqual(self.metadata1, metadata1_copy)

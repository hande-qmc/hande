"""Tests for dmqmc.py."""
import unittest
import copy
import pandas as pd
import numpy as np
import pyhande.dmqmc as dmqmc
import tests.create_mock_df as create_mock_df


class TestAnalyseObservables(unittest.TestCase):
    """Test dmqmc.analyse_observables."""

    def setUp(self):
        # Create mock pandas dataframe containing means of various
        # DMQMC observables as well as their convariances as a function
        # of beta.
        rng = np.random.default_rng(947)
        # Use create_mock_df.create_qmc_frame even though it is not
        # technically a qmc dataframe with qm time series data here.
        # See testsuite/idmqmc/np1/ueg_n14_ec4_rs1 for some inspiration
        # of these values - they are fairly random though.
        # Don't cut down the number of columns in this test as a few
        # are explicitly named in this function.
        means = [
            -0.3, 1243.0, 21434.0, 23111.0, 23111.0, 21234.0, -1534.0, -938.9,
            125.0, 1.3, -345.0, 0.1, 43.0, 421.0, 1000.0, 643.0, 1.2, -0.8,
            0.4, 100.34, 1.45, -0.185
        ]
        self.cols = [
            'Instant shift', 'Trace', r'\sum\rho_{ij}H_{ji}',
            r'\sum\rho_{ij}T_{ji}', r'\sum\rho_{ij}H0{ji}',
            r'\sum\rho_{ij}HI{ji}', r'\sum\rho_{ij}U_{ji}',
            r'\sum\rho_{ij}VI_{ji}', r'\sum\rho_{ij}S_{ji}',
            r'\sum\rho_{ij}M2{ji}', r'\sum\rho_{ij}H2{ji}', 'n_34', 'n_98',
            'S_34', 'S_98', 'Suu_34', 'Suu_98', 'Sud_34', 'Sud_98', 'Sud_123',
            '# particles', 'alt1'
        ]
        sine_periods = [
            34.0, 1.0, 1.1234, 7.3, 2.3, 8.3, 3.0, 21.9,
            4.0, 4.3, 2.2, 0.1, 1.2, 4.2, 21.0, 0.2, 1.4, 4.3, 2.3,
            6.0, 2.3, 4.4
        ]
        noise_facs = [0.1*mean for mean in means]  # being lazy
        num_mc_its = 3
        self.data = create_mock_df.create_qmc_frame(
            rng, self.cols, means, sine_periods, noise_facs, num_mc_its=num_mc_its)
        self.data.index.name = 'Beta'
        # Only scale a fraction of means here to that cov does not get
        # too big (causing NaN otherwise later). A bit of a hack so that
        # the error of a fraction can be evaluated.
        means[1] = 0.01*means[1]
        means[2] = 0.4*means[1]
        means[14] = 0.5*means[14]
        self.cov = create_mock_df.create_cov_frame(
            rng, self.cols, means, num_mc_its, index_name='Beta')
        self.nsamples = pd.Series([4000.0, 10000.0, 2000.0])
        self.nsamples.index.name = 'Beta'
        # Some shared exp results
        self.cols_exp = [
            'Tr[Hp]/Tr[p]', 'Tr[Hp]/Tr[p]_error', 'Tr[H2p]/Tr[p]',
            'Tr[H2p]/Tr[p]_error', 'Tr[Sp]/Tr[p]', 'Tr[Sp]/Tr[p]_error',
            'Tr[Mp]/Tr[p]', 'Tr[Mp]/Tr[p]_error', 'Tr[Tp]/Tr[p]',
            'Tr[Tp]/Tr[p]_error', 'Tr[Up]/Tr[p]', 'Tr[Up]/Tr[p]_error',
            'Tr[H0p]/Tr[p]', 'Tr[H0p]/Tr[p]_error', 'Tr[HIp]/Tr[p]',
            'Tr[HIp]/Tr[p]_error', 'VI', 'VI_error', 'n_34', 'n_34_error',
            'n_98', 'n_98_error', r'\sum\rho_{ij}S_{ji}',
            r'\sum\rho_{ij}S_{ji}_error', 'S_34', 'S_34_error', 'S_98',
            'S_98_error', 'Suu_34', 'Suu_34_error', 'Suu_98', 'Suu_98_error',
            'Sud_34', 'Sud_34_error', 'Sud_98', 'Sud_98_error', 'Sud_123',
            'Sud_123_error'
        ]
        self.result_exp = pd.DataFrame([[
            1.83810838e+01, 7.80937590e-04, -1.95040347e-01, 3.31649238e-03,
            1.05572002e-01, 1.02076313e-03, 1.14046306e-03, 1.04944497e-05,
            1.85468649e+01, 1.16551505e-01, -1.21724797e+00, 1.26930348e-02,
            1.82442637e+01, 2.36289777e-01, 1.64156605e+01, 1.34621439e-01,
            -7.91410966e-01, 6.68855619e-03, 7.86629266e-05, 5.53141981e-07,
            3.70743265e-02, 4.03306784e-04, 1.05572002e-01, 1.02076313e-03,
            3.38243779e-01, 4.38274774e-03, 8.66704774e-01, 2.09509983e-03,
            5.65387573e-01, 3.25985608e-03, 8.78750880e-04, 5.75310439e-06,
            -7.27350210e-04, 9.98046319e-06, 3.65009178e-04, 3.73698471e-06,
            6.66571424e-02, 8.49830261e-04
        ], [
            2.07673527e+01, 1.88450113e-04, -2.01688310e-01,
            2.79207996e-03, 1.18878041e-01, 5.01585276e-04, 1.04475680e-03,
            9.14491249e-06, 1.75648459e+01, 1.77431607e-01,
            -1.24390484e+00, 6.46502950e-03, 2.29668399e+01,
            1.85617194e-01, 1.52025390e+01, 4.35540141e-02,
            -8.19932002e-01, 5.93783493e-03, 8.71851089e-05,
            3.75070737e-07, 4.23113161e-02, 1.58752358e-04, 1.18878041e-01,
            5.01585276e-04, 3.98921268e-01, 2.48228859e-03,
            7.62241208e-01, 3.94163902e-03, 5.69420856e-01, 3.37993732e-03,
            1.14802236e-03, 6.33104260e-06, -7.12324321e-04,
            6.30500182e-06, 3.56932675e-04, 2.25559444e-06, 8.36670781e-02,
            6.87738977e-04
        ], [
            1.78165282e+01, 4.07376230e-03, -3.05565989e-01,
            2.84972287e-03, 1.19663741e-01, 2.33131895e-03, 1.31495009e-03,
            2.47621880e-05, 2.44372903e+01, 3.17639195e-01,
            -1.46645202e+00, 2.72964312e-02, 2.43270977e+01,
            3.96348095e-01, 2.32256391e+01, 4.54600317e-01,
            -8.27709491e-01, 2.03522394e-02, 9.14081367e-05,
            2.12671253e-06, 3.91268398e-02, 8.24337154e-04, 1.19663741e-01,
            2.33131895e-03, 3.91845115e-01, 8.55357649e-03, 9.70475476e-01,
            6.96118155e-03, 6.55752616e-01, 6.31743894e-03, 9.33489265e-04,
            2.14844361e-05, -7.65062946e-04, 7.94781737e-06,
            4.17602440e-04, 1.17137522e-06, 1.12226696e-01, 1.39555634e-03
        ]], columns=self.cols_exp)

    def test_basic_input(self):
        """Basic input - Don't include 'alt1' column and index."""
        result = dmqmc.analyse_observables(
            self.data.drop(columns=['alt1']),
            self.cov.drop(
                columns='alt1').drop(index='alt1', level=1), self.nsamples)
        pd.testing.assert_frame_equal(
            result, self.result_exp, check_exact=False)

    def test_ignore_col(self):
        """Input with a column which should be ignored."""
        result = dmqmc.analyse_observables(self.data, self.cov, self.nsamples)
        pd.testing.assert_frame_equal(
            result, self.result_exp, check_exact=False)

    def test_only_two_cols(self):
        """Input with only two columns, 'Trace' and 'Suu_34'.  Have to
        pass 'Trace'.
        """
        self.cols.remove('Trace')
        self.cols.remove('Suu_34')
        result = dmqmc.analyse_observables(
            self.data.drop(columns=self.cols),
            self.cov.drop(
                columns=self.cols).drop(index=self.cols, level=1),
            self.nsamples)
        self.cols_exp.remove('Suu_34')
        self.cols_exp.remove('Suu_34_error')
        self.result_exp.drop(columns=self.cols_exp, inplace=True)
        pd.testing.assert_frame_equal(
            result, self.result_exp, check_exact=False)

    def test_unchanged_mutable(self):
        """Check that mutable objects, such as pd DataFrames, don't
        change when they shouldn't.
        """
        data_copy = self.data.copy()
        cov_copy = self.cov.copy()
        nsamples_copy = self.nsamples.copy()
        _ = dmqmc.analyse_observables(self.data, self.cov, self.nsamples)
        pd.testing.assert_frame_equal(self.data, data_copy, check_exact=True)
        pd.testing.assert_frame_equal(self.cov, cov_copy, check_exact=True)
        pd.testing.assert_series_equal(
            self.nsamples, nsamples_copy, check_exact=True)


class TestAddObservableToDict(unittest.TestCase):
    """Test dmqmc.add_observable_to_dict."""

    def setUp(self):
        self.observables = {'obs1': 'test', 'alt2': 'test', 'alt3': 'dd'}
        self.columns = ['obs2', 'let4', 'let5', 'ltrp']

    def test_basic_input(self):
        """Test basic input."""
        dmqmc.add_observable_to_dict(self.observables, self.columns, 'le')
        result_exp = {
            'obs1': 'test', 'alt2': 'test', 'alt3': 'dd', 'let4': 'let4',
            'let5': 'let5'
        }
        self.assertDictEqual(self.observables, result_exp)

    def test_not_matching_label(self):
        """Test label which is not present in columns."""
        result_exp = self.observables.copy()
        dmqmc.add_observable_to_dict(self.observables, self.columns, 'x')
        self.assertDictEqual(self.observables, result_exp)

    def test_empty_observables(self):
        """Test passing an empty dictionary for observables."""
        self.observables = {}
        dmqmc.add_observable_to_dict(self.observables, self.columns, 'le')
        result_exp = {'let4': 'let4', 'let5': 'let5'}
        self.assertDictEqual(self.observables, result_exp)

    def test_unchanged_mutable(self):
        """Check that mutable objects, such as pd DataFrames, don't
        change when they shouldn't.
        """
        columns_copy = self.columns.copy()
        dmqmc.add_observable_to_dict(self.observables, self.columns, 'x')
        self.assertListEqual(self.columns, columns_copy)


class FreeEnergyErrorAnalysis(unittest.TestCase):
    """Test dmqmc.free_energy_error_analysis."""

    def setUp(self):
        # Create mock pandas dataframe with "raw" DMQMC data as well
        # as one that contains the means as a function of Beta loops.
        rng = np.random.default_rng(132)
        # Use create_mock_df.create_qmc_frame also for the means df even
        # though it is not technically a qmc dataframe with qm time
        # series data here.
        # See testsuite/idmqmc/np1/ueg_n14_ec4_rs1 for some inspiration
        # of these values - they are fairly random though.
        # r'\sum\rho_{ij}VI_{ji}' and 'Trace' are explicity mentioned
        # in the function so have to be here.
        self.means = [-0.3, 1200.0, 21434.0, -9879, 4.5, 100.34, -0.185]
        self.cols_data = [
            'Instant shift', 'Trace', r'\sum\rho_{ij}H_{ji}',
            r'\sum\rho_{ij}VI_{ji}', 'n_34', '# particles', 'alt1'
        ]
        sine_periods = [34.0, 1.0, 1.1234, 3.7, 7.3, 2.3, 8.3]
        noise_facs = [0.1*mean for mean in self.means]  # being lazy
        num_mc_its = 30
        self.data = create_mock_df.create_qmc_frame(
            rng, self.cols_data, self.means, sine_periods, noise_facs,
            num_mc_its=num_mc_its)
        betas = pd.DataFrame(
            [float(i) for i in range(num_mc_its//2)]
            + [float(i) for i in range(num_mc_its - (num_mc_its//2))],
            columns=['Beta'])
        self.data = pd.concat([betas, self.data], axis=1)
        self.cols_results = self.cols_data.copy()
        # Some columns were renamed (in analyse_observables) before
        # they were passed.
        self.cols_results[2] = 'Tr[Hp]/Tr[p]'
        self.cols_results[3] = 'VI'
        num_mc_its = 2
        self.results = create_mock_df.create_qmc_frame(
            rng, self.cols_results, self.means, sine_periods, noise_facs,
            num_mc_its=num_mc_its)
        betas = pd.DataFrame([0.0, 1.0], columns=['Beta'])
        self.results = pd.concat([betas, self.results], axis=1)

    def test_basic_input(self):
        """Test basic input."""
        dtau = 0.1
        dmqmc.free_energy_error_analysis(self.data, self.results, dtau)
        results_exp = pd.DataFrame([[
            0.00000000e+00, -3.41640932e-01, 1.19243274e+03, 2.16706665e+04,
            -1.01442958e+04, 4.64565821e+00, 1.00524636e+02, -1.66484484e-01,
            0.00000000e+00, 0.0
        ], [
            1.00000000e+00, -2.92399492e-01, 1.07397341e+03,
            1.65884483e+04, -9.93263836e+03, 4.75548513e+00,
            1.02023051e+02, -1.52752934e-01, -1.00384671e+03,
            4.68667904e-03
        ]], columns=[
            'Beta', 'Instant shift', 'Trace', 'Tr[Hp]/Tr[p]', 'VI', 'n_34',
            '# particles', 'alt1', 'f_xc', 'f_xc_error'
        ])
        pd.testing.assert_frame_equal(
            self.results, results_exp, check_exact=False)
        self.assertWarnsRegex(
            UserWarning, r'Ratio error of 992.434317 is not insignificant - '
            'check results.', dmqmc.free_energy_error_analysis, self.data,
            self.results, dtau)

    # def test_further(self):
        # [todo] add more tests, possibly one without a warning!

    def test_unchanged_mutable(self):
        """Check that mutable objects, such as pd DataFrames, don't
        change when they shouldn't.  self.results is supposed to
        change, it is the result.
        """
        data_copy = self.data.copy()
        dtau = 0.1
        dmqmc.free_energy_error_analysis(self.data, self.results, dtau)
        pd.testing.assert_frame_equal(
            self.data, data_copy, check_exact=True)


class AnalyseRenyiEntropy(unittest.TestCase):
    """Test dmqmc.analyse_renyi_entropy."""

    def setUp(self):
        # Not too dissimilar to TestAnalyseObservables.setUp().
        # Create mock pandas dataframe containing means of various
        # DMQMC observables as well as their convariances as a function
        # of beta.
        rng = np.random.default_rng(5922)
        # Use create_mock_df.create_qmc_frame even though it is not
        # technically a qmc dataframe with qm time series data here.
        # See testsuite/idmqmc/np1/ueg_n14_ec4_rs1 for some inspiration
        # of these values - they are fairly random though.
        # Explicity include some columns that are explicitly named in
        # the function.
        means = [
            -0.3, 1222.0, 21434.0, 2311.0, 81.0, 34.0, 154.0, 1111.9, 125.0,
            87.0, 31.0, 159.0,
        ]
        self.cols = [
            'Instant shift', 'Trace', r'\sum\rho_{ij}H_{ji}',
            'Trace 2', 'RDM1 trace 1', 'RDM1 trace 2', 'RDM1 S2',
            'Full S2', 'alt1', 'RDM2 trace 1', 'RDM2 trace 2', 'RDM2 S2'
        ]
        sine_periods = [
            34.0, 1.0, 1.1234, 7.3, 2.3, 8.3, 3.0, 21.9, 4.0, 2.3, 18.3, 3.0
        ]
        noise_facs = [0.1*mean for mean in means]  # being lazy
        num_mc_its = 3
        self.data = create_mock_df.create_qmc_frame(
            rng, self.cols, means, sine_periods, noise_facs,
            num_mc_its=num_mc_its)
        self.data.index.name = 'Beta'
        # Fix up means slightly (they are used for scaling for cov)
        means[-6] = 0.1*means[-6]
        means[-5] = 0.1*means[-5]
        means[-1] = 0.1*means[-1]
        self.cov = create_mock_df.create_cov_frame(
            rng, self.cols, means, num_mc_its, index_name='Beta')
        self.nsamples = pd.Series([4000.0, 1290.0, 3000.0])
        self.nsamples.index.name = 'Beta'

    def test_input_only1(self):
        """Don't pass RDM2 data."""
        cols = ['RDM2 trace 1', 'RDM2 trace 2', 'RDM2 S2']
        self.data.drop(columns=cols, inplace=True)
        self.cov.drop(
            columns=cols, inplace=True)
        self.cov.drop(index=cols, level=1, inplace=True)
        results = dmqmc.analyse_renyi_entropy(
            self.data, self.cov, self.nsamples)
        results_exp = pd.DataFrame([
            [4.159087, 0.03201598, 10.984759, 0.04059646],
            [4.211940, 0.04796966, 10.929251, 0.0648403],
            [3.941436, 0.021987599, 11.168785, 0.030056]
        ], columns=(
            'RDM1 S2', 'RDM1 S2 error', 'Full S2', 'Full S2 error'
        ))
        pd.testing.assert_frame_equal(
            results, results_exp, check_exact=False)

    def test_input_with2(self):
        """Test input with column names ending with 2."""
        results = dmqmc.analyse_renyi_entropy(
            self.data, self.cov, self.nsamples
        )
        results_exp = pd.DataFrame([
            [4.159087, 0.03201598, 10.984759, 0.04059646, 3.895219,
             0.01661310],
            [4.211940, 0.04796966, 10.929251, 0.0648403, 4.031990, 0.0495852],
            [3.941436, 0.021987599, 11.168785, 0.030056, 4.122698, 0.04685934]
        ], columns=(
            'RDM1 S2', 'RDM1 S2 error', 'Full S2', 'Full S2 error',
            'RDM2 S2', 'RDM2 S2 error'
        ))
        pd.testing.assert_frame_equal(
            results, results_exp, check_exact=False)

    def test_unchanged_mutable(self):
        """Check that mutable objects, such as pd DataFrames, don't
        change when they shouldn't.
        """
        data_copy = self.data.copy()
        cov_copy = self.cov.copy()
        nsamples_copy = self.nsamples.copy()
        _ = dmqmc.analyse_renyi_entropy(self.data, self.cov, self.nsamples)
        pd.testing.assert_frame_equal(self.data, data_copy, check_exact=True)
        pd.testing.assert_frame_equal(self.cov, cov_copy, check_exact=True)
        pd.testing.assert_series_equal(
            self.nsamples, nsamples_copy, check_exact=True)


class TestCalcS2(unittest.TestCase):
    """Test dmqmc.calc_s2()"""

    def setUp(self):
        rng = np.random.default_rng(99562)
        # Use create_mock_df.create_qmc_frame even though it is not
        # technically a qmc dataframe with qm time series data here.
        self.cols = ['mean', 'standard error']
        meanA = 98.0
        meanB = 2893.0
        meanC = 487.0
        steA = 0.1
        steB = 98.2
        steC = 64.5
        sine_periods = [34.0, 2.0]
        num_mc_its = 3
        self.stats_A = create_mock_df.create_qmc_frame(
            rng, self.cols, [meanA, steA], sine_periods, [0.1*meanA, 0.1*steA],
            num_mc_its=num_mc_its
        )
        self.stats_B = create_mock_df.create_qmc_frame(
            rng, self.cols, [meanB, steB], sine_periods, [0.1*meanB, 0.1*steB],
            num_mc_its=num_mc_its
        )
        self.stats_C = create_mock_df.create_qmc_frame(
            rng, self.cols, [meanC, steC], sine_periods, [0.1*meanC, 0.1*steC],
            num_mc_its=num_mc_its
        )
        self.cov_AB = pd.Series([129388.0, 384912.2, 192843.1])
        self.cov_BC = pd.Series([-1290088.0, -5829912.0, -5987243.9])
        self.cov_AC = pd.Series([-56839.0, -27483.0, -39482.5])
        self.data_len = pd.Series([3444, 4892, 100239])

    def test_basic_input(self):
        """Test basic input."""
        mean, std_err = dmqmc.calc_S2(
            self.stats_A, self.stats_B, self.stats_C, self.cov_AB, self.cov_AC,
            self.cov_BC, self.data_len
        )
        mean_exp = pd.Series([13.589915, 14.130696, 13.822791], name='mean')
        std_err_exp = pd.Series([0.164870, 0.165638, 0.235193])
        pd.testing.assert_series_equal(mean_exp, mean, check_exact=False)
        pd.testing.assert_series_equal(
            std_err_exp, std_err, check_exact=False)

    def test_neg_mean(self):
        """Test passing a negative mean."""
        self.stats_A.replace(
            self.stats_A['mean'].loc[0], -self.stats_A['mean'].loc[0],
            inplace=True)
        mean, std_err = dmqmc.calc_S2(
            self.stats_A, self.stats_B, self.stats_C, self.cov_AB, self.cov_AC,
            self.cov_BC, self.data_len)
        mean_exp = pd.Series([np.nan, 14.130696, 13.822791], name='mean')
        std_err_exp = pd.Series([0.159725, 0.165638, 0.235193])
        pd.testing.assert_series_equal(mean_exp, mean, check_exact=False)
        pd.testing.assert_series_equal(
            std_err_exp, std_err, check_exact=False)

    def test_unchanged_mutable(self):
        """Check that mutable objects, such as pd DataFrames, don't
        change when they shouldn't.
        """
        self_obs_f = [self.stats_A, self.stats_B, self.stats_C]
        self_obs_s = [self.cov_AB, self.cov_BC, self.cov_AC, self.data_len]
        ob_copies_f = [self_ob_f.copy() for self_ob_f in self_obs_f]
        ob_copies_s = [self_ob_s.copy() for self_ob_s in self_obs_s]
        _ = dmqmc.calc_S2(
            self.stats_A, self.stats_B, self.stats_C, self.cov_AB, self.cov_AC,
            self.cov_BC, self.data_len
        )
        for (self_ob, ob_copy) in zip(self_obs_f, ob_copies_f):
            pd.testing.assert_frame_equal(self_ob, ob_copy, check_exact=True)
        for (self_ob, ob_copy) in zip(self_obs_s, ob_copies_s):
            pd.testing.assert_series_equal(self_ob, ob_copy, check_exact=True)


class TestCalcSplineFit(unittest.TestCase):
    """Test dmqmc.calc_spline_fit()"""

    def setUp(self):
        rng = np.random.default_rng(2842)
        # Use create_mock_df.create_qmc_frame even though it is not
        # technically a qmc dataframe with qm time series data here.
        self.cols = [
            'test', 'test error', 'Shift', 'Shift error', 'alt1', 'alt1 error'
        ]
        col_values = [487.0, 23.0, -1.2, 0.345, 6847.0, 47.0]
        sine_periods = [34.0, 2.0, 5.0, 23.0, 56.0, 2.0]
        num_mc_its = 4
        self.estimates = create_mock_df.create_qmc_frame(
            rng, self.cols, col_values, sine_periods,
            [0.1*colv for colv in col_values], num_mc_its=num_mc_its)

    def test_basic_input(self):
        """Test basic input."""
        result = dmqmc.calc_spline_fit('test', self.estimates)
        result_mock = pd.Series(
            [514.867277, 505.734033, 566.953745, 540.543340])
        pd.testing.assert_series_equal(result, result_mock, check_exact=False)

    def test_different_input(self):
        """Test alternative input."""
        result = dmqmc.calc_spline_fit('Shift', self.estimates)
        result_mock = pd.Series([
            -1.10427403697, -1.179055496740, -0.962132947671, -1.198842742582
        ])
        pd.testing.assert_series_equal(result, result_mock, check_exact=False)

    def test_unchanged_mutable(self):
        """Check that mutable objects, such as pd DataFrames, don't
        change when they shouldn't.
        """
        estimates_copy = self.estimates.copy()
        _ = dmqmc.calc_spline_fit('test', self.estimates)
        pd.testing.assert_frame_equal(
            self.estimates, estimates_copy, check_exact=True)


class TestSortMomentum(unittest.TestCase):
    """Test dmqmc.sort_momentum()."""

    def test_basic_input(self):
        """Test basic input."""
        columns = ['n_47', 'n_2', 'n_4', 'n_6']
        sort = dmqmc.sort_momentum(columns)
        self.assertListEqual(sort, ['n_2', 'n_4', 'n_6', 'n_47'])

    def test_digit_input(self):
        """Test passing post comma digits."""
        columns = ['n_4.7', 'n_2', 'n_4', 'n_6']
        sort = dmqmc.sort_momentum(columns)
        self.assertListEqual(sort, ['n_2', 'n_4', 'n_4.7', 'n_6'])

    def test_neg_input(self):
        """Test passing negative input."""
        columns = ['n_-4', 'n_2', 'n_4', 'n_6']
        sort = dmqmc.sort_momentum(columns)
        self.assertListEqual(sort, ['n_-4', 'n_2', 'n_4', 'n_6'])

    def test_unchanged_mutable(self):
        """Check that mutable objects, such as pd DataFrames, don't
        change when they shouldn't.
        """
        columns = ['n_47', 'n_2', 'n_4', 'n_6']
        columns_copy = columns.copy()
        _ = dmqmc.sort_momentum(columns)
        self.assertListEqual(columns, columns_copy)


class TestExtractMomentumCorrelationFunction(unittest.TestCase):
    """Test dmqmc.extract_momentum_correlation_function()."""

    def setUp(self):
        # Create mock pandas dataframe with "raw" DMQMC data as well
        # as one that contains the means as a function of Beta loops.
        rng = np.random.default_rng(39385)
        # Use create_mock_df.create_qmc_frame also for the means df even
        # though it is not technically a qmc dataframe with qm time
        # series data here.
        # See testsuite/idmqmc/np1/ueg_n14_ec4_rs1 for some inspiration
        # of some values - they are fairly random though.
        # r'\sum\rho_{ij}VI_{ji}' and 'Trace' are explicity mentioned
        # in the function so have to be here.
        # [WARNING] The function to be tested implicitly assumes that
        # [WARNING] an equal amount of colums start with "n_", "Suu_",
        # [WARNING] "S_", "Sud_" and they all have their "_error" cols.
        self.means = [
            -0.3, 1200.0, 2143.0, 320.0, -9879, 783, 4.5, 0.3,
            100.34, 10.4, 7.3, 0.4, 90.2, 10.2, 123.0, 14.8, 5.2, 0.4, -0.185
        ]
        self.cols_data = [
            'Instant shift', 'Trace', 'n_1', 'n_1_error', 'n_2',
            'n_2_error', 'S_1', 'S_1_error', 'S_2', 'S_2_error',
            'Sud_1', 'Sud_1_error', 'Sud_2', 'Sud_2_error', 'Suu_1',
            'Suu_1_error', 'Suu_2', 'Suu_2_error', 'alt1'
        ]
        sine_periods = [
            34.0, 1.0, 1.1234, 0.1, 3.7, 0.2, 7.3, 0.6, 2.3, 0.1, 5, 6,
            2, 4, 6, 8.3, 245, 64, 62
        ]
        noise_facs = [0.1*mean for mean in self.means]  # being lazy
        num_mc_its = 30
        self.data = create_mock_df.create_qmc_frame(
            rng, self.cols_data, self.means, sine_periods, noise_facs,
            num_mc_its=num_mc_its)
        betas = pd.DataFrame(
            [float(i) for i in range(num_mc_its)],
            columns=['Beta']
        )
        self.data = pd.concat([betas, self.data], axis=1)
        full_exp = pd.DataFrame([
            [1.00000000e+00, 2.05989953e+03, 3.38085214e+02, 5.21562026e+00,
             3.10943935e-01, 1.15340672e+02, 1.71763584e+01, 7.39238719e+00,
             4.06649786e-01],
            [2.00000000e+00, -8.99494034e+03, 7.81274739e+02, 8.11826620e+01,
             1.17069276e+01, 6.85358137e+00, 3.47132069e-01, 7.69600784e+01,
             9.97119359e+00]
        ], columns=[
            'k', 'n_k', 'n_k_error', 'S_k', 'S_k_error', 'Suu_k',
            'Suu_k_error', 'Sud_k', 'Sud_k_error'
        ])
        beta = pd.DataFrame([0, 0], columns=['Beta'])
        self.full_exp = pd.concat([beta, full_exp], axis=1)

    def test_basic_input(self):
        """Test basic input."""
        full = dmqmc.extract_momentum_correlation_function(self.data, 0)
        pd.testing.assert_frame_equal(full, self.full_exp)

    def test_diff_numbers_in_col_names(self):
        """Rename col names with "_2" to "_7"."""
        to_rename = ['n_', 'S_', 'Sud_', 'Suu_']
        to_rename_dict = {}
        for to_r in to_rename:
            to_rename_dict[to_r+str(2)] = to_r+str(7)
            to_rename_dict[to_r+str(2)+"_error"] = to_r+str(7)+"_error"
        self.data.rename(columns=to_rename_dict, inplace=True)
        full = dmqmc.extract_momentum_correlation_function(self.data, 0)
        self.full_exp.at[1, 'k'] = 7.0
        pd.testing.assert_frame_equal(full, self.full_exp)

    def test_unchanged_mutable(self):
        """Check that mutable objects, such as pd DataFrames, don't
        change when they shouldn't.
        """
        data_copy = self.data.copy()
        _ = dmqmc.extract_momentum_correlation_function(self.data, 0)
        pd.testing.assert_frame_equal(self.data, data_copy, check_exact=True)


class TestAnalyseData(unittest.TestCase):
    """Test dmqmc.analyse_data()."""

    def setUp(self):
        # Create mock pandas dataframe with "raw" DMQMC data as well
        # as metadata.
        self.rng = np.random.default_rng(331)
        # See testsuite/dmqmc/np2/n2_sto3g for some inspiration.
        means = [-0.5, 1200.0, -12143.0, 1320.0, -9879]
        self.cols_data = [
            'Shift', 'Trace', r'\sum\rho_{ij}H_{ji}', '# H psips', 'alt1'
        ]
        sine_periods = [34.0, 1.0, 0.6, 2.3, 64]
        noise_facs = [0.1*mean for mean in means]  # being lazy
        self.num_mc_its = 6
        self.data = []
        for _ in range(2):
            data = create_mock_df.create_qmc_frame(
                self.rng, self.cols_data, means, sine_periods, noise_facs,
                num_mc_its=self.num_mc_its)
            betas = pd.DataFrame(
                [float(i) for i in range(self.num_mc_its//2)]
                + [float(i) for i in
                   range(self.num_mc_its - (self.num_mc_its//2))],
                columns=['iterations'])
            self.data.append(pd.concat([betas, data], axis=1))
        self.metadata_dmqmc = {
            'calc_type': 'DMQMC', 'dmqmc': {'find_weights': False},
            'ipdmqmc': {'ipdmqmc': False, 'symmetric': True},
            'qmc': {'ncycles': 10, 'tau': 0.1}, 'alt1': 7, 'alt2': True
        }
        self.metadata_ipdmqmc = copy.deepcopy(self.metadata_dmqmc)
        self.metadata_ipdmqmc['ipdmqmc']['ipdmqmc'] = True
        self.basic_data_exp = pd.DataFrame([
            [0.0, -11.374353, 1.325053], [0.1, -11.741264, 2.124244],
            [0.2, -9.134270, 2.904695]
        ], columns=['Beta', 'Tr[Hp]/Tr[p]', 'Tr[Hp]/Tr[p]_error'], index=[
            0.0, 0.1, 0.2
        ])

    def test_basic_input(self):
        """Test basic DMQMC input."""
        meta_exp = copy.deepcopy(self.metadata_dmqmc)
        (meta, data) = dmqmc.analyse_data(
            [(self.metadata_dmqmc, self.data[0])])
        self.assertDictEqual(meta[0], meta_exp)
        pd.testing.assert_frame_equal(
            data, self.basic_data_exp, check_exact=False)

    def test_basic_ipdmqmc_input(self):
        """Test basic IPDMQMC input."""
        meta_exp = copy.deepcopy(self.metadata_ipdmqmc)
        (meta, data) = dmqmc.analyse_data(
            [(self.metadata_ipdmqmc, self.data[0])])
        self.assertDictEqual(meta[0], meta_exp)
        pd.testing.assert_frame_equal(
            data, self.basic_data_exp, check_exact=False)

    def test_multiple_input(self):
        """Test passing more than one calculation."""
        meta_exp = copy.deepcopy(self.metadata_dmqmc)
        (meta, data) = dmqmc.analyse_data(
            zip(2*[self.metadata_dmqmc], self.data))
        data_exp = pd.DataFrame([
            [0.0, -10.924401, 0.785440], [0.1, -10.679439, 0.998526],
            [0.2, -9.814926, 1.373230]
        ], columns=['Beta', 'Tr[Hp]/Tr[p]', 'Tr[Hp]/Tr[p]_error'], index=[
            0.0, 0.1, 0.2
        ])
        self.assertDictEqual(meta[0], meta_exp)
        pd.testing.assert_frame_equal(data, data_exp, check_exact=False)

    def test_shift(self):
        """Test the shift=True setting."""
        meta_exp = copy.deepcopy(self.metadata_dmqmc)
        (meta, data) = dmqmc.analyse_data(
            [(self.metadata_dmqmc, self.data[0])], shift=True
        )
        shift_data = pd.DataFrame([
            [-0.540397, 0.0136476], [-0.497985, 0.0795454],
            [-0.475570, 0.03676148]
        ], columns=['Shift', 'Shift s.d.'], index=[0.0, 0.1, 0.2])
        data_exp = pd.concat([self.basic_data_exp, shift_data], axis=1)
        self.assertDictEqual(meta[0], meta_exp)
        pd.testing.assert_frame_equal(
            data, data_exp, check_exact=False)

    def test_free_energy_dmqmc(self):
        """Test the free_energy=True setting."""
        meta_exp = copy.deepcopy(self.metadata_dmqmc)
        data_extra_exp = pd.DataFrame([
            [-11.374353, 1.325053, 0.0, 0.0],
            [-11.741264, 2.124244, -57.789041, 8.645126],
            [-9.134270, 2.904695, -109.977876, 21.318641]
        ], columns=['VI', 'VI_error', 'f_xc', 'f_xc_error'], index=[
            0.0, 0.1, 0.2
        ])
        data_exp = pd.concat(
            [self.basic_data_exp, data_extra_exp], axis=1)
        # Confusingly, in analysis.free_energy_error_analysis which is
        # called when free_energy=True this is done (see below) which
        # seems inconsistent! Fix? [todo] (index was [0.0, 0.1, 0.2],
        # then becomes [0, 1, 2]).
        data_exp.reset_index(drop=True, inplace=True)
        (meta, data) = dmqmc.analyse_data(
            [(self.metadata_dmqmc, self.data[0])], free_energy=True)
        self.assertDictEqual(meta[0], meta_exp)
        pd.testing.assert_frame_equal(
            data, data_exp, check_exact=False)

    def test_free_energy_ipdmqmc_sym(self):
        """Test free_energy=True for IPDMQMC."""
        # Need extra columns in data, r'\sum\rho_{ij}HI{ji}' and
        # r'\sum\rho_{ij}H0{ji}'.
        # meta['ipdmqmc']['symmetric'] = True.
        extra = create_mock_df.create_qmc_frame(
            self.rng, [r'\sum\rho_{ij}HI{ji}', r'\sum\rho_{ij}H0{ji}'],
            [89.0, 1987.0], [6.7, 2.3], [6.0, 293.0],
            num_mc_its=2*self.num_mc_its//2)
        # These have to inserted between 'Shift' and '# H psips'.
        data_dropped_Hpsips = self.data[0].drop(columns=['# H psips'])
        data0 = pd.concat([
            data_dropped_Hpsips, extra, self.data[0]['# H psips']
        ], axis=1)
        meta_exp = copy.deepcopy(self.metadata_ipdmqmc)
        data_extra_exp = pd.DataFrame([
            [1.568324, 0.406441, 0.079436, 0.005401544, -1.488888, 0.41184216,
             0.0, 0.0],
            [1.775345, 0.039159, 0.087249, 0.006767349, -1.688096, 0.032391712,
             -7.942460, 1.110925],
            [1.95467955e+00, 4.95612310e-01, 7.94175931e-02, 1.50709654e-02,
             -1.87526196e+00, 4.80541345e-01, -1.68508558e+01, 2.40671154e+00]
        ], columns=[
            'Tr[H0p]/Tr[p]', 'Tr[H0p]/Tr[p]_error', 'Tr[HIp]/Tr[p]',
            'Tr[HIp]/Tr[p]_error', 'VI', 'VI_error', 'f_xc', 'f_xc_error'
        ], index=[0.0, 0.1, 0.2])
        data_exp = pd.concat(
            [self.basic_data_exp, data_extra_exp], axis=1)
        # Confusingly, in analysis.free_energy_error_analysis which is
        # called when free_energy=True this is done (see below) which
        # seems inconsistent! Fix? [todo] (index was [0.0, 0.1, 0.2],
        # then becomes [0, 1, 2]).
        data_exp.reset_index(drop=True, inplace=True)
        (meta, data) = dmqmc.analyse_data(
            [(self.metadata_ipdmqmc, data0)], free_energy=True)
        self.assertDictEqual(meta[0], meta_exp)
        pd.testing.assert_frame_equal(data, data_exp, check_exact=False)

    def test_free_energy_ipdmqmc_notsym(self):
        """Test free_energy=True for IPDMQMC with symmetric=False."""
        # Need an extra column in data, r'\sum\rho_{ij}H0{ji}'.
        # meta['ipdmqmc']['symmetric'] = False.
        extra = create_mock_df.create_qmc_frame(
            self.rng, [r'\sum\rho_{ij}H0{ji}'], [1987.0], [2.3], [293.0],
            num_mc_its=2*self.num_mc_its//2)
        # These have to inserted between 'Shift' and '# H psips'.
        data_dropped_Hpsips = self.data[0].drop(columns=['# H psips'])
        data0 = pd.concat([
            data_dropped_Hpsips, extra, self.data[0]['# H psips']
        ], axis=1)
        meta_mock = copy.deepcopy(self.metadata_ipdmqmc)
        meta_mock['ipdmqmc']['symmetric'] = False
        data_extra_exp = pd.DataFrame([
            [1.60219598, 0.01637836, -12.97654889, 1.30867505, 0.0, 0.0],
            [1.9349722, 0.42732628, -13.67623586, 2.55157006, -66.63196188,
             9.67689327],
            [1.74821441, 0.48851518, -10.88248452, 3.39321059, -128.02876281,
             24.65774636]
        ], columns=[
            'Tr[H0p]/Tr[p]', 'Tr[H0p]/Tr[p]_error',
            'VI', 'VI_error', 'f_xc', 'f_xc_error'
        ], index=[0.0, 0.1, 0.2])
        data_exp = pd.concat(
            [self.basic_data_exp, data_extra_exp], axis=1)
        # Confusingly, in analysis.free_energy_error_analysis which is
        # called when free_energy=True this is done (see below) which
        # seems inconsistent! Fix? [todo] (index was [0.0, 0.1, 0.2],
        # then becomes [0, 1, 2]).
        data_exp.reset_index(drop=True, inplace=True)
        (meta, data) = dmqmc.analyse_data(
            [(meta_mock, data0)], free_energy=True)
        self.assertDictEqual(meta[0], meta_mock)
        pd.testing.assert_frame_equal(data, data_exp, check_exact=False)

    def test_spline(self):
        """Test spline=True option.  BROKEN!"""
        # This option is broken so have added assert statements to code.
        self.assertRaises(
            ValueError, dmqmc.analyse_data,
            [(self.metadata_dmqmc, self.data[0])], spline=True)

    def test_trace(self):
        """Test trace=True option."""
        meta_exp = copy.deepcopy(self.metadata_dmqmc)
        (meta, data) = dmqmc.analyse_data(
            [(self.metadata_dmqmc, self.data[0])], trace=True)
        trace_data = pd.DataFrame([
            [1165.920215, 4.415579], [1135.143442, 102.783344],
            [1220.003334, 179.329455]
        ], columns=['Trace', 'Trace s.d.'], index=[0.0, 0.1, 0.2])
        data_exp = pd.concat([self.basic_data_exp, trace_data], axis=1)
        self.assertDictEqual(meta[0], meta_exp)
        pd.testing.assert_frame_equal(
            data, data_exp, check_exact=False)

    def test_calc_number(self):
        """Test calc_number option.  BROKEN!"""
        # This option is broken so have added assert statements to code.
        self.assertRaises(
            ValueError, dmqmc.analyse_data,
            [(self.metadata_dmqmc, self.data[0])], calc_number=1)

    def test_unchanged_mutable(self):
        """Check that mutable objects, such as pd DataFrames, don't
        change when they shouldn't.
        """
        data_copy = self.data.copy()
        meta_exp = copy.deepcopy(self.metadata_dmqmc)
        _ = dmqmc.analyse_data(zip(2*[self.metadata_dmqmc], self.data))
        for i in range(len(self.data)):
            pd.testing.assert_frame_equal(
                self.data[i], data_copy[i], check_exact=True)
        self.assertDictEqual(self.metadata_dmqmc, meta_exp)

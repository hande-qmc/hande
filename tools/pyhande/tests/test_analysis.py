"""Tests for analysis.py."""
import unittest
import warnings
import numpy as np
import pandas as pd
import pyhande.analysis as analysis
import create_mock_df


class TestProjectedEnergy(unittest.TestCase):
    """Test analysis.projected_energy."""

    def setUp(self):
        # Create mock correlated noisy data where means gradually get
        # less noisy, standard errors find a plateau and their errors
        # stay very comparable.
        # Create mock reblock, covariance, data_len dataframes:
        rng = np.random.default_rng(3721)
        cols = [r'\sum H_0j N_j', 'N_0', 'alt1', 'alt2']
        means = [-100.0, 20.0, 1000.0, -1.0]
        it_optbls = [7, 6, 6, 5]
        num_reblock_its = 10
        self.reblock_orig = create_mock_df.create_reblock_frame(
            rng, cols, means, it_optbls, num_reblock_its
        )
        self.cov_orig = create_mock_df.create_cov_frame(
            rng, cols, means, num_reblock_its
        )
        self.data_len_orig = create_mock_df.create_date_length_series(
            rng, num_reblock_its
        )

    def test_basic_input(self):
        """Test basic input."""
        proje = analysis.projected_energy(
            self.reblock_orig, self.cov_orig, self.data_len_orig
        )
        proje_exp_means = [
            -5.70972199, -4.86906699, -5.33661295, -5.45374545, -5.17632254,
            -4.94604358, -5.03227561, -4.89093308, -4.94371799, -4.90833319
        ]
        proje_exp_stes = [
            0.07070339, 0.23508274, 0.38226895, 0.52410126, 0.4916993,
            0.52202386, 1.2421994, 1.82672645, 2.80415653, 3.52565766
        ]
        proje_exp_optb = ['' if i != 7 else '<---' for i in range(10)]
        proje_exp_df = pd.DataFrame.from_dict({
            'mean': proje_exp_means,
            'standard error': proje_exp_stes,
            'optimal block': proje_exp_optb
        })
        proje_exp_df = pd.concat(
            [proje_exp_df], axis=1, keys=['Proj. Energy']
        )
        proje_exp_df.index.name = 'reblock'
        pd.testing.assert_frame_equal(proje, proje_exp_df, check_exact=False)

    def test_custom_cols(self):
        """Test input with different columns and output col_name."""
        proje = analysis.projected_energy(
            self.reblock_orig, self.cov_orig, self.data_len_orig,
            sum_key='alt1', ref_key='alt2', col_name='Proje'
        )
        proje_exp_means = [
            -1276.45810774, -909.18258529, -1061.17358148, -1087.54401339,
            -1108.61543336, -1184.03733717, -1037.86512069, -991.42866714,
            -1024.23417293, -990.97344913
        ]
        proje_exp_stes = [
            14.91741087, 43.4508022, 73.18694845, 77.42004437, 95.91400323,
            143.12724715, 187.12002286, 279.46651658, 102.91404468,
            678.55618338
        ]
        proje_exp_optb = ['' if i != 6 else '<---' for i in range(10)]
        proje_exp_df = pd.DataFrame.from_dict({
            'mean': proje_exp_means,
            'standard error': proje_exp_stes,
            'optimal block': proje_exp_optb
        })
        proje_exp_df = pd.concat(
            [proje_exp_df], axis=1, keys=['Proje']
        )
        proje_exp_df.index.name = 'reblock'
        pd.testing.assert_frame_equal(proje, proje_exp_df, check_exact=False)

    def test_zero_in_ref(self):
        """Test input with a zero entry in the ref_key ('N_0')."""
        reblock_orig_copy = self.reblock_orig.copy()
        reblock_orig_copy.replace(
            reblock_orig_copy['N_0']['mean'].iloc[1], 0.0, inplace=True
        )
        proje = analysis.projected_energy(
            reblock_orig_copy, self.cov_orig, self.data_len_orig
        )
        self.assertEqual(abs(proje['Proj. Energy']['mean'].iloc[1]), np.inf)
        self.assertEqual(
            abs(proje['Proj. Energy']['standard error'].iloc[1]), np.inf
        )

    def test_unchanged_mutable(self):
        """Check that mutable objects, such as pd DataFrames, don't
        change when they shouldn't.
        """
        reblock_orig_copy = self.reblock_orig.copy()
        cov_orig_copy = self.cov_orig.copy()
        data_len_orig_copy = self.data_len_orig.copy()
        _ = analysis.projected_energy(
            self.reblock_orig, self.cov_orig, self.data_len_orig
        )
        pd.testing.assert_frame_equal(
            self.reblock_orig, reblock_orig_copy, check_exact=True
        )
        pd.testing.assert_frame_equal(
            self.cov_orig, cov_orig_copy, check_exact=True
        )
        pd.testing.assert_series_equal(
            self.data_len_orig, data_len_orig_copy, check_exact=True
        )


class TestQMCSummary(unittest.TestCase):
    """Test analysis.qmc_summary"""

    def setUp(self):
        # Create mock reblock dataframe:
        rng = np.random.default_rng(345621)
        cols = [
            r'\sum H_0j N_j', 'N_0', 'Shift', 'Proj. Energy', 'alt1', 'alt2'
        ]
        means = [-231.0, 10.0, -2.3, -2.3, 1000004, 0.0002]
        num_reblock_its = 10
        it_optbls = [8, 6, num_reblock_its, 8, 0, num_reblock_its]
        self.reblock_orig = create_mock_df.create_reblock_frame(
            rng, cols, means, it_optbls, num_reblock_its
        )
        # Also set up some shared exp results.
        self.summary_exp_0 = pd.DataFrame(
            np.asarray([[-243.841, 10.564, 0.258305],
                        [10.0446, 0.347958, 0.0483099],
                        [-2.41226, 0.105098, None]]),
            columns=['mean', 'standard error', 'standard error error'],
            index=[r'\sum H_0j N_j', 'N_0', 'Proj. Energy']
        )
        self.summary_exp_1 = ['Shift']

    def test_basic_input(self):
        """Test basic input."""
        summary = analysis.qmc_summary(self.reblock_orig)
        pd.testing.assert_frame_equal(
            summary[0], self.summary_exp_0, check_exact=False
        )
        self.assertEqual(summary[1], self.summary_exp_1)

    def test_custom_cols(self):
        """Use custom values."""
        summary = analysis.qmc_summary(self.reblock_orig)
        summary_cus = analysis.qmc_summary(
            self.reblock_orig, keys=['alt1', 'alt2'], summary_tuple=summary
        )
        summary_exp_0_add = pd.DataFrame(
            np.asarray([[1.14092e+06, 5010.8, 404.726]]),
            columns=['mean', 'standard error', 'standard error error'],
            index=['alt1']
        )
        summary_exp_0_add = pd.concat([
            self.summary_exp_0, summary_exp_0_add
        ])
        self.summary_exp_1.append('alt2')
        pd.testing.assert_frame_equal(
            summary_cus[0], summary_exp_0_add, check_exact=False
        )
        self.assertEqual(summary_cus[1], self.summary_exp_1)

    def test_unchanged_mutable(self):
        """Check that mutable objects, such as pd DataFrames, don't
        change when they shouldn't.
        """
        reblock_orig_copy = self.reblock_orig.copy()
        _ = analysis.qmc_summary(self.reblock_orig)
        pd.testing.assert_frame_equal(
            self.reblock_orig, reblock_orig_copy, check_exact=True
        )


class TestExtractPopGrowth(unittest.TestCase):
    """Test analysis.extract_pop_growth"""

    def setUp(self):
        # Create copy qmc dataframe
        rng = np.random.default_rng(63105)
        self.cols = ['N_0', 'Shift', 'alt1', 'alt2']
        means = [100.0, -11.0, -0.01, 2.0]
        sine_periods = [5.0, 2.9, 11.0, 8.0]
        noise_facs = [0.1, 0.9, 0.001, 0.2]
        self.data = create_mock_df.create_qmc_frame(
            rng, self.cols, means, sine_periods, noise_facs,
            frac_not_convergeds=[0.5 for _ in range(5)]
        )

    def test_basic_input_shift_immediate_variation(self):
        """Test basic input - shift varied immediately, get only first
        row!
        [todo] - Not sure whether this effect is intended!!!
        """
        grow_data = analysis.extract_pop_growth(self.data, min_ref_pop=-199.0)
        grow_data_exp = pd.DataFrame(
            np.asarray([[0.968889], [-1.059212], [-0.984942], [0.99834]]).T,
            columns=self.cols
        )
        pd.testing.assert_frame_equal(
            grow_data, grow_data_exp, check_exact=False
        )

    def test_basic_input(self):
        """Test basic input 2."""
        self.data.replace(
            self.data['Shift'].loc[0:8].values, 0.0, inplace=True
        )
        grow_data = analysis.extract_pop_growth(self.data)
        grow_data_exp = pd.DataFrame(
            np.asarray([[13.030125], [0.0], [-0.163027], [1.474782]]).T,
            index=[8], columns=self.cols
        )
        pd.testing.assert_frame_equal(
            grow_data, grow_data_exp, check_exact=False
        )

    def test_custom_input(self):
        """Test custom input."""
        self.data.replace(
            self.data['Shift'].loc[0:8].values, 0.0, inplace=True
        )
        self.data.replace(self.data['alt1'].loc[0:8].values, 0.0, inplace=True)
        grow_data = analysis.extract_pop_growth(
            self.data, ref_key='alt2', shift_key='alt1', min_ref_pop=1.42
        )
        grow_data_exp = pd.DataFrame(
            np.asarray([[13.030125], [0.0], [0.0], [1.474782]]).T,
            index=[8], columns=self.cols
        )
        pd.testing.assert_frame_equal(
            grow_data, grow_data_exp, check_exact=False
        )

    def test_unchanged_mutable(self):
        """Check that mutable objects, such as pd DataFrames, don't
        change when they shouldn't.
        """
        data_copy = self.data.copy()
        _ = analysis.extract_pop_growth(self.data)
        pd.testing.assert_frame_equal(self.data, data_copy, check_exact=True)


class TestPlateauEstimator(unittest.TestCase):
    """Test analysis.plateau_estimator"""

    def setUp(self):
        # Create mock qmc dataframe
        rng = np.random.default_rng(94)
        cols = ['N_0', '# H psips', 'Shift', 'alt1', 'alt2']
        means = [1000.0, 20000.0, -0.7, 10000.0, 20000.0]
        sine_periods = [5.0, 2.9, 11.0, 10.1, 3.12]
        noise_facs = [0.1, 0.5, 0.001, 0.1, 0.5]
        self.data = create_mock_df.create_qmc_frame(
            rng, cols, means, sine_periods, noise_facs,
            frac_not_convergeds=[0.6 for _ in range(5)], num_mc_its=70
        )
        self.data.replace(
            self.data['Shift'].loc[:50].values, 0.0, inplace=True
        )

    def test_basic_input(self):
        """Test basic input."""
        shoulder = analysis.plateau_estimator(self.data)
        shoulder_exp = pd.DataFrame(
            np.asarray([[19.278545, 0.928348], [13253.874811, 2095.242403]]),
            columns=['mean', 'standard error'],
            index=['shoulder estimator', 'shoulder height']
        )
        pd.testing.assert_frame_equal(
            shoulder, shoulder_exp, check_exact=False
        )

    def test_custom_input(self):
        """Test custom input."""
        shoulder = analysis.plateau_estimator(
            self.data, total_key='alt2', ref_key='alt1',
            pop_data=self.data.loc[:50]
        )
        # [todo] Note that since pop_data is specified, min_ref_pop and
        # [todo] shift_key have no effect! Is that safe? Not tested for.
        shoulder_exp = pd.DataFrame(
            np.asarray([[12.334566, 0.656519], [11910.146569, 2345.367031]]),
            columns=['mean', 'standard error'],
            index=['shoulder estimator', 'shoulder height']
        )
        pd.testing.assert_frame_equal(
            shoulder, shoulder_exp, check_exact=False
        )

    def test_unchanged_mutable(self):
        """Check that mutable objects, such as pd DataFrames, don't
        change when they shouldn't.
        """
        data_copy = self.data.copy()
        # [todo] Check pop_data once this is fixed.
        pop_data = self.data.loc[:50].copy()  # [todo] check this copy()
        # pop_data_copy = pop_data.copy()
        _ = analysis.plateau_estimator(self.data, pop_data=pop_data)
        pd.testing.assert_frame_equal(self.data, data_copy, check_exact=True)
        # [todo] Is there a better way to do this warning?
        warnings.warn("TestPlateauEstimator.test_unchanged_mutable: "
                      "Mutable pop_data is changed in function! Fix?")


class TestPlateauEstimatorHist(unittest.TestCase):
    """Test analysis.plateau_estimator_hist"""

    def setUp(self):
        # Create mock qmc dataframe
        rng = np.random.default_rng(94)
        cols = ['N_0', '# H psips', 'Shift', 'alt1']
        means = [1000.0, 20000.0, -0.7, 10000.0, 20000.0]
        sine_periods = [5.0, 2.9, 11.0, 10.1]
        noise_facs = [0.1, 0.5, 0.001, 0.1]
        self.data = create_mock_df.create_qmc_frame(
            rng, cols, means, sine_periods, noise_facs,
            frac_not_convergeds=[0.6 for _ in range(4)], num_mc_its=1000
        )
        # Set initial shift iterations to zero.
        self.data.replace(
            self.data['Shift'].loc[0:520].values, 0.0, inplace=True
        )
        # Add an artificial plateau in the growing phase
        self.data.replace(
            self.data['# H psips'].loc[400:450].values,
            self.data['# H psips'].loc[420], inplace=True
        )

    def test_basic_input(self):
        """Test basic input."""
        shoulder = analysis.plateau_estimator_hist(self.data)
        self.assertAlmostEqual(shoulder, 12266.0527266)

    def test_custom_bin_width_fn(self):
        """Test custom bin_width_fn function."""
        shoulder = analysis.plateau_estimator_hist(
            self.data, bin_width_fn=lambda x: 11.0/len(x)
        )
        self.assertAlmostEqual(shoulder, 11953.807495405279)

    def test_custom_part_of_pop_data(self):
        """Test some custom input (only pass a bit of data)."""
        shoulder = analysis.plateau_estimator_hist(
            self.data, total_key='alt1', pop_data=self.data.loc[:150]
        )
        if np.isnan(shoulder):  # Type conversion
            shoulder = None
        self.assertIsNone(shoulder)

    def test_custom_more_pop_data(self):
        """Test some custom input (pass more data)
        (but don't pass all data where shift is zero to actually test
        this custom data option)
        """
        shoulder = analysis.plateau_estimator_hist(
            self.data, total_key='alt1', pop_data=self.data.loc[:518]
        )
        self.assertAlmostEqual(shoulder, 1.1722807735)

    def test_unchanged_mutable(self):
        """Check that mutable objects, such as pd DataFrames, don't
        change when they shouldn't.
        """
        data_copy = self.data.copy()
        _ = analysis.plateau_estimator_hist(self.data)
        pd.testing.assert_frame_equal(self.data, data_copy, check_exact=True)


class TestInefficiency(unittest.TestCase):
    """Test analysis.inefficiency"""

    def setUp(self):
        # Create mock opt_block
        self.opt_block = pd.DataFrame(
            np.asarray([[-243.841, 10.564, 0.258305],
                        [10.0446, 0.347958, 0.0483099],
                        [1023.0, 23.32, 1.34],
                        [-2.41226, 0.105098, None]]),
            columns=['mean', 'standard error', 'standard error error'],
            index=[r'\sum H_0j N_j', 'N_0', '# H psips', 'Proj. Energy']
        )

    def test_basic_input(self):
        """Basic inputs."""
        (dtau, iterations) = (0.1, 10000)
        ineff = analysis.inefficiency(self.opt_block, dtau, iterations)
        ineff_exp = pd.DataFrame(
            np.asarray([[106.299756, 15.034506]]),
            columns=['mean', 'standard error'],
            index=['Inefficiency']
        )
        pd.testing.assert_frame_equal(ineff, ineff_exp, check_exact=False)

    def test_zero_division(self):
        """Basic inputs 2 - What happens with a division by zero?"""
        (dtau, iterations) = (0.1, 10000)
        self.opt_block.replace(
            self.opt_block['standard error']['N_0'], 0.0, inplace=True
        )
        ineff = analysis.inefficiency(self.opt_block, dtau, iterations)
        ineff_exp = pd.DataFrame(
            np.asarray([[106.299756, np.inf]]),
            columns=['mean', 'standard error'],
            index=['Inefficiency']
        )
        pd.testing.assert_frame_equal(ineff, ineff_exp, check_exact=False)

    def test_key_error(self):
        """Key error and associated warning being printed."""
        self.opt_block.drop(index='N_0', inplace=True)
        (dtau, iterations) = (0.1, 10000)
        ineff = analysis.inefficiency(self.opt_block, dtau, iterations)
        self.assertWarnsRegex(
            UserWarning, r'Inefficiency not calculated owing to data '
            'unavailable from \'N_0\'', analysis.inefficiency,
            self.opt_block, dtau, iterations
        )
        self.assertIsNone(ineff)

    def test_unchanged_mutable(self):
        """Check that mutable objects, such as pd DataFrames, don't
        change when they shouldn't.
        """
        opt_block_copy = self.opt_block.copy()
        (dtau, iterations) = (0.1, 10000)
        _ = analysis.inefficiency(self.opt_block, dtau, iterations)
        pd.testing.assert_frame_equal(
            self.opt_block, opt_block_copy, check_exact=True
        )

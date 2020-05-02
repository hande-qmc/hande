"""Tests for lazy.py."""
import unittest
import warnings
import copy
import pandas as pd
import numpy as np
import pyhande.lazy as lazy
import create_mock_df


class TestFindStartingIterationMserMin(unittest.TestCase):
    """Test lazy.find_starting_iteration_mser_min."""

    def setUp(self):
        # Create mock qmc dataframe:
        rng = np.random.default_rng(769)
        cols = [
            r'\sum H_0j N_j', 'N_0', 'Shift', 'Proj. Energy', 'alt1', 'alt2'
        ]
        means = [-231.0, 10.0, -2.3, -2.3, 1000004, 0.0002]
        sine_periods = [4, 2, 6, 2, 6, 7]
        noise_facs = [0.1*mean for mean in means]
        self.data = create_mock_df.create_qmc_frame(
            rng, cols, means, sine_periods, noise_facs,
            frac_not_convergeds=[0.33333 for _ in range(5)], num_mc_its=300
        )
        iterations = pd.DataFrame(list(range(1, 1501, 5)),
                                  columns=['iterations'])
        self.data = pd.concat([iterations, self.data], axis=1)
        self.md = {'sys': {'nbasis': 59}, 'qmc': {'ncycles': 5, 'tau': 0.1}}

    def test_basic_input(self):
        """Test basic input."""
        start_it = lazy.find_starting_iteration_mser_min(self.data, self.md)
        # Note that we might expect the start to lie at 500
        # (0.33333*1500) as it was specified that the data is converged
        # after 0.33333 of the data ("frac_not_converged").  Close!
        self.assertEqual(start_it, 470)

    def test_start_max_frac(self):
        """Test start_max_frac setting."""
        start_it = lazy.find_starting_iteration_mser_min(
            self.data, self.md, start_max_frac=0.2)
        self.assertEqual(start_it, 295)  # Not a good setting...
        self.assertWarnsRegex(
            UserWarning, 'Proj. energy may not be converged. MSER min. may '
            'underestimate the starting iteration. One should check 1:\$3/\$4 '
            'plot.', lazy.find_starting_iteration_mser_min, self.data, self.md,
            start_max_frac=0.2)

    def test_n_blocks(self):
        """Test n_blocks setting."""
        start_it = lazy.find_starting_iteration_mser_min(
            self.data, self.md, n_blocks=30)
        self.assertEqual(start_it, 495)

    def test_verbose(self):
        """Test using a verbose setting (should not do anything)."""
        start_it = lazy.find_starting_iteration_mser_min(
            self.data, self.md, verbose=2)
        self.assertEqual(start_it, 470)

    def test_end(self):
        """Test end setting."""
        start_it = lazy.find_starting_iteration_mser_min(
            self.data, self.md, end=460)
        self.assertEqual(start_it, 405)

    def test_unchanged_mutable(self):
        """Check that mutable objects, such as pd DataFrames, don't
        change when they shouldn't.
        """
        data_copy = self.data.copy()
        md_copy = copy.deepcopy(self.md)
        _ = lazy.find_starting_iteration_mser_min(self.data, self.md)
        pd.testing.assert_frame_equal(self.data, data_copy, check_exact=True)
        self.assertDictEqual(self.md, md_copy)


class TestLazyHybrid(unittest.TestCase):
    """Test lazy.lazy_hybrid()."""

    def setUp(self):
        # Create mock reblock dataframe:
        rng = np.random.default_rng(9843)
        cols = [
            r'\sum H_0j N_j', 'N_0', 'Shift', 'Proj. Energy', 'alt1', 'alt2'
        ]
        means = [-231.0, 10.0, -2.3, -2.3, 1000004, 0.0002]
        sine_periods = [4, 2, 6, 2, 6, 7]
        noise_facs = [0.1*mean for mean in means]
        self.data = create_mock_df.create_qmc_frame(
            rng, cols, means, sine_periods, noise_facs,
            frac_not_convergeds=[0.33333 for _ in range(5)], num_mc_its=300
        )
        iterations = pd.DataFrame(list(range(1, 1501, 5)),
                                  columns=['iterations'])
        self.data = pd.concat([iterations, self.data], axis=1)
        self.md = {'sys': {'nbasis': 59}, 'qmc': {'ncycles': 5, 'tau': 0.1}}

    def test_basic_input(self):
        """Test basic input."""
        info = lazy.lazy_hybrid(self.data, self.md)
        opt_block_exp = pd.DataFrame(
            [[-1.970565, 0.25349, None, '-2.0(3)']], columns=[
                'mean', 'standard error', 'standard error error', 'estimate'
            ], index=['Proj. Energy'])
        pd.testing.assert_frame_equal(opt_block_exp, info.opt_block)

    def test_start(self):
        """Test start parameter.  Set to close to one found above."""
        info = lazy.lazy_hybrid(self.data, self.md, start=500)
        opt_block_exp = pd.DataFrame(
            [[-2.25758, 0.0241037, None, '-2.26(2)']], columns=[
                'mean', 'standard error', 'standard error error', 'estimate'
            ], index=['Proj. Energy'])
        pd.testing.assert_frame_equal(opt_block_exp, info.opt_block)

    def test_end(self):
        """Test end parameter.  Set to close to one found above."""
        info = lazy.lazy_hybrid(self.data, self.md, end=1023)
        opt_block_exp = pd.DataFrame(
            [[-1.81619, 0.314786, None, '-1.8(3)']], columns=[
                'mean', 'standard error', 'standard error error', 'estimate'
            ], index=['Proj. Energy'])
        pd.testing.assert_frame_equal(opt_block_exp, info.opt_block)

    def test_batch_size(self):
        """Test batch_size parameter."""
        info = lazy.lazy_hybrid(self.data, self.md, batch_size=6)
        opt_block_exp = pd.DataFrame(
            [[-1.970565, 0.30242585, None, '-2.0(3)']], columns=[
                'mean', 'standard error', 'standard error error', 'estimate'
            ], index=['Proj. Energy'])
        pd.testing.assert_frame_equal(opt_block_exp, info.opt_block)

    def test_unchanged_mutable(self):
        """Check that mutable objects, such as pd DataFrames, don't
        change when they shouldn't.
        """
        data_copy = self.data.copy()
        md_copy = copy.deepcopy(self.md)
        _ = lazy.find_starting_iteration_mser_min(self.data, self.md)
        pd.testing.assert_frame_equal(self.data, data_copy, check_exact=True)
        self.assertDictEqual(self.md, md_copy)


class TestStdAnalysis(unittest.TestCase):
    """Test lazy.std_analysis().
    This also (indirectly) tests zeroT_qmc, lazy.lazy_block,
    lazy.filter_calcs, lazy.concat_calcs very effectively, so at the
    moment they don't have separate tests.
    Note that for clarity and simplicity, only a part of the (long)
    outputs is checked.  verbosity is set to -1 and is not tested for.
    """

    def test_basic_input(self):
        """Test basic input."""
        infos = lazy.std_analysis(
            ['tests/hande_files/long_calc_ueg.out'], verbosity=-1)
        # Length of list.
        self.assertEqual(len(infos), 1)
        info = infos[0]
        # .metadata
        meta_fciqmc = {
            'select_ref_det_every_nreports': 2147483647,
            'init_spin_inv_D0': False, 'ref_det_factor': 1.5,
            'non_blocking_comm': False, 'doing_load_balancing': False,
            'trial_function': 'single_basis', 'guiding_function': 'none',
            'quadrature_initiator': True, 'replica_tricks': False
        }
        self.assertDictEqual(info.metadata['fciqmc'], meta_fciqmc)
        meta_pyhande = {'reblock_start': 3736}
        self.assertDictEqual(info.metadata['pyhande'], meta_pyhande)
        # .data
        data_it_4340 = pd.Series([
            4.34000000e+03, -2.42034087e-01, -3.88598174e+01, 1.36800000e+02,
            1.18000000e+04, 1.11860000e+04, 1.21700000e+03, 9.33000000e-02,
            7.60000000e-03, -2.84062993e-01
        ], index=[
            'iterations', 'Shift', r'\sum H_0j N_j', 'N_0', '# H psips',
            '# states', '# spawn_events', 'R_spawn', 'time',
            'Proj. Energy'
        ], name=433)
        pd.testing.assert_series_equal(info.data.loc[433], data_it_4340)
        # .data_len
        data_len_exp = pd.Series(
            [2627, 1313, 656, 328, 164, 82, 41, 20, 10, 5, 2],
            name='data length')
        data_len_exp.index.name = 'reblock'
        pd.testing.assert_series_equal(info.data_len, data_len_exp)
        # .reblock
        reblock_shift_loc8 = pd.Series([
            -0.2735138312473437, 0.0034409665708355465, 0.0008110435986913453,
            '<---    '
        ], index=[
            'mean', 'standard error', 'standard error error',
            'optimal block'
        ], name=8)
        pd.testing.assert_series_equal(
            info.reblock['Shift'].loc[8], reblock_shift_loc8)
        # .covariance
        reblock_2_N_0 = pd.Series(
            [-34.893919, 135.700242, -0.238517],
            index=[r'\sum H_0j N_j', 'N_0', 'Shift'], name=(2, 'N_0'))
        pd.testing.assert_series_equal(
            info.covariance.loc[2, 'N_0'], reblock_2_N_0)
        # .opt_block
        opt_block_proje = pd.Series(
            [-0.278837, 0.000568666, None, '-0.2788(6)'], index=[
                'mean', 'standard error', 'standard error error', 'estimate'
            ], name='Proj. Energy')
        pd.testing.assert_series_equal(
            info.opt_block.loc['Proj. Energy'], opt_block_proje)
        # .no_opt_block
        self.assertListEqual(info.no_opt_block, [])

    def test_concat(self):
        """Test concatenation of (restarted) calculations.
        Have run the same calculation twice; once in one file and the
        other with restarting in the middle, consisting of two files.
        The calculation data content should be identical except for
        possibly the time.
        """
        infos_one = lazy.std_analysis(
            ['tests/hande_files/long_calc_ueg.out'], verbosity=-1)
        infos_two = lazy.std_analysis(
            ['tests/hande_files/long_calc_split1_ueg.out',
             'tests/hande_files/long_calc_split2_ueg.out'], verbosity=-1)
        self.assertEqual(len(infos_two), 1)
        data_one_notime = infos_one[0].data.drop(columns=['time'])
        data_two_notime = infos_two[0].data.drop(columns=['time'])
        pd.testing.assert_frame_equal(
            data_one_notime, data_two_notime, check_exact=True)
        pd.testing.assert_series_equal(
            infos_one[0].data_len, infos_two[0].data_len, check_exact=True)
        pd.testing.assert_frame_equal(
            infos_one[0].reblock, infos_two[0].reblock, check_exact=True)
        pd.testing.assert_frame_equal(
            infos_one[0].covariance, infos_two[0].covariance, check_exact=True)
        pd.testing.assert_frame_equal(
            infos_one[0].opt_block, infos_two[0].opt_block, check_exact=True)
        self.assertListEqual(
            infos_one[0].no_opt_block, infos_two[0].no_opt_block)

    def test_multiple(self):
        """Test analysing two different calculations at once."""
        infos = lazy.std_analysis(
            ['tests/hande_files/long_calc_ueg.out',
             'tests/hande_files/long_calc_split2_ueg.out'], verbosity=-1)
        self.assertEqual(len(infos), 2)
        # Sample test the first calculation.
        opt_block_proje = pd.Series(
            [-0.278837168316264, 0.0005686659462504128, None, '-0.2788(6)'],
            index=[
                'mean', 'standard error', 'standard error error', 'estimate'
            ], name='Proj. Energy')
        pd.testing.assert_series_equal(
            infos[0].opt_block.loc['Proj. Energy'], opt_block_proje)
        # Sample test the second calculation.
        opt_block_proje = pd.Series(
            [-0.27942919514874415, 0.0005544343214820229, None, '-0.2794(6)'],
            index=[
                'mean', 'standard error', 'standard error error', 'estimate'
            ], name='Proj. Energy')
        pd.testing.assert_series_equal(
            infos[1].opt_block.loc['Proj. Energy'], opt_block_proje)

    def test_fci(self):
        """Test FCI.  Not to be analysed here!"""
        self.assertRaises(
            ValueError, lazy.std_analysis, ['tests/hande_files/fci_ueg.out'])

    def test_start(self):
        """Test start parameter."""
        infos = lazy.std_analysis(
            ['tests/hande_files/long_calc_ueg.out'], start=12000, verbosity=-1)
        # .metadata
        meta_pyhande = {'reblock_start': 12000}
        self.assertDictEqual(infos[0].metadata['pyhande'], meta_pyhande)
        # .opt_block
        opt_block_proje = pd.Series(
            [-0.27809, 0.000933163, None, '-0.2781(9)'], index=[
                'mean', 'standard error', 'standard error error', 'estimate'
            ], name='Proj. Energy')
        pd.testing.assert_series_equal(
            infos[0].opt_block.loc['Proj. Energy'], opt_block_proje)

    def test_end(self):
        """Test end parameter."""
        infos = lazy.std_analysis(
            ['tests/hande_files/long_calc_ueg.out'], end=4400, verbosity=-1)
        # .metadata
        meta_pyhande = {'reblock_start': 3579}
        self.assertDictEqual(infos[0].metadata['pyhande'], meta_pyhande)
        # .no_opt_block
        self.assertListEqual(infos[0].no_opt_block, ['N_0', 'Proj. Energy'])

    def test_select_function(self):
        """Test select_function parameter.
        Similar to the example in docstring.
        """
        infos = lazy.std_analysis(
            ['tests/hande_files/long_calc_ueg.out'],
            select_function=lambda d: d['iterations'] > 12000, verbosity=-1)
        # .metadata
        warnings.warn(
            "Fix bug: reblock_start in info.metadata not correct with "
            "select_function parameter.")
        # .opt_block
        opt_block_proje = pd.Series(
            [-0.27809, 0.000933163, None, '-0.2781(9)'], index=[
                'mean', 'standard error', 'standard error error', 'estimate'
            ], name='Proj. Energy')
        pd.testing.assert_series_equal(
            infos[0].opt_block.loc['Proj. Energy'], opt_block_proje)

    def test_extract_psips(self):
        """Test extract_psips parameter."""
        infos = lazy.std_analysis(
            ['tests/hande_files/long_calc_ueg.out'], extract_psips=True,
            verbosity=-1)
        # .metadata
        meta_pyhande = {'reblock_start': 3736}
        self.assertDictEqual(infos[0].metadata['pyhande'], meta_pyhande)
        # .opt_block
        opt_block_psips = pd.Series([
            12584.401953125, 85.35691857330413, 20.118818648123774,
            '12580(90)'
        ], index=[
            'mean', 'standard error', 'standard error error', 'estimate'
        ], name='# H psips')
        pd.testing.assert_series_equal(
            infos[0].opt_block.loc['# H psips'], opt_block_psips)

    def test_reweight_history_mean_shift(self):
        """Test reweight_history and mean_shift parameters."""
        infos = lazy.std_analysis(
            ['tests/hande_files/long_calc_ueg.out'], reweight_history=3,
            mean_shift=-0.1, verbosity=-1)
        # .reblock
        reblock_wproje_loc8 = pd.Series(
            [-0.278843, 0.000569041, '<---    '],
            index=['mean', 'standard error', 'optimal block'], name=8)
        pd.testing.assert_series_equal(
            infos[0].reblock['Weighted Proj. E.'].loc[8], reblock_wproje_loc8)
        # .covariance
        reblock_3_W_N_0 = pd.Series(
            [-38.69204990933227, 148.59408342412507, -0.27500051571542944,
             -43.852041898078305, 167.62001388833082], index=[
                 r'\sum H_0j N_j', 'N_0', 'Shift', r'W * \sum H_0j N_j',
                 'W * N_0'
            ], name=(3, 'W * N_0'))
        pd.testing.assert_series_equal(
            infos[0].covariance.loc[3, 'W * N_0'], reblock_3_W_N_0)
        # .opt_block
        opt_block_proje = pd.Series(
            [-0.278842866213868, 0.0005690407973903144, None, '-0.2788(6)'],
            index=[
                'mean', 'standard error', 'standard error error', 'estimate'
            ], name='Weighted Proj. E.')
        # Check more precise (7 digits instead of default 5) here as
        # effects are tiny.
        pd.testing.assert_series_equal(
            infos[0].opt_block.loc['Weighted Proj. E.'], opt_block_proje,
            check_less_precise=7)

    def test_reweight_history_arith_mean(self):
        """Test reweight_history and arith_mean parameters."""
        infos = lazy.std_analysis(
            ['tests/hande_files/long_calc_ueg.out'], reweight_history=3,
            arith_mean=True, verbosity=-1)
        # .reblock
        reblock_wproje_loc8 = pd.Series(
            [-0.278843, 0.000569041, '<---    '],
            index=['mean', 'standard error', 'optimal block'], name=8)
        pd.testing.assert_series_equal(
            infos[0].reblock['Weighted Proj. E.'].loc[8], reblock_wproje_loc8)
        # .covariance
        cov_3_W_N_0 = pd.Series(
            [-39.81896416225519, 152.97450543779814, -0.28322030786266594,
             -46.46286890224775, 177.66923750972612], index=[
                 r'\sum H_0j N_j', 'N_0', 'Shift', r'W * \sum H_0j N_j',
                 'W * N_0'
            ], name=(3, 'W * N_0'))
        pd.testing.assert_series_equal(
            infos[0].covariance.loc[3, 'W * N_0'], cov_3_W_N_0)
        # .opt_block
        opt_block_proje = pd.Series(
            [-0.27884176150306245, 0.0005690386058169649, None, '-0.2788(6)'],
            index=[
                'mean', 'standard error', 'standard error error', 'estimate'
            ], name='Weighted Proj. E.')
        # Check more precise (7 digits instead of default 5) here as
        # effects are tiny.
        pd.testing.assert_series_equal(
            infos[0].opt_block.loc['Weighted Proj. E.'], opt_block_proje,
            check_less_precise=7)

    def test_calc_inefficiency(self):
        """Test calc_inefficiency parameter.
        Requires extract_psips.
        """
        infos = lazy.std_analysis(
            ['tests/hande_files/long_calc_ueg.out'], extract_psips=True,
            calc_inefficiency=True, verbosity=-1)
        # .opt_block
        opt_block_ineff = pd.Series(
            [1.0337627369749376, 0.34460541317375093, None, '1.0(3)'],
            index=[
                'mean', 'standard error', 'standard error error', 'estimate'
            ], name='Inefficiency')
        pd.testing.assert_series_equal(
            infos[0].opt_block.loc['Inefficiency'], opt_block_ineff)

    def test_calc_inefficiency_no_psips(self):
        """Test calc_inefficiency parameter.
        Requires extract_psips but don't supply that here.
        """
        self.assertWarnsRegex(
            UserWarning, "Inefficiency not calculated owing to data "
            "unavailable from '# H psips'", lazy.std_analysis,
            ['tests/hande_files/long_calc_ueg.out'], extract_psips=False,
            calc_inefficiency=True, verbosity=-1)

    def test_starts_reweighting(self):
        """Test starts_reweighting parameter."""
        infos = lazy.std_analysis(
            ['tests/hande_files/long_calc_ueg.out'],
            starts_reweighting=[7634, 11000], verbosity=-1)
        # .metadata
        meta_pyhande = {'reblock_start': 7634}
        self.assertDictEqual(infos[0].metadata['pyhande'], meta_pyhande)
        # .opt_block
        opt_block_proje = pd.Series(
            [-0.2789050612630105, 0.0010293205541024083, None, '-0.279(1)'],
            index=[
                'mean', 'standard error', 'standard error error', 'estimate'
            ], name='Proj. Energy')
        pd.testing.assert_series_equal(
            infos[0].opt_block.loc['Proj. Energy'], opt_block_proje)

    def test_extract_rep_loop_time(self):
        """Test extract_rep_loop_time parameter."""
        infos = lazy.std_analysis(
            ['tests/hande_files/long_calc_ueg.out'],
            extract_rep_loop_time=True, verbosity=-1)
        # .metadata
        meta_pyhande = {'reblock_start': 3736}
        self.assertDictEqual(infos[0].metadata['pyhande'], meta_pyhande)
        # .opt_block
        opt_block_psips = pd.Series([
            0.008899375, 0.0002561202535778128, 4.1548191515764066e-05,
            '0.0089(3)'
        ], index=[
            'mean', 'standard error', 'standard error error', 'estimate'
        ], name='time')
        pd.testing.assert_series_equal(
            infos[0].opt_block.loc['time'], opt_block_psips)

    def test_analysis_method(self):
        """Test analysis_method parameter."""
        infos = lazy.std_analysis(
            ['tests/hande_files/long_calc_ueg.out'], analysis_method='hybrid',
            verbosity=-1)
        # .metadata
        meta_pyhande = {'reblock_start': 3736}
        self.assertDictEqual(infos[0].metadata['pyhande'], meta_pyhande)
        # .opt_block
        opt_block_proje = pd.DataFrame(
            [[-0.2790001401084949, 0.00055865144808341, None, '-0.2790(6)']],
            columns=[
                'mean', 'standard error', 'standard error error', 'estimate'
            ], index=['Proj. Energy'])
        pd.testing.assert_frame_equal(
            infos[0].opt_block, opt_block_proje)

    def test_analysis_method_error(self):
        """Test analysis_method parameter.  Wrong input."""
        self.assertRaises(
            ValueError, lazy.std_analysis,
            ['tests/hande_files/long_calc_ueg.out'], analysis_method='test',
            verbosity=-1)

    def test_warmup_detection(self):
        """Test warmup_detection parameter."""
        infos = lazy.std_analysis(
            ['tests/hande_files/long_calc_ueg.out'],
            warmup_detection='mser_min', verbosity=-1)
        # .metadata
        meta_pyhande = {'reblock_start': 1350}
        self.assertDictEqual(infos[0].metadata['pyhande'], meta_pyhande)
        # .opt_block
        opt_block_proje = pd.Series(
            [-0.2789132653871061, 0.0005685202640691406, None, '-0.2789(6)'],
            index=[
                'mean', 'standard error', 'standard error error', 'estimate'
            ], name='Proj. Energy')
        pd.testing.assert_series_equal(
            infos[0].opt_block.loc['Proj. Energy'], opt_block_proje)

    def test_warmup_detection_error(self):
        """Test warmup_detection parameter.  Wrong input."""
        self.assertRaises(
            ValueError, lazy.std_analysis,
            ['tests/hande_files/long_calc_ueg.out'], warmup_detection='test',
            verbosity=-1)


class TestFindStartingIteration(unittest.TestCase):
    """Test lazy.find_starting_iteration().
    'show_graph' and 'verbosity' are not tested.
    """

    def setUp(self):
        # Create mock qmc dataframe:
        rng = np.random.default_rng(33687)
        cols = [
            r'\sum H_0j N_j', 'N_0', 'Shift', 'Proj. Energy', '# H psips',
            'alt1'
        ]
        means = [-231.0, 10.0, -2.3, -2.3, 1000004, 0.0002]
        sine_periods = [4, 2, 6, 2, 6, 7]
        noise_facs = [0.1*mean for mean in means]
        self.data = create_mock_df.create_qmc_frame(
            rng, cols, means, sine_periods, noise_facs,
            frac_not_convergeds=[0.33333 for _ in range(5)], num_mc_its=300
        )
        iterations = pd.DataFrame(list(range(1, 1501, 5)),
                                  columns=['iterations'])
        self.data = pd.concat([iterations, self.data], axis=1)
        self.md = {'sys': {'nbasis': 59}, 'qmc': {'ncycles': 5, 'tau': 0.1}}

    def test_basic_input(self):
        """Test basic input."""
        start_it = lazy.find_starting_iteration(self.data, self.md)
        self.assertEqual(start_it, 898)

    def test_frac_screen_interval(self):
        """Test frac_screen_interval parameter."""
        start_it = lazy.find_starting_iteration(
            self.data, self.md, frac_screen_interval=100)
        self.assertEqual(start_it, 662)

    def test_frac_screen_interval_error(self):
        """Test frac_screen_interval parameter."""
        self.assertRaises(
            ValueError, lazy.find_starting_iteration, self.data, self.md,
            frac_screen_interval=-80)

    def test_number_of_reblockings(self):
        """Test number_of_reblockings parameter.
        [todo] Improve! Result is not distinguishable from default.
        """
        start_it = lazy.find_starting_iteration(
            self.data, self.md, number_of_reblockings=1)
        self.assertEqual(start_it, 898)

    def test_number_of_reblockings_error(self):
        """Test number_of_reblockings parameter."""
        self.assertRaises(
            ValueError, lazy.find_starting_iteration, self.data, self.md,
            number_of_reblockings=-80)

    def test_frac_screen_interval_number_of_reblockings_error(self):
        """Test further error."""
        self.assertRaises(
            ValueError, lazy.find_starting_iteration, self.data, self.md,
            frac_screen_interval=5, number_of_reblockings=67)

    def test_number_of_reblocks_to_cut_off(self):
        """Test number_of_reblocks_to_cut_off parameter."""
        start_it = lazy.find_starting_iteration(
            self.data, self.md, number_of_reblocks_to_cut_off=0)
        self.assertEqual(start_it, 858)

    def test_number_of_reblockings_to_cut_off_error(self):
        """Test number_of_reblocks_to_cut_off parameter."""
        self.assertRaises(
            ValueError, lazy.find_starting_iteration, self.data, self.md,
            number_of_reblocks_to_cut_off=-80)

    def test_pos_min_frac(self):
        """Test pos_min_frac parameter."""
        self.assertRaisesRegex(
            RuntimeError, "Failed to find starting iteration. The calculation "
            "might not be converged.", lazy.find_starting_iteration,
            self.data, self.md, pos_min_frac=0.7)

    def test_pos_min_frac_error(self):
        """Test pos_min_frac parameter."""
        self.assertRaises(
            ValueError, lazy.find_starting_iteration, self.data, self.md,
            pos_min_frac=-80)
        self.assertRaises(
            ValueError, lazy.find_starting_iteration, self.data, self.md,
            pos_min_frac=2.4)

    def test_end(self):
        """Test end parameter."""
        self.assertRaisesRegex(
            RuntimeError, "Failed to find starting iteration. The calculation "
            "might not be converged.", lazy.find_starting_iteration,
            self.data, self.md, end=145)

    def test_unchanged_mutable(self):
        """Check that mutable objects, such as pd DataFrames, don't
        change when they shouldn't.
        """
        data_copy = self.data.copy()
        md_copy = copy.deepcopy(self.md)
        _ = lazy.find_starting_iteration(self.data, self.md)
        pd.testing.assert_frame_equal(self.data, data_copy, check_exact=True)
        self.assertDictEqual(self.md, md_copy)

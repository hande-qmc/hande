"""Test pyhande.error_analysing.find_starting_iteration.py."""
import unittest
import numpy as np
import pandas as pd
import pyhande.error_analysing.find_starting_iteration as find_startit
from tests.create_mock_df import create_qmc_frame


class TestFindStartingIterationsBlocking(unittest.TestCase):
    """Test `find_starting_iteration_blocking()`."""
    # [todo] - test show_graph.

    def setUp(self):
        """Setting up common mock df."""
        rng = np.random.default_rng(6183)
        cols = ['Shift', 'test', 'N_0', r'\sum H_0j N_j', 'alt',
                'Proj. Energy']
        means = [-0.1, 10.0, 11.0, -9.0, 20.0, -0.2]
        sine_periods = list(range(1, len(cols)+1))
        noise_facs = [0.1*mean for mean in means]
        self._df_mock = create_qmc_frame(rng, cols, means, sine_periods,
                                         noise_facs,
                                         frac_not_convergeds=[0.2]*len(cols),
                                         num_mc_its=300)
        iterations = pd.DataFrame(list(range(1, 1501, 5)),
                                  columns=['iterations'])
        self._df_mock = pd.concat([iterations, self._df_mock], axis=1)

    def test_default_usage(self):
        """Test default usage."""
        start_it = find_startit.find_starting_iteration_blocking(
            self._df_mock,  self._df_mock['iterations'].iloc[-1], 'iterations',
            ['Shift', 'test', 'N_0', r'\sum H_0j N_j', 'alt'], None)
        self.assertEqual(start_it, 491)

    def test_specifying_smaller_end_its(self):
        """Specify smaller end iteration."""
        start_it = find_startit.find_starting_iteration_blocking(
            self._df_mock, 800, 'iterations',
            ['Shift', 'test', 'N_0', r'\sum H_0j N_j', 'alt'], None)
        self.assertEqual(start_it, 681)

    def test_different_it_key(self):
        """Rename iterations and check."""
        self._df_mock.rename(columns={'iterations': 'test_it'}, inplace=True)
        start_it = find_startit.find_starting_iteration_blocking(
            self._df_mock,  self._df_mock['test_it'].iloc[-1], 'test_it',
            ['Shift', 'test', 'N_0', r'\sum H_0j N_j', 'alt'], None)
        self.assertEqual(start_it, 491)

    def test_only_one_col(self):
        """Only pass in one column."""
        start_it = find_startit.find_starting_iteration_blocking(
            self._df_mock, self._df_mock['iterations'].iloc[-1], 'iterations',
            ['alt'], None)
        self.assertEqual(start_it, 461)

    def test_pass_hybrid_col(self):
        """Pass in `hybrid_col` (which should just be ignored)."""
        start_it = find_startit.find_starting_iteration_blocking(
            self._df_mock, self._df_mock['iterations'].iloc[-1], 'iterations',
            ['Shift', 'test', 'N_0', r'\sum H_0j N_j', 'alt'],
            ['Proj. Energy'])
        self.assertEqual(start_it, 491)

    def test_start_max_frac(self):
        """Test small, but permissible `start_max_frac`."""
        with self.assertRaisesRegex(RuntimeError, "Failed to find starting "
                                    "iteration. There might not be enough "
                                    "data."):
            _ = find_startit.find_starting_iteration_blocking(
                self._df_mock,  self._df_mock['iterations'].iloc[-1],
                'iterations',
                ['Shift', 'test', 'N_0', r'\sum H_0j N_j', 'alt'], None,
                start_max_frac=0.1)

    def test_too_small_start_max_frac(self):
        """Test small, but permissible `start_max_frac`."""
        with self.assertRaisesRegex(ValueError, "0.00001 < start_max_frac < 1 "
                                    "not satisfied!"):
            _ = find_startit.find_starting_iteration_blocking(
                self._df_mock,  self._df_mock['iterations'].iloc[-1],
                'iterations',
                ['Shift', 'test', 'N_0', r'\sum H_0j N_j', 'alt'], None,
                start_max_frac=0.0000000001)

    def test_too_large_start_max_frac(self):
        """Test small, but permissible `start_max_frac`."""
        with self.assertRaisesRegex(ValueError, "0.00001 < start_max_frac < 1 "
                                    "not satisfied!"):
            _ = find_startit.find_starting_iteration_blocking(
                self._df_mock,  self._df_mock['iterations'].iloc[-1],
                'iterations',
                ['Shift', 'test', 'N_0', r'\sum H_0j N_j', 'alt'], None,
                start_max_frac=2.0)

    def test_grid_size(self):
        """Test `grid_size` parameter."""
        start_it = find_startit.find_starting_iteration_blocking(
            self._df_mock,  self._df_mock['iterations'].iloc[-1], 'iterations',
            ['Shift', 'test', 'N_0', r'\sum H_0j N_j', 'alt'], None,
            grid_size=4)
        self.assertEqual(start_it, 576)

    def test_large_grid_size(self):
        """Test `grid_size` parameter, larger than size of data."""
        start_it = find_startit.find_starting_iteration_blocking(
            self._df_mock,  self._df_mock['iterations'].iloc[-1], 'iterations',
            ['Shift', 'test', 'N_0', r'\sum H_0j N_j', 'alt'], None,
            grid_size=len(self._df_mock)+2)
        self.assertEqual(start_it, 491)

    def test_number_of_reblocks_to_cut_off(self):
        """Test `number_of_reblocks_to_cut_off` parameter."""
        start_it = find_startit.find_starting_iteration_blocking(
            self._df_mock,  self._df_mock['iterations'].iloc[-1], 'iterations',
            ['Shift', 'test', 'N_0', r'\sum H_0j N_j', 'alt'], None,
            number_of_reblocks_to_cut_off=0)
        self.assertEqual(start_it, 451)

    def test_negative_number_of_reblocks_to_cut_off(self):
        """Test negative `number_of_reblocks_to_cut_off` parameter."""
        with self.assertRaisesRegex(ValueError,
                                    "'number_of_reblocks_to_cut_off' can't be "
                                    "negative!"):
            _ = find_startit.find_starting_iteration_blocking(
                self._df_mock,  self._df_mock['iterations'].iloc[-1],
                'iterations',
                ['Shift', 'test', 'N_0', r'\sum H_0j N_j', 'alt'], None,
                number_of_reblocks_to_cut_off=-10)

    def test_col_not_started_varying(self):
        """Pass a column in `cols` that has not starting varying."""
        not_varying_col = pd.DataFrame([3]*len(self._df_mock),
                                       columns=['notvarying'])
        self._df_mock = pd.concat([not_varying_col, self._df_mock], axis=1)
        with self.assertRaisesRegex(RuntimeError, "notvarying has not started "
                                    "varying in considered dataset."):
            _ = find_startit.find_starting_iteration_blocking(
                self._df_mock,  self._df_mock['iterations'].iloc[-1],
                'iterations',
                ['Shift', 'test', 'N_0', r'\sum H_0j N_j', 'notvarying'], None)

    def test_index_indifference(self):
        """Make sure starting iteration found is indep. of index."""
        self._df_mock.rename(
            index={i: i+9 for i in range(len(self._df_mock))}, inplace=True)
        start_it = find_startit.find_starting_iteration_blocking(
            self._df_mock,  self._df_mock['iterations'].iloc[-1], 'iterations',
            ['Shift', 'test', 'N_0', r'\sum H_0j N_j', 'alt'], None)
        self.assertEqual(start_it, 491)

    def test_unchanged_mutable_data(self):
        """Test that `data` passed in does not change."""
        df_mock_copy = self._df_mock.copy()
        _ = find_startit.find_starting_iteration_blocking(
            self._df_mock,  self._df_mock['iterations'].iloc[-1], 'iterations',
            ['Shift', 'test', 'N_0', r'\sum H_0j N_j', 'alt'], None)
        pd.testing.assert_frame_equal(df_mock_copy, self._df_mock,
                                      check_exact=True)

    def test_unchanged_mutable_cols(self):
        """Test that `cols` passed in does not change."""
        cols = ['Shift', 'test', 'N_0', r'\sum H_0j N_j', 'alt']
        cols_copy = ['Shift', 'test', 'N_0', r'\sum H_0j N_j', 'alt']
        _ = find_startit.find_starting_iteration_blocking(
            self._df_mock,  self._df_mock['iterations'].iloc[-1], 'iterations',
            cols, None)
        self.assertListEqual(cols, cols_copy)


class TestFindStartingIterationsMserMin(unittest.TestCase):
    """Test `find_starting_iteration_mser_min()`."""

    def setUp(self):
        """Setting up common mock df."""
        rng = np.random.default_rng(6183)
        cols = ['Shift', 'test', 'N_0', r'\sum H_0j N_j', 'alt',
                'Proj. Energy']
        means = [-0.1, 10.0, 11.0, -9.0, 20.0, -0.2]
        sine_periods = list(range(1, len(cols)+1))
        noise_facs = [0.1*mean for mean in means]
        self._df_mock = create_qmc_frame(rng, cols, means, sine_periods,
                                         noise_facs,
                                         frac_not_convergeds=[0.2]*len(cols),
                                         num_mc_its=300)
        iterations = pd.DataFrame(list(range(1, 1501, 5)),
                                  columns=['iterations'])
        self._df_mock = pd.concat([iterations, self._df_mock], axis=1)

    def test_default_mser_usage(self):
        """Test using defaults."""
        start_it = find_startit.find_starting_iteration_mser_min(
            self._df_mock,  self._df_mock['iterations'].iloc[-1], 'iterations',
            None, 'alt')
        self.assertEqual(start_it, 376)

    def test_specifying_smaller_end_its(self):
        """Specify smaller end iteration."""
        start_it = find_startit.find_starting_iteration_mser_min(
            self._df_mock,  700, 'iterations', None, 'alt')
        self.assertEqual(start_it, 376)

    def test_diff_hybrid_col(self):
        """Specify a different `hybrid_col`."""
        start_it = find_startit.find_starting_iteration_mser_min(
            self._df_mock,  self._df_mock['iterations'].iloc[-1], 'iterations',
            None, 'Shift')
        self.assertEqual(start_it, 451)

    def test_large_start_max_frac(self):
        """Specify a large `start_max_frac`."""
        start_it = find_startit.find_starting_iteration_mser_min(
            self._df_mock,  self._df_mock['iterations'].iloc[-1], 'iterations',
            None, 'alt', start_max_frac=1.0)
        self.assertEqual(start_it, 376)

    def test_short_start_max_frac(self):
        """Specify a short `start_max_frac`."""
        with self.assertWarns(UserWarning):
            start_it = find_startit.find_starting_iteration_mser_min(
                self._df_mock,  self._df_mock['iterations'].iloc[-1],
                'iterations', None, 'alt', start_max_frac=0.2)
            self.assertEqual(start_it, 296)

    def test_small_n_blocks(self):
        """Specify one `n_blocks`."""
        # [todo] - Should this be allowed?
        start_it = find_startit.find_starting_iteration_mser_min(
            self._df_mock,  self._df_mock['iterations'].iloc[-1],
            'iterations', None, 'alt', n_blocks=1)
        self.assertEqual(start_it, 1)

    def test_large_n_blocks(self):
        """Set more `n_blocks` than data points."""
        # Should this be decreased automatically?
        start_it = find_startit.find_starting_iteration_mser_min(
            self._df_mock,  self._df_mock['iterations'].iloc[-1],
            'iterations', None, 'alt', n_blocks=2000)
        self.assertEqual(start_it, 376)

    def test_specify_non_existing_cols(self):
        """Cols should be ignored anyways.  Test."""
        start_it = find_startit.find_starting_iteration_mser_min(
            self._df_mock,  self._df_mock['iterations'].iloc[-1], 'iterations',
            ['testcol'], 'alt')
        self.assertEqual(start_it, 376)

    def test_different_it_key(self):
        """Rename iterations and check."""
        self._df_mock.rename(columns={'iterations': 'test_it'}, inplace=True)
        start_it = find_startit.find_starting_iteration_mser_min(
            self._df_mock,  self._df_mock['test_it'].iloc[-1], 'test_it',
            None, 'alt')
        self.assertEqual(start_it, 376)

    def test_index_indifference(self):
        """Make sure starting iteration found is indep. of index."""
        self._df_mock.rename(
            index={i: i+9 for i in range(len(self._df_mock))}, inplace=True)
        start_it = find_startit.find_starting_iteration_mser_min(
            self._df_mock,  self._df_mock['iterations'].iloc[-1], 'iterations',
            None, 'alt')
        self.assertEqual(start_it, 376)

    def test_unchanged_mutable(self):
        """Test that `data` passed in does not change."""
        df_mock_copy = self._df_mock.copy()
        _ = find_startit.find_starting_iteration_mser_min(
            self._df_mock,  self._df_mock['iterations'].iloc[-1], 'iterations',
            None, 'alt')
        pd.testing.assert_frame_equal(df_mock_copy, self._df_mock,
                                      check_exact=True)

    def test_compare_to_old_version_in_lazy(self):
        """Compare to old version of this function in lazy.py."""
        # [todo] - delete when lazy.py gets deleted.
        import pyhande.lazy as lazy
        lazy_start_it = lazy.find_starting_iteration_mser_min(
            self._df_mock, {'qmc': {'ncycles': 5}}, start_max_frac=0.84)
        start_it = find_startit.find_starting_iteration_mser_min(
            self._df_mock,  self._df_mock['iterations'].iloc[-1], 'iterations',
            None, 'Proj. Energy')
        # The resuls are almost identical except for the fact that after
        # finding the index of the best start iteration, the new
        # implementation directly finds the corresponding iteration at
        # that index whereas the old implementation calculated
        # ncycles*start index which is shifted by the first element of
        # iterations.
        self.assertEqual(
            lazy_start_it+self._df_mock['iterations'].iloc[0], start_it)


class TestSelectFindStart(unittest.TestCase):
    """Test `select_find_start()`."""

    def setUp(self):
        """Setting up common mock df."""
        rng = np.random.default_rng(6183)
        cols = ['Shift', 'test', 'N_0', r'\sum H_0j N_j', 'alt',
                'Proj. Energy']
        means = [-0.1, 10.0, 11.0, -9.0, 20.0, -0.2]
        sine_periods = list(range(1, len(cols)+1))
        noise_facs = [0.1*mean for mean in means]
        self._df_mock = create_qmc_frame(rng, cols, means, sine_periods,
                                         noise_facs,
                                         frac_not_convergeds=[0.2]*len(cols),
                                         num_mc_its=300)
        iterations = pd.DataFrame(list(range(1, 1501, 5)),
                                  columns=['iterations'])
        self._df_mock = pd.concat([iterations, self._df_mock], axis=1)

    def test_select_blocking(self):
        """Test selecting the blocking find starting it function."""
        start_it = find_startit.select_find_start('blocking')(
            self._df_mock,  self._df_mock['iterations'].iloc[-1], 'iterations',
            ['Shift', 'test', 'N_0', r'\sum H_0j N_j', 'alt'], None)
        self.assertEqual(start_it, 491)

    def test_select_mser(self):
        """Test selecting the mser find starting it function."""
        start_it = find_startit.select_find_start('mser')(
            self._df_mock,  self._df_mock['iterations'].iloc[-1], 'iterations',
            None, 'alt')
        self.assertEqual(start_it, 376)

    def test_select_invalid(self):
        """Test selecting a non existent find starting it function."""
        with self.assertRaisesRegex(ValueError, "The find start iteration "
                                    "selected in 'start_its', 'testi', is not "
                                    "available!"):
            _ = find_startit.select_find_start('testi')(
                self._df_mock,  self._df_mock['iterations'].iloc[-1],
                'iterations', None, 'alt')

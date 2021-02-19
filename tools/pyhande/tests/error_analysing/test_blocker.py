"""Test error_analysing.blocker.Blocker."""
import unittest
import copy
import numpy as np
import pandas as pd
from pyhande.error_analysing.blocker import Blocker
from tests.create_mock_df import create_qmc_frame


class TestAccessPropertiesPreExe(unittest.TestCase):
    """Test error_analysing.blocker.Blocker properties."""

    def test_undefined_property_access(self):
        """Testing undefined property access (before .exe())."""
        # These should all raise AttributeError since .exe() has not
        # been executed yet.
        analyser = Blocker.inst_hande_ccmc_fciqmc()
        with self.assertRaises(AttributeError):
            analyser.start_its
        with self.assertRaises(AttributeError):
            analyser.end_its
        with self.assertRaises(AttributeError):
            analyser.opt_block
        with self.assertRaises(AttributeError):
            analyser.no_opt_block
        with self.assertRaises(AttributeError):
            analyser.reblock
        with self.assertRaises(AttributeError):
            analyser.data_len
        with self.assertRaises(AttributeError):
            analyser.covariance


class TestExe(unittest.TestCase):
    """Test .exe"""

    def setUp(self):
        """Setting up."""
        rng = np.random.default_rng(6183)
        cols = ['Shift', 'test', 'N_0', r'\sum H_0j N_j', 'alt', 'testi']
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
        self._observables = {'it_key': 'iterations',
                             'shift_key': 'Shift',
                             'sum_key': r'\sum H_0j N_j',
                             'ref_key': 'N_0', 'total_key': 'test',
                             'replica_key': 'replica id',
                             'proje_key': 'Proj. Energy',
                             'inst_proje_key': 'testi'}

    def test_defaults(self):
        """Testing default behaviour."""
        analyser = Blocker(
            'iterations', ['Shift', 'test', 'N_0', r'\sum H_0j N_j', 'alt'],
            'replica id',
            {'num': r'\sum H_0j N_j', 'denom': 'N_0', 'name': 'Proj. Energy'},
            None)
        analyser.exe([self._df_mock], self._observables)

        # start_its
        self.assertEqual(len(analyser.start_its), 1)
        self.assertEqual(analyser.start_its[0], 491)

        # end_its
        self.assertEqual(len(analyser.end_its), 1)
        self.assertEqual(analyser.end_its[0], 1496)

        # opt_block - test the mean
        self.assertEqual(len(analyser.opt_block), 1)
        pd.testing.assert_series_equal(
            analyser.opt_block[0]['mean'].astype('float'), pd.Series(
                [-0.10076365416711575, 9.909438497615628, 10.972245740153753,
                 -8.98536024489873, 20.150148053542424, -0.818917153123552],
                index=['Shift', 'test', 'N_0', r'\sum H_0j N_j', 'alt',
                       'Proj. Energy'], name='mean'))

        # no_opt_block
        self.assertEqual(len(analyser.no_opt_block), 1)
        self.assertListEqual(analyser.no_opt_block[0], [])

        # reblock - test 'alt' and the later added 'Proj. Energy'
        self.assertEqual(len(analyser.reblock), 1)
        pd.testing.assert_series_equal(
            analyser.reblock[0]['alt']['mean'].astype('float'), pd.Series(
                [20.124655430667065, 20.12465543066707, 20.150148053542424,
                 20.150148053542424, 20.181967635511462, 20.181967635511462,
                 20.181967635511462],
                index=pd.Index(list(range(7)), name='reblock'), name=('mean')))
        pd.testing.assert_series_equal(
            analyser.reblock[0]['Proj. Energy']['mean'].astype('float'),
            pd.Series(
                [-0.8181173916407062, -0.8181173916407061, -0.818917153123552,
                 -0.818917153123552, -0.8170703480189511, -0.8170703480189512,
                 -0.8170703480189512],
                index=pd.Index(list(range(7)), name='reblock'), name=('mean')))

        # covariance - test 'alt'
        self.assertEqual(len(analyser.covariance), 1)
        pd.testing.assert_series_equal(
            analyser.covariance[0].loc[2]['alt'], pd.Series(
                [-0.00016987733411814773, -0.13494641998291607,
                 -0.0073442015388262795, -0.13340899132333114,
                 2.571435360685306],
                index=pd.Index(
                    ['Shift', 'test', 'N_0', r'\sum H_0j N_j', 'alt'],
                    name=''), name='alt'))

    def test_start_its_number(self):
        """Testing passing an explicit start iteration."""
        analyser = Blocker(
            'iterations', ['Shift', 'test', 'N_0', r'\sum H_0j N_j', 'alt'],
            'replica id',
            {'num': r'\sum H_0j N_j', 'denom': 'N_0', 'name': 'Proj. Energy'},
            None, start_its=[491])
        analyser.exe([self._df_mock], self._observables)

        # start_its
        self.assertEqual(len(analyser.start_its), 1)
        self.assertEqual(analyser.start_its[0], 491)

        # end_its
        self.assertEqual(len(analyser.end_its), 1)
        self.assertEqual(analyser.end_its[0], 1496)

        # Only test opt_block to keep things shorter.
        # opt_block - test the mean
        self.assertEqual(len(analyser.opt_block), 1)
        pd.testing.assert_series_equal(
            analyser.opt_block[0]['mean'].astype('float'), pd.Series(
                [-0.10076365416711575, 9.909438497615628, 10.972245740153753,
                 -8.98536024489873, 20.150148053542424, -0.818917153123552],
                index=['Shift', 'test', 'N_0', r'\sum H_0j N_j', 'alt',
                       'Proj. Energy'], name='mean'))

    def test_start_its_mser(self):
        """Testing using the mser starting it finder."""
        # Have to specify `hybrid_col` which is used to find start it.
        analyser = Blocker(
            'iterations', ['Shift', 'test', 'N_0', r'\sum H_0j N_j', 'alt'],
            'replica id',
            {'num': r'\sum H_0j N_j', 'denom': 'N_0', 'name': 'Proj. Energy'},
            'testi', start_its='mser')
        analyser.exe([self._df_mock], self._observables)

        # start_its
        self.assertEqual(len(analyser.start_its), 1)
        self.assertEqual(analyser.start_its[0], 341)

        # end_its
        self.assertEqual(len(analyser.end_its), 1)
        self.assertEqual(analyser.end_its[0], 1496)

        # Only test opt_block to keep things shorter.
        # opt_block - test the mean
        self.assertEqual(len(analyser.opt_block), 1)
        pd.testing.assert_series_equal(
            analyser.opt_block[0]['mean'].astype('float'), pd.Series(
                [-0.10701527003950902, 9.887515934489185, 10.733592492673694,
                 -8.959426970366986, 20.021419902903823, -0.8341988713484589],
                index=['Shift', 'test', 'N_0', r'\sum H_0j N_j', 'alt',
                       'Proj. Energy'], name='mean'))

    def test_start_its_mser_no_hybrid_col(self):
        """Testing using the mser start it finder but no hybrid col."""
        # Have to specify `hybrid_col` which is used to find start it.
        with self.assertRaises(ValueError):
            _ = Blocker(
                'iterations', ['Shift', 'test', 'N_0', r'\sum H_0j N_j',
                               'alt'],
                'replica id', {'num': r'\sum H_0j N_j', 'denom': 'N_0',
                               'name': 'Proj. Energy'}, None, start_its='mser')

    def test_end_its_number(self):
        """Testing passing an explicit end iteration."""
        analyser = Blocker(
            'iterations', ['Shift', 'test', 'N_0', r'\sum H_0j N_j', 'alt'],
            'replica id',
            {'num': r'\sum H_0j N_j', 'denom': 'N_0', 'name': 'Proj. Energy'},
            None, end_its=[1300])
        analyser.exe([self._df_mock], self._observables)

        # start_its
        self.assertEqual(len(analyser.start_its), 1)
        self.assertEqual(analyser.start_its[0], 546)

        # end_its
        self.assertEqual(len(analyser.end_its), 1)
        self.assertEqual(analyser.end_its[0], 1300)

        # Only test opt_block to keep things shorter.
        # opt_block - test the mean
        self.assertEqual(len(analyser.opt_block), 1)
        pd.testing.assert_series_equal(
            analyser.opt_block[0]['mean'].astype('float'), pd.Series(
                [-0.10177814490742886, 9.880630482723749, 10.989156069488276,
                 -8.99157882776527, 20.166519924757836, -0.8182228708836579],
                index=['Shift', 'test', 'N_0', r'\sum H_0j N_j', 'alt',
                       'Proj. Energy'], name='mean'))

    def test_blocking_find_start_kw_args(self):
        """Pass an argument to blocking find starting it function."""
        analyser = Blocker(
            'iterations', ['Shift', 'test', 'N_0', r'\sum H_0j N_j', 'alt'],
            'replica id',
            {'num': r'\sum H_0j N_j', 'denom': 'N_0', 'name': 'Proj. Energy'},
            None, find_start_kw_args={'grid_size': 7})
        analyser.exe([self._df_mock], self._observables)

        # start_its
        self.assertEqual(len(analyser.start_its), 1)
        self.assertEqual(analyser.start_its[0], 526)

        # end_its
        self.assertEqual(len(analyser.end_its), 1)
        self.assertEqual(analyser.end_its[0], 1496)

        # Only test opt_block to keep things shorter.
        # opt_block - test the mean
        self.assertEqual(len(analyser.opt_block), 1)
        pd.testing.assert_series_equal(
            analyser.opt_block[0]['mean'].astype('float'), pd.Series(
                [-0.10074877158424976, 9.903349304555512, 10.984661142067944,
                 -8.986670126393298, 20.165413424228035, -0.8181108192748029],
                index=['Shift', 'test', 'N_0', r'\sum H_0j N_j', 'alt',
                       'Proj. Energy'], name='mean'))

    def test_mser_find_start_kw_args(self):
        """Pass an argument to mser find starting it function."""
        analyser = Blocker(
            'iterations', ['Shift', 'test', 'N_0', r'\sum H_0j N_j', 'alt'],
            'replica id',
            {'num': r'\sum H_0j N_j', 'denom': 'N_0', 'name': 'Proj. Energy'},
            'testi', start_its='mser', find_start_kw_args={'n_blocks': 51})
        analyser.exe([self._df_mock], self._observables)

        # start_its
        self.assertEqual(len(analyser.start_its), 1)
        self.assertEqual(analyser.start_its[0], 346)

        # end_its
        self.assertEqual(len(analyser.end_its), 1)
        self.assertEqual(analyser.end_its[0], 1496)

        # Only test opt_block to keep things shorter.
        # opt_block - test the mean
        self.assertEqual(len(analyser.opt_block), 1)
        pd.testing.assert_series_equal(
            analyser.opt_block[0]['mean'].astype('float'), pd.Series(
                [-0.10647583851219304, 9.892339064902817, 10.74772794441101,
                 -8.954469963178665, 20.06830627791626, -0.8331500396635119],
                index=['Shift', 'test', 'N_0', r'\sum H_0j N_j', 'alt',
                       'Proj. Energy'], name='mean'))

    def test_inst_hande_ccmc_fciqmc(self):
        """Test using class method to get Blocker instance."""
        analyser = Blocker.inst_hande_ccmc_fciqmc()
        analyser.exe([self._df_mock], self._observables)
        analyser2 = Blocker(
            'iterations', ['Shift', r'\sum H_0j N_j', 'N_0', 'test'],
            'replica id',
            {'num': r'\sum H_0j N_j', 'denom': 'N_0', 'name': 'Proj. Energy'},
            None)
        analyser2.exe([self._df_mock], self._observables)

        self.assertListEqual(analyser.start_its, analyser2.start_its)
        self.assertListEqual(analyser.end_its, analyser2.end_its)
        pd.testing.assert_frame_equal(analyser.opt_block[0],
                                      analyser2.opt_block[0])
        pd.testing.assert_frame_equal(analyser.reblock[0],
                                      analyser2.reblock[0])
        pd.testing.assert_frame_equal(analyser.covariance[0],
                                      analyser2.covariance[0])

    def test_multiple_calcs(self):
        """Test analysing two datasets."""
        analyser = Blocker(
            'iterations', ['Shift', r'\sum H_0j N_j', 'N_0', 'test'],
            'replica id',
            {'num': r'\sum H_0j N_j', 'denom': 'N_0', 'name': 'Proj. Energy'},
            None)
        analyser.exe([self._df_mock, self._df_mock], self._observables)
        self.assertEqual(len(analyser.start_its), 2)
        self.assertEqual(analyser.start_its[0], analyser.start_its[1])
        self.assertEqual(analyser.end_its[0], analyser.end_its[1])
        pd.testing.assert_frame_equal(analyser.opt_block[0],
                                      analyser.opt_block[1])
        self.assertListEqual(analyser.no_opt_block[0],
                             analyser.no_opt_block[1])
        pd.testing.assert_frame_equal(analyser.reblock[0],
                                      analyser.reblock[1])
        pd.testing.assert_frame_equal(analyser.covariance[0],
                                      analyser.covariance[1])

    def test_replica_tricks(self):
        """Test analysing dataset with replica tricks."""
        analyser = Blocker(
            'iterations', ['Shift', r'\sum H_0j N_j', 'N_0', 'test'],
            'replica id',
            {'num': r'\sum H_0j N_j', 'denom': 'N_0', 'name': 'Proj. Energy'},
            None)
        df_mock_replica = pd.concat([self._df_mock]*2, keys=[1, 2])
        df_mock_replica.reset_index(inplace=True)
        df_mock_replica.drop(columns=['level_1'], inplace=True)
        df_mock_replica.rename(columns={'level_0': 'replica id'}, inplace=True)
        analyser.exe([df_mock_replica], self._observables)

        self.assertListEqual(analyser.start_its, [491]*2)
        self.assertListEqual(analyser.end_its, [1496]*2)
        self.assertEqual(len(analyser.opt_block), 1)
        pd.testing.assert_series_equal(
            analyser.opt_block[0].loc[2]['mean'].astype('float'), pd.Series(
                [-0.10076365416711575, -8.98536024489873, 10.972245740153753,
                 9.909438497615628, -0.818917153123552],
                index=['Shift', r'\sum H_0j N_j', 'N_0', 'test',
                       'Proj. Energy'], name='mean'))
        pd.testing.assert_frame_equal(analyser.opt_block[0].loc[1],
                                      analyser.opt_block[0].loc[2])
        self.assertEqual(len(analyser.no_opt_block), 1)
        self.assertListEqual(analyser.no_opt_block[0], [[], []])
        self.assertEqual(len(analyser.reblock), 1)
        pd.testing.assert_frame_equal(analyser.reblock[0].loc[1],
                                      analyser.reblock[0].loc[2])
        pd.testing.assert_frame_equal(analyser.covariance[0].loc[1],
                                      analyser.covariance[0].loc[2])

    def test_unchanged_mutable(self):
        """Test that data and observables passed doesn't change."""
        analyser = Blocker(
            'iterations', ['Shift', 'test', 'N_0', r'\sum H_0j N_j', 'alt'],
            'replica id',
            {'num': r'\sum H_0j N_j', 'denom': 'N_0', 'name': 'Proj. Energy'},
            None)
        df_mock_copy = self._df_mock.copy()
        data_list = [self._df_mock]
        data_list_copy = [df_mock_copy]  # overkill with copying?
        observables_copy = copy.deepcopy(self._observables)
        analyser.exe(data_list, self._observables)
        self.assertDictEqual(observables_copy, self._observables)
        self.assertEqual(len(data_list), len(data_list_copy))
        for dat, dat_copy in zip(data_list, data_list_copy):
            pd.testing.assert_frame_equal(dat, dat_copy, check_exact=True)


class TestStaticMethods(unittest.TestCase):
    """Test Blocker's static methods."""

    def test_check_input_basic(self):
        """Test _check_input with basic input."""
        Blocker._check_input('its', ['a', 'b'], None, 'blocking')

    def test_check_input_basic_mser(self):
        """Test _check_input with basic input but use 'mser'."""
        Blocker._check_input('its', ['a', 'b'], 'd', 'mser')

    def test_check_input_not_it_key(self):
        """Test _check_input but don't specify `it_key`."""
        with self.assertRaisesRegex(ValueError, "'it_key' must be specified!"):
            Blocker._check_input(None, ['a', 'b'], None, 'blocking')

    def test_check_input_not_cols(self):
        """Test _check_input but don't specify `cols`."""
        with self.assertRaisesRegex(ValueError, "'cols' has to be specified."):
            Blocker._check_input('it', None, None, None)

    def test_check_input_not_hybrid_col_mser(self):
        """Don't specify `hybrid_col` with 'mser' start_its."""
        with self.assertRaises(ValueError):
            # [todo] add error message
            Blocker._check_input('it', ['a', 'b'], None, 'mser')

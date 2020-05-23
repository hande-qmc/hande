"""Test results_viewer.results.Results."""
import unittest
from unittest.mock import MagicMock, PropertyMock
import pandas as pd
import numpy as np
from tests.create_mock_df import create_qmc_frame
from pyhande.results_viewer.results_ccmc_fciqmc import ResultsCcmcFciqmc


class TestSummary(unittest.TestCase):
    """Test summary."""

    # See documentation https://docs.python.org/3/library/unittest.mock.html
    # for mocking properties.
    def test_add_single_metadata_item_only_extractor(self):
        """Test adding a single piece of metadata to summary."""
        extractor = MagicMock()
        metadata = PropertyMock(
            return_value=[[{'a': {'b': 90}}]])
        type(extractor).metadata = metadata
        results = ResultsCcmcFciqmc(extractor)
        # First summary is empty.
        pd.testing.assert_frame_equal(results.summary, pd.DataFrame())
        # Now add metadata.
        results.add_metadata(['a:b'])
        pd.testing.assert_frame_equal(
            results.summary, pd.DataFrame(
                {'calc id': [0], 'observable': ['b'], 'value/mean': [90]}))

    def test_add_metadata_replica_tricks(self):
        """Test adding metadata to summary with replica tricks."""
        extractor = MagicMock()
        metadata = PropertyMock(return_value=[[{'a': {'b': 7}}]])
        type(extractor).metadata = metadata
        preparator = MagicMock()
        observables = PropertyMock(return_value={'replica_key': 'replica id'})
        data = PropertyMock(return_value=[
            pd.DataFrame({'a': [1, 2, 3], 'replica id': [4, 5, 6]})])
        type(preparator).observables = observables
        type(preparator).data = data
        results = ResultsCcmcFciqmc(extractor, preparator=preparator)

        # Preset summary.
        results.summary = pd.DataFrame(
            {'calc id': [0, 0], 'replica id': [1, 2], 'observable': ['x', 'x'],
             'value/mean': [90, 93], 'standard error': [0.1, 0.7]})
        # Now add metadata.
        results.add_metadata(['a:b'])
        pd.testing.assert_frame_equal(
            results.summary, pd.DataFrame(
                {'calc id': [0, 0, 0, 0], 'replica id': [1, 1, 2, 2],
                 'observable': ['x', 'b', 'x', 'b'],
                 'value/mean': [90, 7, 93, 7],
                 'standard error': [0.1, None, 0.7, None]}))

    def test_init_summary_opt_block(self):
        """Init summary with opt_block."""
        extractor = MagicMock()
        analyser = MagicMock()
        opt_block = PropertyMock(return_value=[
            pd.DataFrame([[10.0, 0.1, 0.01], [2512.2, 384.8, 15.0]],
                         columns=['mean', 'standard error',
                                  'standard error error'],
                         index=['oba', 'obb'])])
        type(analyser).opt_block = opt_block
        results = ResultsCcmcFciqmc(extractor, analyser=analyser)
        pd.testing.assert_frame_equal(
            results.summary, pd.DataFrame(
                {'calc id': [0, 0],
                 'observable': ['oba', 'obb'],
                 'value/mean': [10.0, 2512.2],
                 'standard error': [0.1, 384.8],
                 'standard error error': [0.01, 15.0]}))

    def test_init_summary_opt_block_multi_calc(self):
        """Init summary with opt_block, two calcs."""
        extractor = MagicMock()
        analyser = MagicMock()
        opt_block = PropertyMock(return_value=[
            pd.DataFrame([[10.0, 0.15, 0.01], [2512.2, 385.8, 16.0]],
                         columns=['mean', 'standard error',
                                  'standard error error'],
                         index=['oba', 'obb']),
            pd.DataFrame([[11.0, 0.1, 0.01], [2513.2, 384.8, 15.0]],
                         columns=['mean', 'standard error',
                                  'standard error error'],
                         index=['oba', 'obb'])])
        type(analyser).opt_block = opt_block
        results = ResultsCcmcFciqmc(extractor, analyser=analyser)
        pd.testing.assert_frame_equal(
            results.summary, pd.DataFrame(
                {'calc id': [0, 0, 1, 1],
                 'observable': ['oba', 'obb']*2,
                 'value/mean': [10.0, 2512.2, 11.0, 2513.2],
                 'standard error': [0.15, 385.8, 0.1, 384.8],
                 'standard error error': [0.01, 16.0, 0.01, 15.0]}))

    def test_init_summary_opt_block_replica(self):
        """Init summary with opt_block, replica tricks."""
        extractor = MagicMock()
        analyser = MagicMock()
        opt_block = PropertyMock(return_value=[pd.concat([
            pd.DataFrame([[10.0, 0.15, 0.01], [2512.2, 385.8, 16.0]],
                         columns=['mean', 'standard error',
                                  'standard error error'],
                         index=['oba', 'obb']),
            pd.DataFrame([[11.0, 0.1, 0.01], [2513.2, 384.8, 15.0]],
                         columns=['mean', 'standard error',
                                  'standard error error'],
                         index=['oba', 'obb'])],
            keys=[1, 2], names=['replica id'])])
        type(analyser).opt_block = opt_block
        results = ResultsCcmcFciqmc(extractor, analyser=analyser)
        pd.testing.assert_frame_equal(
            results.summary, pd.DataFrame(
                {'calc id': [0, 0, 0, 0],
                 'replica id': [1, 1, 2, 2],
                 'observable': ['oba', 'obb']*2,
                 'value/mean': [10.0, 2512.2, 11.0, 2513.2],
                 'standard error': [0.15, 385.8, 0.1, 384.8],
                 'standard error error': [0.01, 16.0, 0.01, 15.0]}))

    def test_summary_pretty(self):
        """Show prettified summary with two calcs."""
        extractor = MagicMock()
        preparator = MagicMock()
        observables = PropertyMock(return_value={'replica_key': 'replica id'})
        # not quite realisitc
        data = PropertyMock(return_value=2*[
            pd.DataFrame({'oba': [1, 2, 3], 'obb': [4, 5, 6]})])
        type(preparator).observables = observables
        type(preparator).data = data
        analyser = MagicMock()
        opt_block = PropertyMock(return_value=[
            pd.DataFrame([[10.0, 0.15, 0.01], [2512.2, 385.8, 16.0]],
                         columns=['mean', 'standard error',
                                  'standard error error'],
                         index=['oba', 'obb']),
            pd.DataFrame([[11.0, 0.1, 0.01], [2513.2, 384.8, 15.0]],
                         columns=['mean', 'standard error',
                                  'standard error error'],
                         index=['oba', 'obb'])])
        type(analyser).opt_block = opt_block
        results = ResultsCcmcFciqmc(extractor, preparator=preparator,
                                    analyser=analyser)
        pd.testing.assert_frame_equal(
            results.summary_pretty, pd.DataFrame(
                [['10.0(2)', '2500(400)'], ['11.0(1)', '2500(400)']],
                columns=pd.Index(['oba', 'obb'], name='observable'),
                index=pd.Index(list(range(2)), name='calc id')))

    def test_summary_pretty_replica(self):
        """Show prettified summary with replica tricks."""
        extractor = MagicMock()
        preparator = MagicMock()
        observables = PropertyMock(return_value={'replica_key': 'replica id'})
        # not quite realisitc
        data = PropertyMock(return_value=[
            pd.DataFrame({'oba': [1, 2], 'obb': [4, 5],
                          'replica id': [1, 2]})])
        type(preparator).observables = observables
        type(preparator).data = data
        analyser = MagicMock()
        opt_block = PropertyMock(return_value=[pd.concat([
            pd.DataFrame([[10.0, 0.15, 0.01], [2512.2, 385.8, 16.0]],
                         columns=['mean', 'standard error',
                                  'standard error error'],
                         index=['oba', 'obb']),
            pd.DataFrame([[11.0, 0.1, 0.01], [2513.2, 384.8, 15.0]],
                         columns=['mean', 'standard error',
                                  'standard error error'],
                         index=['oba', 'obb'])],
            keys=[1, 2], names=['replica id'])])
        type(analyser).opt_block = opt_block
        results = ResultsCcmcFciqmc(extractor, preparator=preparator,
                                    analyser=analyser)
        pd.testing.assert_frame_equal(
            results.summary_pretty, pd.DataFrame(
                [['10.0(2)', '2500(400)'], ['11.0(1)', '2500(400)']],
                columns=pd.Index(['oba', 'obb'], name='observable'),
                index=pd.Index([(0, 1), (0, 2)],
                               name=('calc id', 'replica id'))))

    def test_compare_obs(self):
        """Compare observables."""
        extractor = MagicMock()
        metadata = PropertyMock(
            return_value=[[{'a': {'b': 7}}], [{'a': {'b': 6}}]])
        type(extractor).metadata = metadata
        preparator = MagicMock()
        observables = PropertyMock(return_value={'replica_key': 'replica id'})
        # not quite realisitc
        data = PropertyMock(return_value=[
            pd.DataFrame({'oba': [1, 2], 'obb': [4, 5]})])
        type(preparator).observables = observables
        type(preparator).data = data
        analyser = MagicMock()
        opt_block = PropertyMock(return_value=[
            pd.DataFrame([[10.0, 0.15, 0.01], [2512.2, 385.8, 16.0]],
                         columns=['mean', 'standard error',
                                  'standard error error'],
                         index=['oba', 'obb']),
            pd.DataFrame([[11.0, 0.1, 0.01], [2513.2, 384.8, 15.0]],
                         columns=['mean', 'standard error',
                                  'standard error error'],
                         index=['oba', 'obb'])])
        type(analyser).opt_block = opt_block
        results = ResultsCcmcFciqmc(extractor, preparator=preparator,
                                    analyser=analyser)
        results.add_metadata('a:b')
        pd.testing.assert_frame_equal(
            results.compare_obs(['b', 'oba']), pd.DataFrame(
                [[7.0, 10.0, None, 0.15], [6.0, 11.0, None, 0.10]],
                columns=pd.Index([('value/mean', 'b'), ('value/mean', 'oba'),
                                  ('standard error', 'b'),
                                  ('standard error', 'oba')],
                                 name=(None, 'observable')),
                index=pd.Index([0, 1], name='calc id')
            ).astype('float'))

    def test_compare_obs_replica(self):
        """Compare observables with replica tricks."""
        extractor = MagicMock()
        metadata = PropertyMock(return_value=[[{'a': {'b': 7}}]])
        type(extractor).metadata = metadata
        preparator = MagicMock()
        observables = PropertyMock(return_value={'replica_key': 'replica id'})
        # not quite realisitc
        data = PropertyMock(return_value=[
            pd.DataFrame({'oba': [1, 2], 'obb': [4, 5],
                          'replica id': [1, 2]})])
        type(preparator).observables = observables
        type(preparator).data = data
        analyser = MagicMock()
        opt_block = PropertyMock(return_value=[pd.concat([
            pd.DataFrame([[10.0, 0.15, 0.01], [2512.2, 385.8, 16.0]],
                         columns=['mean', 'standard error',
                                  'standard error error'],
                         index=['oba', 'obb']),
            pd.DataFrame([[11.0, 0.1, 0.01], [2513.2, 384.8, 15.0]],
                         columns=['mean', 'standard error',
                                  'standard error error'],
                         index=['oba', 'obb'])],
            keys=[1, 2], names=['replica id'])])
        type(analyser).opt_block = opt_block
        results = ResultsCcmcFciqmc(extractor, preparator=preparator,
                                    analyser=analyser)
        results.add_metadata('a:b')
        pd.testing.assert_frame_equal(
            results.compare_obs(['b', 'oba']), pd.DataFrame(
                [[7.0, 10.0, None, 0.15], [7.0, 11.0, None, 0.10]],
                columns=pd.Index([('value/mean', 'b'), ('value/mean', 'oba'),
                                  ('standard error', 'b'),
                                  ('standard error', 'oba')],
                                 name=(None, 'observable')),
                index=pd.Index([(0, 1), (0, 2)],
                               name=('calc id', 'replica id'))
            ).astype('float'))


class TestShoulder(unittest.TestCase):
    """Calculate shoulder."""

    def setUp(self):
        # from test_analysis.py
        rng = np.random.default_rng(94)
        cols = ['N_0', '# H psips', 'Shift', 'alt1', 'alt2']
        means = [1000.0, 20000.0, -0.7, 10000.0, 20000.0]
        sine_periods = [5.0, 2.9, 11.0, 10.1, 3.12]
        noise_facs = [0.1, 0.5, 0.001, 0.1, 0.5]
        self._data = create_qmc_frame(
            rng, cols, means, sine_periods, noise_facs,
            frac_not_convergeds=[0.6 for _ in range(5)], num_mc_its=70
        )
        self._data.replace(
            self._data['Shift'].loc[:50].values, 0.0, inplace=True
        )

    def test_shoulder(self):
        """Calculate and access shoulder."""
        extractor = MagicMock()
        preparator = MagicMock()
        observables = PropertyMock(return_value={'total_key': '# H psips',
                                                 'ref_key': 'N_0',
                                                 'shift_key': 'Shift',
                                                 'replica_key': 'replica id'})
        data = PropertyMock(return_value=[self._data])
        type(preparator).observables = observables
        type(preparator).data = data
        results = ResultsCcmcFciqmc(extractor, preparator=preparator)
        pd.testing.assert_frame_equal(
            results.shoulder, pd.DataFrame(
                [[0, 'shoulder estimator', 19.2785454706, 0.9283478679488482],
                 [0, 'shoulder height', 13253.874811, 2095.242403]],
                columns=pd.Index(['calc id', 'observable', 'value/mean',
                                  'standard error'], name=None)))
        results.add_shoulder()
        pd.testing.assert_frame_equal(
            results.summary, pd.DataFrame(
                {'calc id': [0, 0],
                 'observable': ['shoulder estimator', 'shoulder height'],
                 'value/mean': [19.278545470627485, 13253.874811],
                 'standard error': [0.9283478679488482, 2095.242403]}))

    def test_shoulder_replica(self):
        """Calculate and access shoulder."""
        extractor = MagicMock()
        preparator = MagicMock()
        observables = PropertyMock(return_value={'total_key': '# H psips',
                                                 'ref_key': 'alt1',
                                                 'shift_key': 'Shift',
                                                 'replica_key': 'replica id'})
        data_copy = self._data.copy()*0.9
        data = pd.concat([self._data, data_copy], ignore_index=True)
        data = pd.concat([pd.DataFrame([[1]]*70 + [[2]]*70,
                                       columns=['replica id']), data], axis=1)
        data = PropertyMock(return_value=[data])
        type(preparator).observables = observables
        type(preparator).data = data
        results = ResultsCcmcFciqmc(extractor, preparator=preparator)
        pd.testing.assert_frame_equal(
            results.shoulder, pd.DataFrame(
                [[0, 1, 'shoulder estimator', 2.264741, 0.040887],
                 [0, 1, 'shoulder height', 13253.874811, 2095.242403],
                 [0, 2, 'shoulder estimator', 2.264741, 0.040887],
                 [0, 2, 'shoulder height', 11928.487330176968, 1885.71816]],
                columns=pd.Index(['calc id', 'replica id', 'observable',
                                  'value/mean', 'standard error'], name=None)))
        results.add_shoulder()
        pd.testing.assert_frame_equal(
            results.summary, pd.DataFrame(
                {'calc id': [0, 0, 0, 0],
                 'replica id': [1, 1, 2, 2],
                 'observable': ['shoulder estimator', 'shoulder height']*2,
                 'value/mean': [2.264741, 13253.874811, 2.264741, 11928.48733],
                 'standard error': [0.040887, 2095.242403, 0.040887,
                                    1885.71816]}))


class TestInefficiency(unittest.TestCase):
    """Test calculating and accessing the inefficiency."""

    def test_inefficiency(self):
        """Test calculating and accessing inefficiency."""
        extractor = MagicMock()
        metadata = [[{'qmc': {'tau': 0.4}}]]
        type(extractor).metadata = metadata
        preparator = MagicMock()
        observables = PropertyMock(return_value={'total_key': 'obb',
                                                 'ref_key': 'oba',
                                                 'sum_key': 'sum',
                                                 'proje_key': 'proje',
                                                 'replica_key': 'replica id'})
        # not quite realisitc
        data = PropertyMock(return_value=[
            pd.DataFrame({'oba': [1, 54], 'obb': [4, 23], 'sum': [34, 43],
                          'proje': [1, 12]})])
        type(preparator).observables = observables
        type(preparator).data = data
        analyser = MagicMock()
        opt_block = PropertyMock(return_value=[
            pd.DataFrame([[10.0, 0.1, 0.01], [2512.2, 384.8, 15.0],
                          [-35.3, 1.3, 0.2], [-3.48, 0.781, 0.23]],
                         columns=['mean', 'standard error',
                                  'standard error error'],
                         index=['oba', 'obb', 'sum', 'proje'])])
        start_its = [98]
        end_its = [293]
        type(analyser).opt_block = opt_block
        type(analyser).start_its = start_its
        type(analyser).end_its = end_its
        results = ResultsCcmcFciqmc(extractor, preparator=preparator,
                                    analyser=analyser)
        pd.testing.assert_frame_equal(
            results.inefficiency, pd.DataFrame(
                [[0, 'Inefficiency', 345.720746, 68.740286]],
                columns=pd.Index(['calc id', 'observable', 'value/mean',
                                  'standard error'], name=None)))
        results.add_inefficiency()
        pd.testing.assert_frame_equal(
            results.summary, pd.DataFrame(
                {'calc id': [0, 0, 0, 0, 0],
                 'observable': ['oba', 'obb', 'sum', 'proje', 'Inefficiency'],
                 'value/mean': [10.0, 2512.2, -35.3, -3.48, 345.720746],
                 'standard error': [0.1, 384.8, 1.3, 0.781, 68.74028633579539],
                 'standard error error': [0.01, 15.0, 0.2, 0.23, None]}))

    def test_inefficiency_replica(self):
        """Test calculating and accessing inefficiency with replicas."""
        extractor = MagicMock()
        metadata = [[{'qmc': {'tau': 0.4}}]]
        type(extractor).metadata = metadata
        preparator = MagicMock()
        observables = PropertyMock(return_value={'total_key': 'obb',
                                                 'ref_key': 'oba',
                                                 'sum_key': 'sum',
                                                 'proje_key': 'proje',
                                                 'replica_key': 'replica id'})
        # not quite realisitc
        data = PropertyMock(return_value=[
            pd.DataFrame({'oba': [1, 54], 'obb': [4, 23], 'sum': [34, 43],
                          'proje': [1, 12], 'replica id': [1, 2]})])
        type(preparator).observables = observables
        type(preparator).data = data
        analyser = MagicMock()
        opt_block = PropertyMock(return_value=[pd.concat([
            pd.DataFrame([[10.0, 0.1, 0.01], [2512.2, 384.8, 15.0],
                          [-35.3, 1.3, 0.2], [-3.48, 0.781, 0.23]],
                         columns=['mean', 'standard error',
                                  'standard error error'],
                         index=['oba', 'obb', 'sum', 'proje']),
            pd.DataFrame([[11.0, 0.1, 0.01], [2513.2, 384.8, 15.0],
                          [-89.3, 1.5, 0.2], [-2.48, 0.781, 0.23]],
                         columns=['mean', 'standard error',
                                  'standard error error'],
                         index=['oba', 'obb', 'sum', 'proje'])],
            keys=[1, 2], names=['replica id'])])
        start_its = [98]
        end_its = [293]
        type(analyser).opt_block = opt_block
        type(analyser).start_its = start_its
        type(analyser).end_its = end_its
        results = ResultsCcmcFciqmc(extractor, preparator=preparator,
                                    analyser=analyser)
        pd.testing.assert_frame_equal(
            results.inefficiency, pd.DataFrame(
                [[0, 1, 'Inefficiency', 345.720746, 68.740286],
                 [0, 2, 'Inefficiency', 345.78954704501984, 63.42063674025]],
                columns=pd.Index(['calc id', 'replica id', 'observable',
                                  'value/mean', 'standard error'], name=None)))
        results.add_inefficiency()
        pd.testing.assert_frame_equal(results.summary, pd.DataFrame(
            {'calc id': [0, 0, 0, 0, 0]*2,
             'replica id': [1]*5 + [2]*5,
             'observable': ['oba', 'obb', 'sum', 'proje', 'Inefficiency']*2,
             'value/mean': [10.0, 2512.2, -35.3, -3.48, 345.720746]
             + [11.0, 2513.2, -89.3, -2.48, 345.78954704501984],
             'standard error': [0.1, 384.8, 1.3, 0.781, 68.74028633579539]
             + [0.1, 384.8, 1.5, 0.781, 63.42063674024586],
             'standard error error': [0.01, 15.0, 0.2, 0.23, None]*2}))

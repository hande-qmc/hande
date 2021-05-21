"""Test results_viewer.results.Results."""
import unittest
from unittest.mock import MagicMock, PropertyMock
import pandas as pd
from pyhande.results_viewer.results import Results


class TestAddMetadata(unittest.TestCase):
    """Test adding metadata to summary table."""

    # See documentation https://docs.python.org/3/library/unittest.mock.html
    # for mocking properties.
    def test_add_single_metadata_item(self):
        """Test adding a single piece of metadata to summary."""
        extractor = MagicMock()
        metadata = PropertyMock(
            return_value=[[{'a': {'b': 90}}]])
        type(extractor).metadata = metadata
        results = Results(extractor)
        # First summary is empty.
        pd.testing.assert_frame_equal(results.summary, pd.DataFrame())
        # Now add metadata.
        results.add_metadata(['a:b'])
        pd.testing.assert_frame_equal(
            results.summary, pd.DataFrame(
                {'calc id': [0], 'observable': ['b'], 'value/mean': [90]}))

    def test_add_tuple_metadata(self):
        """Calculation was merged from three calcs, differing meta."""
        extractor = MagicMock()
        metadata = PropertyMock(
            return_value=[[{'a': {'b': 90}}, {'a': {'b': 90}},
                           {'a': {'b': 8}}]])
        type(extractor).metadata = metadata
        results = Results(extractor)

        # Now add metadata.
        results.add_metadata(['a:b'])
        pd.testing.assert_frame_equal(
            results.summary, pd.DataFrame(
                {'calc id': [0], 'observable': ['b'],
                 'value/mean': [(90, 90, 8)]}))

    def test_add_same_metadata_multi_calc(self):
        """Calculation was merged from three calcs, same meta."""
        extractor = MagicMock()
        metadata = PropertyMock(
            return_value=[[{'a': {'b': 90}}, {'a': {'b': 90}},
                           {'a': {'b': 90}}]])
        type(extractor).metadata = metadata
        results = Results(extractor)

        # Now add metadata.
        results.add_metadata(['a:b'])
        pd.testing.assert_frame_equal(
            results.summary, pd.DataFrame(
                {'calc id': [0], 'observable': ['b'], 'value/mean': [90]}))

    def test_add_metadata_multi_calc(self):
        """Calculation was merged from three calcs, same meta."""
        extractor = MagicMock()
        metadata = PropertyMock(
            return_value=[[{'a': {'b': 90}}, {'a': {'b': 80}}],
                          [{'a': {'b': 90}}]])
        type(extractor).metadata = metadata
        results = Results(extractor)

        # Now add metadata.
        results.add_metadata(['a:b'])
        pd.testing.assert_frame_equal(
            results.summary, pd.DataFrame(
                {'calc id': [0, 1], 'observable': ['b', 'b'],
                 'value/mean': [(90, 80), 90]}))

    def test_add_further_metadata(self):
        """After adding some metadata, add more."""
        extractor = MagicMock()
        metadata = PropertyMock(
            return_value=[[{'a': {'b': 90, 'c': 'y'}},
                           {'a': {'b': 80, 'c': 'y'}}],
                          [{'a': {'b': 90}}]])
        type(extractor).metadata = metadata
        results = Results(extractor)
        results.add_metadata(['a:b'])
        with self.assertWarnsRegex(UserWarning, "Metadata for calc #1 has no "
                                   "key a:c"):
            results.add_metadata(['a:c'])
            pd.testing.assert_frame_equal(
                results.summary, pd.DataFrame(
                    {'calc id': [0, 0, 1, 1],
                     'observable': ['b', 'c', 'b', 'c'],
                     'value/mean': [(90, 80), 'y', 90, None]}))


class TestSettingSummary(unittest.TestCase):
    """Setting .summary directly."""

    def test_set_to_df(self):
        """Set summary to a dataframe."""
        extractor = MagicMock()
        results = Results(extractor)
        results.summary = pd.DataFrame({'test': [0], 'testi': [1]})
        pd.testing.assert_frame_equal(
            results.summary, pd.DataFrame({'test': [0], 'testi': [1]}))

    def test_set_to_non_df(self):
        """Attempt setting summary to something other than dataframe."""
        extractor = MagicMock()
        results = Results(extractor)
        with self.assertRaisesRegex(TypeError, "Cannot set summary. It has to "
                                    "be a pd.DataFrame."):
            results.summary = 354


class TestStaticMethods(unittest.TestCase):
    """Test Results's static method."""

    def test_access_meta_and_check_tuple(self):
        """Test Results._access_meta_and_check() two merged calcs."""
        metaitem = Results._access_meta_and_check(
            0, [{'a': {'b': 5}}, {'a': {'b': 7}}], 'a:b')
        self.assertTupleEqual(metaitem, (5, 7))

    def test_access_meta_and_check_single(self):
        """Test Results._access_meta_and_check() one calc."""
        metaitem = Results._access_meta_and_check(
            0, [{'a': {'b': 5}}], 'a:b')
        self.assertEqual(metaitem, 5)

    def test_access_meta_and_check_non_existent(self):
        """Test Results._access_meta_and_check(), no `meta_key`."""
        # One of the two calcs (which were later merged to one calc)
        # does not have that metadata item.
        with self.assertWarnsRegex(UserWarning,
                                   "Metadata for calc #0 has no key a:b."):
            metaitem = Results._access_meta_and_check(
                0, [{'a': {'b': 5}}, {'a': {'c': 7}}], 'a:b')
        self.assertTupleEqual(metaitem, (5, None))

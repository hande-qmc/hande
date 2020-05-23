"""Test helpers.simple_callables."""
import unittest
import copy
import pyhande.helpers.simple_callables as simple_callables


class TestDoNothing(unittest.TestCase):
    """Test simple_callables.do_nothing()."""

    def test_pass_list(self):
        """Pass mutable object and check that it is unaffected."""
        test_list = [1, 2, 3, 4, 5, 6]
        test_list_copy = copy.copy(test_list)
        simple_callables.do_nothing(test_list, test_list)
        self.assertListEqual(test_list, test_list_copy)

    def test_pass_nothing(self):
        """Pass nothing."""
        simple_callables.do_nothing()


class TestRaiseValueError(unittest.TestCase):
    """Test simple_callables.RaseValueError."""

    def test_pass_list(self):
        """Pass list and raise value error."""
        valerr = simple_callables.RaiseValueError("test")
        test_list = [1, 2, 3, 4, 5, 6]
        test_list_copy = copy.copy(test_list)
        with self.assertRaises(ValueError):
            valerr(test_list, test_list)
        self.assertListEqual(test_list, test_list_copy)

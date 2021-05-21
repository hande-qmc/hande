"""Test testcode.py."""
import unittest
import pandas as pd
import pyhande.testcode as testcode


class ExtractTestData(unittest.TestCase):
    """Test testcode.extract_test_data().
    Note that most of the tests are already done by the next class.
    """

    def test_basic_input(self):
        """Test basic input."""
        data = testcode.extract_test_data("tests/hande_files/ueg.out")
        data_exp_df = pd.DataFrame([
            [10, 0.0, -0.378518, 2.70000000e+00, 1.42000000e+02, 79, 50,
             0.4879, 0.0004],
            [20, -1.10017479e-01, -8.20941670e-01, 3.00000000e+00,
             1.28200000e+03, 782, 463, 7.62000000e-01, 8.00000000e-04]
        ], columns=['iterations', 'Shift', r'\sum_H_0j_N_j', 'N_0',
                    '#_H_psips', '#_states', '#_spawn_events', 'R_spawn',
                    'time'])
        (key, val) = list(data.items())[0]
        self.assertEqual(key, 'FCIQMC [0]')
        pd.testing.assert_frame_equal(val, data_exp_df, check_exact=False)


class TestTestcodeData(unittest.TestCase):
    """Test testcode.testcode_data()."""

    def test_basic_input(self):
        """Test basic input."""
        data = testcode.testcode_data("tests/hande_files/ueg.out")
        exp_data = {
            'iterations': [10, 20], 'Shift': [0.0, -0.11001747899],
            '\\sum H_0j N_j': [-0.37851801553, -0.82094167005],
            'N_0': [2.7, 3.0], '# H psips': [142.0, 1282.0],
            '# states': [79, 782], '# spawn_events': [50, 463],
            'R_spawn': [0.4879, 0.762], 'time': [0.0004, 0.0008]
        }
        pd.testing.assert_frame_equal(pd.DataFrame.from_dict(data),
                                      pd.DataFrame.from_dict(exp_data))

    def test_shortened_input(self):
        """Two qmc data lines have been artifically removed. Test."""
        data = testcode.testcode_data("tests/hande_files/shorten_ueg.out")
        exp_data = {
            'iterations': [], 'Shift': [],
            '\\sum H_0j N_j': [],
            'N_0': [], '# H psips': [],
            '# states': [], '# spawn_events': [],
            'R_spawn': [], 'time': []
        }
        pd.testing.assert_frame_equal(pd.DataFrame.from_dict(data),
                                      pd.DataFrame.from_dict(exp_data))

    def test_all_shortened_input(self):
        """All QMC table data was artifically removed. Test."""
        data = testcode.testcode_data(
            "tests/hande_files/all_qmc_shorten_ueg.out")
        exp_data = {}
        pd.testing.assert_frame_equal(pd.DataFrame.from_dict(data),
                                      pd.DataFrame.from_dict(exp_data))

    def test_hilbert(self):
        """Test MC Hilbert space size estimation extraction."""
        data = testcode.testcode_data("tests/hande_files/hilbert_ueg.out")
        exp_data = {
            'iterations': [2, 3], 'space size': [20311600.0, 10155800.0],
            'mean': [17772650.0, 15233700.0],
            'std. err.': [2538951.0, 2931728.0]
        }
        pd.testing.assert_frame_equal(pd.DataFrame.from_dict(data),
                                      pd.DataFrame.from_dict(exp_data))

    def test_basic_fci_input(self):
        """Test FCI."""
        data = testcode.testcode_data("tests/hande_files/fci_ueg.out")
        exp_data = {
            'eigv 1': [-0.017888297593], 'eigv 2': [9.451775889133],
            'eigv 3': [9.451775889133], 'eigv 4': [9.525116432464],
            'eigv 5': [9.525116432464]
        }
        pd.testing.assert_frame_equal(pd.DataFrame.from_dict(data),
                                      pd.DataFrame.from_dict(exp_data))

    def test_multiple_fciqmc(self):
        """Have two calculations in this output."""
        data = testcode.testcode_data("tests/hande_files/multi_ueg.out")
        exp_data = {
            'iterations': [10, 20, 10, 20],
            'Shift': [0.0, -0.11001747899, 0.0, -0.11001747899],
            '\\sum H_0j N_j': [
                -0.37851801553000003, -0.82094167005, -0.37851801553000003,
                -0.82094167005
            ], 'N_0': [2.7, 3.0, 2.7, 3.0],
            '# H psips': [142.0, 1282.0, 142.0, 1282.0],
            '# states': [79, 782, 79, 782],
            '# spawn_events': [50, 463, 50, 463],
            'R_spawn': [0.4879, 0.762, 0.4879, 0.762],
            'time': [0.0, 0.0004, 0.0004, 0.0004]
        }
        pd.testing.assert_frame_equal(pd.DataFrame.from_dict(data),
                                      pd.DataFrame.from_dict(exp_data))

    def test_longer_fciqmc(self):
        """Have run the calculation for 11 report loops."""
        data = testcode.testcode_data("tests/hande_files/longer_ueg.out")
        exp_data = {
            'iterations': [10, 20, 50, 80, 110],
            'Shift': [
                0.0, 0.0, -0.27841254151, -0.47541378501000003, -0.46025762956
            ], '\\sum H_0j N_j': [
                -0.11688723712000001, -0.33318324599, -1.2617267184,
                -1.1590407344, -0.94329092759
            ], 'N_0': [2.0, 2.0, 2.5, 3.0, 3.3],
            '# H psips': [8.0, 33.0, 201.0, 442.0, 416.0],
            '# states': [7, 30, 188, 432, 405],
            '# spawn_events': [2, 3, 33, 88, 66],
            'R_spawn': [0.1058, 0.1286, 0.1651, 0.186, 0.1773],
            'time': [0.0, 0.0, 0.0, 0.0004, 0.0004]
        }
        pd.testing.assert_frame_equal(pd.DataFrame.from_dict(data),
                                      pd.DataFrame.from_dict(exp_data))

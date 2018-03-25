import unittest
import pandas as pd
import numpy

import pandas.util.testing as pd_testing

import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))
import pyblock

class ErrorTests(unittest.TestCase):
    def setUp(self):
        rand = numpy.random.RandomState(seed=5)
        ndata = 2**5
        data = pd.DataFrame({1:rand.randn(ndata), 2:rand.randn(ndata)})
        (self.data_len, self.reblock, self.cov_12) = pyblock.pd_utils.reblock(data)
        self.cov_12 = self.cov_12.xs(1, level=1)[2]
        self.ratio = pd.DataFrame({'standard error': {0: 0.66458146620536607, 1: 0.75702195760340563, 2: 0.71587945049229107, 3: 0.68021031918966046, 4: 0.16281610758487003}, 'optimal block': {0: '', 1: '', 2: '', 3: '<---    ', 4: ''}, 'mean': {0: -0.086689548180912318, 1: -0.086689548180912304, 2: -0.086689548180912318, 3: -0.08668954818091229, 4: -0.08668954818091229}})
        self.ratio.index.name = 'reblock'
    def tearDown(self):
        del self.data_len, self.reblock, self.cov_12
    def test_ratio(self):
        ratio = pyblock.error.ratio(self.reblock[1], self.reblock[2], self.cov_12, self.data_len)
        # pain to compare pandas objects so compare numpy arrays directly...
        test = ratio.drop('optimal block', axis=1).values
        benchmark = self.ratio.drop('optimal block', axis=1).values
        numpy.testing.assert_array_almost_equal(benchmark, test, decimal=8)
        self.assertEqual(ratio[ratio['optimal block'] != ''].index[0], 3)
    def test_optimal(self):
        self.reblock[(2,'optimal block')] = ''
        self.reblock.loc[1,(2,'optimal block')] = '<---'
        ratio = pyblock.error.ratio(self.reblock[1], self.reblock[2], self.cov_12, self.data_len)
        self.assertEqual(ratio[ratio['optimal block'] != ''].index[0], 3)
    def test_ratio_single(self):
        ratio = pyblock.error.ratio(self.reblock.ix[4,1], self.reblock.ix[4,2], self.cov_12[4], self.data_len[4])
        for key in ('mean', 'standard error'):
            self.assertAlmostEqual(self.ratio.ix[4,key], ratio[key], 8)
    def test_product_single(self):
        product = pyblock.error.product(self.reblock.ix[4,1], self.reblock.ix[4,2], self.cov_12[4], self.data_len[4])
        benchmark = pd.Series({'mean':-0.00561329, 'standard error':0.01576336})
        for key in ('mean', 'standard error'):
            self.assertAlmostEqual(benchmark[key], product[key], 8)
    def test_subtraction_single(self):
        subtraction = pyblock.error.subtraction(self.reblock.ix[4,1], self.reblock.ix[4,2], self.cov_12[4], self.data_len[4])
        benchmark = pd.Series({'mean':-0.27652270, 'standard error':0.17002344})
        for key in ('mean', 'standard error'):
            self.assertAlmostEqual(benchmark[key], subtraction[key], 8)
    def test_addition_single(self):
        addition = pyblock.error.addition(self.reblock.ix[4,1], self.reblock.ix[4,2], self.cov_12[4], self.data_len[4])
        benchmark = pd.Series({'mean':0.23240407, 'standard error':0.06664526})
        for key in ('mean', 'standard error'):
            self.assertAlmostEqual(benchmark[key], addition[key], 8)

class ErrorFmtTests(unittest.TestCase):
    def test1(self):
        self.assertEqual(pyblock.error.pretty_fmt_err(1.2345, 0.01), '1.23(1)')
    def test2(self):
        self.assertEqual(pyblock.error.pretty_fmt_err(12331, 40), '12330(40)')
    def test3(self):
        self.assertEqual(pyblock.error.pretty_fmt_err(float('inf'), 100), str(float('inf')))
    def test4(self):
        self.assertEqual(pyblock.error.pretty_fmt_err(100.01, float('inf')), '100.01(%s)' % (float('inf')))
    def test5(self):
        self.assertEqual(pyblock.error.pretty_fmt_err(100.01, 0.0), '100.01(0)')

def main():
    unittest.main()

if __name__ == '__main__':

    main()

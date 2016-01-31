import numpy
import unittest

import os
import sys
sys.path.append(
    os.path.abspath(os.path.join(os.path.dirname(__file__), '../'))
)
import pyblock
import pyblock.tests.base as tests_base

class BlockTest(unittest.TestCase):
    def check_stats(self, benchmark, test):
        self.assertEqual(len(benchmark), len(test))
        for i in range(len(benchmark)):
            for (test_val, bench_val) in zip(benchmark[i], test[i]):
                if isinstance(test_val, numpy.ndarray):
                    numpy.testing.assert_array_almost_equal(test_val,
                                                            bench_val,
                                                            decimal=8)
                else:
                    self.assertEqual(test_val, bench_val)

class BlockingTests1D(BlockTest):
    def setUp(self):
        self.data = tests_base.data_1D
        self.weights = None
        self.benchmark = tests_base.reblock_1D
        self.benchmark_opt = tests_base.reblock_1D_opt
        self.benchmark_ddof = tests_base.reblock_1D_ddof
    def tearDown(self):
        del self.data
        del self.weights
        del self.benchmark
        del self.benchmark_opt
    def test_reblock(self):
        stats = pyblock.blocking.reblock(self.data, weights=self.weights)
        self.check_stats(self.benchmark, stats)
    def test_reblock_optimal(self):
        stats = pyblock.blocking.reblock(self.data, weights=self.weights)
        self.assertEqual(self.benchmark_opt,
                         pyblock.blocking.find_optimal_block(len(self.data),
                                                             stats))
    def test_reblock_ddof(self):
        stats = pyblock.blocking.reblock(self.data,
                                         ddof=0,
                                         weights=self.weights)
        self.check_stats(self.benchmark_ddof, stats)


class UnitWeightedBlockingTests1D(BlockingTests1D):
    def setUp(self):
        super(UnitWeightedBlockingTests1D, self).setUp()
        self.weights = numpy.ones((self.data.shape[0]))


class WeightedBlockingTests1D(BlockingTests1D):
    def setUp(self):
        self.data = tests_base.data_1D
        self.weights = tests_base.weights
        self.benchmark = tests_base.weighted_reblock_1D
        self.benchmark_opt = tests_base.weighted_reblock_1D_opt
        self.benchmark_ddof = tests_base.weighted_reblock_1D_ddof
    def test_weight_negative(self):
        with self.assertRaises(RuntimeError):
            pyblock.blocking.reblock(self.data, weights=-self.weights)
    def test_weight_len(self):
        with self.assertRaises(RuntimeError):
            pyblock.blocking.reblock(self.data, weights=self.weights[:-1])


class BlockingTests2D(BlockTest):
    def setUp(self):
        self.data = tests_base.data_2D
        self.weights = None
        self.benchmark = tests_base.reblock_2D
        self.benchmark_opt = tests_base.reblock_2D_opt
    def tearDown(self):
        del self.data
    def test_reblock(self):
        stats = pyblock.blocking.reblock(self.data, weights=self.weights)
        self.check_stats(self.benchmark, stats)
    def test_reblock_row(self):
        stats = pyblock.blocking.reblock(self.data.transpose(),
                                         rowvar=0,
                                         weights=self.weights)
        self.check_stats(self.benchmark, stats)
    def test_reblock_rowvar(self):
        with self.assertRaises(ValueError):
            pyblock.blocking.reblock(self.data, ddof=1.2, weights=self.weights)
    def test_reblock_optimal(self):
        stats = pyblock.blocking.reblock(self.data, weights=self.weights)
        self.assertEqual(self.benchmark_opt,
                         pyblock.blocking.find_optimal_block(len(self.data[0]),
                                                             stats))


class UnitWeightedBlockingTests2D(BlockingTests2D):
    def setUp(self):
        super(UnitWeightedBlockingTests2D, self).setUp()
        self.weights = numpy.ones((self.data.shape[1]))


class WeightedBlockingTests2D(BlockingTests2D):
    def setUp(self):
        self.data = tests_base.data_2D
        self.weights = tests_base.weights
        self.benchmark = tests_base.weighted_reblock_2D
        self.benchmark_opt = tests_base.weighted_reblock_2D_opt
    def test_weight_ndim(self):
        weights_2D = numpy.vstack((self.weights, self.weights))
        with self.assertRaises(RuntimeError):
            pyblock.blocking.reblock(self.data, weights=weights_2D)


class BlockingTests3D(unittest.TestCase):
    def setUp(self):
        self.data = numpy.random.randn(3, 3, 3)
    def tearDown(self):
        del self.data
    def test_reblock(self):
        with self.assertRaises(RuntimeError):
            pyblock.blocking.reblock(self.data)


def main():
    unittest.main()

if __name__ == '__main__':

    main()

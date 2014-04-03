import numpy
import unittest

import base
import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))
import pyblock

class BlockTest(unittest.TestCase):
    def check_stats(self, benchmark, test):
        self.assertEqual(len(benchmark), len(test))
        for i in range(len(benchmark)):
            for (test_val, bench_val) in zip(benchmark[i], test[i]):
                if isinstance(test_val, numpy.ndarray):
                    numpy.testing.assert_array_almost_equal(test_val, bench_val, decimal=8)
                else:
                    self.assertEqual(test_val, bench_val)

class BlockingTests1D(BlockTest):
    def setUp(self):
        self.data = base.corr_data(10, 6, 5)
    def tearDown(self):
        del self.data
    def test_reblock(self):
        stats = pyblock.blocking.reblock(self.data)
        self.check_stats(base.reblock_1D, stats)
    def test_reblock_ddof(self):
        stats = pyblock.blocking.reblock(self.data, ddof=0)
        benchmark = [
            (0, 1024, numpy.array(0.20848489810011175), numpy.array(0.5054792211848923), numpy.array(0.022217831846027897), numpy.array(0.0004909493612996486)),
            (1, 512, numpy.array(0.20848489810011178), numpy.array(0.50060779625705), numpy.array(0.0312689878644089), numpy.array(0.000977155870762778)),
            (2, 256, numpy.array(0.20848489810011184), numpy.array(0.4936879913708464), numpy.array(0.04391433383637248), numpy.array(0.001940757702936802)),
            (3, 128, numpy.array(0.20848489810011173), numpy.array(0.479275025445437), numpy.array(0.061190980841072294), numpy.array(0.0038244363025670184)),
            (4, 64, numpy.array(0.2084848981001117), numpy.array(0.4558013370952951), numpy.array(0.08439132592935121), numpy.array(0.007459209854746044)),
            (5, 32, numpy.array(0.20848489810011178), numpy.array(0.3937176741485477), numpy.array(0.11092194245117652), numpy.array(0.013865242806397065)),
            (6, 16, numpy.array(0.20848489810011178), numpy.array(0.31769785169909), numpy.array(0.1409117302824471), numpy.array(0.02490991000786203)),
            (7, 8, numpy.array(0.20848489810011178), numpy.array(0.20018503798159912), numpy.array(0.15818700878295883), numpy.array(0.03954675219573971)),
            (8, 4, numpy.array(0.20848489810011178), numpy.array(0.1580460711956876), numpy.array(0.19877504319939637), numpy.array(0.07027759048847104)),
            (9, 2, numpy.array(0.20848489810011175), numpy.array(0.005749440656203271), numpy.array(0.05361641845649181), numpy.array(0.026808209228245904))
        ]
        self.check_stats(benchmark, stats)
    def test_reblock_optimal(self):
        stats = pyblock.blocking.reblock(self.data)
        self.assertEqual(base.reblock_1D_opt, pyblock.blocking.find_optimal_block(len(self.data), stats))


class BlockingTests2D(BlockTest):
    def setUp(self):
        self.data = numpy.array((base.corr_data(10, 6, 5), base.corr_data(10, 6, 7)))
    def tearDown(self):
        del self.data
    def test_reblock(self):
        stats = pyblock.blocking.reblock(self.data)
        self.check_stats(base.reblock_2D, stats)
    def test_reblock_row(self):
        stats = pyblock.blocking.reblock(self.data.transpose(), rowvar=0)
        self.check_stats(base.reblock_2D, stats)
    def test_reblock_rowvar(self):
        with self.assertRaises(ValueError):
            pyblock.blocking.reblock(self.data, ddof=1.2)
    def test_reblock_optimal(self):
        stats = pyblock.blocking.reblock(self.data)
        self.assertEqual(base.reblock_2D_opt, pyblock.blocking.find_optimal_block(len(self.data[0]), stats))


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

'''Tests for state_histograms.py.'''


import os
import sys
import pkgutil
import unittest
import pandas as pd
from numpy import (around, floor, log10, arange, zeros, array_equal)
from pyhande.state_histograms import (read_state_histogram_file,
                                      collect_state_histogram_data,
                                      average_histograms)


class GenerateReferenceData:
    ''' A wrapper class to generate the exact data'''
    def __init__(self):
        bpd = 5
        bins_nw_1 = around(arange(0.0, 6.0, 1.0/bpd), 4)
        bins_nw_2 = around(arange(0.0, 7.0, 1.0/bpd), 4)
        dataset_1 = {'bin_edges': bins_nw_1}
        dataset_2 = {'bin_edges': bins_nw_1}
        dataset_3 = {'bin_edges': bins_nw_2}
        dataset_4 = {'bin_edges': bins_nw_2}
        dmqmc_cols = ['bin_edges']
        fciqmc_cols = ['bin_edges']
        for ex1 in range(0, 5):
            dataset_3[f'Ex.Lvl 0 {ex1}'] = zeros(bins_nw_2.shape[0])
            dataset_4[f'Ex.Lvl 0 {ex1}'] = zeros(bins_nw_2.shape[0])
            fciqmc_cols.append(f'Ex.Lvl 0 {ex1}')
            for ex2 in range(0, 5):
                dataset_1[f'Ex.Lvl {ex1} {ex2}'] = zeros(bins_nw_1.shape[0])
                dataset_2[f'Ex.Lvl {ex1} {ex2}'] = zeros(bins_nw_1.shape[0])
                dmqmc_cols.append(f'Ex.Lvl {ex1} {ex2}')

        dataset_1['Ex.Lvl 0 0'][5] = 1E0
        dataset_1['Ex.Lvl 0 2'][0] = 5E1
        dataset_1['Ex.Lvl 0 3'][0] = 4.6E1
        dataset_1['Ex.Lvl 0 4'][0] = 5.3E1
        dataset_1['Ex.Lvl 2 0'][:6] = [4.1E1, 3.4E1, 2.4E1, 3.3E1, 2.1E1, 8E0]
        dataset_1['Ex.Lvl 2 2'][0] = 1.991E3
        dataset_1['Ex.Lvl 2 3'][0] = 2.042E3
        dataset_1['Ex.Lvl 2 4'][0] = 2.057E3
        dataset_1['Ex.Lvl 3 0'][:5] = [9.01E2, 2.62E2, 1.55E2, 7.5E1, 7E0]
        dataset_1['Ex.Lvl 3 2'][0] = 8.606E3
        dataset_1['Ex.Lvl 3 3'][0] = 1.0041E4
        dataset_1['Ex.Lvl 3 4'][0] = 9.921E3
        dataset_1['Ex.Lvl 4 0'][:5] = [3.029E3, 3.21E2, 1.27E2, 2.2E1, 3E0]
        dataset_1['Ex.Lvl 4 2'][0] = 1.676E4
        dataset_1['Ex.Lvl 4 3'][0] = 2.0093E4
        dataset_1['Ex.Lvl 4 4'][0] = 2.2094E4

        dataset_2['Ex.Lvl 0 0'][6] = 1E0
        dataset_2['Ex.Lvl 0 2'][0] = 3.6E1
        dataset_2['Ex.Lvl 0 3'][0] = 3.3E1
        dataset_2['Ex.Lvl 0 4'][0] = 3.0E1
        dataset_2['Ex.Lvl 2 0'][:6] = [5.4E1, 2.4E1, 2.9E1, 3.2E1, 2.3E1, 4E0]
        dataset_2['Ex.Lvl 2 2'][:2] = [2.105E3, 1E0]
        dataset_2['Ex.Lvl 2 3'][0] = 2.129E3
        dataset_2['Ex.Lvl 2 4'][0] = 1.967E3
        dataset_2['Ex.Lvl 3 0'][:5] = [9.26E2, 2.65E2, 1.44E2, 7.7E1, 1E1]
        dataset_2['Ex.Lvl 3 2'][0] = 8.664E3
        dataset_2['Ex.Lvl 3 3'][0] = 1.0062E4
        dataset_2['Ex.Lvl 3 4'][0] = 9.844E3
        dataset_2['Ex.Lvl 4 0'][:4] = [3.064E3, 3.4E2, 1.17E2, 3.5E1]
        dataset_2['Ex.Lvl 4 2'][0] = 1.6672E4
        dataset_2['Ex.Lvl 4 3'][0] = 2.027E4
        dataset_2['Ex.Lvl 4 4'][0] = 2.1764E4

        dataset_3['Ex.Lvl 0 0'][13] = 1E0
        dataset_3['Ex.Lvl 0 2'][:6] = [5.2E1, 1.1E1, 1.3E1, 2.4E1, 2E1, 1.6E1]
        dataset_3['Ex.Lvl 0 2'][6:9] = [2E1, 0, 1.6E1]
        dataset_3['Ex.Lvl 0 3'][:3] = [7.6E2, 8.2E1, 8E1]
        dataset_3['Ex.Lvl 0 4'][:5] = [1.472E3, 0, 2.8E1, 8E0, 4E0]

        dataset_4['Ex.Lvl 0 0'][14] = 1E0
        dataset_4['Ex.Lvl 0 2'][:6] = [5.2E1, 9E0, 1.2E1, 2.4E1, 2E1, 1.6E1]
        dataset_4['Ex.Lvl 0 2'][6:9] = [2E1, 0, 1.6E1]
        dataset_4['Ex.Lvl 0 3'][:3] = [7.75E2, 7.3E1, 8.8E1]
        dataset_4['Ex.Lvl 0 4'][:5] = [1.377E3, 0, 2.8E1, 8E0, 4E0]

        self.dataset_1 = pd.DataFrame(dataset_1)
        self.dataset_2 = pd.DataFrame(dataset_2)
        self.dataset_3 = pd.DataFrame(dataset_3)
        self.dataset_4 = pd.DataFrame(dataset_4)

        self.dmqmc_cols = dmqmc_cols
        self.fciqmc_cols = fciqmc_cols

        self.dmqmc_bins = bins_nw_1
        self.fciqmc_bins = bins_nw_2


class TestReadStateHistogramFile(unittest.TestCase, GenerateReferenceData):
    ''' Test the read_state_histogram_file function for DMQMC and
    FCIQMC. We are just testing whether the read in excitations are
    appropriate and if the data is read in properly. The averaging
    function will be tested as well so we test multiple outputs for
    each method here.
    '''
    def setUp(self):
        ''' Generate the exact data and read in the DMQMC/FCIQMC data'''
        GenerateReferenceData.__init__(self)

        path = 'state_histogram_files'
        dmqmc_1 = f'{path}/EXLEVELPOPS-RNG117-IREPORT00250'
        df_1, mex1_1, mex2_1 = read_state_histogram_file(dmqmc_1)
        df_1['bin_edges'] = around(log10(df_1['bin_edges']), 4)
        self.df_1, self.mex1_1, self.mex2_1 = df_1, mex1_1, mex2_1
        dmqmc_2 = f'{path}/EXLEVELPOPS-RNG121-IREPORT00250'
        df_2, mex1_2, mex2_2 = read_state_histogram_file(dmqmc_2)
        df_2['bin_edges'] = around(log10(df_2['bin_edges']), 4)
        self.df_2, self.mex1_2, self.mex2_2 = df_2, mex1_2, mex2_2
        fciqmc_1 = f'{path}/EXLEVELPOPS-RNG7-IREPORT002250'
        df_3, mex1_3, mex2_3 = read_state_histogram_file(fciqmc_1)
        df_3['bin_edges'] = around(log10(df_3['bin_edges']), 4)
        self.df_3, self.mex1_3, self.mex2_3 = df_3, mex1_3, mex2_3
        fciqmc_2 = f'{path}/EXLEVELPOPS-RNG7-IREPORT002500'
        df_4, mex1_4, mex2_4 = read_state_histogram_file(fciqmc_2)
        df_4['bin_edges'] = around(log10(df_4['bin_edges']), 4)
        self.df_4, self.mex1_4, self.mex2_4 = df_4, mex1_4, mex2_4

    def test_dmqmc_maximum_excitation_one(self):
        ''' Test the DMQMC histograms first maximum excitation.'''
        self.assertEqual(self.mex1_1, 4)
        self.assertEqual(self.mex1_2, 4)

    def test_dmqmc_maximum_excitation_two(self):
        ''' Test the DMQMC histograms second maximum excitation.'''
        self.assertEqual(self.mex2_1, 4)
        self.assertEqual(self.mex2_2, 4)

    def test_dmqmc_state_histogram_bins(self):
        ''' Test the bins for the DMQMC histograms.'''
        self.assertTrue(array_equal(self.df_1['bin_edges'], self.dmqmc_bins))
        self.assertTrue(array_equal(self.df_2['bin_edges'], self.dmqmc_bins))

    def test_dmqmc_state_hsitogram_columns(self):
        ''' Test the column names for the DMQMC histograms.'''
        self.assertTrue(array_equal(self.df_1.columns, self.dmqmc_cols))
        self.assertTrue(array_equal(self.df_2.columns, self.dmqmc_cols))

    def test_dmqmc_state_histogram_data(self):
        ''' Test the data within the histogram for the DMQMC histograms.'''
        self.assertTrue(self.df_1.equals(self.dataset_1))
        self.assertTrue(self.df_2.equals(self.dataset_2))

    def test_fciqmc_maximum_excitation_one(self):
        ''' Test the FCIQMC histograms first maximum excitation.'''
        self.assertEqual(self.mex1_3, 0)
        self.assertEqual(self.mex1_4, 0)

    def test_fciqmc_maximum_excitation_two(self):
        ''' Test the FCIQMC histograms second maximum excitation.'''
        self.assertEqual(self.mex2_3, 4)
        self.assertEqual(self.mex2_4, 4)

    def test_fciqmc_state_histogram_bins(self):
        ''' Test the bins for the FCIQMC histograms.'''
        self.assertTrue(array_equal(self.df_3['bin_edges'], self.fciqmc_bins))
        self.assertTrue(array_equal(self.df_4['bin_edges'], self.fciqmc_bins))

    def test_fciqmc_state_hsitogram_columns(self):
        ''' Test the column names for the FCIQMC histograms.'''
        self.assertTrue(array_equal(self.df_3.columns, self.fciqmc_cols))
        self.assertTrue(array_equal(self.df_4.columns, self.fciqmc_cols))

    def test_fciqmc_state_histogram_data(self):
        ''' Test the data within the histogram for the FCIQMC histograms.'''
        self.assertTrue(self.df_3.equals(self.dataset_3))
        self.assertTrue(self.df_4.equals(self.dataset_4))


class TestCollectStateHistogramData(unittest.TestCase, GenerateReferenceData):
    ''' Test the function to collect data from a state histogram file.
    Because much of this code is inherently tested by the read unit tests
    we only do minimal checks on the returned data and instead focus on if
    the appropriate exceptions are raised and data formats are appropriate.
    '''

    def setUp(self):
        ''' Generate the exact data and call the read in wrapper'''
        GenerateReferenceData.__init__(self)

        path = 'state_histogram_files'
        dmqmc_1 = f'{path}/EXLEVELPOPS-RNG117-IREPORT00250'
        dmqmc_2 = f'{path}/EXLEVELPOPS-RNG121-IREPORT00250'
        fciqmc_1 = f'{path}/EXLEVELPOPS-RNG7-IREPORT002250'
        fciqmc_2 = f'{path}/EXLEVELPOPS-RNG7-IREPORT002500'

        dmqmc_files = [dmqmc_1, dmqmc_2]
        fciqmc_files = [fciqmc_1, fciqmc_2]

        self.dmqmc_group_1 = collect_state_histogram_data(dmqmc_files, False)
        self.dmqmc_group_2 = collect_state_histogram_data(dmqmc_files, True)
        self.fciqmc_group_1 = collect_state_histogram_data(fciqmc_files, True)
        self.fciqmc_group_2 = collect_state_histogram_data(fciqmc_files, False)

        self.mixed_files = [dmqmc_1, dmqmc_2, fciqmc_1, fciqmc_2]

    def test_dmqmc_grouped_keys(self):
        ''' Test the group keys are returned in the same order'''
        self.assertEqual('00250', list(self.dmqmc_group_1[0].keys())[0])
        self.assertEqual('00250', list(self.dmqmc_group_1[1].keys())[0])
        self.assertEqual('00250', list(self.dmqmc_group_1[2].keys())[0])
        self.assertEqual('00250', list(self.dmqmc_group_1[3].keys())[0])

    def test_dmqmc_static_grouped(self):
        ''' Test DMQMC is unaltered by the FCIQMC flag'''
        self.assertEqual(len(self.dmqmc_group_1), len(self.dmqmc_group_2))
        gfiles_1, gdata_1, gmex1_1, gmex2_1 = self.dmqmc_group_1
        gfiles_2, gdata_2, gmex1_2, gmex2_2 = self.dmqmc_group_2
        self.assertEqual(len(gfiles_1), len(gfiles_2))
        self.assertEqual(len(gdata_1), len(gdata_2))
        self.assertEqual(len(gmex1_1), len(gmex1_2))
        self.assertEqual(len(gmex2_1), len(gmex2_2))

    def test_fciqmc_grouped_keys(self):
        ''' Test the group keys are returned in the same order'''
        self.assertEqual('002250', list(self.fciqmc_group_1[0].keys())[0])
        self.assertEqual('002250', list(self.fciqmc_group_1[1].keys())[0])
        self.assertEqual('002250', list(self.fciqmc_group_1[2].keys())[0])
        self.assertEqual('002250', list(self.fciqmc_group_1[3].keys())[0])
        self.assertEqual('002500', list(self.fciqmc_group_2[0].keys())[1])
        self.assertEqual('002500', list(self.fciqmc_group_2[1].keys())[1])
        self.assertEqual('002500', list(self.fciqmc_group_2[2].keys())[1])
        self.assertEqual('002500', list(self.fciqmc_group_2[3].keys())[1])

    def test_fciqmc_static_grouped(self):
        ''' Test DMQMC is unaltered by the FCIQMC flag'''
        self.assertEqual(len(self.fciqmc_group_1), len(self.fciqmc_group_2))
        gfiles_1, gdata_1, gmex1_1, gmex2_1 = self.fciqmc_group_1
        gfiles_2, gdata_2, gmex1_2, gmex2_2 = self.fciqmc_group_2
        self.assertNotEqual(len(gfiles_1), len(gfiles_2))
        self.assertNotEqual(len(gdata_1), len(gdata_2))
        self.assertNotEqual(len(gmex1_1), len(gmex1_2))
        self.assertNotEqual(len(gmex2_1), len(gmex2_2))

    def test_raise_runtime(self):
        ''' Test non-matching data sets throw an exception'''
        self.assertRaises(RuntimeError, collect_state_histogram_data,
                          self.mixed_files, False)
        self.assertRaises(RuntimeError, collect_state_histogram_data,
                          self.mixed_files, True)


class TestAverageHistograms(unittest.TestCase):
    ''' Test the function to average state histogram data'''

    pass

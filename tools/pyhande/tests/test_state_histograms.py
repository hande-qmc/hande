'''Tests for state_histograms.py.'''

import unittest
from pandas import DataFrame
from numpy import (around, log10, arange, zeros, array_equal)
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

        self.dataset_1 = DataFrame(dataset_1)
        self.dataset_2 = DataFrame(dataset_2)
        self.dataset_3 = DataFrame(dataset_3)
        self.dataset_4 = DataFrame(dataset_4)

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

        path = 'tests/state_histogram_files'
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

        path = 'tests/state_histogram_files'
        dmqmc_1 = f'{path}/EXLEVELPOPS-RNG117-IREPORT00000'
        dmqmc_2 = f'{path}/EXLEVELPOPS-RNG121-IREPORT00000'
        dmqmc_3 = f'{path}/EXLEVELPOPS-RNG117-IREPORT00250'
        dmqmc_4 = f'{path}/EXLEVELPOPS-RNG121-IREPORT00250'
        fciqmc_1 = f'{path}/EXLEVELPOPS-RNG7-IREPORT000000'
        fciqmc_2 = f'{path}/EXLEVELPOPS-RNG7-IREPORT002250'
        fciqmc_3 = f'{path}/EXLEVELPOPS-RNG7-IREPORT002500'

        dmqmc_files = [dmqmc_1, dmqmc_2, dmqmc_3, dmqmc_4]
        fciqmc_files = [fciqmc_1, fciqmc_2, fciqmc_3]
        self.mixed_files = dmqmc_files + fciqmc_files

        self.dmqmc_group_1 = collect_state_histogram_data(dmqmc_files, False)
        self.dmqmc_group_2 = collect_state_histogram_data(dmqmc_files, True)
        self.fciqmc_group_1 = collect_state_histogram_data(fciqmc_files, True)
        self.fciqmc_group_2 = collect_state_histogram_data(fciqmc_files, False)

        gfiles_1, gdata_1, gmex1_1, gmex2_1 = self.dmqmc_group_1
        self.gfiles_1 = gfiles_1
        self.gdata_1 = gdata_1
        self.gmex1_1 = gmex1_1
        self.gmex2_1 = gmex2_1
        gfiles_2, gdata_2, gmex1_2, gmex2_2 = self.dmqmc_group_2
        self.gfiles_2 = gfiles_2
        self.gdata_2 = gdata_2
        self.gmex1_2 = gmex1_2
        self.gmex2_2 = gmex2_2
        gfiles_3, gdata_3, gmex1_3, gmex2_3 = self.fciqmc_group_1
        self.gfiles_3 = gfiles_3
        self.gdata_3 = gdata_3
        self.gmex1_3 = gmex1_3
        self.gmex2_3 = gmex2_3
        gfiles_4, gdata_4, gmex1_4, gmex2_4 = self.fciqmc_group_2
        self.gfiles_4 = gfiles_4
        self.gdata_4 = gdata_4
        self.gmex1_4 = gmex1_4
        self.gmex2_4 = gmex2_4

    def test_dmqmc_grouped_keys(self):
        ''' Test the DMQMC group keys are proper even for the FCIQMC flag.'''
        self.assertEqual('00000', list(self.dmqmc_group_1[0].keys())[0])
        self.assertEqual('00000', list(self.dmqmc_group_1[1].keys())[0])
        self.assertEqual('00000', list(self.dmqmc_group_1[2].keys())[0])
        self.assertEqual('00000', list(self.dmqmc_group_1[3].keys())[0])
        self.assertEqual('00250', list(self.dmqmc_group_1[0].keys())[1])
        self.assertEqual('00250', list(self.dmqmc_group_1[1].keys())[1])
        self.assertEqual('00250', list(self.dmqmc_group_1[2].keys())[1])
        self.assertEqual('00250', list(self.dmqmc_group_1[3].keys())[1])
        self.assertEqual('00000', list(self.dmqmc_group_2[0].keys())[0])
        self.assertEqual('00000', list(self.dmqmc_group_2[1].keys())[0])
        self.assertEqual('00000', list(self.dmqmc_group_2[2].keys())[0])
        self.assertEqual('00000', list(self.dmqmc_group_2[3].keys())[0])

    def test_dmqmc_grouped_shapes(self):
        ''' Test DMQMC is altered by the FCIQMC flag in the expected way.'''
        self.assertEqual(len(self.dmqmc_group_1), len(self.dmqmc_group_2))
        self.assertNotEqual(len(self.gfiles_1), len(self.gfiles_2))
        self.assertNotEqual(len(self.gdata_1), len(self.gdata_2))
        self.assertNotEqual(len(self.gmex1_1), len(self.gmex1_2))
        self.assertNotEqual(len(self.gmex2_1), len(self.gmex2_2))
        self.assertEqual(int(0.5*len(self.gfiles_1)), len(self.gfiles_2))
        self.assertEqual(int(0.5*len(self.gdata_1)), len(self.gdata_2))
        self.assertEqual(int(0.5*len(self.gmex1_1)), len(self.gmex1_2))
        self.assertEqual(int(0.5*len(self.gmex2_1)), len(self.gmex2_2))
        self.assertEqual(len(self.gfiles_1), 2)
        self.assertEqual(len(self.gdata_1), 2)
        self.assertEqual(len(self.gmex1_1), 2)
        self.assertEqual(len(self.gmex2_1), 2)
        self.assertEqual(len(list(self.gfiles_1.values())[0]), 2)
        self.assertEqual(len(list(self.gdata_1.values())[0]), 2)
        self.assertEqual(len(list(self.gmex1_1.values())[0]), 2)
        self.assertEqual(len(list(self.gmex2_1.values())[0]), 2)

    def test_dmqmc_returned_data(self):
        ''' Test the DMQMC data is returned correctly.'''
        files_1 = list(self.gfiles_1.values())
        data_1 = list(self.gdata_1.values())
        mex1_1 = list(self.gmex1_1.values())
        mex2_1 = list(self.gmex2_1.values())
        files_2 = list(self.gfiles_2.values())
        data_2 = list(self.gdata_2.values())
        mex1_2 = list(self.gmex1_2.values())
        mex2_2 = list(self.gmex2_2.values())
        df_1a, df_1b = data_1[1]
        _, _, df_2a, df_2b = data_2[0]
        df_1a['bin_edges'] = around(log10(df_1a['bin_edges']), 4)
        df_1b['bin_edges'] = around(log10(df_1b['bin_edges']), 4)
        df_2a['bin_edges'] = around(log10(df_2a['bin_edges']), 4)
        df_2b['bin_edges'] = around(log10(df_2b['bin_edges']), 4)
        self.assertTrue(array_equal(files_1[0] + files_1[1], files_2[0]))
        self.assertTrue(array_equal(data_1[0] + data_1[1], data_2[0]))
        self.assertTrue(array_equal(mex1_1[0] + mex1_1[1], mex1_2[0]))
        self.assertTrue(array_equal(mex2_1[0] + mex2_1[1], mex2_2[0]))
        self.assertTrue(df_1a.equals(self.dataset_1))
        self.assertTrue(df_1b.equals(self.dataset_2))
        self.assertTrue(df_2a.equals(self.dataset_1))
        self.assertTrue(df_2b.equals(self.dataset_2))

    def test_fciqmc_grouped_keys(self):
        ''' Test the FCIQMC group keys are returned appropriately.'''
        self.assertEqual('000000', list(self.fciqmc_group_1[0].keys())[0])
        self.assertEqual('000000', list(self.fciqmc_group_1[1].keys())[0])
        self.assertEqual('000000', list(self.fciqmc_group_1[2].keys())[0])
        self.assertEqual('000000', list(self.fciqmc_group_1[3].keys())[0])
        self.assertEqual('000000', list(self.fciqmc_group_2[0].keys())[0])
        self.assertEqual('000000', list(self.fciqmc_group_2[1].keys())[0])
        self.assertEqual('000000', list(self.fciqmc_group_2[2].keys())[0])
        self.assertEqual('000000', list(self.fciqmc_group_2[3].keys())[0])
        self.assertEqual('002250', list(self.fciqmc_group_2[0].keys())[1])
        self.assertEqual('002250', list(self.fciqmc_group_2[1].keys())[1])
        self.assertEqual('002250', list(self.fciqmc_group_2[2].keys())[1])
        self.assertEqual('002250', list(self.fciqmc_group_2[3].keys())[1])
        self.assertEqual('002500', list(self.fciqmc_group_2[0].keys())[2])
        self.assertEqual('002500', list(self.fciqmc_group_2[1].keys())[2])
        self.assertEqual('002500', list(self.fciqmc_group_2[2].keys())[2])
        self.assertEqual('002500', list(self.fciqmc_group_2[3].keys())[2])

    def test_fciqmc_grouped_shapes(self):
        ''' Test FCIQMC is altered by the FCIQMC flag in the expected way.'''
        self.assertEqual(len(self.fciqmc_group_1), len(self.fciqmc_group_2))
        self.assertNotEqual(len(self.gfiles_3), len(self.gfiles_4))
        self.assertNotEqual(len(self.gdata_3), len(self.gdata_4))
        self.assertNotEqual(len(self.gmex1_3), len(self.gmex1_4))
        self.assertNotEqual(len(self.gmex2_3), len(self.gmex2_4))
        self.assertEqual(len(self.gfiles_3), int(0.34*len(self.gfiles_4)))
        self.assertEqual(len(self.gdata_3), int(0.34*len(self.gdata_4)))
        self.assertEqual(len(self.gmex1_3), int(0.34*len(self.gmex1_4)))
        self.assertEqual(len(self.gmex2_3), int(0.34*len(self.gmex2_4)))
        self.assertEqual(len(self.gfiles_3), 1)
        self.assertEqual(len(self.gdata_3), 1)
        self.assertEqual(len(self.gmex1_3), 1)
        self.assertEqual(len(self.gmex2_3), 1)
        self.assertEqual(len(list(self.gfiles_3.values())[0]), 3)
        self.assertEqual(len(list(self.gdata_3.values())[0]), 3)
        self.assertEqual(len(list(self.gmex1_3.values())[0]), 3)
        self.assertEqual(len(list(self.gmex2_3.values())[0]), 3)

    def test_fciqmc_returned_data(self):
        ''' Test the FCIQMC data is returned correctly.'''
        files_3 = list(self.gfiles_3.values())
        data_3 = list(self.gdata_3.values())
        mex1_3 = list(self.gmex1_3.values())
        mex2_3 = list(self.gmex2_3.values())
        files_4 = list(self.gfiles_4.values())
        data_4 = list(self.gdata_4.values())
        mex1_4 = list(self.gmex1_4.values())
        mex2_4 = list(self.gmex2_4.values())
        files_4 = files_4[0] + files_4[1] + files_4[2]
        data_4 = data_4[0] + data_4[1] + data_4[2]
        mex1_4 = mex1_4[0] + mex1_4[1] + mex1_4[2]
        mex2_4 = mex2_4[0] + mex2_4[1] + mex2_4[2]
        _, df_3a, df_3b = data_3[0]
        _, df_4a, df_4b = data_4
        df_3a['bin_edges'] = around(log10(df_3a['bin_edges']), 4)
        df_3b['bin_edges'] = around(log10(df_3b['bin_edges']), 4)
        df_4a['bin_edges'] = around(log10(df_4a['bin_edges']), 4)
        df_4b['bin_edges'] = around(log10(df_4b['bin_edges']), 4)
        self.assertTrue(array_equal(files_3[0], files_4))
        self.assertTrue(array_equal(data_3[0], data_4))
        self.assertTrue(array_equal(mex1_3[0], mex1_4))
        self.assertTrue(array_equal(mex2_3[0], mex2_4))
        self.assertTrue(df_3a.equals(self.dataset_3))
        self.assertTrue(df_3b.equals(self.dataset_4))
        self.assertTrue(df_4a.equals(self.dataset_3))
        self.assertTrue(df_4b.equals(self.dataset_4))

    def test_mismatched_data_raises_runtime(self):
        ''' Test non-matching data sets throw an exception.'''
        self.assertRaises(RuntimeError, collect_state_histogram_data,
                          self.mixed_files, False)
        self.assertRaises(RuntimeError, collect_state_histogram_data,
                          self.mixed_files, True)


class TestAverageHistograms(unittest.TestCase):
    ''' Test the function to average state histogram data by comparing
    the data to the exactly calculated data external to the python code.
    '''
    def setUp(self):
        self.dmqmc_exact_mean = zeros(30)
        self.dmqmc_exact_mean[0] = 9.877000000000E+04
        self.dmqmc_exact_mean[1] = 1.097500000000E+03
        self.dmqmc_exact_mean[2] = 4.740000000000E+02
        self.dmqmc_exact_mean[3] = 1.760000000000E+02
        self.dmqmc_exact_mean[4] = 3.900000000000E+01
        self.dmqmc_exact_mean[5] = 7.000000000000E+00
        self.dmqmc_exact_mean[6] = 5.000000000000E-01

        self.dmqmc_exact_sem = zeros(30)
        self.dmqmc_exact_sem[0] = 2.180687964840E+02
        self.dmqmc_exact_sem[1] = 1.525614630239E+01
        self.dmqmc_exact_sem[2] = 1.072380529476E+01
        self.dmqmc_exact_sem[3] = 7.314369419164E+00
        self.dmqmc_exact_sem[4] = 3.162277660168E+00
        self.dmqmc_exact_sem[5] = 2.121320343560E+00
        self.dmqmc_exact_sem[6] = 5.000000000000E-01

        self.fciqmc_exact_mean = zeros(35)
        self.fciqmc_exact_mean[0] = 2.565000000000E+03
        self.fciqmc_exact_mean[1] = 3.210000000000E+02
        self.fciqmc_exact_mean[2] = 2.335000000000E+02
        self.fciqmc_exact_mean[3] = 1.090000000000E+02
        self.fciqmc_exact_mean[4] = 7.700000000000E+01
        self.fciqmc_exact_mean[5] = 5.300000000000E+01
        self.fciqmc_exact_mean[6] = 3.700000000000E+01
        self.fciqmc_exact_mean[7] = 1.700000000000E+01
        self.fciqmc_exact_mean[8] = 1.700000000000E+01
        self.fciqmc_exact_mean[9] = 1.000000000000E+00
        self.fciqmc_exact_mean[10] = 1.000000000000E+00
        self.fciqmc_exact_mean[11] = 1.000000000000E+00
        self.fciqmc_exact_mean[12] = 1.000000000000E+00
        self.fciqmc_exact_mean[13] = 1.000000000000E+00
        self.fciqmc_exact_mean[14] = 5.000000000000E-01

        self.fciqmc_exact_sem = zeros(35)
        self.fciqmc_exact_sem[0] = 4.848195540611E+01
        self.fciqmc_exact_sem[1] = 6.164414002969E+00
        self.fciqmc_exact_sem[2] = 4.092676385936E+00
        self.fciqmc_exact_sem[3] = 7.071067811865E-01
        self.fciqmc_exact_sem[4] = 7.071067811865E-01
        self.fciqmc_exact_sem[5] = 7.071067811865E-01
        self.fciqmc_exact_sem[6] = 7.071067811865E-01
        self.fciqmc_exact_sem[7] = 7.071067811865E-01
        self.fciqmc_exact_sem[8] = 7.071067811865E-01
        self.fciqmc_exact_sem[9] = 7.071067811865E-01
        self.fciqmc_exact_sem[10] = 7.071067811865E-01
        self.fciqmc_exact_sem[11] = 7.071067811865E-01
        self.fciqmc_exact_sem[12] = 7.071067811865E-01
        self.fciqmc_exact_sem[13] = 7.071067811865E-01
        self.fciqmc_exact_sem[14] = 5.000000000000E-01

        path = 'tests/state_histogram_files'
        dmqmc_1 = f'{path}/EXLEVELPOPS-RNG117-IREPORT00250'
        dmqmc_2 = f'{path}/EXLEVELPOPS-RNG121-IREPORT00250'
        fciqmc_1 = f'{path}/EXLEVELPOPS-RNG7-IREPORT002250'
        fciqmc_2 = f'{path}/EXLEVELPOPS-RNG7-IREPORT002500'

        dmqmc_files = [dmqmc_1, dmqmc_2]
        fciqmc_files = [fciqmc_1, fciqmc_2]

        dmqmc_group_1 = collect_state_histogram_data(dmqmc_files, False)
        fciqmc_group_1 = collect_state_histogram_data(fciqmc_files, True)

        dmqmc_analysis = average_histograms(dmqmc_group_1)
        self.dmqmc_mean = dmqmc_analysis['ndets_00250']
        self.dmqmc_sem = dmqmc_analysis['ndets_00250_error']
        fciqmc_analysis = average_histograms(fciqmc_group_1)
        self.fciqmc_mean = fciqmc_analysis['ndets_002250']
        self.fciqmc_sem = fciqmc_analysis['ndets_002250_error']

        self.dmqmc_sem = around(self.dmqmc_sem, 10)
        self.dmqmc_exact_sem = around(self.dmqmc_exact_sem, 10)
        self.fciqmc_sem = around(self.fciqmc_sem, 10)
        self.fciqmc_exact_sem = around(self.fciqmc_exact_sem, 10)

        dmqmc_data = dmqmc_group_1[1]['00250'][0]
        fciqmc_data = fciqmc_group_1[1]['002250'][0]
        mixed_files = {0: [dmqmc_1, fciqmc_1]}
        mixed_data = {0: [dmqmc_data, fciqmc_data]}
        mixed_mex1 = {0: [4, 0]}
        mixed_mex2 = {0: [4, 4]}
        self.mixed_group = (mixed_files, mixed_data, mixed_mex1, mixed_mex2)

    def test_dmqmc_average_histogram(self):
        ''' Test the DMQMC mean and error works out appropriately.'''
        self.assertTrue(array_equal(self.dmqmc_mean, self.dmqmc_exact_mean))
        self.assertTrue(array_equal(self.dmqmc_sem, self.dmqmc_exact_sem))

    def test_fciqmc_average_histogram(self):
        ''' Test the FCIQMC mean and error works out appropriately.'''
        self.assertTrue(array_equal(self.fciqmc_mean, self.fciqmc_exact_mean))
        self.assertTrue(array_equal(self.fciqmc_sem, self.fciqmc_exact_sem))

    def test_mismatched_data_raises_runtime(self):
        ''' Test non-matching data sets throw an exception.'''
        self.assertRaises(RuntimeError, average_histograms, self.mixed_group)

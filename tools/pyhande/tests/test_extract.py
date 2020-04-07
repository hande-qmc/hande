"""Test extract.py."""
import unittest
import copy
import warnings
import pandas as pd
import numpy as np
import pyhande.extract as extract


class ExtendedTestSetUp():
    """Helper for the following classes, initialises common data."""

    def __init__(self):
        self.filename1 = "hande_files/ueg.out"
        self.filename2 = "hande_files/ueg2.out"
        self.columns = [
            'iterations', 'Shift', r'\sum H_0j N_j', 'N_0', '# H psips',
            '# states', '# spawn_events', 'R_spawn', 'time'
            ]
        self.dummy_data1 = pd.DataFrame([
            [10, 0.00000000e+00, -3.78518016e-01, 2.70000000e+00,
             1.42000000e+02, 79, 50, 4.87900000e-01, 4.00000000e-04],
            [20, -1.10017479e-01, -8.20941670e-01, 3.00000000e+00,
             1.28200000e+03, 782, 463, 7.62000000e-01, 8.0000000e-04]
            ], columns=self.columns)
        self.dummy_data2 = pd.DataFrame([
            [10, 0., -0.08311986, 2., 8.11834868, 7, 3, 0.2765, 0.],
            [20, 0., -0.5709596, 2., 73.67254134, 66, 33, 0.3288, 0.],
            [30, -1.40493959e-01, -9.03119184e-01, 2.00000000e+00,
             1.71159164e+02, 162, 104, 3.77400000e-01, 0.0]
            ], columns=self.columns)
        self.dummy_metadata = {
            'system': {
                'nbasis': 38, 'nel': 14, 'nvirt': 24, 'Ms': 0, 'nalpha': 7,
                'nbeta': 7, 'nvirt_alpha': 12, 'nvirt_beta': 12, 'nsym': 19,
                'sym0': 1, 'sym_max': 19, 'nsym_tot': 19, 'sym0_tot': 1,
                'sym_max_tot': 19, 'symmetry': 1, 'tot_sym': False,
                'aufbau_sym': True, 'max_number_excitations': 14,
                'ueg': {
                    'r_s': 1.0, 'ecutoff': 1.0, 'k_fermi': 1.91915829,
                    'E_fermi': 1.84158428, 'ktwist': [0.0, 0.0, 0.0],
                    'L': [3.88512994, 3.88512994, 3.88512994]
                    }
                },
            'qmc': {
                'rng_seed': 1472, 'real_amplitudes': False,
                'real_amplitude_force_32': False, 'spawn_cutoff': 0.01,
                'excit_gen': 'renorm', 'pattempt_update': False,
                'pattempt_zero_accum_data': False, 'pattempt_single': 0.0,
                'pattempt_double': 1.0, 'pattempt_parallel': 0.0, 'tau': 0.1,
                'tau_search': False, 'vary_shift_from': 0.0,
                'vary_shift_from_proje': False, 'initial_shift': 0.0,
                'shift_damping': 0.05, 'walker_length': 50000,
                'spawned_walker_length': 5000, 'D0_population': 2.0,
                'target_particles': 60.0, 'target_reference': False,
                'initiator_approx': False, 'initiator_pop': 3.0, 'ncycles': 10,
                'nreport': 3, 'power_pitzer_min_weight': 0.01,
                'quasi_newton': False, 'quasi_newton_threshold': 1e-05,
                'quasi_newton_value': 1.0, 'use_mpi_barriers': False
                },
            'fciqmc': {
                'select_ref_det_every_nreports': 2147483647,
                'init_spin_inv_D0': False, 'ref_det_factor': 1.5,
                'non_blocking_comm': False, 'doing_load_balancing': False,
                'trial_function': 'single_basis', 'guiding_function': 'none',
                'quadrature_initiator': True, 'replica_tricks': False
                },
            'semi_stoch': {
                'start_iter': 1, 'shift_iter': -1, 'space_type': 'none',
                'target_size': 0, 'write_determ_space': False,
                'projection_mode': 'separate', 'read_id': 2147483647,
                'write_id': 2147483647, 'ci_space': {'ex_level': -1}
                },
            'restart': {
                'read_restart': False, 'read_id': 2147483647,
                'write_restart': False, 'write_id': 2147483647,
                'write_freq': 2147483647, 'write_restart_shift': False,
                'write_shift_id': 2147483647, 'restart_rng': True
                },
            'blocking': {
                'blocking_on_the_fly': False, 'start_save_frequency': -1,
                'start_point_number': -1, 'filename': 'BLOCKING',
                'start_point': -1, 'error_limit': 0.0,
                'blocks_used': 2147483647, 'min_blocks_used': 10,
                'auto_shift_damping': False, 'shift_damping_precision': 2.0,
                'force_shift_damping_opt': False
                },
            'load balancing': {
                'nslots': 1, 'pop': 1000, 'percent': 0.05, 'max_attempts': 2,
                'write_info': False
                },
            'reference': {
                'det': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14],
                'det_ms': 0, 'det_symmetry': 1, 'H00': 13.60355734,
                'F0': 15.69278015, 'hilbert_space_det': [
                    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14
                    ], 'hilbert_space_det_ms': 0,
                'hilbert_space_det_symmetry': 1, 'ex_level': 14
                },
            'logging_in': {
                'calc': 0, 'calc_file': 'CALC', 'spawn': 0,
                'spawn_file': 'SPAWN', 'death': 0, 'death_file': 'DEATH',
                'start_iter': 0, 'end_iter': 9223372036854775807
                },
            'logging': {
                'write_highlevel_values': False, 'calc_unit': 2147483647,
                'write_successful_spawn': False, 'write_failed_spawn': False,
                'spawn_unit': 2147483647, 'write_successful_death': False,
                'write_failed_death': False, 'death_unit': 2147483647
                },
            'calc_type': 'FCIQMC', 'nblooms': 6253.0, 'max_bloom': 9.0,
            'mean_bloom': 3.93, 'input': [
                '', '-- Create output with:',
                '-- $[HANDE DIR]/bin/hande.x ueg.lua > ueg.out 2> ueg.err',
                '-- Note that these settings are just for testing...',
                'sys = ueg {', 'dim = 3,', 'nel = 14,', 'ms = 0,',
                'cutoff = 1,', '}', '', 'fciqmc {', 'sys = sys,',
                'qmc = {', 'tau = 0.1,', 'rng_seed = 1472,', 'init_pop = 2,',
                'mc_cycles = 10,', 'nreports = 3,', 'target_population = 60,',
                'state_size = 50000,', 'spawned_state_size = 5000,', '},', '}',
                ''
                ], 'UUID': 'c04c1500-cfea-4cc5-a90b-f525b4b36ec5',
            'wall_time': 0.07, 'cpu_time': 0.07, 'calculation_time': 0.07
            }
        # metadata2 is similar:
        self.dummy_metadata2 = copy.deepcopy(self.dummy_metadata)
        self.dummy_metadata2['qmc']['rng_seed'] = 389
        self.dummy_metadata2['qmc']['real_amplitudes'] = True
        self.dummy_metadata2['qmc']['tau'] = 0.03
        self.dummy_metadata2['UUID'] = '0a5946bd-e0f6-4971-9d3a-6ae056ea9ca4'
        self.dummy_metadata2['wall_time'] = 0.0
        self.dummy_metadata2['cpu_time'] = 0.0
        self.dummy_metadata2['calculation_time'] = 0.0
        self.dummy_metadata2['input'][2] = \
            '-- $[HANDE DIR]/bin/hande.x ueg2.lua > ueg2.out 2> ueg2.err'
        self.dummy_metadata2['input'][8] = 'cutoff = 1.0,'
        self.dummy_metadata2['input'][14] = 'tau = 0.03,'
        self.dummy_metadata2['input'][15] = 'rng_seed = 389,'
        self.dummy_metadata2['input'][22] = 'real_amplitudes = true,'
        self.dummy_metadata2['input'][23] = '},'
        self.dummy_metadata2['input'][24] = '}'
        self.dummy_metadata2['input'].append('')
        del self.dummy_metadata2['nblooms']
        del self.dummy_metadata2['max_bloom']
        del self.dummy_metadata2['mean_bloom']


class TestExtractDataSets(unittest.TestCase, ExtendedTestSetUp):
    """Test extract.extract_data_sets().
    Don't test "private" functions at the moment (i.e. functions with
    names starting with "_").  [todo] add tests possibly
    """

    def setUp(self):
        ExtendedTestSetUp.__init__(self)

    def test_basic_input(self):
        """Test basic input."""
        data = extract.extract_data_sets([self.filename1, self.filename2])
        pd.testing.assert_frame_equal(data[0][1], self.dummy_data1)
        pd.testing.assert_frame_equal(data[1][1], self.dummy_data2)
        self.assertDictEqual(data[0][0], self.dummy_metadata)
        self.assertDictEqual(data[1][0], self.dummy_metadata2)

    def test_only_one(self):
        """Only pass one file (in list)."""
        data = extract.extract_data_sets([self.filename1])
        pd.testing.assert_frame_equal(data[0][1], self.dummy_data1)
        self.assertDictEqual(data[0][0], self.dummy_metadata)

    def test_compressed(self):
        """Pass a list of a file, differently compressed."""
        data = extract.extract_data_sets([
            self.filename1+".bz2", self.filename1+".gz", self.filename1+".xz"
            ])
        for i in range(3):
            pd.testing.assert_frame_equal(data[i][1], self.dummy_data1)
            self.assertDictEqual(data[i][0], self.dummy_metadata)

    def test_unchanged_mutable(self):
        """Check that mutable objects, such as pd DataFrames, don't
        change when they shouldn't.
        """
        list_filename = copy.copy([self.filename1])
        _ = extract.extract_data_sets([self.filename1])
        self.assertListEqual(list_filename, [self.filename1])


class TestExtractData(unittest.TestCase, ExtendedTestSetUp):
    """Test extract.extract_data()."""

    def setUp(self):
        ExtendedTestSetUp.__init__(self)

    def test_basic_fciqmc_input(self):
        """Test basic input."""
        data = extract.extract_data(self.filename1)
        pd.testing.assert_frame_equal(data[0][1], self.dummy_data1)
        self.assertDictEqual(data[0][0], self.dummy_metadata)

    def test_bz2_fciqmc(self):
        """Extract a compressed file - .bz2."""
        data = extract.extract_data(self.filename1+".bz2")
        pd.testing.assert_frame_equal(data[0][1], self.dummy_data1)
        self.assertDictEqual(data[0][0], self.dummy_metadata)

    def test_gz_fciqmc(self):
        """Extract a compressed file - .gz."""
        data = extract.extract_data(self.filename1+".gz")
        pd.testing.assert_frame_equal(data[0][1], self.dummy_data1)
        self.assertDictEqual(data[0][0], self.dummy_metadata)

    def test_xz_fciqmc(self):
        """Extract a compressed file - .xz."""
        data = extract.extract_data(self.filename1+".xz")
        pd.testing.assert_frame_equal(data[0][1], self.dummy_data1)
        self.assertDictEqual(data[0][0], self.dummy_metadata)

    def test_multiple_fciqmc(self):
        """Have two calculations in this output."""
        data = extract.extract_data("hande_files/multi_ueg.out")
        # First calculation.
        # Modify some previous data for comparison
        self.dummy_data1.at[0, 'time'] = 0.0
        self.dummy_data1.at[1, 'time'] = 0.0004
        self.dummy_metadata['input'][2] = (
            '-- $[HANDE DIR]/bin/hande.x multi.ueg.lua > multi.ueg.out 2> '
            'multi.ueg.err'
            )
        self.dummy_metadata['input'].append('')
        self.dummy_metadata['input'][12:] = self.dummy_metadata['input'][11:-1]
        self.dummy_metadata['input'][11] = 'for i=1,2 do'
        self.dummy_metadata['input'][-1] = 'end'
        self.dummy_metadata['input'].append('')
        self.dummy_metadata['UUID'] = 'f11d432f-2eec-45e6-801b-56e458b8f3c0'
        self.dummy_metadata['wall_time'] = 0.12
        self.dummy_metadata['cpu_time'] = 0.12
        self.dummy_metadata['calculation_time'] = 0.06
        # Test.
        pd.testing.assert_frame_equal(data[0][1], self.dummy_data1)
        self.assertDictEqual(data[0][0], self.dummy_metadata)
        # Second calculation.
        # Modify some previous data for comparison
        self.dummy_data1.at[0, 'time'] = 0.0004
        # [todo] This needs to be investigated!!! Why not 0 as before?
        self.dummy_metadata['qmc']['pattempt_parallel'] = 5.28562263e+180
        # Test.
        pd.testing.assert_frame_equal(data[1][1], self.dummy_data1)
        self.assertDictEqual(data[1][0], self.dummy_metadata)

    def test_replica_fciqmc(self):
        """Placeholder warning for replica test."""
        # [todo] Once fixed, there should be a test for replica tricks.
        warnings.warn("Once fixed, add test for replica tricks for extract.py")

    def test_basic_ccmc_input(self):
        """Test CCMC."""
        data = extract.extract_data("hande_files/ccmc_ueg.out")
        self.columns.append(self.columns[-1])
        self.columns[-3:-1] = ['# attempts', 'R_spawn']
        dummy_data = pd.DataFrame([
            [10, 0., -0.59349668, 1.86285714, 27., 12, 10, 18, 0.5119, 0.],
            [20, 0., -4.40680439e-01, 2.02475127e+00, 3.90000000e+01, 26,
             19, 42, 4.79900000e-01, 4.00000000e-04],
            [30, -1.53365134e-02, -5.83660200e-01, 2.05614830e+00,
             5.30000000e+01, 34, 12, 55, 3.91400000e-01, 0.0]
            ], columns=self.columns)
        pd.testing.assert_frame_equal(data[0][1], dummy_data)
        self.dummy_metadata['UUID'] = 'acc004f5-bbc6-4b5c-a831-406b90239e98'
        self.dummy_metadata['calc_type'] = 'CCMC'
        self.dummy_metadata['calculation_time'] = 0.0
        self.dummy_metadata['cpu_time'] = 0.01
        self.dummy_metadata['wall_time'] = 0.01
        self.dummy_metadata['reference']['ex_level'] = 3
        self.dummy_metadata['qmc']['target_particles'] = 30.0
        del self.dummy_metadata['fciqmc']
        del self.dummy_metadata['load balancing']
        self.dummy_metadata['ccmc'] = {
            'cluster_multispawn_threshold': 1.79769313e+308,
            'density_matrices': False,
            'density_matrix_file': 'RDM', 'even_selection': False,
            'full_nc': False, 'linked': False, 'move_freq': 5,
            'multiref': False, 'vary_shift_reference': False
            }
        self.dummy_metadata['input'] = [
            '', '-- Create output with:',
            '-- $[HANDE DIR]/bin/hande.x ccmc_ueg.lua > ccmc_ueg.out 2> '
            'ccmc_ueg.err',
            '-- Note that these settings are just for testing...',
            'sys = ueg {', 'dim = 3,', 'nel = 14,', 'ms = 0,', 'cutoff = 1,',
            '}', '', 'ccmc {', 'sys = sys,', 'qmc = {', 'tau = 0.1,',
            'rng_seed = 1472,', 'init_pop = 2,', 'mc_cycles = 10,',
            'nreports = 3,', 'target_population = 30,', 'state_size = 50000,',
            'spawned_state_size = 5000,', '},', 'reference = {',
            'ex_level = 3,', '},', '}', ''
            ]
        self.dummy_metadata['max_bloom'] = 457.0
        self.dummy_metadata['mean_bloom'] = 20.79
        self.dummy_metadata['nblooms'] = 227.0
        self.assertDictEqual(data[0][0], self.dummy_metadata)

    def test_basic_dmqmc_input(self):
        """Test DMQMC."""
        data = extract.extract_data("hande_files/dmqmc_ueg.out")
        dummy_data = pd.DataFrame([
            [0, 3.41598425e-01, 0.00000000e+00, 2.00000000e+02, 100, 48,
             2.49500000e-01, 6.00000000e-03],
            [2, 6.10389961e-01, 0.00000000e+00, 1.01000000e+02, 58, 33,
             2.59800000e-01, 0.00000000e+00],
            [4, 7.80323874e-01, 0.00000000e+00, 5.90000000e+01, 41, 26,
             2.74300000e-01, 0.00000000e+00],
            [0, 3.41598425e-01, 0.00000000e+00, 2.00000000e+02, 97, 45,
             2.51500000e-01, 6.00000000e-03],
            [2, 3.77511292e-01, 0.00000000e+00, 1.01000000e+02, 93, 46,
             2.41400000e-01, 0.00000000e+00],
            [4, 4.64434757e-01, 0.00000000e+00, 9.40000000e+01, 79, 36,
             2.50600000e-01, 0.00000000e+00]
            ], columns=[
                'iterations', 'Shift', 'Trace', r'# H psips',
                r'# states', r'# spawn_events', 'R_spawn', 'time'
                ])
        dummy_metadata = {
            'system': {
                'nbasis': 38, 'nel': 14, 'nvirt': 24, 'Ms': 0, 'nalpha': 7,
                'nbeta': 7, 'nvirt_alpha': 12, 'nvirt_beta': 12, 'nsym': 19,
                'sym0': 1, 'sym_max': 19, 'nsym_tot': 19, 'sym0_tot': 1,
                'sym_max_tot': 19, 'symmetry': 1, 'tot_sym': False,
                'aufbau_sym': True, 'max_number_excitations': 14, 'ueg': {
                    'r_s': 1.0, 'ecutoff': 1.0, 'k_fermi': 1.91915829,
                    'E_fermi': 1.84158428, 'ktwist': [0.0, 0.0, 0.0],
                    'L': [3.88512994, 3.88512994, 3.88512994]
                }},
            'qmc': {
                'rng_seed': 1472, 'real_amplitudes': False,
                'real_amplitude_force_32': False, 'spawn_cutoff': 0.01,
                'excit_gen': 'renorm', 'pattempt_update': False,
                'pattempt_zero_accum_data': False, 'pattempt_single': 0.0,
                'pattempt_double': 1.0, 'pattempt_parallel': 0.0, 'tau': 0.05,
                'tau_search': False, 'vary_shift_from': 0.0,
                'vary_shift_from_proje': False, 'initial_shift': 0.0,
                'shift_damping': 0.05, 'walker_length': 50000,
                'spawned_walker_length': 5000, 'D0_population': 200.0,
                'target_particles': 100.0, 'target_reference': False,
                'initiator_approx': False, 'initiator_pop': 3.0, 'ncycles': 2,
                'nreport': 2, 'power_pitzer_min_weight': 0.01,
                'quasi_newton': False, 'quasi_newton_threshold': 1e-05,
                'quasi_newton_value': 1.0, 'use_mpi_barriers': False
                },
            'dmqmc': {
                'beta_loops': 2, 'replica_tricks': False, 'start_av_rdm': 0,
                'weighted_sampling': False, 'vary_weights': False,
                'find_weights': False, 'find_weights_start': 0,
                'calc_excit_dist': False, 'all_sym_sectors': False,
                'all_spin_sectors': False, 'initiator_level': -1,
                'sampling_probs': '[]', 'finish_varying_weights': 0,
                'fermi_temperature': False, 'target_beta': 1.0,
                'mom_dist_kmax': 0.0, 'struc_fac_qmax': 0.0
                },
            'ipdmqmc': {
                'ipdmqmc': False, 'initial_matrix': 'hartree_fock',
                'grand_canonical_initialisation': False, 'symmetric': True,
                'chem_pot': 0.0, 'metropolis_attempts': 0
                },
            'rdm': {
                'nrdms': 0, 'spawned_length': 0, 'doing_rdm': False,
                'calc_ground_rdm': False, 'calc_inst_rdm': False,
                'doing_concurrence': False, 'doing_vn_entropy': False,
                'output_rdm': False
                },
            'operators': {
                'energy': False, 'energy_squared': False,
                'kinetic_energy': False, 'potential_energy': False,
                'H0_energy': False, 'HI_energy': False,
                'correlation_fn': False, 'staggered_mad_ind': False,
                'rdm_r2': False, 'full_r2': False, 'mom_dist': False
                },
            'restart': {
                'read_restart': False, 'read_id': 2147483647,
                'write_restart': False, 'write_id': 2147483647,
                'write_freq': 2147483647, 'write_restart_shift': False,
                'write_shift_id': 2147483647, 'restart_rng': True
                },
            'load balancing': {
                'nslots': 1, 'pop': 1000, 'percent': 0.05, 'max_attempts': 2,
                'write_info': False
                },
            'reference': {
                'det': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14],
                'det_ms': 0, 'det_symmetry': 1, 'H00': 13.60355734,
                'F0': 15.69278015, 'hilbert_space_det': [
                    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14
                    ], 'hilbert_space_det_ms': 0,
                'hilbert_space_det_symmetry': 1, 'ex_level': 14
                },
            'calc_type': 'DMQMC', 'input': [
                '', '-- Create output with:', '-- $[HANDE DIR]/bin/hande.x '
                'dmqmc_ueg.lua > dmqmc_ueg.out 2> dmqmc_ueg.err',
                '-- Note that these settings are just for testing...',
                'sys = ueg {', 'dim = 3,', 'nel = 14,', 'ms = 0,',
                'cutoff = 1,', '}', '', 'dmqmc {', 'sys = sys,', 'qmc = {',
                'tau = 0.05,', 'rng_seed = 1472,', 'init_pop = 200,',
                'mc_cycles = 2,', 'nreports = 2,', 'target_population = 100,',
                'state_size = 50000,', 'spawned_state_size = 5000,', '},',
                'dmqmc = {', 'beta_loops = 2,', '}', '}', ''
                ],
            'UUID': '8226e3f0-ee64-4128-9569-225da1f8b913', 'wall_time': 0.03,
            'cpu_time': 0.03, 'calculation_time': 0.03
            }
        pd.testing.assert_frame_equal(data[0][1], dummy_data)
        self.assertDictEqual(data[0][0], dummy_metadata)

    def test_basic_fci_input(self):
        """Test FCI."""
        data = extract.extract_data("hande_files/fci_ueg.out")
        dummy_data = pd.Series([
            -1.78882976e-02, 9.45177589e+00, 9.45177589e+00, 9.52511643e+00,
            9.52511643e+00, 9.52511643e+00, 9.88957877e+00, 1.90174844e+01,
            1.90174844e+01, 1.90174844e+01, 1.90320038e+01, 1.90320038e+01,
            1.90827874e+01, 1.90827874e+01, 1.90827874e+01, 1.92329356e+01,
            1.92329356e+01, 1.92329356e+01, 1.97091785e+01
            ], index=list(range(1, 20)), name='FCI (LAPACK)')
        dummy_data.index.name = 'Eigenvalue'
        dummy_metadata = {
            'system': {
                'nbasis': 38, 'nel': 2, 'nvirt': 36, 'Ms': 0, 'nalpha': 1,
                'nbeta': 1, 'nvirt_alpha': 18, 'nvirt_beta': 18, 'nsym': 19,
                'sym0': 1, 'sym_max': 19, 'nsym_tot': 19, 'sym0_tot': 1,
                'sym_max_tot': 19, 'symmetry': 1, 'tot_sym': False,
                'aufbau_sym': True, 'max_number_excitations': 2, 'ueg': {
                    'r_s': 1.0, 'ecutoff': 1.0, 'k_fermi': 1.91915829,
                    'E_fermi': 1.84158428, 'ktwist': [0.0, 0.0, 0.0],
                    'L': [2.0309826, 2.0309826, 2.0309826]
                }},
            'reference': {'ex_level': 2}, 'calc_type': 'FCI', 'input': [
                '', '-- Create output with:', '-- $[HANDE DIR]/bin/hande.x '
                'fci_ueg.lua > fci_ueg.out 2> fci_ueg.err',
                '-- Note that these settings are just for testing...',
                'sys = ueg {', 'dim = 3,', 'nel = 2,', 'ms = 0,',
                'cutoff = 1,', '}', '', 'fci {', 'sys = sys,', '}', ''
                ],
            'UUID': '7f1b774f-975b-4ec0-a626-185879d0e7a7', 'wall_time': 0.0,
            'cpu_time': 0.0, 'calculation_time': 0.0, 'fci_in': {
                'analyse_fci_wfn': 0, 'block_size': 64,
                'determinant_file': 'DETS', 'direct_lanczos': False,
                'hamiltonian_file': 'HAMIL', 'lanczos_string_len': 40,
                'nlanczos_eigv': 5, 'print_fci_wfn': 0,
                'print_fci_wfn_file': 'FCI_WFN', 'write_determinants': False,
                'write_hamiltonian': False
                }
            }
        pd.testing.assert_series_equal(data[0][1], dummy_data)
        self.assertDictEqual(data[0][0], dummy_metadata)

    def test_basic_simple_fciqmc_input(self):
        """Test Simple FCIQMC."""
        data = extract.extract_data("hande_files/simple_fciqmc_ueg.out")
        dummy_data = pd.DataFrame([
            [10, 0., -4.70181114e-02, 8.00000000e+00, 9.00000000e+00, 0, 0,
             5.00000000e-02, 0.0],
            [20, -8.78004297e-03, -2.29343899e-01, 8.00000000e+00,
             1.00000000e+01, 0, 0, 6.67000000e-02, 0.00000000e+00],
            [30, -0.04794701, -0.33774677, 8., 16., 0, 0, 0.08, 0.]
            ], columns=self.columns)
        del_keys = [
            'load balancing', 'blocking', 'fciqmc', 'logging', 'logging_in',
            'max_bloom', 'mean_bloom', 'nblooms', 'semi_stoch'
            ]
        for key in del_keys:
            del self.dummy_metadata[key]
        dummy_metadata = {
            'system': {
                'nbasis': 38, 'nel': 2, 'nvirt': 36, 'Ms': 0, 'nalpha': 1,
                'nbeta': 1, 'nvirt_alpha': 18, 'nvirt_beta': 18, 'nsym': 19,
                'sym0': 1, 'sym_max': 19, 'nsym_tot': 19, 'sym0_tot': 1,
                'sym_max_tot': 19, 'symmetry': 1, 'tot_sym': False,
                'aufbau_sym': True, 'max_number_excitations': 2,
                'ueg': {
                    'r_s': 1.0, 'ecutoff': 1.0, 'k_fermi': 1.91915829,
                    'E_fermi': 1.84158428, 'ktwist': [0.0, 0.0, 0.0],
                    'L': [2.0309826, 2.0309826, 2.0309826]
                }},
            'qmc': {
                'rng_seed': 1472, 'real_amplitudes': False,
                'real_amplitude_force_32': False, 'spawn_cutoff': 0.01,
                'excit_gen': 'renorm', 'pattempt_update': False,
                'pattempt_zero_accum_data': False, 'pattempt_single': -1.0,
                'pattempt_double': -1.0, 'pattempt_parallel': -1.0,
                'tau': 0.06, 'tau_search': False, 'vary_shift_from': 0.0,
                'vary_shift_from_proje': False, 'initial_shift': 0.0,
                'shift_damping': 1.79769313e+308, 'walker_length': 0,
                'spawned_walker_length': 0, 'D0_population': 8.0,
                'target_particles': 8.0, 'target_reference': False,
                'initiator_approx': False, 'initiator_pop': 3.0, 'ncycles': 10,
                'nreport': 3, 'power_pitzer_min_weight': 0.01,
                'quasi_newton': False, 'quasi_newton_threshold': 1e-05,
                'quasi_newton_value': 1.0, 'use_mpi_barriers': False
                },
            'restart': {
                'read_restart': False, 'read_id': 2147483647,
                'write_restart': False, 'write_id': 2147483647,
                'write_freq': 2147483647, 'write_restart_shift': False,
                'write_shift_id': 2147483647, 'restart_rng': True
                },
            'reference': {
                'det': [4, 13], 'det_ms': 0, 'det_symmetry': 1, 'H00': 9.57078,
                'F0': 0.0, 'ex_level': 2
                },
            'sparse_hamil': True, 'calc_type': 'Simple FCIQMC',
            'input': [
                '', '-- Create output with:', '-- $[HANDE DIR]/bin/hande.x '
                'simple_fciqmc_ueg.lua > simple_fciqmc_ueg.out 2> '
                'simple_fciqmc_ueg.err', '-- Note that these settings are '
                'just for testing...', 'sys = ueg {', 'dim = 3,', 'nel = 2,',
                'ms = 0,', 'cutoff = 1,', '}', '', 'simple_fciqmc {',
                'sys = sys,', 'sparse = true,', 'qmc = {', 'tau = 0.06,',
                'rng_seed = 1472,', 'init_pop = 8,', 'mc_cycles = 10,',
                'nreports = 3,', 'target_population = 8,',
                'state_size = 50000,', 'spawned_state_size = 5000,', '},', '}',
                ''
                ], 'UUID': 'c8132a6b-bc15-4586-adf7-80d7c3119e3c',
            'wall_time': 0.0, 'cpu_time': 0.0, 'calculation_time': 0.0
            }
        pd.testing.assert_frame_equal(data[0][1], dummy_data)
        self.assertDictEqual(data[0][0], dummy_metadata)

    def test_hilbert(self):
        """Test MC Hilbert space size estimation extraction."""
        data = extract.extract_data("hande_files/hilbert_ueg.out")
        dummy_data = pd.DataFrame([
            [1, 15233700.0, 15233700.0, np.nan],
            [2, 20311600.0, 17772650.0, 2538951.0],
            [3, 10155800.0, 15233700.0, 2931728.0]
            ], columns=['iterations', 'space size', 'mean', 'std. err.'])
        dummy_metadata = {
            'system': {
                'nbasis': 38, 'nel': 14, 'nvirt': 24, 'Ms': 0, 'nalpha': 7,
                'nbeta': 7, 'nvirt_alpha': 12, 'nvirt_beta': 12, 'nsym': 19,
                'sym0': 1, 'sym_max': 19, 'nsym_tot': 19, 'sym0_tot': 1,
                'sym_max_tot': 19, 'symmetry': 2147483647, 'tot_sym': False,
                'aufbau_sym': True, 'max_number_excitations': 14, 'ueg': {
                    'r_s': 1.0, 'ecutoff': 1.0, 'k_fermi': 1.91915829,
                    'E_fermi': 1.84158428, 'ktwist': [0.0, 0.0, 0.0],
                    'L': [3.88512994, 3.88512994, 3.88512994]
                }},
            'ex_level': 14, 'nattempts': 1000, 'ncycles': 3, 'occ_list': [],
            'rng_seed': -563090706, 'calc_type': 'Hilbert space',
            'input': ['', '-- Create output with:',
                      '-- $[HANDE DIR]/bin/hande.x hilbert_ueg.lua > '
                      'hilbert_ueg.out 2> hilbert_ueg.err',
                      '-- Note that these settings are just for testing...',
                      'sys = ueg {', 'dim = 3,', 'nel = 14,', 'ms = 0,',
                      'cutoff = 1,', '}', '', 'hilbert_space {', 'sys = sys,',
                      'hilbert = {', 'rng_seed = -563090706,',
                      'nattempts = 1000,', 'ncycles = 3,', '},', '}', ''
                      ],
            'UUID': '6c228dba-68c9-4051-b47a-5704fe261ad8', 'wall_time': 0.0,
            'cpu_time': 0.0, 'calculation_time': 0.0
            }
        pd.testing.assert_frame_equal(data[0][1], dummy_data)
        self.assertDictEqual(data[0][0], dummy_metadata)

    def test_canonical_estimates(self):
        """Test canonical estimates extraction."""
        data = extract.extract_data("hande_files/cano_ueg.out")
        dummy_data = pd.DataFrame([
            [1, 2.7943370015E+01, -1.5740506968E+00, 2.5170736424E+01,
             -1.4249954801E+00, 9.0253160272E-01, 3.4800000000E-02],
            [2, 2.7726172803E+01, -1.5817153286E+00, 2.5010826838E+01,
             -1.4341367788E+00, 9.0391912264E-01, 3.4700000000E-02]
            ], columns=[
                'iterations', '<T>_0', '<V>_0', r'Tr(T\rho_HF)',
                r'Tr(V\rho_HF)', r'Tr(\rho_HF)', 'N_ACC/N_ATT'
                ])
        dummy_metadata = {
            'system': {
                'nbasis': 38, 'nel': 14, 'nvirt': 24, 'Ms': 0, 'nalpha': 7,
                'nbeta': 7, 'nvirt_alpha': 12, 'nvirt_beta': 12, 'nsym': 19,
                'sym0': 1, 'sym_max': 19, 'nsym_tot': 19, 'sym0_tot': 1,
                'sym_max_tot': 19, 'symmetry': 2147483647, 'tot_sym': False,
                'aufbau_sym': True, 'max_number_excitations': 14, 'ueg': {
                    'r_s': 1.0, 'ecutoff': 1.0, 'k_fermi': 1.91915829,
                    'E_fermi': 1.84158428, 'ktwist': [0.0, 0.0, 0.0],
                    'L': [3.88512994, 3.88512994, 3.88512994]
                }}, 'all_spin_sectors': False, 'beta': 0.2,
            'fermi_temperature': False, 'nattempts': 10000,
            'free_energy_corr': -96.66701818, 'ncycles': 2,
            'chem_pot': -0.64446167, 'rng_seed': 748,
            'calc_type': 'Canonical energy',
            'input': [
                '', '-- Create output with:',
                '-- $[HANDE DIR]/bin/hande.x cano_ueg.lua > cano_ueg.out 2> '
                'cano_ueg.err', '-- Note that these settings are just for '
                'testing...', 'sys = ueg {', 'dim = 3,', 'nel = 14,',
                'ms = 0,', 'cutoff = 1,', '}', '', 'canonical_estimates {',
                'sys = sys,', 'canonical_estimates = {', 'ncycles = 2,',
                'nattempts = 10000,', 'beta = 0.2,', 'rng_seed = 748', '},',
                '}', ''
            ], 'UUID': '98db19fc-8ab2-4f5a-82d0-5c4c10340a1c',
            'wall_time': 0.0, 'cpu_time': 0.0, 'calculation_time': 0.0
            }
        pd.testing.assert_frame_equal(data[0][1], dummy_data)
        self.assertDictEqual(data[0][0], dummy_metadata)

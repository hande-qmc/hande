'''Obtain metadata for legacy output files (i.e. not containing JSON blocks).'''
import re

def extract_metadata(fh):
    '''Extract metadata from a legacy output file.

Parameters
----------
fh : file
    File handle to (open) file containing HANDE output, positioned at the
    beginning of the file (or at least at the start of the input section).

Returns
-------
metadata : dict
    Metadata extracted from analysing the (plain text) output and input options.
'''

    # metadata from ...
    # ... input (echoed in output)
    md_input = dict(
        sym = r'\bsym\b +\d+',
        ms = r'\bms\b +-*\d+',
        nel = r'\bnel\b|\belectrons\b',
        tau = r'\btau\b',
        truncation = 'truncation_level',
        target = 'varyshift_target',
        shift_damping = 'shift_damping',
        mc_cycles = 'mc_cycles',
    )
    # ... main body of output (ie after input but before QMC data table)
    md_body = dict(
        nbasis = 'Number of basis functions:',
        seed = 'random number generator with a seed of',
        bit_length = 'Bit-length',
        ref = 'Reference determinant, |D0> = ',
        ref_energy = r'E0 = <D0|H|D0>',
        psingle = 'Probability of attempting a single',
        pdouble = 'Probability of attempting a double',
        init_pop = 'Initial population on',
    )

    md_int = ['sym', 'ms', 'nel', 'nbasis', 'truncation', 'seed', 'mc_cycles',
         'bit_length', 'mc_cycles']
    md_float = ['tau', 'ref_energy', 'psingle', 'pdouble', 'init_pop',
         'shift_damping', 'target']

    metadata = {}

    input_pattern = 'Input options'
    underline_regex = re.compile('----+')

    for line in fh:
        if input_pattern in line:
            next(fh)
            break
    # Parse metadata from input block.
    for line in fh:
        if underline_regex.search(line):
            break
        for (key, pattern) in md_input.items():
            if re.search(pattern, line):
                val = line.split()[-1]
                if val[-1] == ',':
                    val = val[:-1]
                if key in md_int:
                    metadata[key] = int(float(val))
                elif key in md_float:
                    metadata[key] = float(val)
                else:
                    metadata[key] = val
    # Parse metadata from body.
    for line in fh:
        if '# iterations' in line:
            # Finished with body! (Only worry about stopping early if there's
            # potentially lots of output left...)
            break
        for (key, pattern) in md_body.items():
            if pattern in line:
                if key in md_int:
                    metadata[key] = int(float(line.split()[-1]))
                elif key in md_float:
                    metadata[key] = float(line.split()[-1])
                elif key == 'ref':
                    metadata[key] = ' '.join(line.split(pattern)).strip()
                else:
                    metadata[key] = line.split()[-1]

    return convert_metadata(metadata)

def convert_metadata(legacy_metadata):
    '''Convert metadata from original/legacy format to the JSON format.

Parameters
---------
legacy_metadata : dict
    Metadata in legacy format (i.e. a single dict).

Returns
-------
metadata : dict
    Metadata in JSON format (i.e. a nested dict).
'''

    # Handle changed names of metadata.
    if "mc_cycles" in legacy_metadata:
        legacy_metadata["ncycles"] = legacy_metadata["mc_cycles"]

    # Convert single dict to nested.
    json_keys = {
        "system": ["nbasis", "nel", "nvirt", "Ms", "nalpha", "nbeta",
            "nvirt_alpha", "nvirt_beta", "nsym", "sym0", "sym_max", "nsym_tot",
            "sym0_tot", "sym_max_tot", "symmetry", "max_number_excitations",
            "lattice", "read_in", "hubbard", "heisenberg", "ueg"],
        "qmc": ["rng_seed", "real_amplitudes", "real_amplitude_force_32",
            "spawn_cutoff", "excit_gen", "pattempt_single", "pattempt_double",
            "tau", "tau_search", "vary_shift_from", "vary_shift_from_proje",
            "initial_shift", "shift_damping", "walker_length",
            "spawned_walker_length", "D0_population", "target_particles",
            "initiator_approx", "initiator_pop", "ncycles", "nreport",
            "use_mpi_barriers"],
        "ccmc": ["move_freq", "cluster_multispawn_threshold", "full_nc",
            "linked"],
        "dmqmc": ["beta_loops", "replica_tricks", "start_av_rdm",
            "weighted_sampling", "vary_weights", "find_weights",
            "find_weights_start", "calc_excit_dist", "all_sym_sectors",
            "all_spin_sectors", "initiator_level"],
        "ipdmqmc": ["propagate_to_beta", "initial_matrix",
            "grand_canonical_initialisation", "symmetric", "chem_pot",
            "metropolis_attempts"],
        "operators": ["energy", "energy_squared", "kinetic_energy",
            "potential_energy", "H0_energy", "HI_energy", "correlation_fn",
            "staggered_magnetisation", "rdm_r2", "full_r2"],
        "rdm": ["nrdms", "spawned_length", "doing_rdm", "calc_ground_rdm",
            "calc_inst_rdm", "doing_concurrence", "doing_vn_entropy",
            "output_rdm"],
        "fciqmc": ["select_ref_det_every_nreports", "init_spin_inv_D0",
             "ref_det_factor", "non_blocking_comm", "doing_load_balancing",
             "trial_function", "guiding_function"],
        "semi_stoch": ["start_iter", "shift_iter", "space_type", "target_size",
             "write_determ_space", "projection_mode", "read_id", "write_id",
             "ci_space"],
        "restart": ["read_restart", "read_id", "write_restart", "write_id",
             "write_freq", "write_restart_shift", "write_shift_id"],
        "load balancing": ["nslots", "pop", "percent", "max_attempts",
             "write_info"],
        "reference": ["det", "det_ms", "det_symmetry", "H00",
             "hilbert_space_det", "hilbert_space_det_ms",
             "hilbert_space_det_symmetry", "ex_level"],
        "fci": ["write_hamiltonian", "hamiltonian_file", "write_determinants",
            "determinant_file", "print_fci_wfn", "print_fci_wfn_file",
            "analyse_fci_wfn", "block_size", "nlanczos_eigv",
            "lanczos_string_len", "direct_lanczos"],
         }

    metadata = {key:{} for key in json_keys.keys()}
    for main_key, keys in metadata:
        for key in keys:
            if key in legacy_metadata:
                metadata[main_key][key] = legacy_metadata[key]

    return legacy_metadata

def extract_input(metadata, variable):
    '''Extract input information about specific variable.

Parameters
----------
metadata : :class:`pandas.DataFrame`
    metadata (i.e. calculation information, parameters and settings) extracted
    from output files.
variable : string
    input variable whose value is desired.

Returns
-------
value : string
    value associated with variable in input file.
'''

    value = [x.split()[2].split(',')[0] for x in metadata['input'] if variable in x]
    return (value)

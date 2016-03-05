sys = ueg {
    nel = 7,
    ms = 7,
    sym = 1,
    dim = 3,
    cutoff = 10,
    rs = 1,
}

dmqmc {
    sys = sys,
    qmc = {
        tau = 0.001,
        rng_seed = 7,
        init_pop = 10000,
        mc_cycles = 10,
        nreports = 100,
        target_population = 10000,
        state_size = -200,
        spawned_state_size = -100,
    },
    dmqmc = {
        fermi_temperature = true,
        all_sym_sectors = true,
        beta_loops = 100,
    },
    ipdmqmc = {
        target_beta = 1.0,
        initial_matrix = 'free_electron',
        grand_canonical_initialisation = true,
        symmetric = false,
    },
    operators = {
        energy = true,
    },
}

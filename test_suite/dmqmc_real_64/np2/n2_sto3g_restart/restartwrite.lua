dmqmc {
    sys = read_in {
        int_file = 'FCIDUMP',
        nel = 10,
        ms = 0,
        sym = 0,
    },
    reference = {
        det = {1, 3, 5, 7, 11, 2, 4, 6, 8, 12},
    },
    qmc = {
        tau = 0.001,
        init_pop = 1E4,
        mc_cycles = 10,
        nreports = 1000,
        target_population = 1E4,
        state_size = -100,
        spawned_state_size = -150,
        real_amplitudes = true,
        rng_seed = 1972011788,
    },
    dmqmc = {
        beta_loops = 3,
        final_beta = 1.0,
        symmetric_bloch = false,
    },
    ipdmqmc = {
        target_beta = 0.5,
        grand_canonical_initialisation = true,
        symmetric_interaction_picture = false,
    },
    operators = {
        energy = true,
    },
    restart = {
        write = 0,
        write_frequency = 10,
    },
}

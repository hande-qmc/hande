dmqmc {
    sys = read_in {
        int_file = 'FCIDUMP',
        nel = 6,
        ms = 0,
        sym = 0,
    },
    reference = {
        det = {1, 2, 3, 4, 5, 6},
    },
    qmc = {
        tau = 0.001,
        init_pop = 2E4,
        mc_cycles = 10,
        nreports = 1,
        target_population = 2E4,
        state_size = -100,
        spawned_state_size = -150,
        real_amplitudes = true,
        rng_seed = 537858000,
    },
    dmqmc = {
        beta_loops = 3,
        final_beta = 2.0,
    },
    operators = {
        energy = true,
    },
    restart = {
        write = 0,
        write_frequency = 20,
    },
}

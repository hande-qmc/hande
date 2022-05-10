dmqmc {
    sys = ueg {
        ms = 0,
        nel = 4,
        cutoff = 2.5,
        rs = 5.0,
        dim = 3,
    },
    qmc = {
        tau = 0.001,
        init_pop = 1E5,
        mc_cycles = 10,
        nreports = 250,
        target_population = 1E2,
        state_size = -400,
        spawned_state_size = -500,
        real_amplitudes = true,
        rng_seed = 117,
        state_histograms = true,
    },
    dmqmc = {
        beta_loops = 2,
    },
    operators = {
        energy = true,
    },
}

fciqmc {
    sys = ueg {
        ms = 0,
        nel = 4,
        cutoff = 2.5,
        rs = 5.0,
        dim = 3,
    },
    qmc = {
        tau = 0.001,
        init_pop = 5E3,
        mc_cycles = 10,
        nreports = 2500,
        target_population = 5E3,
        state_size = -400,
        spawned_state_size = -500,
        real_amplitudes = true,
        rng_seed = 7,
        state_histograms = true,
        state_histograms_nreport = 2250,
    },
}

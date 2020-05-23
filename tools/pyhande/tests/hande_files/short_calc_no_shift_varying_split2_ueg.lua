-- Note that these settings are just for testing...
sys = ueg {
    dim = 3,
    nel = 14,
    ms = 0,
    cutoff = 1,
}

fciqmc {
    sys = sys,
    qmc = {
        tau = 0.001,
        rng_seed = 1472,
        init_pop = 2,
        mc_cycles = 10,
        nreports = 2,
        target_population = 6000,
        state_size = 50000,
        spawned_state_size = 5000,
    },
    restart = {
        read = true,
    },
}

sys = read_in {
    int_file = "INTDUMP",
    nel = 4,
    ms = 0,
    sym = 0,
}

fciqmc {
    sys = sys,
    qmc = {
        tau = 0.01,
        rng_seed = 7,
        init_pop = 100,
        mc_cycles = 10,
        nreports = 900,
        target_population = 1000,
        state_size = -5,
        spawned_state_size = -1,
        real_amplitudes = false,
        spawn_cutoff = 0.1,
    },
    restart = {
        write = 0,
    },
}

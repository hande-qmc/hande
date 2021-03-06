sys = read_in {
    int_file = "INTDUMP",
    nel = 10,
    ms = 0,
    sym = 0,
}

ccmc {
    sys = sys,
    qmc = {
        tau = 0.01,
        init_pop = 1000,
        rng_seed = 1660032958,
        mc_cycles = 10,
        nreports = 175,
        target_population = 50000,
        state_size = -100,
        spawned_state_size = -50,
    },
    reference = {
        ex_level = 2,
    },
    restart = { write = 0, },
}

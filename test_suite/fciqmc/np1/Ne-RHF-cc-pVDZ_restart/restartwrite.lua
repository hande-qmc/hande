sys = read_in {
    int_file = "INTDUMP",
    nel = 10,
    ms = 0,
    sym = 0,
}

fciqmc {
    sys = sys,
    qmc = {
        tau = 0.005,
        rng_seed = 5691,
        init_pop = 10,
        mc_cycles = 10,
        nreports = 1000,
        target_population = 20000,
        state_size = -100,
        spawned_state_size = -50,
    },
    restart = { write= 0}
}

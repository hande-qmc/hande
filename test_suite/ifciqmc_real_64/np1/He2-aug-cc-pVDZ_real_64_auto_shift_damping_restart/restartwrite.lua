sys = read_in {
    int_file = "INTDUMP",
    nel = 4,
    ms = 0,
    sym = 0,
}

fciqmc {
    sys = sys,
    qmc = {
        tau = 0.05,
        rng_seed = 7,
        mc_cycles = 5,
        nreports = 300,
        target_population = 10000,
        state_size = 20000,
        spawned_state_size = 10000,
        real_amplitudes = true,
        spawn_cutoff = 0.1,
        init_pop = 40,
        initiator = true,
        shift_damping = 0.1,
    },
    blocking = {
        blocking_on_the_fly = true,
        auto_shift_damping = true,
        shift_damping_precision = 1.5,
    },
    restart = {
        write = 0,
    },    
}

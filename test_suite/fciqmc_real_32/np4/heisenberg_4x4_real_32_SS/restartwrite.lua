sys = heisenberg {
    lattice = { {4, 0}, {0, 4} },
    ms = 0,
    J = -0.0625,
}

fciqmc {
    sys = sys,
    qmc = {
        tau = 0.05,
        rng_seed = 7,
        init_pop = 10,
        mc_cycles = 1,
        nreports = 2500,
        target_population = 10000,
        real_amplitudes = true,
        real_amplitude_force_32 = true,
        spawn_cutoff = 0.01,
        state_size = 40000,
        spawned_state_size = 40000,
    },
    semi_stoch = {
        separate_annihilation = false,
        start_iteration = 4000,
        space = "read",
    },
    restart = {
        write = 0,
    },
}

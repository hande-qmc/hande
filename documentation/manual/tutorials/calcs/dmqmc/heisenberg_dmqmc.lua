sys = heisenberg {
    lattice = {
        {6, 0},
        {0, 6},
    },
    J = -1.0,
    ms = 0,
}

dmqmc {
    sys = sys,
    qmc = {
        tau = 0.001,
        init_pop = 10^6,
        rng_seed = 19838,
        mc_cycles = 10,
        nreports = 400,
        target_population = 10^6,
        state_size = -400,
        spawned_state_size = -400,
    },
    dmqmc = {
        beta_loops = 100,
    },
    operators = {
        energy = true,
    },
}

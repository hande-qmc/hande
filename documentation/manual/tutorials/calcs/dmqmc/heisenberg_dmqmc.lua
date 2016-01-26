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
        tau = 0.0025,
        init_pop = 10000,
        mc_cycles = 10,
        nreports = 40,
        target_population = 10000,
        state_size = -100,
        spawned_state_size = -100,
    },
    dmqmc = {
        beta_loops = 50,
    },
    operators = {
        energy = true,
        staggered_magnetisation = true,
    },
}

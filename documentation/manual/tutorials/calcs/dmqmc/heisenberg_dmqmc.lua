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
        init_pop = 5e6,
        rng_seed = 19838,
        mc_cycles = 10,
        nreports = 500,
        shift_damping = 0.5,
        target_population = 5e6,
        state_size = -400,
        spawned_state_size = -400,
    },
    dmqmc = {
        beta_loops = 1,
    },
    operators = {
        energy = true,
        excit_dist = true,
    },
}

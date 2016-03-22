sys = heisenberg {
    lattice = { {6, 0}, {0, 6} },
    ms = 0,
    J = -1,
}

qmc_table = {
    tau = 0.001,
    init_pop = 5e6,
    mc_cycles = 10,
    nreports = 500,
    target_population = 5e6,
    shift_damping = 0.1,
    initial_shift = -1.0,
    state_size = -400,
    spawned_state_size = -400,
}

qmc_state, weights = dmqmc {
    sys = sys,
    qmc = qmc_table,
    dmqmc = {
        beta_loops = 10,
        find_weights = true,
        find_weights_start = 3000,
    },
    operators = {
        energy = true,
    }
}

qmc_state:free()

dmqmc {
    sys = sys,
    qmc = qmc_table,
    dmqmc = {
        beta_loops = 100,
        sampling_weights = weights,
    },
    operators = {
        energy = true,
    }
}

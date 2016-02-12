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
        mc_cycles = 10,
        nreports = 1000,
        shift_damping = 0.5,
        target_population = 5e6,
        state_size = -400,
        spawned_state_size = -400,
    },
    dmqmc = {
        beta_loops = 1,
        sampling_weights = {1.9001E+02,4.8842E+01,2.0936E+01,1.1026E+01,6.4820E+00,4.0698E+00,2.6569E+00,1.7802E+00,1.2086E+00,8.2743E-01,5.6173E-01,3.7639E-01,2.4571E-01,1.5427E-01,9.0698E-02,4.7764E-02,2.0474E-02,5.2630E-03},
    },
    operators = {
        energy = true,
        excit_dist = true,
    },
}

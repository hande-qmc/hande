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

dmqmc {
    sys = sys,
    qmc = qmc_table,
    dmqmc = {
        beta_loops = 1,
        sampling_weights = {  1.7903E+02,    4.7470E+01,    2.0594E+01,    1.0904E+01,    6.4310E+00,    4.0442E+00,    2.6456E+00,    1.7707E+00,    1.2021E+00,    8.3191E-01,    5.6474E-01,    3.7799E-01,    2.4727E-01,    1.5550E-01,    9.1708E-02,    4.8557E-02,    2.1066E-02,    5.5856E-03,  },
    },
    operators = {
        energy = true,
        excit_dist = true,
    }
}

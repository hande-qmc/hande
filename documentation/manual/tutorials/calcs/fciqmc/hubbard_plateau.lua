hubbard = hubbard_k {
    lattice = {
        { 3,  3 },
        { 3, -3 },
    },
    electrons = 18,
    ms = 0,
    U = 1.3,
    t = 1,
    sym = 0,
}

fciqmc {
    sys = hubbard,
    qmc = {
        tau = 0.002,
        mc_cycles = 20,
        nreports = 500,
        init_pop = 100,
        target_population = 10^10,
        state_size = -1000,
        spawned_state_size = -100,
    },
}

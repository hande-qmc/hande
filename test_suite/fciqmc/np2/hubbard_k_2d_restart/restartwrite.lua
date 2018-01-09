sys = hubbard_k {
    electrons = 18,
    ms = 0,
    sym = 'tot_sym',
    U = 0.5,
    lattice = { {3,3},{3,-3} },
}

fciqmc {
    sys = sys,
    qmc = {
        tau = 0.01,
        rng_seed = 7,
        init_pop = 10,
        mc_cycles = 25,
        nreports = 1000,
        target_population = 5500,
        state_size = -500,
        spawned_state_size = -200,
    },
    restart = { write = 0 },
}

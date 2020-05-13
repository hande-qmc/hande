-- Create output with:
-- $[HANDE DIR]/bin/hande.x 0.2_ccsd.lua > 0.2_ccsd.out 2> 0.2_ccsd.err
-- Note that these settings are just for testing...
sys = ueg {
    dim = 3,
    nel = 14,
    ms = 0,
    cutoff = 1,
}

ccmc {
    sys = sys,
    qmc = {
        tau = 0.01,
        rng_seed = 438,
        init_pop = 10,
        mc_cycles = 10,
        nreports = 1500,
        target_population = 7000,
        state_size = 50000,
        spawned_state_size = 5000,
    },
    reference = {
        ex_level = 2,
    },
}

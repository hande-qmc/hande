-- Create output with:
-- $[HANDE DIR]/bin/hande.x simple_fciqmc_ueg.lua > simple_fciqmc_ueg.out 2> simple_fciqmc_ueg.err
-- Note that these settings are just for testing...
sys = ueg {
    dim = 3,
    nel = 2,
    ms = 0,
    cutoff = 1,
}

simple_fciqmc {
    sys = sys,
    sparse = true,
    qmc = {
        tau = 0.06,
        rng_seed = 1472,
        init_pop = 8,
        mc_cycles = 10,
        nreports = 3,
        target_population = 8,
        state_size = 50000,
        spawned_state_size = 5000,
    },
}

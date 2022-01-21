-- Create output with:
-- $[HANDE DIR]/bin/hande.x long_calc_ueg.lua > long_calc_ueg.out 2> long_calc_ueg.err
-- Note that these settings are just for testing...
sys = ueg {
    dim = 3,
    nel = 14,
    ms = 0,
    cutoff = 1,
}

fciqmc {
    sys = sys,
    qmc = {
        tau = 0.01,
        rng_seed = 438,
        init_pop = 10,
        mc_cycles = 10,
        nreports = 3000,
        target_population = 7000,
        state_size = 50000,
        spawned_state_size = 5000,
    },
}

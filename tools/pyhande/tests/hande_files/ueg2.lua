-- Create output with:
-- $[HANDE DIR]/bin/hande.x ueg2.lua > ueg2.out 2> ueg2.err
-- Note that these settings are just for testing...
sys = ueg {
    dim = 3,
    nel = 14,
    ms = 0,
    cutoff = 1.0,
}

fciqmc {
    sys = sys,
    qmc = {
        tau = 0.03,
        rng_seed = 389,
        init_pop = 2,
        mc_cycles = 10,
        nreports = 3,
        target_population = 60,
        state_size = 50000,
        spawned_state_size = 5000,
        real_amplitudes = true,
    },
}

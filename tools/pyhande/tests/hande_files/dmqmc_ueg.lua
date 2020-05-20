-- Create output with:
-- $[HANDE DIR]/bin/hande.x dmqmc_ueg.lua > dmqmc_ueg.out 2> dmqmc_ueg.err
-- Note that these settings are just for testing...
sys = ueg {
    dim = 3,
    nel = 14,
    ms = 0,
    cutoff = 1,
}

dmqmc {
    sys = sys,
    qmc = {
        tau = 0.05,
        rng_seed = 1472,
        init_pop = 200,
        mc_cycles = 2,
        nreports = 2,
        target_population = 100,
        state_size = 50000,
        spawned_state_size = 5000,
    },
    dmqmc = {
        beta_loops = 2,
    }
}

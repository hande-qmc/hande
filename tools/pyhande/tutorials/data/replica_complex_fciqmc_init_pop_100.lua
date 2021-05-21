--Adapted from test_suite/fciqmc/np2/polyyne_complex_replica.
-- Create output with:
-- $[HANDE DIR]/bin/hande.x replica_complex_fciqmc_init_pop_100.lua > replica_complex_fciqmc_init_pop_100.out
-- Note that these settings are just for demonstration...
sys = read_in{
    int_file = "../../../../test_suite/fciqmc/np2/polyyne_complex_replica/FCIDUMP",
    complex = true,
    CAS = {8,8},
}

fciqmc {
    sys = sys,
    qmc = {
        tau = 0.001,
        rng_seed = 793,
        init_pop = 100,
        mc_cycles = 10,
        nreports = 7000,
        target_population = 2500,
        state_size = 600000,
        spawned_state_size = 400000,
        shift_damping = 0.01,
    },
    fciqmc = {
        replica_tricks = true,
    },
}

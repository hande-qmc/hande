-- Create output with:
-- $[HANDE DIR]/bin/hande.x cano_ueg.lua > cano_ueg.out 2> cano_ueg.err
-- Note that these settings are just for testing...
sys = ueg {
    dim = 3,
    nel = 14,
    ms = 0,
    cutoff = 1,
}

canonical_estimates {
    sys = sys,
    canonical_estimates = {
        ncycles = 2,
        nattempts = 10000,
        beta = 0.2,
        rng_seed = 748
    },
}

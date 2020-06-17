-- Create output with:
-- $[HANDE DIR]/bin/hande.x hilbert_ueg.lua > hilbert_ueg.out 2> hilbert_ueg.err
-- Note that these settings are just for testing...
sys = ueg {
    dim = 3,
    nel = 14,
    ms = 0,
    cutoff = 1,
}

hilbert_space {
    sys = sys,
    hilbert = {
        rng_seed = -563090706,
        nattempts = 1000,
        ncycles = 3,
    },
}

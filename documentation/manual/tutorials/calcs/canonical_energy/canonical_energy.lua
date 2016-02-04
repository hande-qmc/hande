sys = ueg {
    nel = 7,
    ms = 7,
    dim = 3,
    cutoff = 10,
    rs = 1,
    chem_pot = 2.91602409282,
}

canonical_energy {
    sys = sys,
    canonical_energy = {
        beta = 1,
        nattempts = 10000,
        ncycles = 100,
        rng_seed = 7,
        fermi_temperature = true,
    },
}

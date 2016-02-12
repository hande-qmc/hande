sys = ueg {
    nel = 7,
    ms = 7,
    dim = 3,
    cutoff = 10,
    rs = 1,
}

canonical_energy {
    sys = sys,
    canonical_energy = {
        beta = 1,
        nattempts = 10000,
        ncycles = 1000,
        fermi_temperature = true,
    },
}

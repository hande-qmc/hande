hubbard = hubbard_k {
    lattice = {
        { 3,  3 },
        { 3, -3 },
    },
    electrons = 18,
    ms = 0,
    U = 1.3,
    t = 1,
    sym = 0,
}

hilbert_space {
    sys = hubbard,
    hilbert = {
        nattempts = 100000,
        ncycles = 30,
    }
}

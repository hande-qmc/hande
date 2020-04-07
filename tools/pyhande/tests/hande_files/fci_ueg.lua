-- Create output with:
-- $[HANDE DIR]/bin/hande.x fci_ueg.lua > fci_ueg.out 2> fci_ueg.err
-- Note that these settings are just for testing...
sys = ueg {
    dim = 3,
    nel = 2,
    ms = 0,
    cutoff = 1,
}

fci {
    sys = sys,
}

sys = read_in {
    int_file = "CO.CCPVDZ.FCIDUMP",
    nel = 14,
    ms = 0,
}

ccmc {
    sys = sys,
    qmc = {
        tau = 1e-3,
        mc_cycles = 10,
        nreports = 1e5,
        state_size = -500,
        spawned_state_size = -200,
        init_pop = 1e4,
        target_population = 1e6,
	real_amplitudes = true,
    },
    reference = {
        ex_level = 3,
    },
}

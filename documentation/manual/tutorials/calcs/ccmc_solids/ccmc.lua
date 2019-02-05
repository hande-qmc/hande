sys = read_in {
    int_file = "fcidumpfile",
    complex = true,
    ex_int_file = "fcidumpfile_X",
}
ex_l=2

ccmc {
    sys = sys,
    qmc = {
        tau = 0.02,
        rng_seed = 13086,
        mc_cycles = 10,
        init_pop = 200,
        nreports = 50000,
        target_population = 1e4,
        state_size = -800,
        spawned_state_size = -500,
        excit_gen = "heat_bath",
        real_amplitudes = true,
    },
    ccmc = {
	    even_selection = true,
            full_non_composite=true,
           },
    reference = {
        ex_level = ex_l,
    },
    blocking = {
        blocking_on_the_fly = true,
        auto_shift_damping = true,
    },
    restart = {
        write = true,
    },
-- [review] - AJWT: Can these comments be removed?
   -- logging = {
   --     kpoint = 2,
  --  },
}


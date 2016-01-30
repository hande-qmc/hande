hubbard = hubbard_k {
    lattice = {
        { 3,  3 },
        { 3, -3 },
    },
    electrons = 18,
    ms = 0,
    U = 1.3,
    t = 1,
    sym = 1,
}

targets = {2.5*10^3, 5*10^3, 7.5*10^3, 10^4, 2.5*10^4, 5*10^4, 1*10^5}
for i,target in ipairs(targets) do
    qmc_state = fciqmc {
        sys = hubbard,
        qmc = {
            tau = 0.002,
            mc_cycles = 20,
            nreports = 10000,
            init_pop = 100,
            target_population = target,
            state_size = -1000,
            spawned_state_size = -100,
            initiator = true,
            real_amplitudes = true,
        },
    }
    -- For memory efficiency, explicitly free qmc_state after each calculation.
    qmc_state:free()
end

-- Function to generate a random twist.
function get_twist()
    ks = {}
    -- 3D UEG.
    for i = 0, 2 do
        -- For the UEG, we only need to generate a twist vectors whose components lie in
        -- the range [-pi/L, pi/L). In HANDE we interpret the input ks as being in terms
        -- of 2pi/L, so we need to randomly pick components in the range [-0.5, 0.5).
        sign = math.pow(-1, math.random(0, 1))
        ks[i] = 0.5*sign*math.random()
    end
    return ks
end

-- The number of simulations to average over.
ntwists = 3000
math.randomseed( os.time() )

for i = 1, ntwists do
    ks = get_twist()
    sys = ueg {
        nel = 19,
        ms = 19,
        dim = 3,
        cutoff = 20,
        rs = 0.5,
        twist = ks,
        verbose = false,
    }
    mc_state = canonical_estimates {
        sys = sys,
        canonical_estimates = {
            beta = 16,
            nattempts = 10000,
            ncycles = 10,
            fermi_temperature = true,
        },
    }
    sys:free() -- Free up memory.
end

nel = 8
ms = nel

fci {
    sys = ueg {
        nel = nel,
        ms = ms,
        sym = isym,
        dim = 3,
        cutoff = 1.01,
        rs = 1.0,
        use_mom_vec = true,
    },
    fci = {
        determinant_table_only = true,
    },
}

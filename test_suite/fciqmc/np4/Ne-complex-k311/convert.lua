sys = read_in {
    int_file = "fcidumpfile",
    complex = true,
    ex_int_file = "fcidumpfile_X",
    CAS = {6, 7},
}
hdf5_file = write_read_in_system {
    sys = sys,
    }

#include "dSFMT.h"

// Wrap around the required dSFMT functions so that they're accessible from
// fortran.  We only expose functions as needed.

// We can't bind directly to dSFMT functions as they either:
//
// a) are inline functions and I am reluctant to change the dSFMT source code
//    directly
// b) require a dsfmt_t state and dsfmt_t contains a union struct which is not
//    interoperable.  I did try to make a compatible derived type but had
//    issues with different data alignment in Fortran and C.
//
// There is an overhead associated with calling a wrapping routine, but this is
// minimal if an array of random numbers is produced rather than obtaining one
// random number at a time.

// extern C as we already use a C++ compiler (which we need for the hashing
// routines) and don't want to use a C compiler as well.
extern "C"
{
    void global_init_gen_rand(uint32_t seed)
    {
        // Initialise random number generator.
        // See also the main dSFMT code for the ability to initialise using an
        // array of seeds.
        dsfmt_gv_init_gen_rand(seed);
    }

    double global_genrand_close_open(void)
    {
        // Return a random number in the interval [0,1).
        //
        // This uses the global state function.
        // dSFMT also has the ability to have "local" states---useful for
        // threading?
        return dsfmt_gv_genrand_close_open();
    }

    void global_fill_array_close_open(double array[], int size)
    {
        // Fill an array of length size with random numbers in the interval
        // [0,1).
        // This is much faster than than repeated calls to genrand_close_open_.
        //
        // This uses the global state function.
        // dSFMT also has the ability to have "local" states---useful for
        // threading?
        dsfmt_gv_fill_array_close_open(array, size);
    }

}

#include "dSFMT.h"
#include <stdlib.h>

/* copyright James Spencer 2012.
 * New BSD License, see License.txt for details.
 */

/* Utility (memory-access) functions to enable use of dSFMT from Fortran. */

void* malloc_dsfmt_t(void) {
    /* Allocate sufficient memory for a dSFMT state (ie a variable of type dsfmt_t). */
    return malloc(sizeof(dsfmt_t));
}

void free_dsfmt_t(dsfmt_t* ptr) {
    /* Free memory associated with a dSFMT state (ie a variable of type dsfmt_t). */
    free(ptr);
}

// Wrapper to detect if an error occurred. This allows for simple error checking in
// Fortran rather than calling c_associated on a C pointer, which can trigger a bug in
// gfortran 7.1 (see https://gcc.gnu.org/bugzilla/show_bug.cgi?id=82869).
char *dsfmt_str_to_state_wrapper(dsfmt_t *dsfmt, char *str, char *prefix, int32_t *error) {
    char* err_str = dsfmt_str_to_state(dsfmt, str, prefix);
    *error = 0;
    if (err_str)
        *error = 1;
    return err_str;
}

#include <stdint.h> 
#include <iostream> 
#include "../src/cdefs.h"

using namespace std;

// Implement the FNV-1 and FNV-1a hash for 32-bit and 64-bit integers.
// See http://www.isthe.com/chongo/tech/comp/fnv/index.html.
// A consistent interface is provided to both:
//
//   fnv1_hash_(const void *key, int &len);
//   fnv1a_hash_(const void *key, int &len);
//
// In:
//    key: integer array of kind i0 to be hashed.
//    len: number of elements in key.
// Returns:
//    hash: hash of array key using the FNV-1 or FNV-1a algorithm.
//    The returned hash is of the same type as the key (i.e. 32-bit or 64-bit
//    integer).
//
// Similarly there are procedures which are always available for 32-bit
// integers, no matter how DET_SIZE is defined:
//
//   fnv1_hash32_(const void *key, int &len);
//   fnv1a_hash32_(const void *key, int &len);
//
// Neither key nor length are altered in the hashing.

extern "C"
{

#if DET_SIZE == 32

    uint32_t fnv1_hash_(const void *key, int &len)
    {
        const unsigned char *p = (const unsigned char *)key;
        uint32_t h = 2166136261;
        int i;

        for (i = 0; i < 4*len; i++)
        {
            h = ( h * 16777619 ) ^ p[i];
        }

        return h;
    }

#elif DET_SIZE == 64

    uint64_t fnv1_hash_(const void *key, int &len)
    {
        const unsigned char *p = (const unsigned char *)key;
        uint64_t h = 14695981039346656037;
        int i;

        for (i = 0; i < 8*len; i++)
        {
            h = ( h * 1099511628211) ^ p[i];
        }

        return h;
    }

#else

    unsigned int fnv1_hash( const void * key, int &len)
    {
        cout << "fnv1_hash only works with 32 or 64 bit integers." << endl;
        exit(1);
    }

#endif

#if DET_SIZE == 32

    uint32_t fnv1a_hash_(const void *key, int &len)
    {
        const unsigned char *p = (const unsigned char *)key;
        uint32_t h = 2166136261;
        int i;

        for (i = 0; i < 4*len; i++)
        {
            h = ( h ^ p[i] ) * 16777619 ;
        }

        return h;
    }

#elif DET_SIZE == 64

    uint64_t fnv1a_hash_(const void *key, int &len)
    {
        const unsigned char *p = (const unsigned char *)key;
        uint64_t h = 14695981039346656037;
        int i;

        for (i = 0; i < 8*len; i++)
        {
            h = ( h ^ p[i] ) * 1099511628211;
        }

        return h;
    }

#else

    unsigned int fnv1a_hash( const void * key, int &len)
    {
        cout << "fnv1_hash only works with 32 or 64 bit integers." << endl;
        exit(1);
    }

#endif

    uint32_t fnv1_hash32_(const void *key, int &len)
    {
        const unsigned char *p = (const unsigned char *)key;
        uint32_t h = 2166136261;
        int i;

        for (i = 0; i < 4*len; i++)
        {
            h = ( h * 16777619 ) ^ p[i];
        }

        return h;
    }

    uint32_t fnv1a_hash32_(const void *key, int &len)
    {
        const unsigned char *p = (const unsigned char *)key;
        uint32_t h = 2166136261;
        int i;

        for (i = 0; i < 4*len; i++)
        {
            h = ( h ^ p[i] ) * 16777619 ;
        }

        return h;
    }

} // extern "C"

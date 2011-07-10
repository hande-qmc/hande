#include <stdint.h> 
#include <iostream>
#include "../../src/cdefs.h"

using namespace std;

// murmurhash2 has a consistent interface:
//  murmurhash2( const void * key, int &len, unsigned int &seed )
// note that len is destroyed by the process.
//
// murmurhash2 is only defined for 32 and 64 bit integers.

extern "C"
{

    #if DET_SIZE == 32

    // 32-bit integers used as bit strings for storing determinants.

    //-----------------------------------------------------------------------------
    // murmurhash2, by Austin Appleby

    // Note - This code makes a few assumptions about how your machine behaves -

    // 1. We can read a 4-byte value from any address without crashing
    // 2. sizeof(int) == 4

    // And it has a few limitations -

    // 1. It will not work incrementally.
    // 2. It will not produce the same results on little-endian and big-endian
    //    machines.

    unsigned int murmurhash2( const void * key, int &len, unsigned int &seed )
    {
        // 'm' and 'r' are mixing constants generated offline.
        // They're not really 'magic', they just happen to work well.

        const unsigned int m = 0x5bd1e995;
        const int r = 24;

        // Initialize the hash to a 'random' value

        unsigned int h = seed ^ len;

        // Mix 4 bytes at a time into the hash

        const unsigned char * data = (const unsigned char *)key;

        while(len >= 4)
        {
            unsigned int k = *(unsigned int *)data;

            k *= m; 
            k ^= k >> r; 
            k *= m; 
            
            h *= m; 
            h ^= k;

            data += 4;
            len -= 4;
        }
        
        // Handle the last few bytes of the input array

        switch(len)
        {
        case 3: h ^= data[2] << 16;
        case 2: h ^= data[1] << 8;
        case 1: h ^= data[0];
                h *= m;
        };

        // Do a few final mixes of the hash to ensure the last few
        // bytes are well-incorporated.

        h ^= h >> 13;
        h *= m;
        h ^= h >> 15;

        return h;
    } 

    #elif DET_SIZE == 64

    // 64-bit integers used as bit strings for storing determinants.

    //-----------------------------------------------------------------------------
    // murmurhash2, 64-bit versions, by Austin Appleby

    // The same caveats as 32-bit murmurhash2 apply here - beware of alignment 
    // and endian-ness issues if used across multiple platforms.

    #ifdef BIT32

    // 64-bit hash for 32-bit platforms

    uint64_t murmurhash2( const void * key, int &len, unsigned int &seed )
    {
        const unsigned int m = 0x5bd1e995;
        const int r = 24;

        unsigned int h1 = seed ^ len;
        unsigned int h2 = 0;

        const unsigned int * data = (const unsigned int *)key;

        while(len >= 8)
        {
            unsigned int k1 = *data++;
            k1 *= m; k1 ^= k1 >> r; k1 *= m;
            h1 *= m; h1 ^= k1;
            len -= 4;

            unsigned int k2 = *data++;
            k2 *= m; k2 ^= k2 >> r; k2 *= m;
            h2 *= m; h2 ^= k2;
            len -= 4;
        }

        if(len >= 4)
        {
            unsigned int k1 = *data++;
            k1 *= m; k1 ^= k1 >> r; k1 *= m;
            h1 *= m; h1 ^= k1;
            len -= 4;
        }

        switch(len)
        {
        case 3: h2 ^= ((unsigned char*)data)[2] << 16;
        case 2: h2 ^= ((unsigned char*)data)[1] << 8;
        case 1: h2 ^= ((unsigned char*)data)[0];
                h2 *= m;
        };

        h1 ^= h2 >> 18; h1 *= m;
        h2 ^= h1 >> 22; h2 *= m;
        h1 ^= h2 >> 17; h1 *= m;
        h2 ^= h1 >> 19; h2 *= m;

        uint64_t h = h1;

        h = (h << 32) | h2;

        return h;
    } 

    #else

    // 64-bit hash for 64-bit platforms

    uint64_t murmurhash2( const void * key, int &len, unsigned int &seed )
    {
        const uint64_t m = 0xc6a4a7935bd1e995;
        const int r = 47;

        uint64_t h = seed ^ (len * m);

        const uint64_t * data = (const uint64_t *)key;
        const uint64_t * end = data + (len/8);

        while(data != end)
        {
            uint64_t k = *data++;

            k *= m; 
            k ^= k >> r; 
            k *= m; 
            
            h ^= k;
            h *= m; 
        }

        const unsigned char * data2 = (const unsigned char*)data;

        switch(len & 7)
        {
        case 7: h ^= uint64_t(data2[6]) << 48;
        case 6: h ^= uint64_t(data2[5]) << 40;
        case 5: h ^= uint64_t(data2[4]) << 32;
        case 4: h ^= uint64_t(data2[3]) << 24;
        case 3: h ^= uint64_t(data2[2]) << 16;
        case 2: h ^= uint64_t(data2[1]) << 8;
        case 1: h ^= uint64_t(data2[0]);
                h *= m;
        };
     
        h ^= h >> r;
        h *= m;
        h ^= h >> r;

        return h;
    } 

    #endif // 32-bit or 64-bit platform using 64-bit integers for bit strings.

    #else

    unsigned int murmurhash2( const void * key, int len, unsigned int seed )
    {
        cout << "MurmurHash only works with 32 or 64 bit integers." << endl;
        exit(1);
    }

    #endif

} // extern "C"

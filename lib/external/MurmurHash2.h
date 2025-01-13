/* ANSI C port from https://github.com/abrandoned/murmur2 */

#ifndef MURMUR_PLATFORM_H
#define MURMUR_PLATFORM_H

  void SetAffinity ( int cpu );
  #include <stdint.h>
  #include <strings.h>

  #define	FORCE_INLINE __attribute__((always_inline))

  FORCE_INLINE uint32_t rotl32 ( uint32_t x, int8_t r ){ return (x << r) | (x >> (32 - r)); }
  FORCE_INLINE uint64_t rotl64 ( uint64_t x, int8_t r ){ return (x << r) | (x >> (64 - r)); }
  FORCE_INLINE uint32_t rotr32 ( uint32_t x, int8_t r ){ return (x >> r) | (x << (32 - r)); }
  FORCE_INLINE uint64_t rotr64 ( uint64_t x, int8_t r ){ return (x >> r) | (x << (64 - r)); }
  // This isn't referenced anywhere, and breaks for non x86 architectures.
  //FORCE_INLINE unsigned long long int rdtsc(){ unsigned long long int x; __asm__ volatile ("rdtsc" : "=A" (x)); return x; }

  #define	ROTL32(x,y)	rotl32(x,y)
  #define ROTL64(x,y)	rotl64(x,y)
  #define	ROTR32(x,y)	rotr32(x,y)
  #define ROTR64(x,y)	rotr64(x,y)
  #define BIG_CONSTANT(x) (x##LLU)
  #define _stricmp strcasecmp

#endif 
/* MURMUR_PLATFORM_H */

#ifdef __cplusplus
  extern "C" {
#endif

uint32_t MurmurHash2        ( const void * key, int len, uint32_t seed );
uint32_t MurmurHash2A       ( const void * key, int len, uint32_t seed );
uint32_t MurmurHashNeutral2 ( const void * key, int len, uint32_t seed );
uint32_t MurmurHashAligned2 ( const void * key, int len, uint32_t seed );
uint64_t MurmurHash64A      ( const void * key, int len, uint64_t seed );
uint64_t MurmurHash64B      ( const void * key, int len, uint64_t seed );

#ifdef __cplusplus
  }
#endif


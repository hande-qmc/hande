! Common file to set defaults for definitions not set.

! Default to using 32 bit integers for bit strings for storing determinants.
#if !defined(DET_SIZE)
#warning "DET_SIZE is not defined.  I assumed DET_SIZE is 32."
#define DET_SIZE 32
#endif

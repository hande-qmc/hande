add_library(external OBJECT
  MurmurHash2.c
  MurmurHash2.h
  )

add_library(HANDE::external ALIAS external)

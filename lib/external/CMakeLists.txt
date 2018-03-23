add_library(external OBJECT
  MurmurHash2.c
  MurmurHash2.h
  )

target_sources(external
  PUBLIC
    $<$<BOOL:${USE_LANCZOS}>:trl_info.f90>
  )

add_library(HANDE::external ALIAS external)

add_library(external
  OBJECT
    MurmurHash2.c
    MurmurHash2.h
  )

target_sources(external
  PUBLIC
    $<$<BOOL:${USE_LANCZOS}>:trl_info.f90>
  )

set_target_properties(external
  PROPERTIES
    POSITION_INDEPENDENT_CODE 1
  )

add_library(HANDE::external ALIAS external)

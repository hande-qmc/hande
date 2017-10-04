add_library(external
  OBJECT
    MurmurHash2.c
    MurmurHash2.h
  )
if(NOT USE_LANCZOS)
  target_sources(external
    PUBLIC
      trl_info.f90
    )
endif()
set_target_properties(external
  PROPERTIES
    POSITION_INDEPENDENT_CODE 1
  )
add_library(HANDE::external ALIAS external)

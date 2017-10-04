add_subdirectory(aotus)
add_subdirectory(cephes)
add_subdirectory(external)
add_subdirectory(local)

# Compile dSFMT and dSFMT_F03_interface
list(APPEND dSFMT_sources
  dSFMT-src-2.2/dSFMT-common.h
  dSFMT-src-2.2/dSFMT-params.h
  dSFMT-src-2.2/dSFMT-params11213.h
  dSFMT-src-2.2/dSFMT-params1279.h
  dSFMT-src-2.2/dSFMT-params132049.h
  dSFMT-src-2.2/dSFMT-params19937.h
  dSFMT-src-2.2/dSFMT-params216091.h
  dSFMT-src-2.2/dSFMT-params2203.h
  dSFMT-src-2.2/dSFMT-params4253.h
  dSFMT-src-2.2/dSFMT-params44497.h
  dSFMT-src-2.2/dSFMT-params521.h
  dSFMT-src-2.2/dSFMT-params86243.h
  dSFMT-src-2.2/dSFMT.c
  dSFMT-src-2.2/dSFMT.h
  )
list(APPEND dSFMT_F03_interface_sources
  dSFMT_F03_interface/dSFMT_interface.F90
  dSFMT_F03_interface/dSFMT_utils.c
  )
add_library(dSFMT
  OBJECT
    ${dSFMT_sources}
    ${dSFMT_F03_interface_sources}
  )
target_compile_definitions(dSFMT
  PRIVATE
    -DDSFMT_MEXP=${DSFMT_MEXP}
  )
target_include_directories(dSFMT
  PRIVATE
    dSFMT-src-2.2
  )
add_library(HANDE::dSFMT ALIAS dSFMT)

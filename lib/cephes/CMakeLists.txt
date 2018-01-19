add_library(cephes
  OBJECT
    chbevl.c 
    mconf.h
    mtherr.c
    polevl.c
    protos.h
    psi.c
  )
set_target_properties(cephes
  PROPERTIES
    POSITION_INDEPENDENT_CODE 1
  )
add_library(HANDE::cephes ALIAS cephes)

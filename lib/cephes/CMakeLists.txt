add_library(cephes OBJECT
  chbevl.c
  mconf.h
  mtherr.c
  polevl.c
  protos.h
  psi.c
  )
add_library(HANDE::cephes ALIAS cephes)

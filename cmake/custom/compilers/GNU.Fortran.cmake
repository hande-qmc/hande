# Overrides contents of all variables previously set by CMake
if(NOT DEFINED ENV{FCFLAGS})
  if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fimplicit-none -fautomatic -fmax-errors=5")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -fbacktrace -Wall")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -funroll-all-loops -ftree-vectorize")
  endif()
endif()

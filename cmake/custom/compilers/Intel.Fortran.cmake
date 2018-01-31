if(NOT DEFINED ENV{FCFLAGS})
  if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
    add_definitions(-DVAR_IFORT)
    set(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -w -fpp -assume byterecl -traceback -nosave")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -warn all")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -ip")
  endif()
endif()
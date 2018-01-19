# Checks Fortran and C/C++ compilers compatibility, relying on the
# FortranCInterface CMake module.

# 12.5.0 broken (LAB), above should trap pre-Mavericks
set(BROKEN_MACOSX FALSE)
if((${CMAKE_SYSTEM_NAME} MATCHES "Darwin") AND (CMAKE_SYSTEM_VERSION VERSION_LESS 13))
    set(BROKEN_MACOSX TRUE)
endif()

if(BROKEN_MACOSX)
    set(CMAKE_C_FLAGS_HOLD ${CMAKE_C_FLAGS})
    set(CMAKE_CXX_FLAGS_HOLD ${CMAKE_CXX_FLAGS})
    message(STATUS "Trying Fortran name mangling appending -m32 flag")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -m32")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m32")
endif()

include(FortranCInterface)
FortranCInterface_VERIFY()

if(BROKEN_MACOSX)
    set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS_HOLD})
    set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS_HOLD})
endif()

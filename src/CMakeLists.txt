add_library(hande-shared SHARED
            $<TARGET_OBJECTS:HANDE::cephes>
            $<TARGET_OBJECTS:HANDE::coreutils>
            $<TARGET_OBJECTS:HANDE::dSFMT>
            $<TARGET_OBJECTS:HANDE::external>
            )

target_link_libraries(hande-shared
  PUBLIC
    HANDE::aotus
  )

add_subdirectory(ccmc)
add_subdirectory(core)
add_subdirectory(determinants)
add_subdirectory(dmqmc)
add_subdirectory(fci)
add_subdirectory(hellmann_feynman_sampling)
add_subdirectory(importance_sampling)
add_subdirectory(lua_bindings)
add_subdirectory(qmc)
add_subdirectory(system)

set_target_properties(hande-shared PROPERTIES SOVERSION 1 #${PROJECT_VERSION_MAJOR}
                                            OUTPUT_NAME "hande"
                                            EXPORT_NAME "hande")

# TODO Make the HANDE library discoverable
#install(TARGETS hande-shared
#        EXPORT "${PROJECT_NAME}Targets-shared"
#        RUNTIME DESTINATION ${CMAKE_INSTALL_LIBDIR}
#        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

# TODO This needs to be rewritten for MKL and ScaLAPACK
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

add_executable(hande.x hande.f90)
target_link_libraries(hande.x
  PRIVATE
    hande-shared
    HANDE::aotus
  PUBLIC
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
  )

if(USE_LIBUUID)
  target_link_libraries(hande.x
    PUBLIC
      ${UUID_LIBRARIES}
    )
endif()

if(USE_HDF5)
  target_link_libraries(hande.x
    PUBLIC
      ${HDF5_Fortran_LIBRARIES}
    )
endif()

target_sources(hande-shared
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/ct_fciqmc.F90
    ${CMAKE_CURRENT_LIST_DIR}/fci_lanczos.F90
    ${CMAKE_CURRENT_LIST_DIR}/fci_lapack.F90
    ${CMAKE_CURRENT_LIST_DIR}/fci_utils.F90
    ${CMAKE_CURRENT_LIST_DIR}/fciqmc.f90
    ${CMAKE_CURRENT_LIST_DIR}/ifciqmc.f90
    ${CMAKE_CURRENT_LIST_DIR}/semi_stoch.F90
    ${CMAKE_CURRENT_LIST_DIR}/simple_fciqmc.f90
)

if(NOT USE_TRLAN)
  target_compile_definitions(hande-shared
    PUBLIC
      DISABLE_LANCZOS
    )
else()
  target_sources(hande-shared
    PUBLIC
      ${CMAKE_CURRENT_LIST_DIR}/fci_lanczos.F90
  )
endif()

if(NOT USE_HDF5)
  target_compile_definitions(hande-shared
    PUBLIC
      DISABLE_HDF5
    )
endif()
